//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Implementation of Flutcurve: tracing flutter curves using a
// pseudo-arclength continuation method (pac.h/pac.cpp)
// Public:
//   trace_curves (Fstack& curves, int threadno)

#include <cmath>
#include <mutex>
#include <numeric>
#include <thread>

#include "config.h"
#include "concurrent.h"
#include "Exception.h"
#include "exim.h"
#include "fio.h"
#include "flutcurve.h"
#include "fma.h"
#include "lapack.h"
#include "matrix.h"
#include "nrc.h"
#include "pac.h"
#include "process.h"
#include "qr.h"
#include "settings.h"
#include "specs.h"
#include "target.h"
#include "trace.h"

using namespace std;

static double
precession(pset& plt, Matrix* gyro);

ostream&
operator<<(ostream& s, const Ivar& t) {
	s << t.par->name << " = " << t.par->value();
	if (!t.par->conv.eunits.empty())
		s << " (" << t.par->conv.eunits << ")";
	s << "*" << t.cu2eu();
	return s;
}

double
Cfo::
operator()(std::vector<double> const& x, std::vector<double>& c) {
// Cfo: Constraint function object, for setting Pac::confcn
		return fcv->constraintfcn(target, x, c);
}

double
Flutcurve::
constraintfcn(Target* t, const vector<double> & x, vector<double>& c) {
	// constraint function for holding "par" to "value"
	T_(Trace trc(2,"constraintfcn");)
	T_(trc.dprint("target ",*t);)
	// XXX update from x?
	std::fill(c.begin(), c.end(), 0.0);
	Par* par = findp(t->parname);
	for (size_t i=0; i<ad2ivar.size(); i++) {
		int idx = ad2ivar[i];	// 0b element # in x, tan
		if (idx < 0) continue;
		double dt = par->deriv(i);
		c[idx] = dt*ivar[idx].cu2eu(); // convert dt from EU to CU
	}
	// return p - \hat{p} = f[nf]
	double rval = par->value() - t->value;
	T_(trc.dprint("returning ",rval," = ",par->value()," - ",t->value);)
	T_(trc.dprint("returning c = ",c);)
	return rval;
}

// XXX make this inherit from Cfo
class TanConstraintfcn {
	private:
		vector<double> tan;
	public:
		TanConstraintfcn(vector<double> const& t) : tan(t) {};
		double operator()(vector<double> const& x, vector<double>& c) {
			// constraint function for forcing interates to be orthogonal to tan
			T_(Trace trc(2,"TanConstraintfcn");)
			c = tan;
			double rval = 0.0;
			T_(trc.dprint("returning ",rval);)
			return rval;
		}
};

Target*
find_target(vector<Target*> targets, Par* par) {
// given a set of Targets and a parameter, search the
// targets for a limit parameter at the current par value
	T_(Trace trc(2, "find_target");)
	double val = par->value();
	T_(trc.dprint("searching for a target for ",par->name," at ",val);)
	for (auto t : targets) {
		if (t->parname == par->name && is_equal(val, t->value, 6)) {
			T_(trc.dprint("returning target: ",*t);)
			return t;
		}
	}
	T_(trc.dprint("no target found");)
	return nullptr;
}

Flutcurve::
Flutcurve(const string& aid, const string& cid, const string& vzid) : Curve(aid,cid,vzid) {
// Flutcurve constructor: the Curve constructor copies gpset to "params",
// we create the origin from that after creating ivar_
	T_(Trace trc(2,"Flutcurve constructor");)

	T_(trc.dprint("aid <",aid,"> cid<",cid,"> vzid<",vzid,">");)

	// clear the solns arrays
	this->clear_solns();

	// set default and user-specified options
	this->specs = flutspecs();

	// set up the ivar_ vector in the order the independent variables
	// will be in. First the eigenvector components...
	evstart_ = 0;
	vector<Eigpar*>& eigpar = params.eigv();
	for (auto ep : eigpar)
		ivar.push_back({ep,1.0});	// cu2eu=1 for eigenvector elements
	ncev = ivar.size()/2;  // number of complex eigenvector elements

	// ... then any indep that are not eigenvector components.
	vector<Par*> indep = params.get_indep();
	for (auto pp : indep)
		ivar.push_back({pp,cu2eu(pp)});
	T_(trc.dprint("ivar: ",ivar);)

	// create workspace for use in pac_fjac
	adev = vector<complex<Ad>>(ncev);
	dynmat = vector<complex<Ad>>(ncev*ncev);
	work = vector<complex<Ad>>(ncev*ncev);
	Dy = vector<complex<Ad>>(ncev);

	int nx = get_nx();
	int nf = 2*get_order() + 2;

	// create the start point
	vector<double> start;
	for (auto& iv : ivar)
		start.push_back(iv.par->value()/iv.cu2eu());

	// set "phase_index" - the 0b index in origin.x of the component of the
	// eigenvector that is to be held to zero, i.e. the imaginary part
	// of the largest component of the complex eigenvector
	int indexreal;  // 0b index of the largest component of the complex vector
	normalize(ncev, (complex<double>*)&start[evstart_], indexreal);
	this->phase_index_ = evstart_ + indexreal*2 + 1;

	// the tangent only needs to point in the approximate direction desired;
	// check independents for limits, point away from them
	vector<double> tan(nx, 0.0);
	string atlim;
	Constraintfcn con{nullptr};
	for (auto pp : indep) {
#ifdef NEVER // use sp.targets
		double val = pp->value();
		if (pp->has_min() && is_equal(val, pp->min(),4)) {
			int idx = index(pp->name);
			tan[idx] = 1.0;
			cfo = new Cfo(Target(pp->name,pp->min()), this);
			con = *cfo;
		}
		if (pp->has_max() && is_equal(val, pp->max(),4)) {
			int idx = index(pp->name);
			tan[idx] = -1.0;
			cfo = new Cfo(Target(pp->name,pp->max()), this);
			con = *cfo;
		}
#else // NEVER // use sp.targets
		Target* t = find_target(this->specs.targets, pp);
		if (t != nullptr) {
			cfo = Cfo(t, this);
			con = cfo;
			// set tan to point towards opposite limit
			int idx = index(pp->name);
			if (t->is_lowerlimit)
				tan[idx] = 1.0;
			else
				tan[idx] = -1.0;
		}
#endif // NEVER // use sp.targets
	}
	// if no independents are an obvious choice for tan, set
	// it to a random vector
	if (blas::snrm2(nx, tan.data(), 1) == 0.0) {
		for (int i=0; i<nx; i++)
			tan[i] = flaps::xrand();
		blas::normalize2(nx, tan.data(), 1, true);
		T_(trc.dprint("random start tan: ",tan);)
	}

	// create the origin with a lambda fjac which calls
	// pac_fjac and possibly a constraint function (con)
	// Note: using a lambda is necessary because I can't use
	// member function pac_fjac directly and pac_fjac has to
	// be a member of Flutcurve
	origin = Pac(start, tan, nf,  [&](const vect& x, vect& f, vect& jac) {
			return pac_fjac(x, f, jac);
		});
	// give the origin a zero stepsize so it corrects it and
	// computes the tangent
	origin.step.hk = 0.0;
	origin.confcn = con;

	// Create projectee from project, projdir
	if (nx-nf > 1) {
		Specs& sp = this->specs;
		sp.projectee = vect(nx, 0.0);
		if (sp.project.empty())
			throw runtime_error("an optimization parameter must be specified "
				"when there are more than 3 independent parameters");
		for (size_t j=0; j<sp.project.size(); j++) {
			int i = index(sp.project[j]);
			sp.projectee[i] = sp.projdir[j];
		}
		T_(trc.dprint("projectee: ",sp.projectee);)
	}

	// create vector ad2ivar: adnames[i] -> x[ad2ivar[i]]
	ad2ivar.clear();
	for (auto& nm : Ad::adnames())
		ad2ivar.push_back(index(nm));
	T_(trc.dprint("ad2ivar: ",ad2ivar);)

	this->vzid(this->specs.vzid);	// set Curve::vzid for viz

}  // Flutcurve(aid,cid) constructor

Flutcurve::
Flutcurve(const string& aid, const string& cid, const string& vzid, Pac& orig) :
	Curve(aid,cid,vzid), origin(orig) {
// Flutcurve origin constructor: the Curve constructor copies gpset to
// "params", "orig" is the origin 
	T_(Trace trc(2,"Flutcurve origin constructor");)

	T_(trc.dprint("aid <",aid,"> cid<",cid,"> vzid<",vzid,">");)

	// clear the solns arrays
	this->clear_solns();

	// set default and user-specified options
	this->specs = flutspecs();

	// set up the ivar vector in the order the independent variables
	// will be in. First the eigenvector components...
	evstart_ = 0;
	vector<Eigpar*>& eigpar = params.eigv();
	for (auto ep : eigpar)
		ivar.push_back({ep,1.0});	// cu2eu=1 for eigenvector elements
	ncev = ivar.size()/2;  // number of complex eigenvector elements

	// ... then any indep that are not eigenvector components.
	vector<Par*> indep = params.get_indep();
	for (auto pp : indep)
		ivar.push_back({pp,cu2eu(pp)});
	T_(trc.dprint("ivar: ",ivar);)

	if (ivar.size() != origin.x.size())
		throw runtime_error(vastr("attempt to construct a Flutcurve with nx(origin) ",
			origin.x.size(), ", ivar size ",ivar.size()));

	// create workspace for use in pac_fjac
	adev = vector<complex<Ad>>(ncev);
	dynmat = vector<complex<Ad>>(ncev*ncev);
	work = vector<complex<Ad>>(ncev*ncev);
	Dy = vector<complex<Ad>>(ncev);

	// set "phase_index" - the 0b index in origin.x of the component of the
	// eigenvector that is to be held to zero, i.e. the imaginary part
	// of the largest component of the complex eigenvector
	int indexreal;  // 0b index of the largest component of the complex vector
	normalize(ncev, (complex<double>*)&origin.x[evstart_], indexreal);
	this->phase_index_ = evstart_ + indexreal*2 + 1;

	// create vector ad2ivar: adnames[i] -> x[ad2ivar[i]]
	ad2ivar.clear();
	for (auto& nm : Ad::adnames())
		ad2ivar.push_back(index(nm));
	T_(trc.dprint("ad2ivar: ",ad2ivar);)

	this->vzid(this->specs.vzid);

}  // Flutcurve(aid,cid,origin) constructor

void
Flutcurve::
update(const Pac& v) {
// update this->params, from v.x, and e9n data from v.step if requested
	T_(Trace trc(2,"update Pac");)

	if (v.x.size() != ivar.size())
		throw runtime_error(vastr("attempt to update with ",v.x.size()," variables"));
	
	// update each independent parameter with v.x converted from CU to EU
	//!! double* xi = v.x.data();
	//!! for (auto& iv : ivar)
		//!! iv.cuvalue(*xi++);
	this->update(v.x);

	// auxilliary parameters: coord, 
	Par* coord = this->params.findp("coord");
	if (coord != nullptr)
		coord->value(v.coord);
	// .. from e9n data if requested
	if (specs.e9n)
		e9n_update(v);
}

void
Flutcurve::
update(const vect& x) {
// update this->params, from x...

	const double* xi = x.data();
	for (auto& iv : ivar)
		iv.cuvalue(*xi++);

	// mark all parameters out-of-date...
	this->params.touch();
	// ... then set the autodiff derivative parameters derivative spot to 1,
	// and evaluate all derived parameters
	this->params.eval();
}

void
Flutcurve::
e9n_update(Pac const& pac) {
// Update my exploration (e9n) parameters
	T_(Trace trc(2,"Flutcurve::e9n_update");)
	const vector<pair<string,string>>& par = pac.step.parnames;
	double const* data = pac.step.cbegin();
	for (auto pi : par) {
		Par* pp = this->params.findp(std::get<0>(pi));
		if (pp != nullptr) {
			pp->value(*data++);
			T_(trc.dprint(pp->name," = ",pp->value());)
		}
	}
	// also update from Pac::stepsize...
	Par* pp = this->params.findp("stepsize");
	if (pp != nullptr)
		pp->value(pac.stepsize);
	// ... from Pac::projnorm
	pp = this->params.findp("projnorm");
	if (pp != nullptr)
		pp->value(pac.projnorm);
	// ... from Pac::specs.minf
	pp = this->params.findp("ecc");
	if (pp != nullptr) {
		double v = pac.specs.minf.empty() ? -1.0 : 1.0;
		pp->value(v);
	}
}

/*------------------------------------------------------------------
 * Local (static) routines
 *------------------------------------------------------------------*/

double
vuelta (const vector<double>& a, const vector<double>& b,
	const vector<double>& x, double& ratio);

int
Flutcurve::
processfcn(Pac& from, Pac& to, bool& looped, Issue* issue) {
// function to process results from Pac::trace
	T_(Trace trc(1,"processfcn coord ",to.coord);)
	Specs& sp = this->specs;
	int rval{1};		// default: continue

	bool was_constrained{false};
	if (to.confcn != nullptr) {
		was_constrained = true;
		Target* tp = cfo.target;
		T_(trc.dprint("solution was constrained: ",*tp);)
		if (from.coord != 0) {
			if (tp->is_lowerlimit || tp->is_upperlimit) {
				Par* pp = this->findp(tp->parname);
				assert (pp != nullptr);
				if (is_equal(pp->value(), tp->value, 3)) {
					string msg{vastr(tp->parname," reached ",tp->value)};
					this->atlimit = msg;
					add2msg(this->finished, msg);
					T_(trc.dprint("returning 0: ",this->atlimit);)
					rval = 0;	// quit
				}
			}
		} else {
			T_(trc.dprint("ignoring constraint: coord 0");)
		}
		to.confcn = nullptr;
	}

	// deal with issues from pac - most issues are dealt with in pac
	if (issue != nullptr) {
		T_(trc.dprint("got issue \"",issue,"\"");)
		if (issue->is_angle()) {
#ifdef NEVER // repeat step
			flaps::warning("possible discontinuity in ",cid(),
					" encountered at coordinate ",to.coord,", angle=",to.step.alphak);
			return 1;	// continue with large angle - a discontinuity?
#else // NEVER // repeat step
			if (discon == nullptr) {
				discon = new Discontinuity(to.coord,to.step.alphak);
			} else if (discon->angles.size() > 4) {
				flaps::warning("possible discontinuity in ",cid(),
					" encountered at coordinate ",to.coord,", angles=",discon->angles);
				delete discon;
				discon = nullptr;
				return 1;	// continue with large angle - a discontinuity?
			} else {
				// check that coords are the same?
				discon->angles.push_back(to.step.alphak);
			}
			// reduce the stepsize and repeat the step
			from.step.hk /= 3.0;
			from.step.red_angle += 1;	// XXX this won't get into the plotfile
			return -1;	// repeat step
#endif // NEVER // repeat step
		} else if (issue->is_noconv() || issue->is_divergence() || issue->is_slow()) {
			flaps::warning(cid()," coord ",from.coord,": ",issue->msg);
			return 0;
		} else if (issue->is_zerotan()) {
			// ignore a small tan for coord==0 unless this is optimization
			if (!sp.projectee.empty())
				this->atlimit = issue->msg;
			else if (from.coord == 0)
					return 1;
			return 0;		// can't continue
		} else {
			throw runtime_error(vastr("curve ",cid(),", coord ",from.coord,
				", pac issue \"",issue->msg,"\" not treated yet"));
		}
	}

	// was this step successful?
	if (to.step.conv_crit < 0) {
		flaps::warning(cid(),": no convergence at ",params.summary(sp.toprint));
		return 0;
	}
	// check for bifurcation...
	if (sp.bifurcation) {
		if (checkdsc(from, to) < 0) {
			T_(trc.dprint("retry step: bif near");)
			return -1;
		}
	}
	// ... compute lcostab
	if (sp.lcostab)
		lcostability(to);

	// update params: updating from a Pac also updates e9n parameters
	// this must be done prior to saving solns
	this->update(to);
	// save the solution on my "solns" array
	if (to.coord >= 0.0)
		back_solns();
	else
		front_solns();
	T_(trc.dprint("now have ",nsolns()," solutions");)

	// if this is coord 0 skip the rest & return 1 (sucess)
	if (to.coord == 0.0) {
		T_(trc.dprint("return 1: coord 0");)
		return 1;
	}
	// check that *all* parameters (except rf) are in range
	vector<string> oor = this->inrange();
	if (!oor.empty()) {
		for (auto& pi : oor) {
			if (pi == "rf" && !sp.rflimits)
				continue;
			this->atlimit = vastr(oor," is out of limits");
			add2msg(this->finished, this->atlimit);
		}
		if (!this->atlimit.empty()) {
			T_(trc.dprint("returning: ",this->atlimit,", error<",this->error);)
			return 0;	// atlimit
		}
	}
	// check for limits/targets needing constraint on "to" (to.confcn)
	// unless "to" was constrained
	if (!was_constrained) {
		string atlim;			// XXX why not this->atlimit?
		this->cfo = pac_constraint(to, atlim);
		if (this->cfo.fcv != nullptr)
			to.confcn = this->cfo;
		if (!atlim.empty()) {
			T_(trc.dprint("returning 0: ",atlim);)
			this->atlimit = atlim;
			add2msg(this->finished, atlim);
			return 0;
		}
	}
	// check for looping
	if (this->specs.checklooping && abs(to.coord) > 4) {
		double ratio{0.0};
		double t = vuelta(from.x, to.x, origin.x, ratio);
		if (0.0 <= t && t <= 1.0 && ratio < 0.01) {
			looped = true;
			add2msg(this->finished, "curve is looping");
			// add the origin as the next and final point
			this->origin.coord = to.coord;
			update(this->origin);
			if (to.coord >= 0.0)
				back_solns();
			else
				front_solns();
			T_(trc.dprint("closed the loop: now have ",nsolns()," solutions");)
			return 0;
		}
	}
	// check Jacobian?
	if (sp.checkjacobian > 0) {
		if (to.coord != 0 && abs(to.coord)%sp.checkjacobian == 0) {
			Issue reason(vastr("user-requested Jacobian check at coord ",
					to.coord), Kind::user_request);
			Issue* res = to.check_jac(nullptr, nullptr, &reason);
			if (res != nullptr)
				flaps::info(res->msg);
		}
	}
	return rval;
}	// processfcn

void
trace_curves (Fstack* curves, int threadno) {
// Given a stack of Flutcurves, get the next available Flutcurve by
// calling curves->pop(), then extend the curve starting from the
// Curve's origin
	T_(Trace trc(1,"trace_curves ",threadno);)
	int ntrace{0};

	try {
		// pop Flutcurves until nullptr returned
		Flutcurve* curve;
		while ((curve = curves->pop()) != nullptr) {
			bool looped{false};
			string curveid = curve->cid();
			Specs& sp = curve->specs;
			ntrace++;

			T_(trc.dprint("thread ",threadno," tracing ",curveid);)
			curve->thread = threadno;

			// if the origin is at a limit (as evidenced by a constraint)
			// we don't need to trace in reverse
			bool reverse{true};
			if (curve->origin.confcn != nullptr)
				reverse = false;

			if (curve->specs.pacspecs != nullptr)
				curve->origin.specs = *curve->specs.pacspecs;
			// start at coord 0
			curve->origin.coord = 0;

			vect* projectee{nullptr};
			if (!sp.projectee.empty())
				projectee = &sp.projectee;
			// Trace the curve from "origin" until Processfcn returns false.
			// Starting with the predicted stepsize (step.hk) =0 causes Pac::trace
			// to correct origin and compute an initial tangent; then if processfcn
			// returns true it will continue tracing the curve, otherwise it returns
			// here and we trace in the opposite direction
			// Use a lambda as the process function, which just calls processfcn
			curve->origin.step.hk = 0.0;
			curve->origin.trace([&](Pac& from, Pac& to, Issue* issue) {
				return curve->processfcn(from, to, looped, issue);}, projectee);

			if (looped) {
				T_(trc.dprint(curve->cid(),": not tracking in reverse: first direction looped");)
				reverse = false;
			}

			// reverse the tangent & projectee & track again
			if (!reverse) {
				if (curve->params.nsolns() == 0)
					flaps::warning(vastr(curveid,": first tracking direction yielded ",
						"no results: check limits"));
			} else {
				blas::scal(curve->origin.tan.size(), -1.0, &curve->origin.tan[0], 1);
				if (!curve->specs.projectee.empty())
					blas::scal(curve->specs.projectee.size(), -1.0,
						&curve->specs.projectee[0], 1);
				curve->origin.step.hk = curve->origin.specs.initial_stepsize;
				//!! curve->finished = "";	// XXX no: accumulate both directions
				curve->error = "";
				curve->atlimit = "";
				// coords go negative starting at origin.coord (0)
				int direction = -1;
				// trace the curve from "origin" but this time with a non-zero
				// stepsize since "origin" has already been corrected and given a tan
				curve->origin.trace([&](Pac& from, Pac& to, Issue* issue) {
						return curve->processfcn(from, to, looped, issue);
					}, projectee, direction);
			}
			if (!curve->origin.specs.minf.empty())
				cerr << curve->cid() << " used ECC\n";
		}
	} catch (runtime_error& s) {
		flaps::error(s.what());
	}

	// T_(trc.dprint("thread ",this_thread::get_id()," traced ",ntrace," curves");)
	T_(trc.dprint("thread ",threadno," traced ",ntrace," curves");)
	return;
}  // trace_curves

double
vuelta (const vector<double>& a, const vector<double>& b,
	const vector<double>& x, double& ratio) {
// Test to see if going from a to b passes x, a known solution
// on the curve (e.g. the origin).
// Let S = b - a (secant), then find a point c = a + tS on the secant
// such that the distance x-c is minimal. Let
//   f(t) = (x-c)'(x-c) = x'x - 2x'c + c'c
// then the min is where
//   df/dt = 0 = -2x'dc/dt + 2(dc/dt)'c = -2x'S + 2(a+tS)'S
// so
//   t = (x - a)'S/S'S
	T_(Trace trc(2,"vuelta");)

	int nx = x.size();
	vector<double> S{b};
	blas::axpy(nx, -1.0, a.data(), 1, S.data(), 1);
	vector<double> xma{x};
	blas::axpy(nx, -1.0, a.data(), 1, xma.data(), 1);
	double SS;
	blas::dot(nx, S.data(), 1, S.data(), 1, SS);
	double xmaS;
	blas::dot(nx, xma.data(), 1, S.data(), 1, xmaS);
	double t = xmaS/SS;
	// if the minimizing point (t) is on the secant check for small f(t)
	if (0.0 < t && t <= 1.0) {
		vector<double> xmc{x};
		blas::axpy(nx, -1.0, a.data(), 1, xmc.data(), 1);
		blas::axpy(nx, -t, S.data(), 1, xmc.data(), 1);
		double fnorm = blas::snrm2(nx, xmc.data(), 1);
		double xnorm = blas::snrm2(nx, x.data(), 1);
		ratio = fnorm/xnorm;
	}
	T_(trc.dprint("returning t ",t,", ratio ",ratio);)
	return t;
}

void
Flutcurve::
normalize (int n, complex<double>* cv, int& indexReal) {
/*----------------------------------------------------------------
 * Scales a complex vector so that the 2-norm is
 * 1.0 (i.e. x(*)x = 1.0) and element "indexReal"
 * (0-based) is positive and real.
 * The condition imag(cv[indexReal]) = 0.0 is required
 * for a unique normalization.
 *
 * input:
 *  n         number of elements in vector cv
 *  cv        (n complex) complex vector
 * output:
 *  cv        on output cv has been scaled so that cv(*)*cv = 1.0
 *            and cv(indexReal) is real.
 *  indexReal   index (0-based) of the largest element in the vector
 *----------------------------------------------------------------*/
	T_(Trace trc(2,"normalize");)
	double eps = std::numeric_limits<double>::epsilon();

	indexReal = blas::icamax (n, cv, 1) - 1;  // 0b

	// Compute the 2-norm of the vector
	double norm = blas::scnrm2 (n, cv, 1);

	T_(trc.dprint("norm ",norm,", indexReal ",indexReal);)

	// with cmax = cv[indexReal],
	// scale factor = abs(cmax)/(norm*cmax)
	//              = abs(cmax)*conjg(cmax)/(norm*abs(cmax)**2)
	//              = conjg(cmax)/(norm*abs(cmax))

	complex<double> cmax = cv[indexReal];

	double t = norm * abs(cmax);
	if (abs(t) < eps) {
		throw runtime_error("normalize() passed zero vector");
	}

	complex<double> scale = conj(cmax)/t;

	T_(trc.dprint("scale ",scale,", cmax ",cmax);)

	blas::scal (n, scale, cv, 1);
}

void
equalize_det(double& d1,double& exp1,double& d2,double& exp2) {
// check for different exponents on d1, d2
// replace d2 <- d2*10^(exp2-exp1) so det2 = d2*10^exp1
	T_(Trace trc(1,"equalize_det");)
	T_(trc.dprint("d1 = ",d1,"*10^",exp1,", d2 = ",d2,"*10^",exp2);)
	if (exp1 != exp2) {
		d2 *= pow(10.0, exp2-exp1);
		exp2 = exp1;
		T_(trc.dprint("adjusted d1 = ",d1,"*10^",exp1,", d2 = ",d2,"*10^",exp2);)
	}
}

double
Flutcurve::
cu2eu (Par* pp) {
// Compute a scale factor to go from continuation units (CU) to
// equation units (EU):
//   EU = CU*scale
// to give this parameter a range of approximately 2 in the continuation
// solution. If limits have not been specified for this parameter the
// scale factor will be 1.0, otherwise the scale factor will be the
// closest value of 2^n. The factor may also be specified in flutspecs.cu2eu.
	T_(Trace trc(2,"cu2eu");)
	double scale{1.0};
	Specs& sp = flutspecs();

	// check if the cu2eu is set as a spec
	for (auto& ci : sp.cu2eu) {
		if (ci.name == pp->name) {
			T_(trc.dprint("returning spec: ",pp->name," = ",ci.value);)
			return ci.value;
		}
	}
	// no limits?
	if (!(pp->has_max() && pp->has_min())) {
		static bool visit{false};
		if (!visit)
			flaps::warning("no limits given for ",pp->name," so it will not be "
				"scaled, which may result in less than ideal tracking");
		visit = true;
		return scale;
	}

	scale = std::abs(pp->max() - pp->min())/2.0;
#ifdef NEVER // no frexp
	// std::frexp returns scale = s*2^exp
	int exp;
	double s = std::frexp(scale, &exp);
	scale = pow(2.0, exp);
#endif // NEVER // no frexp
	T_(trc.dprint("scale factor for ",pp->name," based on min,max: ",scale);)
	return scale;
}

int
Flutcurve::
index(const string& parname) const {
// Returns the index (zero-based) in an Ivar of the
// parameter with the name "parname" or -1 if it is not present
	T_(Trace trc(2,"Flutcurve::index ",parname);)
	int rval{-1};

	int i{0};
	for (auto& iv : ivar) {
		if (iv.par->name == parname) {
			rval = i;
			break;
		}
		i++;
	}
	T_(trc.dprint("returning ",rval);)
	return rval;
}

// functions for printing Flutcurves
string
Flutcurve::
current_values() {
// returns a one-line summary of the current parameter values
// specified by policy "toprint"
	return this->params.summary(this->specs.toprint);
}

string
Flutcurve::
print_solns(const vector<string>& leaders, bool firstlast) {
// returns a summary of all the parameter values, one line
// per value, specified by policy "toprint"
// Note: this function cannot be const because summary_solns is not
	string header{make_header()};
	return this->params.summary_solns(this->specs.toprint, header, leaders, firstlast);
}

string
Flutcurve::
make_header() const {
// returns a string of parameter names to be printed with
// e.g. pset::summary() which prints each value with a field width of 14
	ostringstream os;
	for (auto s : this->specs.toprint) {
		const Par* pp = params.findp(s);
		// s may not be a parameter name, e.g. "eigenvector"
		if (pp != nullptr && !pp->conv.punits.empty())
			s += " (" + pp->conv.punits + ")";
		os << std::setw(14) << s;
	}
	return os.str();
}

string
Flutcurve::
summarize(vector<double> x) {
// write a summary of nx-vector x to a string
	if (x.size() != ivar.size()) {
		throw runtime_error(vastr("attempt to summarize a ",x.size(),"-vector",
					" in a curve with ",ivar.size()," indep variables"));
	}
	ostringstream os;
	for (size_t i=0; i<ivar.size(); i++)
		os << setw(10) << ivar[i].par->name << "  " << x[i] << endl;
	return os.str();
}

void
Flutcurve::
plot(const string& aid,const string& pf, const vector<Flutcurve*>& fcurves) {
// create 2 flaps plotfiles:
// 1) a compressed tarfile of "fcurves" plus data needed for amviz
// 2) an ASCII file with just "fcurves"
	T_(Trace trc(1,"Flutcurve::plot");)
	Specs& sp = fcurves[0]->specs;

	vector<string> toplot = sp.toplot;

	string basename = pf;
	auto dot = pf.rfind('.');
	if (dot != string::npos) {
		string ext = pf.substr(dot);
		if (ext == ".pf" || ext == ".apf")
			basename = pf.substr(0,dot);
	}
	T_(trc.dprint("plotting aid<",aid,">, base plotfile ",basename);)

	// 1) the ASCII plotfile is basename+.apf
	//    It will not have amviz data
	//    Include only parameters with valid State (not nostate)
	for (auto& pp : fcurves[0]->params.pmap())
		if (!pp.second->is_nostate() && !pp.second->is_eigv())
			toplot.push_back(pp.first);

	// Create the ASCII plotfile using Apf::exporter; first
	// convert fcurves to Curve*
	vector<Curve*> curves;
	for (auto i : fcurves)
		curves.push_back(i);
	Apf::exporter(curves, toplot, basename+".apf", false);

	// 2) the .pf file is just a Flaps savefile with a list
	//    of matrices in "tosave"
	//    if we are to exclude nostate parameters we would have to create
	//    new curves without them OR delete them in Flutcurve constructor
	vector<string> tosave;

	// XXX allow for new "period separators"
	tosave.push_back(vastr("curve[,.]",aid,"[,.].*"));

	// If a vzid was given in specs or gotten from the source run,
	// read the fma, check that ngc (# columns in gct) is equal to either
	// basic or n
	string vzid{sp.vzid};
	if (vzid.empty()) {
		flaps::info("amviz data will not be included in the plotfile: no vzid given");
	} else {
		try {
			Fma* fma = Fma::fetch(vzid);
			size_t basic;
			size_t n = get_order(&basic);
			size_t ngc = fma->gct.size()/fma->coords.size();
			if (ngc == n || ngc == basic) {
				tosave.push_back(Fma::mid(vzid));
			} else {
				string msg = vastr("vzid \"",vzid,"\" has ",ngc,
					" columns in gct, should have ", basic);
				if (n != basic)
					msg += vastr(" or ",n);
				else
					flaps::info("not including amviz kit \"",sp.vzid,"\" in the plot file: ",msg);
			}
		} catch (const std::exception& s) {
			flaps::info("amviz data will not be included in the plot file: ",s.what());
		}
	}

	// the plotfile is just basename with a .pf extension
	string file = basename+".pf";
	vector<string> saved = fio::save(tosave, file);
	T_(trc.dprint("created ",file," with:",saved);)
}

Fstack::
Fstack(vector<Flutcurve*>& curves) {
// Fstack constructor
	next = 0;
	for (auto fcp : curves) {
		fcp->params.clear_solns();
		this->push_back(fcp);
	}
}

Flutcurve*
Fstack::
pop() {
// return the next Flutcurve on the stack or nullptr if empty
	std::lock_guard<std::mutex> guard(this->mtx);
	if (next < size()) {
		return (*this)[next++];
	}
	return nullptr;
}

Cfo
Flutcurve::
pac_constraint(Pac& to, string& atlim) {
/*------------------------------------------------------------------
 * If there are targets within [x:x+h], where h = stepsize*tan,
 * decrease stepsize to hit the target and addd a Constraintfcn
 * to "to" for the corrector to ensure the target is hit.
 *------------------------------------------------------------------*/
	T_(Trace trc(2,"pac_constraint");)
	double meps(std::numeric_limits<double>::epsilon());
	int nx = to.x.size();
	Cfo rval{};	// default: no constraint

	// first check for user-requested constraint: force iterates to be
	// orthogonal to tangent
#ifdef NEVER // convert TanCon... to Cfo
	Specs& sp = this->specs;
	if (sp.constrained) {
		//!! to.confcn = TanConstraintfcn(to.tan);
		return TanConstraintfcn(to.tan);
	}
#endif // NEVER // convert TanCon... to Cfo

	if (this->specs.targets.empty()) {
		T_(trc.dprint("quick return: no targets");)
		return rval;
	}

	// Check each target's parameter to see if x+h exceeds
	// the target value (g); if so compute a constraint
	//  T = \Vector{\left\{\frac{\partial p}{\partial x_i} \right\}}
	// and a stepsize stepj = (g-p)/T^t h where p is the current
	// value of the parameter, such that x + sh yields the target.
	// Take the smallest such scale factor and the corresponding T
	double newstep{to.step.hk};  // stepsize needed to hit target
	Target* the_target{nullptr};	// closest target approaching

	for (auto& t : this->specs.targets) {
		// check Target.windows (EU)
		// XXX should update first but this will work *most* of the time
		bool oor{false};
		for (auto& wi : t->windows) {
			if (!wi.inrange(this->params)) {
				oor = true;
				break;
			}
		}
		if (oor) {
			T_(trc.dprint("skip target ",*t,", out of window");)
			continue;
		}
		vector<double> c(nx,0.0);
		double cfval = constraintfcn(t, to.x, c);
		double dcf;
		blas::dot(nx, c.data(), 1, to.tan.data(), 1, dcf);
		// cfval approx 0? (i.e. at limit) check which direction tan points
		if (is_equal(cfval, 0.0, 5)) {
			if ((t->is_lowerlimit && dcf <= meps) ||
				(t->is_upperlimit && dcf >= -meps)) {
				atlim = vastr(t->parname," is at a limit (",t->value,")");
				return rval;
			}
		}
		if (abs(dcf) < meps) {
			T_(trc.dprint("skip target ",*t,", dcf = ",dcf);)
			continue;
		}
		// note: cfval is p - \hat{p} but the stepsize is (\hat{p}-p)/dcf
		// see sect. 2.2.4 in the manual
		double stepj = -cfval/dcf;
		T_(trc.dprint("stepj = ",stepj," = -(",cfval,")/",dcf);)
		// is this a smaller stepsize?
		T_(trc.dprint("stepj ",stepj," < newstep? ",newstep);)
		if (stepj > 0.0 && stepj < newstep) {
			newstep = stepj;
   		the_target = t;
		}
		// ... or are we approaching this target?
		if (stepj > newstep && stepj-newstep < 0.1*newstep) {
			newstep /= 2.0;
			to.step.hk = newstep;
			T_(trc.dprint("approaching target ",*t," decreasing stepsize to ",newstep);)
		}
	}

	if (the_target == nullptr) {
		to.confcn = nullptr;
		T_(trc.dprint("no target found");)
		return rval;
	}

	// If newstep is between tol and 0.99*stepsize give "to" the Constraintfcn
	// functor and the new predicted stepsize (hk)
	// XXX if it is < 1.5 just scale h by 1/2
	double tol{to.specs.minstepsize};
	if (newstep > tol) {
		T_(trc.dprint("decreased stepsize from ",to.step.hk," to ",newstep," to hit ",*the_target);)
		to.step.hk = newstep;
		//!! to.confcn = Cfo(the_target, this);
		return Cfo(the_target, this);
	} else {
		// newstep small => we are at the target
		T_(trc.dprint("newstep ",newstep," too small");)
	}
	return rval;
} // pac_constraint

int
Flutcurve::
pac_fjac (vector<double> const& x, vector<double>& f, vector<double>& jac) {
// compute the function and associated Jacobian
// This is a non-static member function so it cannot be used as an Fjacfcn (pac.h)
// so instead we call this function from a lambda in calls to pac::trace
	T_(Trace trc(2,"pac_fjac ");)
	int rval{0};
	int nx = x.size();
	int nf = f.size();

	// initialize
	jac.assign(jac.size(), double(0.0));
	f.assign(f.size(), double(0.0));

	// Update my parameters from the input x
	this->update(x);

	// get the eigenvector as a complex AD array, put it in Flutcurve::adev
	this->params.getadev(this->adev);

	// sanity check on adev
	if (ncev != (int)adev.size())
		throw runtime_error(vastr("ncev from get_order (",
				ncev, ") does not agree with getadev (", adev.size(),")"));


	// ncev is the length of the complex eigenvector and the size of the
	// dynamic matrix; nrev is the length of the equivalent real vector (2*ncgc)
	int nrev = ncev*2;
	T_(trc.dprint("AD eigenvector:\n", adev);)

	// f should be nrev+2
	if (nf != nrev+2) {
		if (this->specs.constrained && nf != nrev+3)
			throw runtime_error(vastr("f is ",nf,", should be ",nrev+2));
	}
	// the equations must be underdetermined or square:
	if (nf > nx)
		throw runtime_error(vastr("nf(",nf,") > nx(",nx,"): too many constraints?"));

	// compute the dynamic matrix as a complex AD matrix
	dmatrix(this->params, this->dynmat, this->work);

	// multiply Dy
	blas::gemv("n", ncev, ncev, 1.0, &dynmat[0], ncev, &adev[0], 1,
		0.0, &Dy[0], 1);

	T_(trc.dprint("Dy",Dy);)

	// function f(x): copy the "value" part of each ad
	for (size_t i=0; i<Dy.size(); i++) {
		f[2*i] = Ad::real(Dy[i]).data()[0];
		f[2*i+1] = Ad::imag(Dy[i]).data()[0];
	}

	// Extract the dynamic matrix as a double vector, treat dynmatval as
	// a (nrev,ncev) real double matrix; this will be inserted into the
	// Jacobian as derivatives of f wrt linear eigenvector components
#ifdef NEVER // use extract
	std::vector<double> dynmatval(nrev*ncev);
	for (size_t i=0; i<dynmat.size(); i++) {
		dynmatval[2*i] = Ad::real(dynmat[i]).data()[0];
		dynmatval[2*i+1] = Ad::imag(dynmat[i]).data()[0];
	}
#else // NEVER // use extract
	std::vector<double> dynmatval(nrev*ncev);
	extract(dynmat, "", dynmatval.data());
#endif // NEVER // use extract

	// insert the dynamic matrix (values, not derivatives) into the
	// Jacobian in the columns corresponding to the eigenvector, i.e.
	// beginning with column evstart(); nonlinear eigenvector components
	// will replace columns inserted here with the correct derivatives
	blas::real_rep (ncev, ncev, (complex<double>*)&dynmatval[0], ncev,
		&jac[IJ(0,this->evstart(),nf)], nf);

	// Columns of the Jacobian (j = 0:nx-1).
	// Each column has a parameter associated with it: ivar[j].par
	// Two cases to consider depending on whether ivar[j] is
	// an AD parameter:
	// 1) Find which ad derivative parameter ivar[j] is, extract that from
	//    Dy, the ad-valued residual and insert it into column j of the Jacobian
	// 2) If ivar[j] is not an AD parameter it must be the real or imaginary
	//    part of a linear eigenvector element (nonlinear eigenvector components
	//    are AD parameters) and it was accounted for in real_rep().
	std::complex<double> ci(double(0.0), double(1.0));
	std::vector<double> df(nrev);
	vector<string> const& adpar = Ad::adnames();
	for (size_t j=0; j<adpar.size(); j++) {
		int idx = index(adpar[j]);
		if (idx >= 0) {
			for (size_t i=0; i<Dy.size(); i++) {
				df[2*i] = Ad::real(Dy[i]).data()[j+1];
				df[2*i+1] = Ad::imag(Dy[i]).data()[j+1];
			}
			blas::copy (nrev, &df[0], 1, &jac[IJ(0,idx,nf)], 1);
			T_(trc.dprint("copied partial Dy wrt ",adpar[j]," to column ",idx," 0b");)
		}
	}
	// The first nrev rows are now in equation units (EU); multiply them
	// by ivar.cu2eu to put them in continuation units (CU):
	//   p_{eu} = \alpha p_{cu} where \alpha = ivar.cu2eu
	// so
	//   \partial{f}/partial{p_{cu}} = \alpha \partial{f}/\partial{p_{eu}}
	for (int j=0; j<nx; j++) {
		double sf = this->ivar[j].cu2eu();
		if (sf != 1.0)
			blas::scal(nrev, sf, &jac[IJ(0,j,nf)], 1);
	}

	// next row (0b) of f and jac to insert stuff
	size_t nextrow = nrev;

	// The last two eqns are eigenvector normalization conditions.
	// Since x is in CU, whereas adev and rev are in EU, use eigenvector
	// components directly from x
	// 1) phase normalization: \Im{\hat{y_k}} = 0
	int start = this->evstart();
	const double* y = &x[start];
	int k = this->phase_index();  // index of the component of x to be zero
	T_(trc.dprint("phase normalization: f[",nextrow,"] = x[",k,"] = ",x[k]);)
	f[nextrow] = x[k];
	jac[IJ(nextrow,k,nf)] = 1.0;
	nextrow++;

	// 2) next equation is the 2-norm squared of the eigenvector y = x[start:]
	//      f = y^t y - 1.0
	//      df/dy = 2y^t
	double evnorm(1.0);	// eigenvector norm is always 1.0 in CU
	double evcurrent;
	blas::dot(nrev, y, 1, y, 1, evcurrent);
	if (evcurrent == double(0.0))
		throw runtime_error("zero eigenvector encountered");
	f[nextrow] = evcurrent - evnorm;     // 2-norm squared
	T_(trc.dprint("normalization: f[",nextrow,"] = ",f[nextrow], " = ", evcurrent," - ",evnorm);)
	// the eigenvector starts in column "start"
	double* jacev = &jac[IJ(nextrow,start,nf)];
	blas::copy(nrev, y, 1, jacev, nf);
	blas::scal(nrev, 2.0, jacev, nf);
	nextrow++;

	static int visit{0};
	if (visit++ == 0)
		Dprint::dprintm(nf,nx,nf,&jac[0],"Jacobian");

	return rval;
}  // pac_fjac

void
Flutcurve::
dmatrix (pset& plt, vector<complex<Ad>>& result,
		vector<complex<Ad>>& work) {
/*------------------------------------------------------------------
 * Computes the dynamic matrix at the current values of the
 * parameters in plt
 * Returns: an (n,n) complex matrix, the dynamic matrix
 * Input:
 *   plt      a pset of parameters to use in evaluating D
 *   result   (n,n) complex Ad matrix
 *   work     (n,n) complex Ad matrix
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"dmatrix");)
	complex<Ad> ctmp;
	size_t n = sqrt(result.size());
	int nsq = n*n;
	Specs& sp = flutspecs();

	// required parameters...
	Ad freq = plt.parval("freq");
	Ad sigma = plt.parval("sigma");
	Ad sdamp = plt.parval("sdamp");
	Ad dpress = plt.parval("dpress");

	complex<Ad> s(sigma, freq);

	complex<Ad>* wp = &work[0];
	complex<Ad>* rp = &result[0];

	// initialize result
	for (size_t i=0; i<result.size(); i++) {
		Ad::real(result[i]).zero();
		Ad::imag(result[i]).zero();
	}

	// Mass matrix: s^2M
	Matrix* cp = Matrix::find_desc("mass");
	if (cp) {
		ctmp = s*s;
		cp->eval(plt, work);
		blas::axpy(nsq, ctmp, wp, 1, rp, 1);
	}

	// Viscous damping matrix: s*V
	cp = Matrix::find_desc("vdamp");
	if (cp) {
		cp->eval(plt, work);
		// result.caxpy (s, work);
		blas::axpy(nsq, s, wp, 1, rp, 1);
	}

	// Gyro matrix: s*G
	cp = Matrix::find_desc("gyro");
	if (cp) {
		cp->eval(plt, work);
		blas::axpy(nsq, s, wp, 1, rp, 1);
		plt.setpar(precession(plt,cp), "precession");
	}

	// Stiffness matrix
	cp = Matrix::find_desc("stif");
	if (cp) {
		cp->eval(plt, work);
		complex<Ad> cdamp(1.0, sdamp);
		blas::axpy(nsq, cdamp, wp, 1, rp, 1);
	} else {
		throw runtime_error("no stiffness matrix");
	}

	// Aero matrix
	if (dpress.value() > 0.0) {
		cp = Matrix::find_desc("gaf");
		if (cp) {
			cp->eval(plt, work);
			complex<Ad> scale(-dpress);
			blas::axpy(nsq, scale, wp, 1, rp, 1);
		}
	}

	// Control-law matrix
	cp = Matrix::find_desc("controls");
	if (cp) {
		cp->eval(plt, work);
		complex<Ad> scale(1.0);
		blas::axpy(nsq, scale, wp, 1, rp, 1);
	}

	// print the dynamic matrix once if settings printmatrices is true
	static bool visit{false};
	if (!visit && sp.printmatrices) {
		vector<complex<double>> values(nsq);
		extract(result, "",&values[0]);
		string comment = vastr("sigma ",sigma.value(),", freq ",freq.value(),
			", dpress ",dpress.value());
		MM::exporter("dynamic_matrix.mm",comment,&values[0], n, n);
		visit = true;
	}
	T_(trc.dprintm(n,n,n,result,"dynamic matrix");)
}

vector<string>
Flutcurve::
dependson(pset& plt) {
// get a list of *all* parameters the dmatrix depends on
// and mark them active
	T_(Trace trc(1,"Flutcurve::dependson");)
	Specs& sp = flutspecs();
	int n = get_order();
	// set independents to their midpoints
	vector<Par*> indeps = plt.get_indep();
	for (auto pp : indeps)
		if (pp->has_min() && pp->has_max())
			pp->value((pp->min()+pp->max())/2.0);
	plt.eval();

	vector<complex<Ad>> result(n*n);
	vector<complex<Ad>> work(n*n);
	// turn on monitoring and evaluate the dmatrix...
	plt.monitor();
	Flutcurve::dmatrix (plt, result, work);
	// ... also get target, print, and plot parameters
	for (auto& ti : sp.targets)
		gpset::find(ti->parname);
	for (auto& pi : sp.toprint)
		gpset::find(pi);
	for (auto& pi : sp.toplot)
		gpset::find(pi);
	// rotinom returns the immediate dependents 
	vector<string> dep = plt.rotinom();
	// now deep-dive each immediate dependent
	vector<string> rval;
	for (auto ni : dep) {
      Par* pi = plt.findp(ni);
      if (pi == nullptr) continue;
		if (std::find(rval.begin(), rval.end(), ni) == rval.end())
			rval.push_back(ni);
      vector<string> depi = pi->dependson(plt);
      for (auto nj : depi) {
         Par* pj = plt.findp(nj);
         if (pj == nullptr)
				continue;
			if (std::find(rval.begin(), rval.end(), nj) == rval.end())
				rval.push_back(nj);
      }
   }
	T_(trc.dprint("dmatrix(",plt.desc(),") deps:",rval);)
	return rval;
}

static double
precession(pset& plt, Matrix* gyro) {
// Check the direction of precession of a gyroscopic system
//    p = y_i^t G y_r
// Returns +1 if it is forward precession, -1 if backward
	double rval{0.0};
	vector<double> y = plt.getev();
	// y is complex, split into yr, yi
	size_t n = y.size()/2;
	assert(n == gyro->rsize());
	vector<double> yr(n, 0.0);
	vector<double> yi(n, 0.0);
	blas::copy(n, &y[0], 2, &yr[0], 1);
	blas::copy(n, &y[1], 2, &yi[0], 1);
	vector<double> gy(n, 0.0);
	blas::gemv("n", n, n, 1.0, gyro->elem(), n, &yr[0], 1, 0.0, &gy[0], 1);
	blas::dot(n, &yi[0], 1, &gy[0], 1, rval);
	return rval;
}



static double
lcosvd(int nf, int nx, vector<double>& jac, int jeta, vector<double> dfdsigma) {
// compute lco stability by replacing a column (preferably vtas)
// of the Jacobian, SVD'ing it and computing dsigma/deta from the
// null vector
	T_(Trace trc(2,"lcosvd");)

	T_(trc.dprint("nx = ",nx,", nf = ",nf,", jeta = ",jeta, ", size(dfdsigma) = ",dfdsigma.size());)
	vector<double> jacobian(nf*(nx+1), 0.0);
	for (int j=0; j<nx; j++) {
		blas::copy(nf, &jac[IJ(0,j,nf)], 1, &jacobian[IJ(0,j,nf)], 1);
	}
	blas::copy (dfdsigma.size(), &dfdsigma[0], 1, &jacobian[IJ(0,nx,nf)], 1);

	SVD jf(nf,nx+1,&jacobian[0]);
	// the last 2 rows of (nx+1,nx+1) vt are the null vectors
	// project [0 ... 1 ...0] onto them
	double deta1 = jf.vt[IJ(nx-1,jeta,nx+1)];
	double deta2 = jf.vt[IJ(nx,jeta,nx+1)];
	vector<double> tan(nx+1, 0.0);
	blas::axpy(nx+1, deta1, &jf.vt[nx-1], nx+1, &tan[0], 1);
	blas::axpy(nx+1, deta2, &jf.vt[nx], nx+1, &tan[0], 1);
	
	double dsigma = tan[nx];
	double deta = tan[jeta];
	double rval{0.0};
	if (abs(deta) > numeric_limits<double>::epsilon())
		rval = dsigma/deta;
	T_(trc.dprint("returning ",dsigma,"/",deta," = ",rval);)
	return rval;
}

void
Flutcurve::
lcostability(Pac& p) {
// evaluate LCO stability by computing \frac{\partial{\sigma}}{\partial{\tau}}
//   = (\partial{sigma}/\partial{tau})/(\partial{eta}/\partial{tau})
// where tau is the arclength along a tangent to the curve
// Using an expensive technique for now: evaluate the Jacobian, add
// a column (dsigma), svd it, and project onto the nullspace
	T_(Trace trc(2,"lcostability");)

	// parameter "lcostab" must be in my params
	Par* lcopar = params.findp("lcostab");
	if (lcopar == nullptr)
		return;
	// get dsigma index into Ad::data
	int idx = Ad::find("sigma");
	if (idx < 0)
		return;
	// which column is eta?
	int jeta = index("gcnorm");
	if (jeta < 0)
		return;

	// evalutate the jacobian
	this->pac_fjac(p.x, p.f, p.jac);

	int nf = p.f.size();
	int nx = p.x.size();

	vector<double> dfdsigma(nf, 0.0);
	for (size_t i=0; i<Dy.size(); i++) {
		dfdsigma[2*i] = Ad::real(Dy[i]).data()[idx];
		dfdsigma[2*i+1] = Ad::imag(Dy[i]).data()[idx];
	}
	// compute lcostab by svd XXX this is the *expensive* way
	double lcostab = lcosvd(nf, nx, p.jac, jeta, dfdsigma);
	// update lcostab in my params
	lcopar->value(lcostab);
}
