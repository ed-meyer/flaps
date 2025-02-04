//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// start points for the continuation process are created in
// this hierarchy and returned as a vector of Flutcurve pointers
//
//                        startpts
//                           |
//                           |
//        ---------------------------------------
//        |                                     |
//    low_speed --------------            fetch_curves
//        |                  |                   |
//    eigensolution     homotopy           search_curve
//                                              |
//                                         interp_curve
//
// Curves are identified by
// 1) analysis id (aid): set with the "id" or "aid" option in flut
// 2) curve id (cid):
//        a string of digits separated by periods '.' with the first digit
//        the mode number assigned by increasing frequency in the
//        free-vibration solution in this analysis if the "source" option
//        was not included, or in the non-source analysis this analysis is
//        source'd from. The next digits are included if there are multiple-
//        fixed-valued (mfv) parameters, giving the ordinal number of these
//        multiple values.

//   the basic name for all curves is mode_n where "n" is
//   the frequency-ascending mode number. In addition there may be
//   tags for multiple-fixed-valued (mfv) parameters and start regions:
//   - srtag is a suffix (a,b,..) for the start-region number
//   - mfvtag is a suffix (1,2,..) for the combination of mfv parameters
//   so for example a curve named mode_3.4.2 in a VSO run is mode_3 with
//   the 1st mfv parameter set to its 4th value, and the 2nd mfv
//   parameter set to its 2nd value.

#include <cassert>

#include "config.h"
#include "Ad.h"
#include "concurrent.h"
#include "conv.h"
#include "Exception.h"
#include "exim.h"
#include "flutcurve.h"
#include "lapack.h"  // polyeig()
#include "pac.h"
#include "process.h"
#include "specs.h"
#include "startpts.h"
#include "trace.h"

using namespace std;

/*------------------------------------------------------------------
 * Private functions
 *------------------------------------------------------------------*/

void
fetch_curves (const string& sid,
		vector<Flutcurve*>& startpts, Region& sr);
void
search_curve (Curve& source_curve, vector<string>& messages,
		vector<Flutcurve*>& rval, Region& sr);

void
low_speed(vector<Flutcurve*>& rval, Region& sr);

void
eigensolution(vector<Flutcurve*>& rval, Region& startregion);

vector<Flutcurve*>
interp_curve(Curve* curve, string const& name, const Region& rl,
		double val, bool closest, vector<string>& messages);

static ComplexMatrix*
getComplexMatrix (const string& desc, int n);

/*------------------------------------------------------------------*/

// Public functions

vector<Region>&
startregions(string const& options) {
// add to, or get a reference to *the* list of Regions where start
// points are to be searched for
	T_(Trace trc(1,"startregions");)
	// static: not threadsafe - should be ok since this is just
	// used in parsing, before any tracking
	static vector<Region> thelist;

	// If we were called without an arg just return thelist
	if (options.empty()) {
		T_(trc.dprint("quick return: no options");)
		return thelist;
	}

	// split the options string into Tok*
	vector<Tok*> toks = flaps::lexer(options);

	// look for mode(s) = (m : s : n), convert to mode_m, mode_(m+s), ... mode_n
	// Also look for multi-valued-fixed (mvf) parameters
	vector<size_t> ld;	// leading-dimension for each mvf par
	size_t nmvf{1};		// number of mvf combinations: product of ld's
	vector<Tok*> mvftoks;
	for (auto tok : toks) {
		size_t sz = tok->svec.size();
		Par* pp = gpset::find(tok->lhs);
		if (pp != nullptr && sz > 0) {
			ld.push_back(sz);
			nmvf *= sz;
			mvftoks.push_back(tok);
		}
	}
	// for each combination of the multi-valued parameters create
	// a Region (this also works if there are no mvf toks)
	for (size_t i=0; i<nmvf; i++) {
		vector<size_t> mvfidx = vec2mdim(i, ld);
		for (size_t j=0; j<mvfidx.size(); j++) {
			mvftoks[j]->srhs = mvftoks[j]->svec[mvfidx[j]];
		}
		Region region;
		for (auto tok : toks) {
			region.push_back(Bndry(*tok));
		}
		thelist.push_back(region);
		T_(trc.dprint("got startregion<",region,">, now have ",thelist.size());)
	}

	return thelist;
}

vector<Flutcurve*>
startpts() {
// Compute all start points, return them as Flutcurve*.
// "curveid" will be something like "mode_5.aec.db" where
// "aec" means there are 3 mfv parameters and their altval
// are set to the 1st, 5th, and 3rd, respectively, and "db"
// means there are 2 fixed parameters in the startregions, and
// for this startregion they are set to the 4th and 2nd altval,
// respectively.
	T_(Trace trc(1,"startpts");)
	ostringstream os;
	vector<Flutcurve*> rval;
	vector<Region>& srlist = startregions();
	Specs& sp = flutspecs();

	if (srlist.empty()) {
		T_(trc.dprint("no start regions: adding \"anyregion\"");)
		srlist.push_back(Region());
	}

	// for each start region and each combination of mfv
	// parameters, compute start points either by a normal-mode
	// calculation or interpolating source curves
	//   srtag is a suffix (a,b,..) for the start-region number
	//   mfvtag is a suffix (1,2,..) for the combination of mfv parameters
	// Interpolated curves can have another integer for the ordinal
	// number of the target, so a curveid could be
	//   mode_15.2.a.3
	// meaning this came from a source curve "mode_15", this is the
	// second target found (either fixed parameter or point start region),
	// it is in the first start region (a), and it is the 3rd mfv
	// parameter value.
	string srtag;
	size_t nsr = srlist.size();
	T_(trc.dprint(nsr," start regions");)
	for (size_t i=0; i<nsr; i++) {
		setPointBndrys(srlist[i]);
		if (nsr > 1)
			srtag = vastr(char('a' + i));

		// do this for each mv combination - first initialize
		string mvftag = incr_mvf(0);
		do {
			vector<Flutcurve*> startpts;
			try {
				// interpolate curves from previous runs if source option
				// was included, otherwise solve the normal-modes problem
				// If a vzid was not specified take it from a source curve
				if (!sp.sid.empty()) {
					fetch_curves (sp.sid, startpts, srlist[i]);
				} else {
					low_speed (startpts, srlist[i]);
					if (startpts.empty()) {
						os.str("");
						os << "no normal modes found in specified start region";
						T_(trc.dprint("throwing exception: ",os.str());)
						throw runtime_error(os.str());
					}
				}

				// Add mvf and startregion tags to the curve id and
				// give the same curve id to the params desc member
				for (auto& curve : startpts) {
					string curveid = curve->cid();
					if (!mvftag.empty()) {
						curveid += "." + mvftag;
					}
					if (!srtag.empty()) {
						curveid += "." + srtag;
					}
					curve->cid(curveid);
					curve->params.desc(curveid);
					rval.push_back(curve);
				}
			} catch (runtime_error& s) {
				flaps::error(s.what());
			}
			if (!mvftag.empty())
				mvftag = incr_mvf();
		} while (!mvftag.empty());
	}

	T_(trc.dprint("returning ",rval.size()," start points");)
	return rval;
}

static
string
mindpress() {
// Find an indep parameter which when set to it's min or max
// results in zero (or very small) dynamic pressure. For example
// vtas, veas, vcas, alt. Note that alt is used in a constant-mach-
// variable-alt analysis and it must be set to its maximum to minimize
// dynamic pressure.
	T_(Trace trc(1,"mindpress");)
	string rval;
	vector<Par*> indeps = gpset::get().get_indep();

	// a GAF matrix is required
	Matrix* gaf = Matrix::find_desc("gaf");
	if (!gaf) {
		throw runtime_error("a GAF matrix was not specified");
	}

	// Candidate parameters at their upper bounds...
	vector<string> hi;
	hi.push_back("alt");
	hi.push_back("rf");
	// ... and candidates at their lower bounds
	// lo.push_back("cdpress");
	vector<string> lo;
	lo.push_back("dpress");
	lo.push_back("mach");
	lo.push_back("rho");
	lo.push_back("tpress");
	lo.push_back("vcas");
	lo.push_back("veas");
	lo.push_back("vtas");

	for (auto pp : indeps) {
		if (std::find(hi.begin(),hi.end(),pp->name) != hi.end()) {
			if (!pp->has_max())
				throw runtime_error(vastr("no max given for ", pp->name));
			pp->value(pp->max());
			T_(trc.dprint("set ",pp->name," to its max: ",pp->max());)
			return pp->name;
		}
		if (std::find(lo.begin(),lo.end(),pp->name) != lo.end()) {
			if (!pp->has_min())
				throw runtime_error(vastr("no min given for ", pp->name));
			pp->value(pp->min());
			T_(trc.dprint("set ",pp->name," to its min: ",pp->min());)
			return pp->name;
		}
	}

	throw runtime_error("none of the independent parameters "
			"can minimize dynamic pressure");
}

void
setPointBndrys(Region const& rl) {
// Set the values of all parameters in "rl" which are
// in point regions: regions with a parameter set to
// a value instead of an interval, e.g.
//   startregion{vtas=20, ...}
// Set the value in the global pl
// XXX should call make_equations() after setting a par to Fixed
	T_(Trace trc(1,"setPointBndrys");)

	for (size_t j=0; j<rl.size(); j++) {
		if (rl[j].isPoint()) {
			Par* pp = gpset::find(rl[j].parname);
			if (pp != nullptr) {
				T_(trc.dprint("setting ",pp->name," to ",rl[j].min);)
				pp->valuef(rl[j].min);
			} else {
				string exc = vastr("parameter \"",rl[j].parname,
						"\" has not been defined");
				throw runtime_error(exc);
			}
		}
	}
}

void
low_speed(vector<Flutcurve*>& rval, Region& startrl) {
//------------------------------------------------------------------
// Solves the flutter equations at a very low speed so that
// unsteady aerodynamics are small or non-existant.
// Requires that sigma and freq are indep parameters
// Returns: new Flutcurves are added to "rval": one for each
// freq within the start region.
//------------------------------------------------------------------
	T_(Trace trc(1,"low_speed");)
	ostringstream os;

	// find an indep parameter which minimizes dynamic
	// pressure - set it to that value
	string minpar = mindpress();
	// ... reset point regions: mindpress might have changed them
	setPointBndrys(startrl);
	// evaluate all parameters
	gpset::get().touch();
	gpset::get().eval();

	vector<Flutcurve*> modes;
	eigensolution(modes, startrl);

	// Under some circumstances the eigensolution will not be correct:
	// 1) if there are nonlinear gc and gcnorm > 0.
	// 2) if there is a "T" matrix, e.g. with controls equations
	// 3) there are terms not included in the eigensolution, e.g. gyro, vdamp
	// Correct this with a homotopy holding minpar constant
	// add these to rval
	for (auto mp : modes)
		rval.push_back(mp);

	T_(trc.dprint("returning ",modes.size()," modes, now have ",rval.size());)
}

void
getFreqLimits(Region& startregion, double& freqMin, double& freqMax) {
// See if the user has specified a freq startregion - set
// freqMin, freqMax if he has; otherwise they default to
// freq[min:max] if available, otherwise to some resonable default?
	T_(Trace trc(1,"getFreqLimits");)
	Par* freqp = gpset::get().findg("freq");
	double fconv = freqp->convFactor();

	freqMin = 0;
	freqMax = 50.0/fconv; // default limits
	if (freqp) {
		if (freqp->has_min())
			freqMin = freqp->min();
		if (freqp->has_max())
			freqMax = freqp->max();
	}
	// see if the startregion has frequency limits
	Bndry* rp = startregion.get("freq");
	if (rp != nullptr) {
		T_(trc.dprint("check startregion freq limits [",rp->min,":",rp->max,"]");)
		if (rp->min >= freqMin) {
			freqMin = rp->min;
		} else {
			flaps::info("requested min start frequency (",rp->min*fconv,
					") is less than the freq min (",freqMin,")");
		}
		if (rp->max <= freqMax) {
			freqMax = rp->max;
		} else {
			flaps::info("requested max start frequency (",rp->max*fconv,
					") is greater than the freq max (",freqMax,")");
		}
	}
	T_(trc.dprint("returning [",freqMin,":",freqMax,"]");)
}

// Eigmatrices evaluates and casts to complex all matrices
// needed for the low speed eigensolution
// An example of RAII
class Eigmatrices {
public:
	Matrix* gaf{nullptr};
	Matrix* mass{nullptr};
	Matrix* stif{nullptr};
	Matrix* vdamp{nullptr};
	Matrix* gyro{nullptr};
	Matrix* controls{nullptr};

	Eigmatrices() {
		T_(Trace trc(1,"Eigmatrices constructor");)
		vector<string> names;
		int n = get_order();
		gaf = getComplexMatrix("gaf", n);
		if (gaf != nullptr) names.push_back("gaf");
		mass = getComplexMatrix("mass", n);
		if (mass != nullptr) names.push_back("mass");
		stif = getComplexMatrix("stif", n);
		if (stif != nullptr) names.push_back("stif");
		vdamp = getComplexMatrix("vdamp", n);
		if (vdamp != nullptr) names.push_back("vdamp");
		controls = getComplexMatrix("controls", n);
		gyro = getComplexMatrix("gyro", n);
		if (gyro != nullptr)
			names.push_back("gyro");

		T_(trc.dprint("returning: ",names);)
	}
	~Eigmatrices() {
		if (gaf)
			delete gaf;
		if (mass)
			delete mass;
		if (stif)
			delete stif;
		if (vdamp)
			delete vdamp;
		if (gyro)
			delete gyro;
		if (controls)
			delete controls;
	}
}; // class Eigmatrices

void
free_vibration() {
// Compute the free-vibration frequencies and eigenvectors:
// Kx = \lambda Mx
	// Eigmatrices evaluates & casts all matrices; it's destructor
	// will delete all matrices: RAII
	Eigmatrices matrices;
	// put the stiffness (A_0) and mass (A_2) on a vector
	// for polyeig which solves
	//   (A_0 + \lambda A_1 + ... + \lambda^r A_r) x = 0
	vector<complex<double>*> ml;
	// stif: A_0
	ml.push_back(matrices.stif->celem());
	if (ml.size() == 1)  // neither gyro nor vdamp
		ml.push_back(nullptr);
	// mass: A_2
	ml.push_back(matrices.mass->celem());

	int n = matrices.mass->rsize();
	int maxeig{n};
	int nrefine{0};
	vector<complex<double> > w(n);
	vector<complex<double> > vl(n*n);
	vector<complex<double> > vr(n*n);
	vector<int> pil(n);
	vector<int> pir(n);

	// Solve the quadratic eigenvalue problem
	lapack::polyeig (n, maxeig, ml, &w[0], &vl[0], &vr[0],
			&pil[0], &pir[0], nrefine);

	// print the eigensolution (free-vibration)
	string head{" Free-vibration solution"};
	cout << stringCenter(page_width(), head) << endl;
	cout << stringCenter(page_width(), string(head.size(), '-')) << endl;
	cout << "          sigma           freq (Hz)        eigenvector\n";
	double fconv{flaps::radps2Hz};
	for (size_t i=0; i<w.size(); i++) {
		complex<double> wi(w[i].real(), fconv*w[i].imag());
		cout << std::setw(6) << i+1 << ")  "
			<< lapack::eigstring(wi,  n, &vr[IJ(0,i,n)]) << endl;
	}
	cout << endl;
}

#ifdef NEVER // not yet
static void
refine() {
		// refine the solution with nleig() using a lambda to evaluate D & derivatives
		// iff freq and sigma are independents
		Par* freqp = gpset::get().findp("freq");
		Par* sigmap = gpset::get().findp("sigma");
		vector<double> results;
		vector<double>* plot_results{nullptr};
		vector<double> stepdata;
		vector<double>* plot_stepdata{nullptr};
		vector<pair<string,string>> params;
		vector<pair<string,string>>* plot_params{nullptr};
		if (trc()) {plot_results = &results; plot_stepdata = &stepdata; plot_params = &params;}

		vector<complex<double>> x(n, 0.0);
		std::copy_n(vr.begin()+i*n, n, x.begin());
		// scale the frequency (imag part of w) for continuation
		double freqscale{100.0};
		complex<double> s = complex<double>(w[i].real(), w[i].imag()/freqscale);
		{
			nleig(s, x, [&freqscale](complex<double> s,
					const vector<complex<double>>& x,
					vector<double>& Dx,			// D(s,x)*x (2n)
					vector<double>& Dxsigma,	// d(Dx)/dsigma (2n)
					vector<double>& Dxfreq,		// d(Dx)/dfreq (2n)
					vector<double>& Dxx) {		// d(Dx)/dx (2n,2n)
					// compute d(Dx)/ds and d(Dx)/dx for nleig
						T_(Trace trc(1,"eigensolution dmat");)
						//!! Par* freqp = gpset::get().findg("freq");
						//!! Par* sigmap = gpset::get().findg("sigma");
						size_t n = x.size();
						sigmap->value(s.real());
						freqp->value(s.imag()*freqscale);
						// set the eigenvector in gpset
						gpset::get().updateEv(2*n, (double*)&x[0]);
						// eval all parameters with new freq, sigma, ev
						gpset::get().eval();
						// get the eigenvector as a complex Ad
						vector<complex<Ad>> adev;
						gpset::get().getadev(adev);
						// evaluate the dynamic matrix as a complex Ad (n,n) matrix
						vector<complex<Ad>> Dad(n*n);
						vector<complex<Ad>> work(n*n);
						Flutcurve::dmatrix(gpset::get(), Dad, work);
						// multiply Dad*adev -> work
						work.resize(n);
						blas::gemv("n",n,n,1.0, Dad.data(), n, adev.data(), 1, 0.0, work.data(), 1);
						// extract from work: Dx (double 2n-vector)
						extract(work, "", (complex<double>*)&Dx[0]);
						// ... now extract Dsigmax, Dfreqx: double 2n-vectors 
						extract(work, "sigma", (complex<double>*)Dxsigma.data());
						extract(work, "freq", (complex<double>*)Dxfreq.data());
						// scale Dxfreq
						for (auto& di : Dxfreq)
							di *= freqscale;
						// ... first part of Dxx = D + dD/dx*x: complex (n,n) matrices
						// extract into a complex (n,n) matrix...
						vector<complex<double>> cdxx(n*n, 0.0);
						extract(Dad, "", cdxx.data());
						// ... then insert the real representation in Dxx...
						blas::real_rep(n,n,cdxx.data(), n, Dxx.data(), 2*n);
						// .. and any Ad eigenvector components: these *replace* columns of Dxx
						vector<string> const& adpar = Ad::adnames();
						for (size_t j=0; j<adpar.size(); j++) {
							int column = evparse(adpar[j]);		// 0b component of real-equiv ev
							if (column >= 0)
								extract(work, adpar[j], (complex<double>*)&Dxx[column*2*n]);
						}
						return 0;
					}, plot_results, plot_params, plot_stepdata);
		}

		// plot results
		if (!results.empty()) {
			int np = params.size();
			params[np-3] = pair<string,string>("sigma","stability factor");
			params[np-2] = pair<string,string>("freq","Frequency Hz");
			// scale frequencies
			double t = freqscale*flaps::radps2Hz;
			int nc = results.size()/np;			// results is (np,nc)
			blas::scal(nc, t, &results[np-2], np);
			Apf::exporter(vastr("mode_",mode),params, results, "nleig.apf", append);
		}
		if (!stepdata.empty())
			Apf::exporter(vastr("mode_",mode), Step::parnames, stepdata, "stepdata.apf", append);
		append = true;
}
#endif // NEVER // not yet

void
eigensolution(vector<Flutcurve*>& rval, Region& startregion) {
	T_(Trace trc(1,"eigensolution");)
	double freq{0.0};
	double sdamp{0.0};
	Par* freqp = gpset::get().findg("freq");
	double freqMin = 0.0;
	double freqMax = 0.0;
	Par* sigmap = gpset::find("sigma");
	Par* growthp = gpset::find("growth");
	Par* sdampp = gpset::find("sdamp");
	Specs& sp = flutspecs();

	// do a free-vibration eigensolution if requested
	if (sp.free_vibration)
		free_vibration();

	// Get frequency limits from start regions or gpset
	getFreqLimits(startregion, freqMin, freqMax);
	
	// frequency *must* be indep
	if (!freqp->is_indep())
		throw runtime_error("frequency is not an independent parameter");

	// Get an indep "damping" parameter: either sigma,
	// growth, or sdamp must be indep.
	bool hasDamping = false;
	if (sigmap && sigmap->is_indep())
		hasDamping = true;
	if (growthp && growthp->is_indep())
		hasDamping = true;
	if (sdampp && sdampp->is_indep())
		hasDamping = true;
	if (!hasDamping)
		throw runtime_error("either sigma, growth, or sdamp must be independent");

	if (sdampp && !sdampp->is_indep()) {
		sdampp->eval(gpset::get());
		sdamp = sdampp->value();
	}

	// Before creating Eigmatrices set rf to it's max - this means
	// that the gaf matrix is an approximation and a subsequent
	// homotopy may be needed to correct this
	Par* rfp = gpset::find("rf");
	if (!rfp->is_fixed())
		rfp->value(rfp->max());
	// likewise, set freq to it's midpoint in case there is a T-matrix
	freqp->value((freqMin+freqMax)/2.0);
	T_(trc.dprint("eval T at freq = ", freqp->value());)

	// and set the eigenvector components to 1 in case there are
	// freeplay nonlinearities
	vector<double> ones(2*get_order(), 1.0);
	for (auto& i : ones)
		i = flaps::xrand();
	gpset::get().updateEv(ones);

	// evaluate the gpset before evaluating the Eigmatrices
	gpset::get().eval();

	// Eigmatrices evaluates & casts all matrices; it's destructor
	// will delete all matrices: RAII
	Eigmatrices matrices;

	// the size of these matrices may be smaller than that
	// returned by get_order() if a control-law matrix
	// was included
	assert(matrices.mass != nullptr);
	int n = matrices.mass->rsize();
	int nsq = n*n;

	// add structural damping to stiffness matrix
	if (sdamp != 0.0) {
		complex<double> ctmp(1.0, sdamp);
		blas::scal(nsq, ctmp, matrices.stif->celem(), 1);
	}

	// add T to the stiffness matrix
	if (matrices.controls != nullptr) {
		assert((int)matrices.controls->rsize() == n);
		T_(trc.dprint("adding T to stif");)
		blas::axpy (2*nsq, 1.0, matrices.controls->elem(), 1,
				matrices.stif->elem(), 1);
		T_(trc.dprint("gpset for eigensolution\n",gpset::get());)
	}

	// these are the results we need to create Flutcurves
	vector<complex<double>> w;
	vector<complex<double>> vr;

	// simpler algorithms may be available
	if (matrices.gyro != nullptr && matrices.vdamp == nullptr) {
		vector<int> pi;
		vector<double> m(nsq, 0.0);
		blas::copy(nsq, matrices.mass->elem(), 2, &m[0], 1);
		// the gyro matrix was scaled by spin in the Eigmatrices constructor
		vector<double> g(nsq, 0.0);
		blas::copy(nsq, matrices.gyro->elem(), 2, &g[0], 1);
		vector<double> k(nsq, 0.0);
		blas::copy(nsq, matrices.stif->elem(), 2, &k[0], 1);
		vector<double> iw;
		lapack::gyroeig(n, m, g, k, iw, vr, pi);		// eigensolution
		for (auto wi : iw)
			w.push_back(complex<double>(0.0,wi));
	} else {
		// put the stiffness (A_0), vdamp (A_1) and mass (A_2)
		// on a vector for polyeig which solves
		//   (A_0 + \lambda A_1 + ... + \lambda^r A_r) x = 0
		// add the gyro and vdamp matrices together, only use vdamp hereafter...
		if (matrices.gyro != nullptr) {
			assert((int)matrices.gyro->rsize() == n);
			if (matrices.vdamp != nullptr) {
				assert((int)matrices.vdamp->rsize() == n);
				// gyro was scaled by spin in Eigmatrices
				blas::axpy (2*nsq, 1.0, matrices.gyro->elem(), 1,
						matrices.vdamp->elem(), 1);
			} else {
				matrices.vdamp = matrices.gyro;
			}
			matrices.gyro = nullptr;
		}
		vector<complex<double>*> ml;
		ml.push_back(matrices.stif->celem());
		if (matrices.vdamp != nullptr)
			ml.push_back(matrices.vdamp->celem());
		else
			ml.push_back(nullptr);
		ml.push_back(matrices.mass->celem());

		// compute all 2n eigenvalue & vectors
		int maxeig{-1};
		int nrefine{0};
		int n2 = 2*n;
		w = vector<complex<double>>(n2,0.0);
		vector<complex<double> > vl(n*n2);
		vr = vector<complex<double>>(n*n2, 0.0);
		vector<int> pil(n2);
		vector<int> pir(n2);

		// Solve the quadratic eigenvalue problem
		lapack::polyeig (n, maxeig, ml, &w[0], &vl[0], &vr[0],		// eigensolution
				&pil[0], &pir[0], nrefine);

		// extract w with imag >= 0 and their right vectors
		vector<complex<double>> posw;
		vector<complex<double>> posvl;
		vector<complex<double>> posvr;
		for (int i=0; i<n2; i++) {
			if (w[i].imag() >= 0.0) {
				posw.push_back(w[i]);
				vector<complex<double>> tmp(n);
				blas::copy(n, &vr[IJ(0,i,n)], 1, tmp.data(), 1);
				blas::append(tmp, posvr);
				blas::copy(n, &vl[IJ(0,i,n)], 1, tmp.data(), 1);
				blas::append(tmp, posvl);
			}
			if ((int)posw.size() == n)
				break;
		}
		T_(trc.dprint("positive eigenvalues: ",posw);)
		T_(size_t m = posw.size();)
		T_(trc.dprintm(n,m,n,posvr.data(), "pos vectors");)
		// sort by increasing imag part
		lapack::sortEigen(n, posw.size(), posw.data(), posvl.data(), posvr.data());
		// copy these back to w, vr, vl
		w = posw;
		vr = posvr;
		vl = posvl;
		T_(trc.dprint("sorted eigenvalues: ",w);)
		T_(trc.dprintm(n,m,n,vr.data(), "sorted vectors");)
	}

	// create a new Flutcurve for each frequency and curveid requested
	double eps{std::numeric_limits<double>::epsilon()};
	int mode{0};  // 1b ordinal mode no. in freq range
	vector<int> outofrange;
	bool append{false};
	for (size_t i=0; i<w.size(); i++) {
		// skip negative frequencies - don't count them in ordinal mode #s
		// check if this curveid was requested
		mode++;
#ifdef NEVER // new cid format
		string curveid{vastr("mode_", mode)};
#else // NEVER // new cid format
		//!! string curveid{vastr(sp.aid, ".", mode)};
		string curveid{std::to_string(mode)};
#endif // NEVER // new cid format
		if (!startregion.contains_mode(mode) && !startregion.curveidOk(curveid)) {
			T_(trc.dprint("reject mode ",mode,": not requested");)
			continue;
		}
		// freq in range?
		freq = w[i].imag();
		bool inrange = false;
		if (freqMin <= freq && freq <= freqMax) {
			inrange = true;
		}
		if (!inrange) {
			outofrange.push_back(mode);
			T_(trc.dprint("ignoring mode ",mode," freq (",freq, ") out of range [", freqMin,':',freqMax,']');)
			continue;
		}
		// watch for zero vector XXX why does polyeig return zero eigenvector for i>n?
		double vnorm = blas::snrm2(2*n, (double*)&vr[IJ(0,i,n)], 1);
		if (vnorm < eps)
			continue;

		// update sigma, freq, and eigenvector in gpset
		vector<complex<double>> x(n, 0.0);
		std::copy_n(vr.begin()+i*n, n, x.begin());
		complex<double> s(w[i].real(), w[i].imag());
		sigmap->value(s.real());
		freqp->value(s.imag());
		gpset::get().updateEv(2*n, (double*)&x[0]);
		// eval all parameters with new freq, sigma, ev
		gpset::get().eval();

		if (sigmap->is_indep()) {
			// refine the solution with nleig() using a lambda to evaluate D & derivatives
			// iff freq and sigma are independents
			vector<double> results;
			vector<double>* plot_results{nullptr};
			vector<double> stepdata;
			vector<double>* plot_stepdata{nullptr};
			vector<pair<string,string>> params;
			vector<pair<string,string>>* plot_params{nullptr};
			T_(if (trc()) {plot_results = &results; plot_stepdata = &stepdata; plot_params = &params;})

			// scale the frequency (imag part of w) for continuation
			double freqscale{100.0};
			s = complex<double>(s.real(), s.imag()/freqscale);
			// nonlinear eigenvalue function, part of pac.cpp/h
			nleig(s, x, [&freqscale,&freqp,&sigmap](complex<double> s,
					const vector<complex<double>>& x,
					vector<double>& Dx,			// D(s,x)*x (2n)
					vector<double>& Dxsigma,	// d(Dx)/dsigma (2n)
					vector<double>& Dxfreq,		// d(Dx)/dfreq (2n)
					vector<double>& Dxx) {		// d(Dx)/dx (2n,2n)
					// compute d(Dx)/ds and d(Dx)/dx for nleig
						T_(Trace trc(1,"eigensolution dmat");)
						size_t n = x.size();
						sigmap->value(s.real());
						freqp->value(s.imag()*freqscale);
						// set the eigenvector in gpset
						gpset::get().updateEv(2*n, (double*)&x[0]);
						// eval all parameters with new freq, sigma, ev
						gpset::get().eval();
						// get the eigenvector as a complex Ad
						vector<complex<Ad>> adev;
						gpset::get().getadev(adev);
						// evaluate the dynamic matrix as a complex Ad (n,n) matrix
						vector<complex<Ad>> Dad(n*n);
						vector<complex<Ad>> work(n*n);
						Flutcurve::dmatrix(gpset::get(), Dad, work);
						// multiply Dad*adev -> work
						work.resize(n);
						blas::gemv("n",n,n,1.0, Dad.data(), n, adev.data(), 1, 0.0, work.data(), 1);
						// extract from work: Dx (double 2n-vector)
						extract(work, "", (complex<double>*)&Dx[0]);
						// ... now extract Dsigmax, Dfreqx: double 2n-vectors 
						extract(work, "sigma", (complex<double>*)Dxsigma.data());
						extract(work, "freq", (complex<double>*)Dxfreq.data());
						// scale Dxfreq
						for (auto& di : Dxfreq)
							di *= freqscale;
						// ... first part of Dxx = D + dD/dx*x: complex (n,n) matrices
						// extract into a complex (n,n) matrix...
						vector<complex<double>> cdxx(n*n, 0.0);
						extract(Dad, "", cdxx.data());
						// ... then insert the real representation in Dxx...
						blas::real_rep(n,n,cdxx.data(), n, Dxx.data(), 2*n);
						// .. and any Ad eigenvector components: these *replace* columns of Dxx
						vector<string> const& adpar = Ad::adnames();
						for (size_t j=0; j<adpar.size(); j++) {
							int column = evparse(adpar[j]);		// 0b component of real-equiv ev
							if (column >= 0)
								extract(work, adpar[j], (complex<double>*)&Dxx[column*2*n]);
						}
						return 0;
					}, plot_results, plot_params, plot_stepdata);

			// plot results
			if (!results.empty()) {
				int np = params.size();
				params[np-3] = pair<string,string>("sigma","stability factor");
				params[np-2] = pair<string,string>("freq","Frequency Hz");
				// scale frequencies
				double t = freqscale*flaps::radps2Hz;
				int nc = results.size()/np;			// results is (np,nc)
				blas::scal(nc, t, &results[np-2], np);
				Apf::exporter(vastr("mode_",mode),params, results, "nleig.apf", append);
			}
#ifdef NEVER // new Fv
			if (!stepdata.empty())
				Apf::exporter(vastr("mode_",mode), Step::parnames,
						stepdata, "stepdata.apf", append);
#endif // NEVER // new Fv
			append = true;
		}

		// Create a new Flutcurve based on gpset
		Flutcurve* newcurve = new Flutcurve(sp.aid, curveid, sp.vzid);
		// add to return list
		rval.push_back(newcurve);
		T_(trc.dprint("new start guess ",curveid, ": ",newcurve->current_values());)
	}

	if (!outofrange.empty())
		flaps::warning(outofrange.size()," modes had frequencies out of range");

	T_(trc.dprint("returning ",rval.size()," roots");)
	return;
}  // eigensolution


class IgnoreMsg {
	// An IgnoreMsg contains the message and a description of
	// the aeroelastic mode that is being ignored (e.g. "mode 3")
public:
	string msg;
	vector<string> modes;

	IgnoreMsg(const string& m, const string& cid) : msg(m) { modes.push_back(cid); }
};

int
consolidateIgnoreMsg(vector<IgnoreMsg>& ignoreMsg) {
// consolidate multiple IgnoreMsg's into one per unique msg,
// return the number of modes in the collection
	vector<IgnoreMsg> con;
	size_t i, j, k;
	int rval = 0;

	if (ignoreMsg.size() < 2) {
		for (i=0; i<ignoreMsg.size(); i++)
			rval += ignoreMsg[i].modes.size();
		return rval;
	}

	for (i=0; i<ignoreMsg.size(); i++) {
		// see if this IgnoreMsg::msg has been encountered; if
		// so consolidate it's modes
		for (j=0; j<con.size(); j++) {
			if (ignoreMsg[i].msg == con[j].msg) {
				for (k=0; k<ignoreMsg[i].modes.size(); k++) {
					con[j].modes.push_back(ignoreMsg[i].modes[k]);
				}
				break;
			}
		}
		// if not encountered add it to the return set
		if (j == con.size()) {
			con.push_back(ignoreMsg[i]);
		}
	}
	ignoreMsg = con;
	for (i=0; i<ignoreMsg.size(); i++)
		rval += ignoreMsg[i].modes.size();
	return rval;
}


void
fetch_curves (const string& sid, vector<Flutcurve*>& rval,
		Region& startregion) {
	T_(Trace trc(1,"fetch_curves ",sid);)
	size_t i, j;
	ostringstream os;
	vector<IgnoreMsg> ignoreMsg;
	Specs& sp = flutspecs();

	// Eigmatrices evaluates & casts all matrices; it's destructor
	// will delete all matrices: RAII
	// XXX Do this here so that if the user does static Interps
	// they will be created here instead of in a thread
	Eigmatrices matrices;

	// get a vector of all the source curve names
	string mid = Curve::midrx(sid, ".*");
	vector<string> mids = fio::catalog(mid);

	// Read each curve for each requested mode. Quit searching
	// when a start is found (only track each mode once) XXX ordinal deprecated?
	// OR when we exhaust all avail curves for the mode (fetch returns nullptr)
	for (auto& mi : mids) {
		// before fetching this curve, check that it's cid matches
		// any in "startregion"
		Curve* ci{nullptr};
		try {
			string aid, cid;
			Curve::mid2comp(mi, aid, cid);
			if (!startregion.curveidOk(cid)) {
				T_(trc.dprint("rejecting this region: curveid does not match");)
				ignoreMsg.push_back(IgnoreMsg("does not match requested id", cid));
				continue;
			}
			ci = Curve::fetch(mi);
			if (ci->nsolns() == 0) {
				T_(trc.dprint(ci->cid()," is an empty curve");)
				continue;
			}
			// if a vzid was not specified take it from a source curve
			if (sp.vzid.empty())
				sp.vzid = ci->vzid();
			// XXX need RAII to delete ci
		} catch (runtime_error const& s) {
			T_(trc.dprint("no more results available for analysis id \"", sid,"\": ",s.what());)
			break;
		}
		T_(trc.dprint("working on curve \"",ci->cid(),"\"");)

		try {
			vector<string> messages;
			// get all points matching startregion, interpolating if necessary
			search_curve (*ci, messages, rval, startregion);

			for (auto& msg : messages)
				ignoreMsg.push_back(IgnoreMsg(msg, ci->cid()));

			delete ci;
		} catch (StartSearch& s) {
			ignoreMsg.push_back(IgnoreMsg(s.what(), ci->cid()));
			T_(trc.dprint("caught StartSearch exc: ",s," now have ",ignoreMsg.size()," messages");)
		}
	}

	// consolidate messages
	int nignored = consolidateIgnoreMsg(ignoreMsg);

	// Print messages regarding ignored modes
	size_t width = page_width() - 10;
	if (!ignoreMsg.empty()) {
		cout << endl;
		string blanks(20,' ');
		os.str("");
		os << "Ignoring " << nignored << " start point";
		if (nignored > 1)
			os << "s";
		os << ':';
		cout << blanks << os.str() << endl;
		cout << blanks << string(os.str().size(), '-') << endl;
		for (j=0; j<ignoreMsg.size(); j++) {
			if (ignoreMsg[j].modes.size() > 1) {
				string::size_type idx = ignoreMsg[j].msg.find("does", 0);
				if (idx != string::npos)
					ignoreMsg[j].msg.replace(idx, 4, "do");
			}
			// message: modes line...
			os.str("");
			os << "   " << ignoreMsg[j].modes.size() << " "
				<< ignoreMsg[j].msg << ':';
			i = 0;
			while(os.str().size() < width && i < ignoreMsg[j].modes.size()) {
				os << "  " << ignoreMsg[j].modes[i++];
			}
			cout << os.str() << endl;
			// ...remaining modes on new lines
			while(i < ignoreMsg[j].modes.size()) {
				os.str("");
				os << "       ";
				while(os.str().size() < width && i < ignoreMsg[j].modes.size()) {
					os << "  " << ignoreMsg[j].modes[i++];
				}
				cout << os.str() << endl;
			}
		}
		cout << separator() << endl;
	}

	T_(trc.dprint("returning ",rval.size()," start points");)
}

void
search_curve (Curve& source_curve, vector<string>& messages,
		vector<Flutcurve*>& rval, Region& rl) {
/*------------------------------------------------------------------
 * Given a curve ("solns" array of "source_curve" parameters),
 * interpolate to find parameter values which match the values of:
 * 1) point Bndrys given in the Region (rl), e.g. sigma=0, or
 * 2) the Fixed parameters that were given as options.
 * Returns:
 *   a new parameter list, or throws exceptions:
 *     OutOfRange    the entire curve has parameters out of range
 *     DiffValue     the value of a Fixed parameter was different
 *                   in the source analysis
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"search_curve ",source_curve.cid());)
	ostringstream os;
	string interpPar;
	string curveid = source_curve.cid();
	int nnew{0};

	T_(trc.dprint(" start region: ",rl,", already have ",rval.size()," start pts");)

	// 1) does this Region (rl) have a point Bndry (e.g. vtas=10)?
	//    if so interpolate to that point
	Bndry* pointBndry = rl.pointBndry();
	if (pointBndry != nullptr) {
		// interpolate to pointBndry->parname=min
		vector<Flutcurve*> newcurves = interp_curve(&source_curve,
				pointBndry->parname, rl, pointBndry->min,
				pointBndry->closest, messages);

		T_(trc.dprint("got ",newcurves.size(),", ",pointBndry->parname, " = ",pointBndry->min," Solutions, ", messages.size()," messages");)

		// Create a new Flutcurve for each interpolated pset
		// give each interpolated value a description which is
		// the existing description appended with the ordinal
		// number of this interpolation
		for (size_t k=0; k<newcurves.size(); k++) {
			int ordinal = k+1;
			if (!rl.contains_ordinal(ordinal)) {
				delete newcurves[k];
				continue;
			}
			// check startregion limits
			string oormsg;
			const Par* oorpar{nullptr};
			oorpar = rl.outOfRange(newcurves[k]->params, oormsg);
			if (oorpar == nullptr) {
				T_(trc.dprint("ok, adding to rval");)
				rval.push_back(newcurves[k]);
				nnew++;
				T_(trc.dprint("added ",newcurves[k]->cid(),", now have ",rval.size());)
			}
		}
	}

	// if we added start points with pointregrions and
	// curveid's we are done
	if (nnew > 0) {
		T_(trc.dprint("returning ",nnew," new start points");)
		return;
	}

	// Quick return if only one point
	if (source_curve.nsolns() < 2) {
		T_(trc.dprint("returning empty pset: only ", source_curve.nsolns()," points in source");)
	}

	// 2) interpolate on a Fixed or multiple-fixed-value (mvf) parameter
	//    For each parameter that was declared Fixed as an option
	//    (so Par::pref=1) in the current analysis:
	//       a) if it was constant in the source analysis and the values
	//          in the two analyses are different, reject this curve.
	//       b) if it was varied in the source analysis interpolate
	//          to get the source solution at the desired Fixed value.
	// Check that values of all the parameters are in range
	// and check that the values are in the startregions.
	for (auto& par : gpset::get().pmap()) {
		Par* si = par.second;
		Par* fromi = source_curve.params.findp(si->name);
		if (fromi == nullptr)  // not in source
			continue;
		// pmin, pmax: limits of the source curve
		double pmin = fromi->min_solns();
		double pmax = fromi->max_solns();

		// is the entire curve out of range?
		// Ignore limits on rf
		if (si->name != "rf") {
			if (si->has_min() && is_lessthan(pmax, si->min(), 5)) {
				os.str("");
				os << si->name << " is below its minimum (" << si->min() << ")"
					<< " throughout the curve (max is " << pmax << ')';
				T_(trc.dprint("throwing exception: ",os.str());)
				throw OutOfRange(os.str());
			}
			if (si->has_max() && is_greaterthan(pmin, si->max(), 5)) {
				os.str("");
				os << si->name << " is above its maximum (" << si->max() << ")"
					<< " throughout the curve (min is " << pmin << ')';
				T_(trc.dprint("throwing exception: ",os.str());)
				throw OutOfRange(os.str());
			}
		}
		if (si->pref && si->is_fixed()) {
			T_(trc.dprint("trying ",si->name," as interpolation parameter");)
			// was it constant? if so check to see if it matches input value...
			if (is_equal(pmin, pmax, 5)) {
				if (!is_equal(si->value(), pmin, 5)) {
					os.str("");
					os << fromi->name << " was constant in the source anaysis ("
						<< pmin << ") but not " << si->prevalue();
					flaps::warning(os.str());
				} else {
					T_(trc.dprint("value (",si->prevalue(),") matches source");)
				}
			} else {
				// ... otherwise interpolate to get a param list at the
				// desired value(s) and return these
				// It is an error if more than one Fixed parameter
				// requires interpolation
				// interp_curve() will throw an exception if the curve
				// does not contain the value - don't bother catching it
				if (!interpPar.empty())
					throw runtime_error(vastr("curves require interpolation in both ",
						interpPar, " and ", si->name));

				// interpolate to si->value()
				vector<Flutcurve*> newcurves =
					interp_curve(&source_curve, si->name,
							rl, si->value(), false, messages);

				// check each returned Flutcurve against start region
				// ordinals and limits
				for (size_t j=0; j<newcurves.size(); j++) {
					int ordinal = j+1;
					if (!rl.contains_ordinal(ordinal)) {
						delete newcurves[j];
						continue;
					}
					// check start region limits
					string oormsg;
					const Par* oorpar{nullptr};
					oorpar = rl.outOfRange(newcurves[j]->params, oormsg);
					if (oorpar == nullptr) {
						rval.push_back(newcurves[j]);
						T_(trc.dprint("added ",newcurves[j]->cid(), ", now have ",rval.size());)
					} else {
						messages.push_back(oormsg);
						delete newcurves[j];
					}
				}
				// save the name of the interpolation parameter to check
				// for multiple interpolations required
				interpPar = si->name;
			}
		}
	}
	T_(trc.dprint("returning ",rval.size()," solutions");)
}

vector<Flutcurve*>
interp_curve(Curve* curve, string const& name, const Region& rl,
		double val, bool closest, vector<string>& messages) {
// find all points in "curve" where parameter "name" equals "val".
// Returns a vector of new Flutcurves at the interpolated values
	T_(Trace trc(1,"interp_curve ",curve->cid());)
	ostringstream os;
	Specs& sp = flutspecs();

	T_(trc.dprint("searching for ",name,closest?" closest to ":" = ",val);)
	T_(trc.dprint("curve parameters:",curve->params);)
	T_(trc.dprint("already have ",messages.size()," ignore-messages");)

	Par* interp_par = curve->params.findp(name);
	if (interp_par == nullptr) {
		string exc = vastr("\"",name,"\" is not a valid parameter");
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}
	
	// find all places where parameter "name" equals val
	vector<pair<size_t,double> > loc;
	if (closest)
		loc = interp_par->find_closest(val);
	else
		loc = interp_par->find_solns(val);

	if (loc.empty()) {
		double convfactor = interp_par->convFactor();
		string exc = vastr("does not have ", interp_par->name, " = ", val*convfactor);
		T_(trc.dprint("throwing exception: ",exc);)
		throw OutOfRange(exc);
	}

	T_(trc.dprint("startregions: ",rl);)
	T_(trc.dprint("found ",loc.size()," values where ",interp_par->name," = ",val);)

	// for each loc create a Flutcurve with this value
	vector<Flutcurve*> curves;
	for (size_t j=0; j<loc.size(); j++) {
		// interpolate each parameter in "curve": sets the value
		T_(trc.dprint("interpolating parameters, start ",loc[j].first,", t=",loc[j].second);)
		curve->params.interp(loc[j]);
		// re-set the interp_par to "exactly" val - may be important for
		// downstream flut processes that use this as source
		interp_par->valuef(val);

		vector<string> names = curve->inrange();
		if (!names.empty()) {
			string msg;
			for (auto& pi : names) {
				if (pi == "rf" && !sp.rflimits)
					continue;
				add2msg(msg, vastr(pi," is out of range"));
			}
			if (!msg.empty()) {
				T_(trc.dprint("ignoring crossing: ",msg);)
				continue;
			}
		}

		// ... and in the list of Regions
		T_(trc.dprint("checking target limits");)
		const Par* oorpar{nullptr};
		string oormsg;
		oorpar = rl.outOfRange(curve->params, oormsg);
		if (oorpar != nullptr) {
			messages.push_back(oormsg);
			continue;
		}

		// all within limits - create new curve id
		string curveid{curve->cid()};

		// update gpset with the interpolated curve values,
		// except for fixed parameters
		gpset::get().update(curve->params);

		// create a Flutcurve from gpset, and put it
		// on the return vector
		Flutcurve* newcurve = new Flutcurve(sp.aid, curveid, sp.vzid);
		// add to return list
		curves.push_back(newcurve);
		T_(trc.dprint("new interpolated start guess ",curveid, ": ",newcurve->current_values());)
	}

	// check for ordinals
	vector<Flutcurve*> rval;
	for (size_t j=0; j<curves.size(); j++) {
		Flutcurve* curve = curves[j];
		int ordinal = j+1;
		if (rl.contains_ordinal(ordinal)) {
			// multiple locations? add .int to curveid
			//!! if (loc.size() > 1)
			if (curves.size() > 1) {
				curve->cid(vastr(curve->cid(),".", ordinal));
				T_(trc.dprint("adding ordinal to cid: ",curve->cid(),", ",loc.size()," locations");)
			}
			rval.push_back(curve);
		} else {
			delete curve;
		}
	}

	T_(trc.dprint("returning ",rval.size()," start points");)
	return rval;
}  // interp_curve

static
Matrix*
getComplexMatrix (const string& desc, int n) {
// Returns a *copy* of the "desc" matrix, evaluated using the current
// values of parameters in the gpset, and inserted into an (n,n)
// complex Matrix. It is the caller's responsibility to delete the
// returned matrix. Returns nullptr if the matrix has not been declared.
	T_(Trace trc(1,"getComplexMatrix ", desc);)
	ostringstream os;

	Matrix* mp = Matrix::find_desc(desc);

	if (mp == nullptr) {
		T_(trc.dprint("returning nullptr: ",desc," is missing");)
		return nullptr;
	}

	// evaluate the matrix at current values of gpset, insert
	// into an (n,n) complex AD array
	vector<complex<Ad>> work(n*n);
	mp->eval(gpset::get(), work);

	// create a new (n,n) complex matrix to return
	Matrix* rval = new Matrix(desc, desc, n, n, true);

	// extract the values from the AD array, put them into rval
	extract (work, 0, rval->celem());

	T_(trc.dprintm(n,n,n,rval->celem(), rval->mid());)

	return rval;
}
