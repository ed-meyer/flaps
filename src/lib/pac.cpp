//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "exim.h"
#include "lapack.h"
#include "matrix.h" // flaps::summarize()
#include "pac.h"
#include "qr.h"
#include "trace.h"

using namespace std;

static double
vector_distance (int n, const double* a, const double* b);

#ifdef NEVER // no Fv
ostream&
operator<<(ostream& s, const Fv& t) {
	s << "  " << t.name << " (" << t.desc << ") = " << t.value << endl;
	return s;
}

ostream&
operator<<(ostream& s, const Iv& t) {
	s << "  " << t.name << " (" << t.desc << ") = " << t.value << endl;
	return s;
}
#endif // NEVER // no Fv

std::ostream&
operator<<(std::ostream& s, const Issue* t) {
	if (t != nullptr) {
		s << t->msg;
	}
	return s;
}

void
Det::
equalize(Det const& p) {
// adjust my det and exp so that exp == p.exp
	Trace trc(1,"Det::equalize");
	trc.dprint("this = ",*this,", p = ",p);
	if (this->exp != p.exp) {
		this->coef *= pow(10.0, this->exp-p.exp);
		this->exp = p.exp;
	}
	trc.dprint("equalized this = ",*this,", p = ",p);
}

std::ostream&
operator<<(std::ostream& s, const Det& t) {
	s << t.coef << "*10^" << t.exp;
	return s;
}


Pac::
Pac(int nf, int nx, Fjacfcn fj) : fjac(fj) {
// Pac constructor: allocate x, f, jac and tan
	Trace trc(1,"Pac constructor");
	f = vect(nf, 0.0);
	x = vect(nx, 0.0);
	jac = vect(nf*nx, 0.0);
	tan = vect(nx, 0.0);
}

Pac::
Pac(const vect& start, const vect& t, int nf, Fjacfcn fj) :
	x(start), tan(t), fjac(fj) {
// Pac constructor: allocate x, f, jac and tan and set starting x and tan
	Trace trc(1,"Pac x constructor");
	f = vect(nf, 0.0);
	int nx = x.size();
	jac = vect(nf*nx, 0.0);
	trc.dprint("nx = ",x.size(),", nf = ",nf);
}

std::ostream&
operator<<(std::ostream& s, const Pac& t) {
	s << "coord " << t.coord << " x = " << t.x;
	return s;
}

Issue*
Pac::
sigue (Pac& to, vector<double>* projectee) {
// Extend this Pac to "to" by continuation using the predicted stepsize (hk)
	Trace trc(1,"sigue from coord ",coord," to ",to.coord,", stepsize ",step.hk);
	//  int nf = f.size();
	int nx = x.size();
	Issue* rval{nullptr};

	// give "to" either my tan or a projectee to use in corrector to compute
	// new tangent by projecting onto the Jacobian nullspace
	if (projectee != nullptr)
		to.tan = *projectee;
	else
		to.tan = tan;
	// initialize to.stepsize with my predicted stepsize -
	// it may be reduced below
	to.stepsize = step.hk;

	trc.dprint("from tan: ",tan);

	// Loop until we get an acceptable solution...
	while(true) {
		// predictor: to = this + to.stepsize*tan
		to.x = x;
		blas_axpy(nx, to.stepsize, tan.data(), 1, to.x.data(), 1);
		// ...corrector
		if (confcn == nullptr)
			rval = to.corrector();
		else
			rval = to.corrector(confcn);

		// converged? get a stepsize (to.hk) for the next step
		if (rval == nullptr || rval->converged()) {
			Issue* iss = to.rheinboldt(*this);
			// did rheinboldt return an issue (other than large angle)?
			if (iss == nullptr)
				break;
			if (rval == nullptr)
				rval = iss;
			else
				add2msg(rval->msg, iss->msg);
			if (iss->is_angle())
				break;
		}

		// non-converged Issue: re-try with a smaller stepsize,
		// unconstrained, unless the stepsize is too small
		to.stepsize /= specs.reductionfac*specs.reductionfac;
		//!! to.step.red_angle += 1.0;
		if (rval->is_angle())
			to.step.red_angle += 1.0;
		else
			to.step.red_conv += 1.0;
		confcn = nullptr;
		//!! if (stepsize < specs.minstepsize*8.0)
		if (to.stepsize < specs.minstepsize*pow(specs.reductionfac,4.0))
			return rval;
		int red = step.red_angle + step.red_conv + step.red_dsc;
		trc.dprint("reduction ",red,", retry with stepsize ",to.stepsize,", issue: ",rval);
	}

	return rval;
} // sigue

Det
Pac::
determinant() {
// compute the determinant of the Jacobian augmented by
// the tangent in the last row
	Trace trc(2,"Pac::determinant");
	Det rval;
	int nx = x.size();
	int nf = f.size();
	// only for nx == nf+1
	if (nx != nf+1) {
		trc.dprint("quick return: ",nf,", ",nx);
		return rval;
	}
	fjac (x, f, jac);
	vector<double> A(nx*nx, 0.0);
	for (int j=0; j<nx; j++)
		blas_copy(nf, &jac[j*nf], 1, &A[j*nx], 1);
	blas_copy(nx, tan.data(), 1, &A[nf], nx);
	vector<int> ipiv(nx, 0);
	lapack::dgetrf(nx,nx,A.data(), nx, ipiv.data());
	double d{1.0};
	double exp{0.0};
	// keep 0 <= abs(d) < 10
	for (int j=0; j<nx; j++) {
		d *= A[j*(1+nx)];
		if (d == 0.0) break;
		if (ipiv[j] != j+1) d = -d;
		while (abs(d) < 1.0) {
			d *= 10.0;
			exp -= 1.0;
		}
		while(abs(d) > 10.0) {
			d /= 10.0;
			exp += 1.0;
		}
	}
	rval = Det(d,exp);
	trc.dprint("returning ",rval);
	return rval;
}

Issue*
Pac::
corrector() {
//------------------------------------------------------------------
// Correct a predicted solution to
//    y = f(x) = 0
// using Gauss-Newton iteration.
// Returns the convergence criteria (int):
//  Negative values are errors and this->error gives the reason:
//  -3    some parameter is out of range that prevents further processing
//  -2    converging too slowly - more than maxiter iterations
//  -1    residual is not decreasing
//  Positive values indicate success:
//   1    strong acceptance:
//   2    weak acceptance: small residual or small iterate
//   3    weak acceptance: two consecutive small residuals
//   4    weak acceptance: 2 consecutive small iterates
//------------------------------------------------------------------
	Trace trc(2,"Pac::corrector");
	Issue* rval{nullptr};
	int nx = x.size();
	int nf = f.size();
	double meps = std::numeric_limits<double>::epsilon();
	// convergence tolerance for norm(f):
	QR* jf{nullptr};

	trc.dprint ("nx ", nx, ", nf ", nf,", x: ", x);

	// Workspace
	vect h(nx, 0.0);		// Newton correction
	vect hold(nx, 0.0);	// previous correction
	vect x0(x);				// initial x: copy constructor

	// Newton-iteration loop - evalutate jacobian each time through.
	double prev_hnorm{0.0};
	double prev_fxnorm{0.0};

	for (step.niter=0; step.niter < specs.maxiter; step.niter++) {
		// evaluate the function & jacobian
		fjac (x, f, jac);

		// compute abs convergence tolerance by estimating the smallest
		// possible residual norm
		if (specs.epsabs == 0.0)
			specs.epsabs = ecc();

		// compute function norm, save previous
		if (step.niter > 0)
			prev_fxnorm = step.fxnorm;
		step.fxnorm = blas_snrm2 (nf, &f[0], 1);
		step.xnorm = blas_snrm2 (nx, &x[0], 1);
		trc.dprint("iter ",step.niter,", residual ",step.fxnorm, ", norm(x) ",step.xnorm,", f: ",flaps::summarize(f,80));

		// check for convergence XXX move to after x+h
		if (specs.minf.empty())
			rval = converged(prev_fxnorm, prev_hnorm);
		else
			rval = convergev(prev_fxnorm, prev_hnorm);

		// sanity check - return Issue (rval) must be set if conv_crit < 0
		if (step.conv_crit < 0 && rval == nullptr)
			throw runtime_error(vastr("conv crit ",step.conv_crit," with no Issue"));

		// return if converged returned an Issue
		if (rval != nullptr) {
			trc.dprint("returning: ",rval->msg);
			return rval;
		}
		// if it is converged break out of the loop only if we have
		// been through at least once so that jf != nullptr
		if (step.conv_crit > 0 && jf != nullptr)
			break;

		// Factor the jacobian: may throw exceptions?
		if (jf != nullptr) delete jf;
		jf = new QR(nf, nx, &jac[0]);

		step.rcond = jf->rcond;
		step.svd = (jf->svd == nullptr ? 0 : 1);

		// solve for the Newton correction (h)
		jf->solve(f, h);

		// norm of the correction (h)
		if (step.niter > 0)
			prev_hnorm = step.hnorm;
		step.hnorm = blas_snrm2 (nx, &h[0], 1);
		trc.correction(nx,&this->x[0],&h[0],"correction ",step.niter);

		// subtract h from x...
		blas_axpy(nx, -1.0, &h[0], 1, &x[0], 1);

		// check for convergence again in case jf was nullptr
		if (step.conv_crit > 0)
			break;
	}  // end of Newton iteration loop

	// Compute last and total corrector distances
	if (step.niter > 0) {
		step.deltak = vector_distance(nx, &x[0], &x0[0]);
	} else {
		step.hnorm = 0.0;
		step.deltak = 0.0;
	}

	// compute the tangent by projecting the input tan onto the
	// Jacobian nullspace; watch out for small projection
	double scale{0.0};
	//!! double tol{0.01};
	double tol = sqrt(meps);
	tan = jf->nullProj(tan, scale);
	// normalize
	projnorm = blas_snrm2(nx, &tan[0], 1);
	if (projnorm > meps)
		blas_scal(nx, 1.0/projnorm, &tan[0], 1);
	if (projnorm < tol)
		rval = new Issue("small tangent encountered", Kind::zerotan);
	trc.dprint("tangent: ",tan);

	// delete the factorization
	if (jf != nullptr)
		delete jf;

	trc.dprint("returning convergence crit: ", step.conv_crit, ", rcond ",step.rcond);
	return rval;
}  // corrector

Issue*
Pac::
corrector(Constraintfcn constraint) {
//------------------------------------------------------------------
// Correct a predicted solution to
//    y = f(x) = 0
// using Gauss-Newton iteration with a constraint.
// Corrections to x will be orthogonal to the constraint vector.
// This function simply calls the unconstrained corrector() with
// a special Fjacfcn function which adds the constraint.
//------------------------------------------------------------------
	Trace trc(2,"Pac::corrector constrained");
	size_t nx = x.size();
	size_t nf = f.size();
	int nf_c = nf + 1;
	Issue* rval{nullptr};

	// lambda to call the input Fjacfcn and add the constraint
	Fjacfcn fjac_c = [&](const vect& x_c, vect& f_c, vect& jac_c) {
		Trace trc(1,"fjac_c");
		this->fjac(x_c, this->f, this->jac);
		std::copy_n(this->f.begin(), nf, f_c.begin());
		for (size_t j=0; j<nx; j++)
			std::copy_n(this->jac.begin()+j*nf, nf, jac_c.begin()+j*nf_c);
		// add the constraint to the last row of f and the Jacobian
		vect c(nx, 0.0);		// constraint vector
		f_c[nf] = constraint(x, c);	// constraint returns p-\hat{p}
		blas_copy(nx, &c[0], 1, &jac_c[nf], nf_c);
		return 0;
	};
		
	// create new "to" and "from" Pacs with the constraint
	Pac to_c(this->x, this->tan, nf_c, fjac_c);
	to_c.coord = this->coord;
	// also give to_c this->step in case it has info from
	// previous failures - copy it back below
	to_c.step = this->step;

	// correct it using the "unconstrained" corrector
	rval = to_c.corrector();
	// copy info from to_c to "this"
	this->x = to_c.x;
	this->step = to_c.step;

	if (rval != nullptr)
		return rval;

	// the tangent in to_c is not correct - compute it using the
	// last this->jac XXX why is it not correct?
	this->fjac(this->x, this->f, this->jac);
	QR jf(nf,nx,this->jac.data());
	double scale{0.0};
	this->tan = jf.nullProj(this->tan, scale);
	// normalize
	projnorm = blas_snrm2(nx, this->tan.data(), 1);
	blas_scal(nx, 1.0/projnorm, this->tan.data(), 1);
	trc.dprint("tangent: ",this->tan);

	return rval;
}  // corrector constrained


Issue*
Pac::
converged(double prev_fxnorm, double prev_hnorm) {
// Check for convergence or non-convergence of this point, and set
// step.conv_crit to one of the following:
//    0  neither convergence or non-convergence
//   -1 residual is not decreasing and small hnorm
//   -2  no convergence after "maxiter" iterations
//    1  strong acceptance
//    2  weak acceptance: small h
//    3  weak acceptance: two consecutive small residuals
//    4  weak acceptance: 2 consecutive small f and h
// Returns a message if step.conv_crit != 0 telling why it
// did (>0) or did not (<0) converge
// References
//   Rheinboldt, Werner C., "Numerical Analysis of Parameterized
//   Nonlinear Equations", John Wiley & Sons, New York, 1986
//   ISBN 0-471-88814-1 (book)
	Trace trc(1,"converged");
	Issue* rval{nullptr};
	double meps = std::numeric_limits<double>::epsilon();
	double relerr = sqrt(meps);
	double weakEps = 8.0*meps;
	double fmp{2.0};
	// relative error measure, p. 122 ref
	double tlstep = specs.epsabs + relerr*step.xnorm;
	Form sci1;
	sci1.scientific().precision(1);

	trc.dprint("hnorm ",step.hnorm,"(prev ",prev_hnorm, ") fxnorm ",step.fxnorm, "(prev ",prev_fxnorm,")");

	step.conv_crit = 0;

	// (1) strong acceptance: the only check when niter==0 XXX hnorm is
	//     always 0 when niter==0
	bool fx1 = (step.fxnorm < specs.epsabs);
	if (fx1 && step.hnorm <= tlstep) {
		step.conv_crit = 1;
		trc.dprint("strong acceptance (1): ",step.fxnorm,
				" <= ",specs.epsabs, " and ", step.hnorm, " <= ", tlstep);
		return rval;
	}
	// if this is before the first iteration just check for strong acceptance
	if (step.niter == 0) {
		trc.dprint("returning 0: niter=0");
		return rval;
	}
	// (2) weak acceptance: small residual or small iterate
	if (fx1 || step.hnorm <= weakEps*step.xnorm) {
		step.conv_crit = 2;
		trc.dprint("weak acceptance (2): norm(f)(", sci1(step.fxnorm),
			") < ", sci1(specs.epsabs), " or hnorm (", sci1(step.hnorm),
			") < ", sci1(weakEps*step.xnorm));
		return rval;
	}

	if (step.niter > 1) {
		// (3) weak acceptance: two consecutive small residuals
		if ((step.fxnorm + prev_fxnorm) <= specs.epsabs &&
				step.hnorm <= 8.0*prev_hnorm) {
			step.conv_crit = 3;
			trc.dprint("weak acceptance (3): ", sci1(step.fxnorm), "+",
				sci1(prev_fxnorm), " < ", sci1(specs.epsabs), " or hnorm(",
				sci1(step.hnorm), " <= 8*prev hnorm(",sci1(prev_hnorm),")");
			return rval;
		}
		// (4) weak acceptance: 2 consecutive small iterates
		//   norm(f_i) < 8*epsabs and
		//   norm(h_i) + norm(h_{i-1}) < epsabs + relerr*norm(x);
		bool fx8 = (step.fxnorm < 8.0*specs.epsabs);
		if (fx8 && (step.hnorm + prev_hnorm) <= tlstep) {
			trc.dprint("weak acceptance (4): ",sci1(step.fxnorm)," <= 8*",
					sci1(specs.epsabs)," and ",sci1(step.hnorm),
					"+",sci1(prev_hnorm), " <= ",sci1(tlstep));
			step.conv_crit = 4;
			return rval;
		}
	}

	// Error conditions: all get conv_crit -1 and error set:
	//
	// (-1) residual not decreasing, small hnorm
	if ((abs(step.fxnorm-prev_fxnorm) < relerr*step.fxnorm) &&
			step.hnorm < relerr) {
		step.conv_crit = -1;
		string error = vastr("residual is not decreasing (",
			prev_fxnorm," -> ",step.fxnorm,") and small h: ",step.hnorm);
		rval = new Issue(error, Kind::noconv);
		trc.dprint("returning ",step.conv_crit,": ",error);
		return rval;
	}
	// Divergence or too many iterations.
	// NOTE: in PITCON the criteria was (STEPX.GT.(FMP*STEPXL+ABSERR))
	// which I think should not have the ABSERR, see ref pg. 122
	if (step.niter > 1)
		fmp = 1.05;
	if (step.niter >= 3 && step.hnorm > fmp*prev_hnorm && step.hnorm > tlstep) {
		// ignore if the residual is decreasing:
		if (step.fxnorm > prev_fxnorm) {
			string error = vastr("divergence after ",step.niter," iterations, h: ",
				step.hnorm," > ",fmp, " * ",prev_hnorm);
			step.conv_crit = -1;
			trc.dprint(error);
			rval = new Issue(error, Kind::divergence);
			return rval;
		}
	}
	if (step.niter > 1
			&& step.fxnorm > (fmp*prev_fxnorm + specs.epsabs)) {
		string error = vastr("divergence: norm(f) ",sci1(step.fxnorm),
				", prev norm(f) ",sci1(prev_fxnorm)," + epsabs ",sci1(specs.epsabs));
		step.conv_crit = -1;
		trc.dprint("returning -1: ",error);
		rval = new Issue(error, Kind::divergence);
		return rval;
	}
	// Too many iterations?
	if (step.niter >= specs.maxiter) {
		string error = vastr("converging too slowly (", specs.maxiter, " iterations)");
		step.conv_crit = -2;
		trc.dprint("returning -2: ",error);
		rval = new Issue(error, Kind::slow);
		return rval;
	}
	trc.dprint("returning ",step.conv_crit);
	return rval;
} // Pac::converged

static bool
testf(const vect& f, const vect& eps) {
// test every abs(f[i]) against eps[i], return false the first
// time abs(f[i]) > eps[i]. The sizes of f and eps may not be the same
// due to constraints.
	Trace trc(2,"testf");
	size_t n = std::min(f.size(), eps.size());
	for (size_t i=0; i<n; i++) {
		if (abs(f[i]) > eps[i]) {
			trc.dprint("returning false: f[",i,"] ",abs(f[i])," > ",eps[i]);
			return false;
		}
	}
	trc.dprint("returning true");
	return true;
}

Issue*
Pac::
convergev(double prev_fxnorm, double prev_hnorm) {
// Check for convergence or non-convergence of this point, and set
// step.conv_crit to one of the following:
//    0  neither convergence or non-convergence
//   -1 residual is not decreasing and small hnorm
//   -2  no convergence after "maxiter" iterations
//    1  strong acceptance
//    2  weak acceptance: small h
//    3  weak acceptance: two consecutive small residuals
//    4  weak acceptance: 2 consecutive small f and h
// Returns a message if step.conv_crit != 0 telling why it
// did (>0) or did not (<0) converge
// References
//   Rheinboldt, Werner C., "Numerical Analysis of Parameterized
//   Nonlinear Equations", John Wiley & Sons, New York, 1986
//   ISBN 0-471-88814-1 (book)
	Trace trc(1,"convergev");
	Issue* rval{nullptr};
	double meps = std::numeric_limits<double>::epsilon();
	double relerr = sqrt(meps);
	double weakEps = 8.0*meps;
	double fmp{2.0};
	// relative error measure, p. 122 ref
	double tlstep = specs.epsabs + relerr*step.xnorm;
	Form sci1;
	sci1.scientific().precision(1);

	trc.dprint("hnorm ",step.hnorm,"(prev ",prev_hnorm, ") fxnorm ",step.fxnorm, "(prev ",prev_fxnorm,")");

	step.conv_crit = 0;

	// (1) strong acceptance: the only check when niter==0 XXX hnorm is
	//     always 0 when niter==0
	bool fx1{false};
	if (specs.minf.empty())
		fx1 = step.fxnorm < specs.epsabs;
	else
		fx1 = testf(f, specs.minf);
	if (fx1 && step.hnorm <= tlstep) {
		step.conv_crit = 1;
		trc.dprint("strong acceptance (1): ",step.fxnorm,
				" <= ",specs.epsabs, " and ", step.hnorm, " <= ", tlstep);
		return rval;
	}
	// if this is before the first iteration just check for strong acceptance
	if (step.niter == 0) {
		trc.dprint("returning 0: niter=0");
		return rval;
	}
	// (2) weak acceptance: small residual or small iterate
	if (fx1 || step.hnorm <= weakEps*step.xnorm) {
		step.conv_crit = 2;
		trc.dprint("weak acceptance (2): norm(f)(", sci1(step.fxnorm),
			") < ", sci1(specs.epsabs), " or hnorm (", sci1(step.hnorm),
			") < ", sci1(weakEps*step.xnorm));
		return rval;
	}

	if (step.niter > 1) {
		// (3) weak acceptance: two consecutive small residuals
		if ((step.fxnorm + prev_fxnorm) <= specs.epsabs &&
				step.hnorm <= 8.0*prev_hnorm) {
			step.conv_crit = 3;
			trc.dprint("weak acceptance (3): ", sci1(step.fxnorm), "+",
				sci1(prev_fxnorm), " < ", sci1(specs.epsabs), " and hnorm(",
				sci1(step.hnorm), " <= 8*prev hnorm(",sci1(prev_hnorm),")");
			return rval;
		}
		// (4) weak acceptance: 2 consecutive small iterates
		//   norm(f_i) < 8*epsabs and
		//   norm(h_i) + norm(h_{i-1}) < epsabs + relerr*norm(x);
		bool fx8 = (step.fxnorm < 8.0*specs.epsabs);
		if (fx8 && (step.hnorm + prev_hnorm) <= tlstep) {
			trc.dprint("weak acceptance (4): ",sci1(step.fxnorm)," <= 8*",
					sci1(specs.epsabs)," and ",sci1(step.hnorm),
					"+",sci1(prev_hnorm), " <= ",sci1(tlstep));
			step.conv_crit = 4;
			return rval;
		}
	}

	// Error conditions: all get conv_crit -1 and error set:
	//
	// (-1) residual not decreasing, small hnorm
	if ((abs(step.fxnorm-prev_fxnorm) < relerr*step.fxnorm) &&
			step.hnorm < relerr) {
		step.conv_crit = -1;
		string error = vastr("residual is not decreasing (",
			prev_fxnorm," -> ",step.fxnorm,") and small h: ",step.hnorm);
		rval = new Issue(error, Kind::noconv);
		trc.dprint("returning ",step.conv_crit,": ",error);
		return rval;
	}
	// Divergence or too many iterations.
	// NOTE: in PITCON the criteria was (STEPX.GT.(FMP*STEPXL+ABSERR))
	// which I think should not have the ABSERR, see ref pg. 122
	if (step.niter > 1)
		fmp = 1.05;
	if (step.niter >= 3 && step.hnorm > fmp*prev_hnorm && step.hnorm > tlstep) {
		// ignore if the residual is decreasing:
		if (step.fxnorm > prev_fxnorm) {
			string error = vastr("divergence after ",step.niter," iterations, h: ",
				step.hnorm," > ",fmp, " * ",prev_hnorm);
			step.conv_crit = -1;
			trc.dprint(error);
			rval = new Issue(error, Kind::divergence);
			return rval;
		}
	}
	if (step.niter > 1
			&& step.fxnorm > (fmp*prev_fxnorm + specs.epsabs)) {
		string error = vastr("divergence: norm(f) ",sci1(step.fxnorm),
				", prev norm(f) ",sci1(prev_fxnorm)," + epsabs ",sci1(specs.epsabs));
		step.conv_crit = -1;
		trc.dprint("returning -1: ",error);
		rval = new Issue(error, Kind::divergence);
		return rval;
	}
	// Too many iterations?
	if (step.niter >= specs.maxiter) {
		string error = vastr("converging too slowly (", specs.maxiter, " iterations)");
		step.conv_crit = -2;
		trc.dprint("returning -2: ",error);
		rval = new Issue(error, Kind::slow);
		return rval;
	}
	trc.dprint("returning ",step.conv_crit);
	return rval;
} // Pac::convergev

//------------------------------------------------------------------
//  Compute a reasonable stepsize for the next continuation step
//  based on information at the current solution.
//  This is probably the most important function in the
//  continuation method. Most of the technical stuff is in
//  dtau(). It is based on the work of W.C. Rheinboldt.
//
// Public interface:
//   double rheinboldt (Point const* prev, Point* current)
// 
//  References in reverse chronological order:
//   1) Rheinboldt, W.C., and Burkardt, John, PITCON 6.1,
//      www.netlib.org/contin/dpcon61.f (Fortran source code for
//      PITCON, the method described in [2])
//
//   2) Rheinboldt, W.C., Numerical Analysis of Parameterized
//      Nonlinear Equations, The University of Arkansas Lecture
//      Notes in the Mathematical Sciences Vol. 7,
//      John Wiley & Sons, 1986.
// 
//   3) Rheinboldt, W.C., and Burkardt, J.V., "A Locally Parameterized
//      Continuation Process", ACM Trans Math Software, Vol 9 No. 2,
//      June 1983, pp 215-235.
//
//   4) Den Heijer, C. and Rheinboldt, W.C., "On Steplength Algorithms
//      for a Class of Continuation Methods", SIAM J. Numer. Analysis,
//      Vol. 18 No. 5, Oct 1981, pp 925-948
//
//
//
// Names used in Step for some variables from ref. 3:
//   deltak      (\delta_k eqn. 4.3) total corrector distance
//   omegatilde  (\tilde{\omega} eqn. 4.10) hnorm/deltak
//   thetak      (\theta_k table 1 & eqn 4.5) convergence quality
//   wk          (2-norm of w^k, eqn. 4.16) curvature between previous
//               and current point
//   alphak      (\alpha_k, eqn. 4.16) angle between prev and current tangents
//   gammak      (\gamma_k, eqn 4.17a,b) predicted curvature between current
//               point and next point
//   hk1         (h_k^1, eqn. 4.21) tentative predicted stepsize
//   hk2         (h_k^2, eqn. 4.22, 4.23) predicted stepsize before truncation
//   hk          (h_k, eqn. 4.23) predicted optimal stepsize
// In some places the relevant line number in the dpcon61.f source code
// are given
//------------------------------------------------------------------

Issue*
Pac::
rheinboldt (const Pac& from) {
// Compute a reasonable stepsize for the next continuation step
// from information at the current and previous solutions.
// Based mostly on Ref. 3, and 1; equation numbers refer to ref 3
// Throws:
//   stalled: tracking stalled with small stepsize
//   retry:  large angle between tangents: retry the step
//                with retry.stepsize
	Trace trc(1,"rheinboldt");
	Issue* rval{nullptr};
	size_t nx = x.size();

	// step.coord duplicates Pac::coord for plotting
	step.coord = coord;

	// theta_k: the quality of the correction
	step.quality = step.convergence_quality(specs.maxiter);

	// total # of stepsize reductions getting to this
	int stepred = step.red_conv + step.red_angle + step.red_dsc;
	trc.dprint("stepsize reductions: conv ",step.red_conv,", angle ",step.red_angle,", dsc ",step.red_dsc);

	// 2-norm of the secant
	step.secant = vector_distance(nx, &x[0], &from.x[0]);

	// how does secant compare with previous stepsize? should be close
	double eps = std::numeric_limits<double>::epsilon();
	double ps = from.stepsize;
	if (ps > eps)
		step.secssratio = step.secant/ps;
	else
		step.secssratio = 1.0;

	trc.dprint("secssratio = ",step.secant," / ", ps," = ",step.secssratio);
	if (abs(coord) > 2 && (step.secant < specs.minstepsize &&
			step.secssratio < 0.0001))
		rval = new Issue(vastr("no progress: secant ",step.secant,
				", step size ",ps), Kind::stalled);

	// alphak: angle between the previous and current tangents (4.16)
	double tt;
	blas_dot(nx, &from.tan[0], 1, &tan[0], 1, tt);
	if (tt >= 1.0)
		step.alphak = 0.0;
	else
		step.alphak = std::acos(tt);
	string tandiff = vector_diff(from.tan,tan);
	trc.dprint("alphak = ",step.alphak);

	// compute the predicted stepsize (step.hk) for the next step
	step.dtau(from, *this, stepred);

	// return an Issue if angle is large
	constexpr double deg2rad{atan(1.0)/45.0};
	if (std::abs(coord) > 3 && step.alphak > specs.maxangle) {
		string exc{vastr("large angle (",step.alphak/deg2rad," deg)")};
		rval = new Issue(exc, Kind::angle);
		trc.dprint("returning issue: ",exc);
		return rval;
	}

	// ratio of the stepsize to the previous stepsize
	if (ps > 0.0)
		step.stepratio = stepsize/ps;
	else
		step.stepratio = 1.0;

	trc.dprint("returning predicted stepsize ",step.hk);
	return rval;
} // rheinboldt

#ifdef NEVER // unused
const vector<pair<string,string>>&
Pac::
parnames() { return parnames_; }
#endif // NEVER // unused


// Parameter names for class Step, used for printing in operator<<
vector<pair<string, string>>
Step::parnames {
	{"deltak","corrector dist"},				// dist between predicted & corrected (4.3)
	{"rcond","recip cond of Jacobian"},		// repciprocal cond of Jacobian
	{"hnorm","last correction norm"},		// norm of last correction
	{"omegatilde","hnorm/deltak"},			// ratio hnorm/deltak
	{"quality","thetak 4.5"},					// thetak (4.5 & table 1)
	{"secant","secant norm"},					// norm of secant between this and prev solution
	{"alphak","tan angle (rad)"},     		// angle (rad) between this and prev tangents
	{"wk","curvature"},							// curvature between current and prev
	{"gammak","predicted curvature"},	  	// predicted curvature on next step
	{"epsk","starting error"},					// truncated ideal starting error
	{"fxnorm","fx norm"},						// 2-norm of f(x)
	{"xnorm","x norm"},							// 2-norm of x
	{"secssratio","secant/stepsize"},		// secant/stepsize (should be near 1)
	{"stepratio","stepsize/prev"},			// stepsize/prev stepsize
	{"hk1","1st est stepsize"},				// first est of stepsize
	{"hk2","2nd est stepsize"},				// second est of stepsize, before truncation
	{"hk","predicted stepsize"},				// predicted optimal stepsize
	{"coord","steps from the origin"},		// duplicates Pac::coord for plotting
	{"niter","# iterations"},					// number of iterations to converge
	{"conv_crit","convergence crit"},		// see Pac::convergence()
   {"red_conv","# reductions: converg."},	// no convergence
   {"red_angle","# reductions: angle"},	// large angle between tangents
   {"red_dsc","# reductions: dsc"},			// Dsc: determinant sign change
	{"svd","qr or svd 0/1"}};					// qr/svd: 0/1

double
Step::
dtau(const Pac& from, Pac& to, int stepred) {
//------------------------------------------------------------------
// Predictor stepsize computation.
// This routine was adapted from the PITCON stepsize routine (setstp)
// Equation numbers and variable names are from ref 3; ref 2 has the
// same equations but are not all numbered.
//------------------------------------------------------------------
	Trace trc(1,"Step::dtau");
	constexpr double meps = std::numeric_limits<double>::epsilon();
	double growthfac = to.specs.growthfac;
	double curvaturefac = to.specs.curvaturefac;
	double reductionfac = to.specs.reductionfac;
	double minstepsize = to.specs.minstepsize;
	double maxstepsize = to.specs.maxstepsize;

	// if the secant is zero, return with stepsize=??
	if (this->secant <= 0.0) {
		*this = from.step;
		trc.dprint("returning prev step ",from.stepsize);
		return from.stepsize;
	}
	
	// eqn 4.16 (ref 3): alphak is the angle between the previous and current
	// tangents; lower limit on angle (dpcon61.f 2628) is arccos(1-meps)
	constexpr double minangle = acos(1.0 - meps);
	if (this->alphak < minangle) {
		trc.dprint("replacing alphak (",this->alphak,") with ",minangle);
		this->alphak = minangle;
	}

	// eqn 4.16 ref 3: wk, the curvature somewhere between this
	// and the previous step
	this->wk = 2.0*curvaturefac*abs(sin(alphak/2.0))/secant;
	double prevwk{0.0};
	prevwk = from.step.wk;
	if (prevwk <= 0.0)
		prevwk = this->wk;
	trc.dprint("curvature (wk) = 2*curvaturefac*sin(alphak/2)/secant = 2*", curvaturefac,"*sin(",alphak,"/2)/", secant," = ",wk,", prev wk ",prevwk);

	// eqn 4.17a ref 3 (dpcon61.f:2636)
	// simple linear extrapolation to predict curvature during
	// next step when info is avail at previous point.
	gammak = wk;
	if (from.step.secant > 0.0)
		gammak += secant*(wk - prevwk)/(secant + from.step.secant);

	// Eqn 4.17b: gammak can be small, zero, or negative;
	// do not let it fall below gammaMin = max(0.001, 0.01/secant)
	// per dpcon61.f:2639
	// double gammaMin = std::max(0.001, 0.01/secant);
	double gammaMin = 0.001;
	if (gammak < gammaMin) {
		trc.dprint("increasing gammak from ",gammak," to ",gammaMin);
		gammak = gammaMin;
	}
	trc.dprint("gammak = ",gammak);

	// Eqn 4.20 (Ref 3): epsk, the desired starting error for the next step.
	// Nominally epsk is deltak_star (Eqn 4.11), the quality times the corrector
	// distance, truncated to lie between 0.01*secant and secant
	double deltak_star = quality*deltak;
	double eps_min{0.01};
	if (deltak_star <= eps_min*secant) {
		this->epsk = eps_min*secant;
	} else if (deltak_star >= secant) {
		this->epsk = secant;
	} else {
		this->epsk = deltak_star;
	}
	trc.dprint("deltak_star = quality*distance = ",deltak_star,", epsk = ",this->epsk);

	// Eqn 4.21 (Ref 3): tentative step h_k^1. If convergence distance was
	// zero, set estimate to growthfactor*secant
	if (deltak == 0.0) {
		hk1 = growthfac*secant;
		trc.dprint("corrector dist = 0 so use growthfactor*secant = ", growthfac,"*",secant," = ",hk1);
	} else {
		hk1 = sqrt(2.0*this->epsk/gammak);
		trc.dprint("hk1 = sqrt(2*epsk<",this->epsk,">/gammak<",gammak,">) = ",hk1);
	}

	// if there were stepsize reductions limit hk1 to
	//  (growthfactor-1)*secant/2  (dpcon61.f:2654)
	// slightly different than Ref 3 comment in paragraph between
	// 4.22 and 4.23 (hk1 <= secant) but the same if growthfactor=3
	if (stepred > 0) {
		double gf = std::max(growthfac-1.0, 1.0);
		hk1 = std::min(hk1, gf*secant/2.0);
		trc.dprint(stepred," step reductions: trunc hk1 to ",hk1);
	}

	// Eqn 4.22: hk2 (dpcon61.f:2657)
	// adjust the step to account for the estimated curvature in
	// the continuation parameter; since we don't use local-parameterization
	// like dpcon61 does, find the largest component of the current
	// tangent and use this. Ref 2 suggests doing this in double precision;
	// using long double seems to make very little diff but do it anyway
	size_t nx = to.x.size();
	int ipc = blas_isamax(nx, &to.tan[0], 1) - 1;  // 0b
	long double tk = to.tan[ipc];
	long double tkm1 = from.tan[ipc];
	trc.dprint(" using tan[",ipc,"] = ",tk,", prev = ",tkm1);
	long double lhk1{hk1};
	long double lhk2{lhk1};
	long double lsecant{secant};
	long double one{1.0};
	long double two{2.0};
	if (tk != 0.0)
		lhk2 = lhk1*(one + lhk1*(one - tkm1/tk)/(two*lsecant));
	hk2 = lhk2;
	trc.dprint("hk2 = ",hk2,", hk1 = ",hk1);

	// ... a couple of lambdas for checking limits on hk2:
	// Eqn 4.23, dpcon61.f:2663
	auto hk2min = [&](double h, double low) {
		Trace trc(1,"hk2min");
		if (h > low) {
			trc.dprint("decreasing hk2 to ",low);
			h = low;
		}
		return h;
	};
	auto hk2max = [&](double h, double hi) {
		Trace trc(1,"hk2max");
		if (h < hi) {
			trc.dprint("increasing hk2 to ",hi);
			h = hi;
		}
		return h;
	};
	trc.dprint("truncate to min/maxstepsize: ",minstepsize, ", ",maxstepsize);
	// hk2 = std::max(hk2, secant/reductionfac);
	hk = hk2max(hk2, secant/reductionfac);
	// hk2 = std::min(hk2, secant*growthfac);
	hk = hk2min(hk, secant*growthfac);
	// hk2 = std::max(hk2, minstepsize);
	hk = hk2max(hk, minstepsize);
	// hk2 = std::min(hk2, maxstepsize);
	hk = hk2min(hk, maxstepsize);
	trc.dprint("returning stepsize ",hk);
	return hk;
} // dtau

double
Step::
convergence_quality (int maxiter) {
// Compute the convergence quality for Newton's method (not Modified).
// Reference 3 gives values for it in Table 1 as a function of
// $\tilde{\omega}$, the ratio of the last correction to total correction
// (hnorm/deltak), Eqn 4.2 ref 3. The Pitcon code does something
// different; the 2 methods are contained in table1() and coqual()
// which one is used is determined by "usecoqual":
	bool usecoqual{true};
// the 2 methods are plotted in main() below
	Trace trc(1,"convergence_quality");
	double rval{0.0};
	double eight{8.0};
	double thgie{0.125};
	int aveIter = maxiter/2 - 1;
	double meps{std::numeric_limits<double>::epsilon()};

	trc.dprint("niter ",niter,", max iter ",maxiter,", hnorm ",hnorm, ", distance ",deltak);

	if (niter <= 1 || deltak <= 8.0*meps)
		rval = eight;
	else if (niter == aveIter)
		rval = 1.0;
	else if (niter >= maxiter)
		rval = thgie;
	if (rval > 0.0) {
		trc.dprint("returning ",rval);
		return rval;
	}

	omegatilde = hnorm/deltak;
	if (usecoqual) {
		rval = coqual(niter, omegatilde);
	} else {
		rval = table1(niter, omegatilde);
	}
	trc.dprint("returning ",rval);
	return rval;
}

std::ostream&
operator<<(std::ostream& s, const Step& t) {
	// enclose the list in {}
	s << "{\n";
#ifdef NEVER
	const double* spar = &t.deltak;	// just point to the first element...
	for (auto& pi : t.parnames)
		s << get<0>(pi) << " (" << get<1>(pi) << "): " << *spar++ << endl;
#endif // NEVER
	s << "}\n";
	return s;
}

double
Step::
coqual (int niter, double omegatilde) {
// translation of subroutine coqual in dpcon61.f
	int aveIter{4};
	double rval{0.125};
	int itop = pow(2, niter - aveIter);
	int ibot = pow(2, niter - 1) - 1;
	double expo = 1.0/(double)ibot;
	double base = pow(omegatilde, expo);
	double top = base + 1.0 + (1.0/base);
	double term = pow(base, itop);
	double bot = term + 1.0 + (1.0/term);
	rval = top/bot;
	if (rval > 8.0) rval = 8.0;
	if (rval < 0.125) rval = 0.125;
	return rval;
}

double
Step::
table1(int niter, double omegatilde) {
	double rval{0.125};
	if (niter == 2) {
		if (omegatilde >= 0.8735115 && omegatilde <= 1.0) {
			rval = 1.0;
		} else if (omegatilde >= 0.1531947 && omegatilde < 0.8735115) {
			rval = 0.9043128 - 0.7075675*log(omegatilde);
		} else if (omegatilde >= 0.03191815 && omegatilde < 0.1531947) {
			rval = -4.667383 - 3.677482*log(omegatilde);
		} else if (omegatilde >= 0.0 && omegatilde < 0.03191815) {
			rval = 8.0;
		}
	} else if (niter == 3) {
		if (omegatilde >= 0.4677788 && omegatilde <= 1.0) {
			rval = 1.0;
		} else if (omegatilde >= 0.6970123e-3 && omegatilde < 0.4677788) {
			rval = 0.8516099 - 0.1953119*log(omegatilde);
		} else if (omegatilde >= 0.1980863e-5 && omegatilde < 0.6970123e-3) {
			rval = -4.830636 - 0.9770528*log(omegatilde);
		} else if (omegatilde >= 0.0 && omegatilde < 0.1980863e-5) {
			rval = 8.0;
		}
	} else if (niter == 4) {
		if (omegatilde >= 0.0 && omegatilde <= 1.0) {
			rval = 1.0;
		}
	} else if (niter == 5) {
		if (omegatilde >= 0.333946e-10 && omegatilde <= 1.0) {
			rval = 1.040061 + 0.03793395*log(omegatilde);
		} else if (omegatilde >= 0 && omegatilde < 0.333946e-10) {
			rval = 0.125;
		}
	} else if (niter == 6) {
		if (omegatilde >= 0.1122789e-8 && omegatilde <= 1.0) {
			rval = 1.042177 + 0.04450706*log(omegatilde);
		} else if (omegatilde >= 0.0 && omegatilde < 0.1122789e-8) {
			rval = 0.125;
		}
	} else if (niter >= 7) {
		if (omegatilde >= 0.0 && omegatilde <= 1.0) {
			rval = 0.125;
		}
	}
	if (rval > 8.0) rval = 8.0;
	if (rval < 0.125) rval = 0.125;
	return rval;
}

static double
vector_distance (int n, const double* a, const double* b) {
// Compute the 2-norm distance between two n-vectors
	double rval = 0.0;
	for (int i=0; i<n; i++) {
		double t = a[i] - b[i];
		rval += t*t;
	}
	rval = sqrt(rval);
	return rval;
}

double
Pac::
ecc() {
/*------------------------------------------------------------------
 * Given a non-linear equation f(x) = 0, and f'(x) (the Jacobian)
 * compute the smallest attainable residual norm by assuming x represents
 * an exact solution, so that f(x) = 0, and estimating the interval
 * the residual will lie in if x lies in an interval whose
 * relative width is eps.
 *
 * Input
 *   m        number of rows in f and jac
 *   n        number of elements in x
 *   eps      error estimate for elements of x; if eps <= 0 on input
 *            it defaults to 128*machine eps
 *   jac      (m,n) double matrix f'(x)
 *   x        n double vector
 * Returns: the 2-norm of the estimated interval f(x) lies in
 *------------------------------------------------------------------*/
	Trace trc(1,"ecc");

	int nx = x.size();
	int nf = f.size();
	//!! double eps = 512.0*nx*std::numeric_limits<double>::epsilon();
	double eps = nx*specs.ecc_eps;
	double tol = sqrt(eps);
	double lo{std::numeric_limits<double>::max()};
	double hi{std::numeric_limits<double>::min()};
	double accum{0.0};

	specs.minf = vector<double>(nf,0.0);
	for (int i=0; i<nf; i++) {
		double ri = 0.0;
		for (int j=0; j<nx; j++) {
			double xj = abs(x[j]);
			if (xj < tol) xj = tol;
			double jacij = abs(jac[IJ(i,j,nf)]);
			if (jacij < tol) jacij = tol;
			ri += jacij*xj;
		}
		if (ri < 1.0)
			ri = 1.0;
		double t = eps*ri;
		accum += t*t;
		lo = std::min(lo, t);
		hi = std::max(hi, t);
		specs.minf[i] = t;
	}
	double rval = sqrt(accum);
	double ratio = lo/hi;
	if (specs.ecc < 0 || (specs.ecc == 0 && ratio > tol)) {
		specs.minf.clear();
		trc.dprint("returning ",rval,", not using minf, ratio = ",ratio,", ecc=",specs.ecc);
	} else {
		trc.dprint("using minf, ratio = ",ratio,", ecc = ",specs.ecc);
	}
	return rval;
}

bool
homotopy(vect& x0, Fjacfcn fjac, vect* results,
	vector<pair<string,string>>* params, vect* stepdata) {
// given an initial guess (x0), compute the nearest
// solution to f(x) = 0 where f and x are the same size
// using homotopy:
//   let g(y) = f(x) + (h-1)*f0
//   where f0 = f(x0), y = [ x h ]',  y0 = [ x0 0 ]'
//   so g(y0) = f(x0) - f0 = 0
//   and using continuation to h=1, g(y1) = f(x1) = 0
//   where y1 = [ x1 1 ]' and x1 is the desired solution
//   dg/dh = f0
	Trace trc(1,"homotopy");
	int n = x0.size();
	int m = n+1;		// x0 + homotopy parameter (h)
	bool rval{true};

	vect f0(n, 0.0);
	vect jac(n*n, 0.0);
	vect tan(n, 0.0);
	// compute f0 = f(x0) ...
	fjac(x0, f0, jac);
	trc.dprint("x0 = ",x0,", f0 = ",f0);
	// ... and return if it is small
	double tol{std::numeric_limits<double>::epsilon()};
	double fnorm = blas_snrm2(n, f0.data(), 1);
	if (fnorm < tol)
		return true;

	// continuation on g(y)
	// lambda to compute g(y), gjac(y)
	auto gj = [&](const vect& y, vect& g, vect& gjac) {
		Trace trc(1,"homotopy gj");
		vect x(y.begin(), y.end()-1);
		vect f(x.size(), 0.0);
		fjac(x, f, jac);
		double t = y[n] - 1.0;
		for (int i=0; i<n; i++) {
			g[i] = f[i] + t*f0[i];
			for (int j=0; j<n; j++)
				gjac[i+j*n] = jac[i+j*n];
			gjac[i+n*n] = f0[i];
		}
		trc.dprint("g(",y,") = ",g);
		trc.dprintm(n,m,n,&gjac[0],"g Jacobian");
		return 0;
	};
		
	// create the starting Pac
	Pac y0(n, m, gj);
	// copy n elements of x0 to y0...
	std::copy_n(x0.begin(), n, y0.x.begin());
	// ... then homotopy parameter h is y[n]
	y0.x[n] = 0.0;
	// initial tan = column n of the identity
	y0.tan[n] = 1.0;
	y0.step.hk = 0.0;		// let trace correct the initial y0
	// lambda to constrain h=0
	y0.confcn = [&](const vector<double>& x, vector<double>& c) {
		Trace trc(1,"h=0 constraint");
		int m = x.size()-1;
		c[m] = 1.0;
		return x[m];
	};
	// trace the curve from 0 -> 1 with a lambda Processfcn
	rval = y0.trace([&n,&results,&stepdata](Pac& y0, Pac& y1, Issue* issue) {
			// process the results from trace() for homotopy
			Trace trc(1,"process (homotopy trace)");
			bool rval{true};
			double h = y1.x[n];
			double dh = y1.tan[n];
			trc.dprint("h = ",h);
			// watch for homotopy at the limit (1)
			y1.confcn = nullptr;
			if (abs(h-1.0) < 0.001)
					rval = false;	// tell trace() to quit
			else if (h + y1.stepsize*dh > 1.0) {
				double newsz = (1.0 - h)/dh;
				trc.dprint("decreasing stepsize from ",y1.step.hk," to ",newsz," to hit 1, h: ",h,", dh: ",dh);
				y1.step.hk = newsz;
				// return a constraint function
				y1.confcn = [&](const vector<double>& x, vector<double>& c) {
					Trace trc(1,"h=1 constraint");
					int m = x.size()-1;
					c[m] = 1.0;
					return x[m] - 1.0;
				};
			}
			if (results != nullptr)
				results->insert(results->cend(), y1.x.begin(), y1.x.end());
			if (stepdata != nullptr)
				stepdata->insert(stepdata->end(), y1.step.cbegin(), y1.step.cend());
			return rval;
		});

	// success? copy the first elements of y1 to x0 to return...
	if (rval) {
		std::copy_n(y0.x.begin(), n, x0.begin());
		// ... and set the parameter name & descriptions
		if (params != nullptr) {
			for (int i=0; i<n; i++) {
				string nm{vastr("y",i+1)};
				params->push_back({nm,nm});
			}
			params->push_back({"homotopy","lambda"});
		}
	}
	return rval;
} // homotopy

bool
Pac::
trace (Processfcn process, vect* projectee, int direction) {
// trace a curve from "this" until "process" returns 0
	Trace trc(1,"Pac::trace" );
	bool rval{true};
	int maxsteps = this->specs.maxsteps;
	int stat{0};

	// correct "this" iff the predicted stepsize hk==0
	if (this->step.hk == 0.0) {
		if (projectee != nullptr)
			this->tan = *projectee;
		Issue* issue{nullptr};
		if (confcn != nullptr)
			issue = this->corrector(confcn);
		else
			issue = this->corrector();
		stat = process(*this,*this,issue);
		// if process returns 0 let the caller reverse the tangent
		// but if start.tan is not zero it will determine the direction
		if (stat == 0) return true;
	}

	Pac x0 = *this;
	if  (x0.step.hk == 0.0)
		x0.step.hk = this->specs.initial_stepsize;
	// extend from x0 -> x1 on each step
	Pac x1{x0};

	// continuation steps; a step needs to be repeated if process
	// returns <0 assuming process has decreased x0 predicted stepsize (step.hk)
	for (int i=0; i<maxsteps; i++) {
		x1.coord = x0.coord + direction;
		x1.step = {};		// set defaults
		while(true) {
			// continuation x0 -> x1
			Issue* issue = x0.sigue(x1,projectee);
			// process results ...
			stat = process(x0,x1,issue);
			// ... if return is +/0/- continue/quit/retry
			if (stat == 0)
				return rval;
			else if (stat > 0)
				break;
		}
		x0 = x1;
	}
	return rval;
}	// trace

int
nleig(complex<double>& s, vector<complex<double>>& x, Dfcn dmat,
	vect* results, vector<pair<string,string>>* params, vect* stepdata) {
// given an estimate of an eigenvalue/eigenvector for the complex nonlinear
// eigenvalue problem D(s)x = 0, refine the estimate using homotopy
// let g(y) = | D(s)x | = 0,  y = [ x s ]'
//            | x[k]-1 |
	Trace trc(1,"nleig");
	size_t n = x.size();
	vect y(2*n+2, 0.0);
	trc.dprint("s = ",s,", x = ",x);
	for (size_t i=0; i<n; i++) {
		y[2*i] = x[i].real();
		y[2*i+1] = x[i].imag();
	}
	// normalize y: y'y = 1
	double ynorm = blas_snrm2(2*n, &y[0], 1);
	if (ynorm == 0.0)
		throw runtime_error("zero x passed to nleig");
	blas_scal(2*n, 1.0/ynorm, &y[0], 1);

	// find the largest component (0 based) of x, set its phase to zero
	int phase_idx = blas_isamax(n, &x[0], 1) - 1;

	// last 2 components: s
	y[2*n] = s.real();
	y[2*n+1] = s.imag();

	// the fjac arg to homotopy is a lambda, defined here with access to dmat
	bool rval = homotopy(y,
			[&dmat,&phase_idx](const vect& y, vect& g, vect& jac) {
				// lambda to evaluate g(y) and g'(y) = jac where
				// y = [x s]' and g = [ f(x) x'x-1 x_k ]'
				Trace trc(1,"nleig fjac");
				size_t m = y.size();		// 2n + 2
				size_t n = (m - 2)/2;
				vector<double> Dx(2*n,0.0);			// D*x
				vector<double> Dxsigma(2*n,0.0);	// d(Dx)/dsigma
				vector<double> Dxfreq(2*n,0.0);		// d(Dx)/dfreq
				vector<double> Dxx(2*n*2*n,0.0);		// d(Dx)/dx
				vector<complex<double>> x(n,0.0);
				for (size_t i=0; i<n; i++)
					x[i] = complex<double>(y[2*i],y[2*i+1]);
				complex<double> s(y[2*n],y[2*n+1]);

				// zero g & jac
				std::fill(g.begin(), g.end(), 0.0);
				std::fill(jac.begin(), jac.end(), 0.0);

				// get Dx, Dxsigma, Dxfreq & Dxx
				dmat(s, x, Dx, Dxsigma, Dxfreq, Dxx);
				// f(x) = Dx
				for (size_t i=0; i<2*n; i++)
					g[i] = Dx[i];

				// last 2 equations are normalization: y'y - 1 & x[phase_idx] real
				blas_dot(2*n, &y[0], 1, &y[0], 1, g[2*n]);
				g[2*n] -= 1.0;
				g[2*n+1] = y[2*phase_idx+1];	// set imag(x[phase_idx]) = 0

				// copy Dxx to first 2n columns of jac
				double* dxx = Dxx.data();
				double* jacd = jac.data();
				for (size_t j=0; j<2*n; j++)
					blas_copy(2*n, &dxx[j*2*n], 1, &jacd[j*m], 1);
				// copy Dxsigma & Dxfreq to last 2 columns XXX just jac[2*n*m] to dmat?
				blas_copy(2*n, Dxsigma.data(), 1, &jac[2*n*m], 1);
				blas_copy(2*n, Dxfreq.data(), 1, &jac[(2*n+1)*m], 1);
				// last 2 eqns: normalization
				blas_copy(2*n, &y[0], 1, &jac[2*n], m);
				blas_scal(2*n, 2.0, &jac[2*n], m);
				jac[(2*n+1)+(2*phase_idx+1)*m] = 1.0;		// phase row 2*n+1, col 2*phase_idx+1

				trc.dprint("x = ",x);
				trc.dprintm(n,n,n,Dx.data(),"Dx");
				trc.dprint("g(",y,") = ",g);
				trc.dprintm(m,m,m,&jac[0],"Jacobian");
				return 0;
			},
			results, params, stepdata);

	if (rval) {
		s = complex<double>(y[2*n], y[2*n+1]);
		for (size_t i=0; i<n; i++)
			x[i] = complex<double>(y[2*i], y[2*i+1]);
	}
	//!! cout << "homotopy returned " << rval << ", s = " << s << ", x = " << x << endl;
	return 0;
}	// nleig

// check the Jacobian by finite differencing using Ridder's technique
#include "nrc.h"

Issue*
Pac::
check_jac(double* xmin, double* xmax,
		Issue* reason, vector<string>* pnames, double tol) {
// pnames:   names of the parameters in x to check
// Just call the vector version
	Trace trc(1,"check_jac(Pac) ", reason);
	return checkjac(x, f, jac, fjac, xmin, xmax, reason, pnames, tol);
}

Issue*
checkjac(vect& x, vect& f, vect& jac, Fjacfcn fjac, double* xmin, double* xmax,
		Issue* reason, vector<string>* pnames, double tol) {
// pnames:   names of the parameters in x to check XXX unused??
// tol:      Largest allowable error in Jacobian
	Trace trc(1,"check_jac ", reason->msg);
	// Jacobian is (nf,nx)
	size_t nx = x.size();
	size_t nf = f.size();
	assert(jac.size() == nx*nf);
	double dfnorm, exnorm;
	double rdiff, maxRdiff = 0.0;
	double reldiff, maxReldiff = 0.0;
	double eps = 128.0*sqrt(std::numeric_limits<double>::epsilon());
	double jij, fdjij;
	double t;
	double xminj, xmaxj;
	size_t i, j;
	int maxj = -1;
	int maxreli = -1, maxrelj = -1;
	std::ostringstream os;
	std::ostringstream erros;

	trc.dprint("nx ",nx,", nf ",nf," <",reason->msg,">");

	// Accumulate messages regarding the check, throw an exception
	// if there is an error or the jacobian has errors
	erros << "checking the Jacobian ";
	if (reason != nullptr)
		erros << '(' << reason << ") at: " << x;
	erros << endl;

	vector<double> anjac(nx*nf, 0.0);
	vector<double> diffJac(nx*nf, 0.0);
	vector<double> relErrJac(nx*nf, 0.0);
	vector<double> relNormErrJac(nx*nf, 0.0);

	// Compute the analytical Jacobian
	if (fjac (x, f, anjac) != 0) {
		erros << "could not evaluate the jacobian for checking";
		trc.dprint("throwing exception \"", erros.str(),"\"");
		throw runtime_error(erros.str());
	}

	// Compute each requested column of the Jacobian using finite differencing
	vector<double> df(nf, 0.0);
	for (size_t j=0; j<nx; j++) {
		for (i=0; i<nf; i++)
			df[i] = 0.0;

		double del = std::max(0.02*abs(x[j]), 1.0);
		if (xmin)
			xminj = xmin[j];
		else
			xminj = x[j] - del;
		if (xmax)
			xmaxj = xmax[j];
		else
			xmaxj = x[j] + del;
		trc.dprint("-----------------> checking column ",j,", [",xminj,":",xmaxj,"]");
		double err;
		if (!NRC::fdapprox (nf, x[j], xminj, xmaxj, &df[0],
			[&x,&f,&jac,&fjac,&j](int n, double t, double* y) {
			// lambda to evaluate y(t) for ridders
				Trace trc(2,"fdfcn col ",j," 0b");
				// save and restore the current value of x
				double save = x[j];
				x[j] = t;
				// compute the jacobian and function
				fjac(x, f, jac);
				// copy f to y
				blas_copy(f.size(), f.data(), 1, y, 1);
				// restore x
				x[j] = save;
				trc.dprint("y(",t,") = ",f);
			}, err)) {
			//!! flaps::info("Jacobian check failed (",reason->msg,") for ",x);
			return new Issue(vastr("Jacobian check failed (",reason->msg,
				") for ",x), reason->kind);
		}
		dfnorm = blas_snrm2 (nf, &df[0], 1);
		exnorm = blas_snrm2 (nf, &anjac[j*nf], 1);
		if (exnorm < eps && dfnorm < eps)
			continue;
		blas_copy(nf, &df[0], 1, &diffJac[j*nf], 1);
	}

	// For each (i,j) compute the difference between the
	// analytical and f-d Jacobians relative to the norms
	// of row i and column j
	for (j=0; j<nx; j++) {
		for (i=0; i<nf; i++) {
			jij = anjac[IJ(i,j,nf)];
			fdjij = diffJac[IJ(i,j,nf)];
			double jn = blas_snrm2(nf, &anjac[IJ(0,j,nf)], 1);
			double fdn = blas_snrm2(nf, &diffJac[IJ(0,j,nf)], 1);
			t =  std::max(jn, fdn);
			rdiff = fabs(fdjij - jij)/std::max(0.1, t);
			reldiff = fabs(fdjij - jij)/std::max(0.1, t); // w/o col norm
			relNormErrJac[IJ(i,j,nf)] = rdiff;
			relErrJac[IJ(i,j,nf)] = reldiff;
			if (rdiff > maxRdiff) {
				maxj = j + 1;		// one-based index
				maxRdiff = rdiff;
			}
			if (reldiff > maxReldiff) {
				maxreli = i + 1;				// one-based index
				maxrelj = j + 1;			// one-based index
				maxReldiff = reldiff;
			}
		}
	}

	string maxrelerror = vastr("max relative error in Jacobian (",reason->msg,") = ",
		maxReldiff, " in the (", maxreli, ",", maxrelj, ") term");
	erros << maxrelerror << std::endl;
	trc.dprint("Result of checking Jacobian: \"",erros.str(),"\"");

	os.str("");
	os << "column " << maxj << " (1b): exact and by diff";
	for (size_t i=0; i<nf; i++)
		os << "   " << anjac[IJ(i,maxj-1,nf)] << "   "
			<< diffJac[IJ(i,maxj-1,nf)] << endl;
	trc.dprint(os.str());

		MM::exporter("jac_exact.mm","analytic",&anjac[0],nf,nx);
		MM::exporter("jac_fd.mm","finite_diff",&diffJac[0],nf,nx);

	if (maxReldiff > tol) {
		os.str("");
		os << "jacobian errors rel to col norm (max " << maxRdiff << ")";
		trc.dprintm(nf,nx,nf,relNormErrJac, "jacobian rel errors");
		os.str("");
		os << "jacobian errors rel (max " << maxReldiff << " (" << maxreli << ','
			<< maxrelj << "))";
		trc.dprintm(nf,nx,nf,relErrJac, "analytical jacobian");
		MM::exporter("jac_rel.mm","rel errors",&relErrJac[0],nf,nx);
	}

	// throw an exception if the largest error is too large
	if (maxReldiff > tol) {
		trc.dprint("throwing exception: large error: ",erros.str());
		throw runtime_error(erros.str());
	}

	trc.dprint("returning (",reason->msg,") max rel err ",maxReldiff);
	//!! return maxReldiff;
	return new Issue(vastr(reason->msg," (",maxReldiff,")"), reason->kind);
}

#ifdef MAIN
#undef MAIN


#include "mcheck.h"

int
compare_cq() {
// plot 2 ways of computing convergence quality:
	int npt = 101;
	double meps{std::numeric_limits<double>::epsilon()};
	double om0 = 16.0*meps;
	double del = (1.0-om0)/(double)(npt-1);
	vector<pair<string,string>> params{
		{"omegatilde","hnorm/deltak 4.10"},
		{"table1", "quality from table1"},
		{"coqual", "quality from coqual"}};
	bool append{false};
	vector<double> point;
	vector<double> data;
	// create curves for niter=2,3,4,5,6,7,8
	for (int niter=2; niter<8; niter++) {
		data.clear();
		for (int i=0; i<npt; i++) {
			point.clear();
			double omega = om0 + i*del;
			point.push_back(omega);
			point.push_back(Step::table1(niter,omega));
			point.push_back(Step::coqual(niter,omega));
			data.insert(data.cend(), point.begin(), point.end());
		}
		Apf::exporter(vastr("niter",niter),params,data,"cq.apf",append);
		append = true;
	}
	return 0;
}

int
test_nleig() {
// simple nonlinear eigenvalue problem
	Trace trc(1,"test_nleig");
	int nf{2};
	complex<double> s(0.1,0.5);
	vector<complex<double>> x(nf, 0.0);
	x[0] = 1.0;
	x[1] = 1.0;
	vector<double> results;
	vector<pair<string,string>> params;
	vector<double> stepdata;
	// diagonal D
	nleig(s, x, [&](complex<double> s,
			const vector<complex<double>>& x,
			vector<double>& Dx,
			vector<double>& Dxsigma,
			vector<double>& Dxfreq,
			vector<double>& Dxx) {
			Trace trc(1,"Dfcn");
			// lambda to compute Dx, d(Dx)/dsigma, etc
			int n = x.size();
			vector<complex<double>> D(n*n, complex<double>(0.0));
			D[0] = 1.0 - s;
			D[3] = 2.0 - s;
			trc.dprintm(n,n,n,D.data(),"D matrix");
			vector<complex<double>> cDx(n);
			blas::gemv("n",n,n,1.0,D.data(),n,x.data(),1,0.0,cDx.data(),1);
			Dx[0] = cDx[0].real();
			Dx[1] = cDx[0].imag();
			Dx[2] = cDx[1].real();
			Dx[3] = cDx[1].imag();
			// dD/dsigma = | -1   0 |   dD/dfreq= | -i  0 |
			//             |  0  -1 |             |  0 -i |
			Dxsigma[0] = -x[0].real();
			Dxsigma[1] = -x[0].imag();
			Dxsigma[2] = -x[1].real();
			Dxsigma[3] = -x[1].imag();
			Dxfreq[0] = x[0].imag();
			Dxfreq[1] = -x[0].real();
			Dxfreq[2] = x[1].imag();
			Dxfreq[3] = -x[1].real();
			// Dxx[0] = complex<double>(1.0) - s;
			// Dxx[3] = complex<double>(2.0) - s;
			blas::real_rep(n,n,D.data(),n,Dxx.data(),2*n);
			trc.dprintm(2*n,2*n,2*n,Dxx.data(),"Dxx");
			return  0;
		}, &results, &params, &stepdata);
	trc.dprint("nleig returned s = ",s,", x = ",x);
	// plot the homotopy and stepdata results
	bool append{false};
	if (!results.empty())
		Apf::exporter("nleig_homotopy", params, results, "nleig_homotopy.apf", append);
	if (!stepdata.empty()) {
		Apf::exporter("nleig_stepdata", Step::parnames, stepdata, "nleig_stepdata.apf", append);
	}

	return 0;
}

int
test_ex24() {
// example 2.4 Seydel
// f(x) = | 1 + p*(x1^2 + x2^2 - 1)          |
//        | 10*x2 - p*x2*(1 + 2*x1^2 + x2^2) | = 0
// where p = \lambda
	Trace trc(1,"test_ex24");

	auto fjac = [&](const vect& x, vect& f, vect& jac) {
		// lambda for computing f(x), f'(x) with x=(x1,x2,p)
		// jac: (2,3)
		Trace trc(1,"fjac");
		double p = x[2];
		f[0] = 1.0 + p*(x[0]*x[0] + x[1]*x[1] - 1.0);
		f[1] = 10.0*x[1] - p*x[1]*(1.0 + 2.0*x[0]*x[0] + x[1]*x[1]);
		jac[0] = 2.0*p*x[0];
		jac[1] = -4.0*p*x[1]*x[0];
		jac[2] = 2.0*p*x[1];
		jac[3] = 10.0 - p - 2.0*p*x[0]*x[0] - 3.0*p*x[1]*x[1];
		jac[4] = x[0]*x[0] + x[1]*x[1] - 1.0;
		jac[5] = -x[1]*(1.0 + 2.0*x[0]*x[0] + x[1]*x[1]);
		trc.dprint("f(",x,") = ",f);
		trc.dprintm(2,2,2,&jac[0],"Jacobian");
		return 0;
	};

	// get a starting point for tracing the u-shaped curve in the y1-p plane
	vector<double> y(2,0.0);
	y[0] = 0.5;
	double p0{2.0};		// \lambda in the book
	homotopy(y, [&p0](const vect& y, vect& f, vect& jac) {
		// lambda for computing a start point with p fixed at p0
		// f(y) = 0, y = (x1,x2), jac: (2,2)
		Trace trc(1,"fjh");
		f[0] = 1.0 + p0*(y[0]*y[0] + y[1]*y[1] - 1.0);
		f[1] = 10.0*y[1] - p0*y[1]*(1.0 + 2.0*y[0]*y[0] + y[1]*y[1]);
		jac[0] = 2.0*p0*y[0];
		jac[1] = -4.0*p0*y[1]*y[0];
		jac[2] = 2.0*p0*y[1];
		jac[3] = 10.0 - p0 - 2.0*p0*y[0]*y[0] - 3.0*p0*y[1]*y[1];
		trc.dprint("f(",y,") = ",f);
		trc.dprintm(2,2,2,&jac[0],"Jacobian");
		return 0;
	});
	trc.dprint("homotopy returned ",y);

	// ... then trace the u-shaped curve in the y1-p plane
	int nf{2};
	int nx{3};
	Pac x0(nf, nx, fjac);
	// start at x = (y, p0)
	x0.x[0] = y[0];
	x0.x[1] = y[1];
	x0.x[2] = p0;
	// initial tan in +p direction
	x0.tan[nf] = 1.0;
	x0.step.hk = 0.001;
	vector<double> plotdata;
	// trace until "process" returns false
	x0.trace([nf,&plotdata](Pac& from, Pac& to, Issue* issue) {
		// process: lambda to process results from trace
			if (issue != nullptr)
				cerr << issue << endl;
			if (to.x[nf] > 5.0 || to.x[0] < 0.0) {
				vector<pair<string,string>> params{{"y1","y1"},{"y2","y2"},{"lambda","lambda"}};
				Apf::exporter("ex24",params,plotdata,"ex24.apf",false);
				return false;	// quit
			}
			if ((from.det.coef*to.det.coef) < 0.0) {
				cerr << "determinant sign change between " << from.x[2] << " and " << to.x[2] << endl;
				to.det.equalize(from.det);
			}
			plotdata.insert(plotdata.cend(), to.x.begin(), to.x.end());
			return true;
		}, nullptr);
	return 0;
}

int
test_fr() {
// Freudenstein-Roth function, 2 equations in 3 unknowns:
//   \Vector{f}_1 = x_1 - x_2^3 + 5x_2^2 - 2x_2 - 13 + 34(x_3-1)
//   \Vector{f}_2 = x_1 + x_2^3 + x_2^2 - 14x_2 - 29 + 10(x_3-1)
// Target: (5,4,1)
	int nf{2};
	int nx{3};
	vector<double> plotdata;
	vector<double> stepdata;

	auto fjac = [&](const vect& x, vect& f, vect& jac) {
		// lambda for computing f(x), f'(x) with x=(x1,x2,p)
		// jac: (2,3)
		Trace trc(1,"fjac");
		double x2 = x[1]*x[1];
		double x3 = x2*x[1];
		f[0] = x[0] - x3 + 5.0*x2 - 2.0*x[1] - 13.0 + 34.0*(x[2]-1.0);
		f[1] = x[0] + x3 + x2 - 14.0*x[1] - 29.0 + 10.0*(x[2]-1.0);
		jac[0] = 1.0;
		jac[1] = 1.0;
		jac[2] = -3.0*x2 + 10.0*x[1] - 2.0;
		jac[3] = 3.0*x2 + 2.0*x[1] - 14.0;
		jac[4] = 34.0;
		jac[5] = 10.0;
		trc.dprint("f(",x,") = ",f);
		trc.dprintm(2,2,2,&jac[0],"Jacobian");
		return 0;
	};

	// lambda to process results from pac::trace
	auto process = [nf,&plotdata](Pac& from, Pac& to, Issue* issue) {
			if (issue != nullptr)
				cerr << issue << endl;
			// stop at x[0]=70, x[0]=4
			if (to.x[0] > 70.0 || to.x[0] < 0.0) {
				vector<pair<string,string>> params{{"x1","x1"},{"x2","x2"},{"x3","x3"}};
				Apf::exporter("fr",params,plotdata,"fr.apf",false);
				return false;	// quit
			}
			if (to.coord >= 0)
				plotdata.insert(plotdata.end(), to.x.cbegin(), to.x.cend());
			else
				plotdata.insert(plotdata.begin(), to.x.cbegin(), to.x.cend());
			return true;
	};

	// initial point
	Pac x0(nf, nx, fjac);
	// set some specs
	x0.specs.maxstepsize = 3.0;
	// start at x = (15, -2, 0)
	x0.x[0] = 15.0;
	x0.x[1] = -2.0;
	x0.x[2] = 0.0;
	// initial tan
	x0.tan[0] = 1.0;
	// step.hk 0 causes trace to correct the initial point
	x0.step.hk = 0.0;
	// trace in positive direction until "process" returns false...
	x0.trace(process);
	// ... then in negative direction
	x0.tan[0] = -1.0;
	x0.step.hk = 0.0001;
	x0.trace(process, nullptr, -1);

	// plot the results
	vector<pair<string,string>> params{{"x1","x1"},{"x2","x2"},{"x3","x3"}};
	Apf::exporter("fr",params,plotdata,"fr.apf",false);
	return 0;
}

static int
test_hopf() {
// simple 2-dof problem which shows Hopf bifurcation
// Eqn 2.17 Seydel
//   \dot{y}_1 = -y_2 + y_1 (p - y_1^2 - y_2^2)
//   \dot{y}_2 = y_1 + y_2 (p - y_1^2 - y_2^2)
	int nf{2};
	int nx{3};

	auto fjac = [&](const vect& x, vect& f, vect& jac) {
		// lambda for computing f(x), f'(x) with x=(x1,x2,p)
		// jac: (2,3)
		Trace trc(1,"fjac");
		double p = x[2];
		f[0] = -x[1] + x[0]*(p - x[0]*x[0] - x[1]*x[1]);
		f[1] = x[0] + x[1]*(p - x[0]*x[0] - x[1]*x[1]);
		jac[0] = p;
		jac[1] = 1.0;
		jac[2] = -1.0;
		jac[3] = p;
		jac[4] = x[0];
		jac[5] = x[1];
		trc.dprint("f(",x,") = ",f);
		trc.dprintm(2,2,2,&jac[0],"Jacobian");
		return 0;
	};
	Pac x0(nf, nx, fjac);
	// the only equilibrium is y = (0,0) for all p
	// start at x = (y, p0)
	vector<double> y(nf, 0.0);	// only equilibrium at any lambda
	double p0{-1.0};
	x0.x[0] = y[0];
	x0.x[1] = y[1];
	x0.x[2] = p0;
	// initial tan in +p direction
	x0.tan[nf] = 1.0;
	x0.step.hk = 0.001;
	vector<double> plotdata;
	// trace until "process" returns false
	x0.trace([nf,&plotdata](Pac& from, Pac& to, Issue* issue) {
		// process: lambda to process results from trace
			if (issue != nullptr)
				cerr << issue << endl;
			if (to.x[nf] > 5.0) {
				vector<pair<string,string>> params{{"y1","y1"},{"y2","y2"},{"lambda","lambda"}};
				Apf::exporter("hopf",params,plotdata,"hopf.apf",false);
				return false;	// quit
			}
			if ((from.det.coef*to.det.coef) < 0.0) {
				cerr << "determinant sign change between " << from.x[2] << " and " << to.x[2] << endl;
				to.det.equalize(from.det);
			}
			plotdata.insert(plotdata.cend(), to.x.begin(), to.x.end());
			return true;
		});
	return 0;
}


static void
usage(char const* prog) {
	cerr << "Usage: " << prog << " option\n";
	cerr << "  options:\n";
	cerr << "   -c:   create an apf plot file comparing 2 methods of computing\n";
	cerr << "         convergence quality\n";
	cerr << "   -e:   trace a curve from example 2.4 in cite{seydel1988equilibrium}\n";
	cerr << "   -h:   Hopf bifurcation: Eqn 2.17 in cite{seydel1988equilibrium}\n";
	cerr << "   -n:   solve a trival nonlinear eigenvalue problem\n";
	exit(1);
}

#include "main.h"
int
main(int argc, char** argv) {
	Trace trc(1,"main");

Pac q;
	mtrace();	// run with MALLOC_TRACE=path

	if (argc < 2 || argv[1][0] != '-')
		usage(argv[0]);

	char opt = argv[1][1];

	if (opt == 'c') {
		return compare_cq();
	} else if (opt == 'n') {
		return test_nleig();
	} else if (opt == 'e') {
		// ex24 with p=2 there is a root at (+/- sqrt(0.5), 0)
		return test_ex24();
	} else if (opt == 'f') {
		// Freudenstein-Roth function
		return test_fr();
	} else if (opt == 'h') {
		// simple 2-dof Hopf bifurcation
		return test_hopf();
	} else {
		cerr << "unrecognized option (" << opt << ")\n";
		usage(argv[0]);
	}
	return 0;
}

#endif // MAIN
