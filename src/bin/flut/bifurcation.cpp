//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

constexpr int tracelvl{1};
// computation of bifurcation curves using the determinants
// computed while tracking Flutcurves
// References:
// 1) Meyer, E.E., "Continuation and Bifurcation in Linear Flutter Equations",
//    AIAA Journal, V.53, No. 10, 2015, pp 3113-3116
// 2) Allgower, E.L. and Georg, K., "Numerical Continuation Methods",
//    Springer Series in Computational Mathematics, Berlin, 1990

#include <string>

#include "config.h"
#include "blas.h"
#include "flutcurve.h"
#include "matrix.h"
#include "message.h"
#include "nrc.h"
#include "svd.h"

using namespace std;

static string bif_tan(Pac& bifp, vector<double>& t1, vector<double>& t2);
static string alpha_beta (Pac& bifp, const vector<double>& um,
		const vector<double>& vm, const vector<double>& vmp1,
		double& alpha1, double& beta1, double& alpha2, double& beta2);

static void plotTangent (Flutcurve& fc, vector<double>& origin,
	vector<double>& tan, vector<double>& t1, vector<double>& t2,
	string const& plotfile);

int
Flutcurve::
checkdsc(Pac& a, Pac& b) {
// Returns int: -1 if the step should be repeated to localize the dsc.
// The predicted stepsize (hk) in Point "a" has been reduced to put
// the next point just before the dsc so that the next step will have
// a dsc but just beyond it.
	Specs& sp = this->specs;
	int rval{0};

	int nf = a.f.size();
	int nx = a.x.size();

	// compute determinant at point "b"...
	b.det = b.determinant();
	// ... and set the "detsign" parameter if available
	// XXX save both d & exp? det = d*10^exp
	Par* p = this->params.findp("detsign");
	if (p != nullptr)
		p->value(b.det.coef);

	// quit reducing stepsize for determinant sign change k
	// reductions before min stepsize or maxdscred reductions
	double adet = a.det.coef;
	double bdet = b.det.coef;
	double dsc = adet*bdet;

	// no sign change? set b.step.red_dsc to a.step.red_dsc (may not be zero)
	if (dsc >= 0.0) {
		b.step.red_dsc = a.step.red_dsc;
		return rval;
	}

	// determinant has changed signs
	if (!this->atlimit.empty()) {
		flaps::info("ignoring determinant sign change: ",this->atlimit);
		return rval;
	}
	// if the stepsize is large, re-do the step with small stepsize,
	// otherwise localize the bifurcation with NRC::root
	b.det.equalize(a.det);

	// debug tracing for when we encountered bif
	Trace trc(tracelvl,"checkdsc");

	double t = -a.det.coef/(b.det.coef-a.det.coef);
	trc.dprint("predicted bif stepsize factor: ",t);
	double dsctol = a.specs.minstepsize*243.0;
	if (a.step.hk > dsctol) {
		if (t > 0.01) {
			double tstep = 0.7*t*a.step.hk;
				a.step.hk = tstep;
				a.step.red_dsc++;
				trc.dprint("returning -1: retry step with hk=",a.step.hk);
				return -1;	// retry
		}
	}

	// lambda to evaluate the determinant at (1-t)*a + t*b
	Pac bifp{a};
	auto rootfcn = [&](double t) {
		Trace trc(tracelvl,"rootfcn");
		for (int i=0; i<nx; i++) {
			bifp.x[i] = (1.0-t)*a.x[i] + t*b.x[i];
			bifp.tan[i] = (1.0-t)*a.tan[i] + t*b.tan[i];
		}
		bifp.corrector();
		Det det = bifp.determinant();
		det.equalize(a.det);
		trc.dprint("returning f(",t,") = ",det.coef);
		return det.coef;
	};

	double bif = t;

	trc.dprint("det sign change localized between ",a.coord," and ",b.coord," with step ",a.stepsize);

	// recompute bifp in case root did not return the last iteration
	double det = rootfcn(bif);
	trc.dprint("converged det = ",det);

	// compute starting tangents for the bifurcations
	vector<double> t1(nx, 0.0);
	vector<double> t2(nx, 0.0);
	string error = bif_tan(bifp, t1, t2);
	if (!error.empty()) {
		trc.dprint("bif_tan failed: ",error);
		flaps::warning("suspected bifurcation between coord ",a.coord,
				" and ",b.coord," failed: ",error);
		return 0;
	}

	trc.dprintvvv(bifp.tan, t1, t2, "tan, t1, t2");
	string plotfile("biftan");
	plotTangent (*this, bifp.x, bifp.tan, t1, t2, plotfile);
	// compare these 2 tangents with bifp's tangent: one
	// should be close, use the other one
	double tp1;
	blas_dot(nx, bifp.tan.data(), 1, t1.data(), 1, tp1);
	double tp2;
	blas_dot(nx, bifp.tan.data(), 1, t2.data(), 1, tp2);
	trc.dprint("p1 * t1, t2 = ", tp1, ", ", tp2);
//!! #ifdef NEVER // both
	if (abs(tp2) < abs(tp1))
		t1 = t2;
//!! #endif // NEVER // both
	// create a Flutcurve origin to return, using the smaller of t1, t2
	Pac b1(bifp.x, t1, nf, bifp.fjac);
	// give the origin a non-zero stepsize so pac does not try
	// to correct the origin
	b1.stepsize = sp.bif_init;
	this->bif_origins.push_back(b1);
#ifdef NEVER // only one
   Pac b2(bifp.x, t2, nf, bifp.fjac);
   b2.stepsize = 0.00001;
   this->bif_origins.push_back(b2);
#endif // NEVER // only one

	// when tracing bifurcation curves it may be necessary
	// to stop when a dsc occurs
	if (sp.stopatdsc)
		this->atlimit = "determinant sign change encountered";
	return 1;
} // checkdsc

static string
bif_tan(Pac& bifp, vector<double>& t1, vector<double>& t2) {
	Trace trc(tracelvl,"bif_tan");

	int nf = bifp.f.size();
	int nx = bifp.x.size();

	// SVD the jacobian to see how many rank def and to get U and V
   SVD jf(nf, nx, bifp.jac.data());
	trc.dprint("rcond: ",jf.s[nf-1]/jf.s[0]);

	// how many rank deficiencies?
	size_t nrdef = 0;
	//!! double eps(sqrt(std::numeric_limits<double>::epsilon()));
	//!! double tol = 20.0*4.0*eps*jf.s[0];
	double eps(std::numeric_limits<double>::epsilon());
	double tol = 8.0*nf*eps*jf.s[0];
	trc.dprint("rank tol = ",tol);
	//!! trc.dprint("singular values: ",jf.s);
 	for (int i=0; i<nf; i++) {
		if (jf.s[nf-1-i] < tol)
			nrdef++;
		else
			break;
	}

	// often we will not have clean rank deficiency: just use the smallest
	if (nrdef == 0) {
		if (jf.s[nf-1] > sqrt(eps)*jf.s[0])
			return vastr("Jacobian is not rank-deficient: rcond = ",
				jf.s[nf-1]/jf.s[0]);
		nrdef = 1;
	} else if (nrdef != 1) {
		//!! flaps::warning("there are ",nrdef," rank deficiencies - ignoring all but one");
		trc.dprint("there are ",nrdef," rank deficiencies - ignoring all but one");
		nrdef = 1;
	}

	// the Jacobian is (m,m+1) where m = nf
	size_t m = nf;

	// extract from jac svd:
	//    um = U_m            last column of U: (m,1) left nullspace of J
	//    vm = [V_m V_{m+1}]  (m+1,nrdef+1) nullspace (last nrdef+1 columns of V)
	vector<double> um(m);
	blas_copy(m, &jf.u[IJ(0,m-1,m)], 1, &um[0], 1);
	vector<double> vm((m+1));
	vector<double> vmp1((m+1));
	vector<double> work(m+1);
	// first vm is the last column of V, then the nrdef
	// penultimate columns
	blas_copy(m+1, &jf.vt[IJ(m,0,m+1)], m+1, &vm[0], 1);
	blas_copy(m+1, &jf.vt[IJ(m-1,0,m+1)], m+1, &vmp1[0], 1);

	// compute alpha_k, beta_k, k=1,2
	double alpha1, beta1;
	double alpha2, beta2;
	string rval = alpha_beta (bifp, um, vm, vmp1, alpha1, beta1, alpha2, beta2);
	if (!rval.empty())
		return rval;

	// mult t_k = alpha_k*v_m + beta_k*v_{m+1}
	t1 = vm;
	blas_scal(nx, alpha1, &t1[0], 1);
	blas_axpy(nx, beta1, &vmp1[0], 1, &t1[0], 1);
	t2 = vm;
	blas_scal(nx, alpha2, &t2[0], 1);
	blas_axpy(nx, beta2, &vmp1[0], 1, &t2[0], 1);

	// normalize new tangents
	double norm = blas_snrm2(nx, &t1[0], 1);
	blas_scal(nx, 1.0/norm, &t1[0], 1);
	norm = blas_snrm2(nx, &t2[0], 1);
	blas_scal(nx, 1.0/norm, &t2[0], 1);
	return rval;
}

static string
alpha_beta (Pac& bifp, const vector<double>& um,
		const vector<double>& vm, const vector<double>& vmp1,
		double& alpha1, double& beta1, double& alpha2, double& beta2) {
// Compute (alphak,betak), k=1,2 as part of computing bifurcation
// tangents. Returns a string: if !empty the computation failed and
// it gives the message
	Trace trc(tracelvl,"alpha_beta");
	string rval;

	int nf = bifp.f.size();
	int nx = bifp.x.size();
	if ((int)vm.size() != nx)
		throw runtime_error(vastr("vm is ",vm.size(),", nx is ",nx));
	if ((int)vmp1.size() != nx)
		throw runtime_error(vastr("vmp1 is ",vm.size(),", nx is ",nx));

	// define a lambda g(xi1,xi2) (eqn 25, ref 1)
	auto g = [&](double xi1, double xi2) {
		vector<double> xp(bifp.x);
		if (xi1 != 0.0)
			blas_axpy(nx, xi1, &vm[0], 1, &xp[0], 1);
		if (xi2 != 0.0)
			blas_axpy(nx, xi2, &vmp1[0], 1, &xp[0], 1);
		bifp.fjac(xp, bifp.f, bifp.jac);
		double rval;
		blas_dot(nf, &um[0], 1, &bifp.f[0], 1, rval);
		return rval;
	};

	double g00 = g(0.0, 0.0);
	trc.dprint("g00 = ",g00);
	double eps = pow(std::numeric_limits<double>::epsilon(), 1.0/3.0);
	double epssq{eps*eps};
	double a11 = (g(eps,0.0) - 2.0*g00 + g(-eps,0))/epssq;
	double a12 = (g(eps,eps) + g(-eps,-eps) - g(eps,-eps) - g(-eps,eps))/(4.0*epssq);
	double a22 = (g(0,eps) - 2.0*g00 + g(0,-eps))/epssq;
	double disc = a12*a12 - a11*a22;
	trc.dprint("a11 ",a11,", a12 ",a12,", a22 ",a22,", disc ",disc);
	if (disc <= 0.0) {
		rval = "not a simple bifurcation: negative discriminant";
		//!! flaps::warning("suspected bifurcation at\n", this->summarize(bifp.x),
			//!! "\nhas a negative discriminant(", disc, "): not a simple bifurcation");
		return rval;
	}

	disc = sqrt(disc);
	double c1 = -a12 + disc;
	double c2 = -a12 - disc;
	// compute the two (alpha,beta) pair
	alpha1 = c1/sqrt(c1*c1 + a11*a11);
	beta1 = c1/sqrt(c1*c1 + a22*a22);
	alpha2 = c2/sqrt(c2*c2 + a11*a11);
	beta2 = c2/sqrt(c2*c2 + a22*a22);
	trc.dprint("alpha1 ",alpha1,", beta1 ",beta1);
	trc.dprint("alpha2 ",alpha2,", beta2 ",beta2);
	return rval;
}	// alpha_beta

static
void
plotTangent (Flutcurve& fc, vector<double>& origin,
vector<double>& tan,
vector<double>& t1,
vector<double>& t2,
string const& plotfile) {
	Trace trc(2,"plotTangent");
	int nx = origin.size();
	string aid = fc.aid();
	// compute 3 tangents: origin->tan, origin->t1, origin->t2
	vector<double> ot{origin};
	blas_axpy(nx, 1.0, tan.data(), 1, ot.data(), 1);
	vector<double> ot1{origin};
	blas_axpy(nx, 1.0, t1.data(), 1, ot1.data(), 1);
	vector<double> ot2{origin};
	blas_axpy(nx, 1.0, t2.data(), 1, ot2.data(), 1);

	Flutcurve ftan(aid, vastr(fc.cid(),".orig"),fc.vzid());
	//!! ftan.clear_solns();
	//!! ftan.cid("biforigin");
	//!! ftan.params.desc("biforigin");
	ftan.update(origin);
	ftan.back_solns();
	ftan.update(ot);
	ftan.back_solns();

	Flutcurve btan1(aid, vastr(fc.cid(),".t1"), fc.vzid());
	btan1.update(origin);
	btan1.back_solns();
	btan1.update(ot1);
	btan1.back_solns();

	Flutcurve btan2(aid, vastr(fc.cid(),".t2"), fc.vzid());
	btan2.update(origin);
	btan2.back_solns();
	btan2.update(ot2);
	btan2.back_solns();

	vector<Flutcurve*> fcurves;
	fcurves.push_back(&ftan);
	fcurves.push_back(&btan1);
	fcurves.push_back(&btan2);
	Flutcurve::plot(ftan.aid(), plotfile, fcurves);
}
