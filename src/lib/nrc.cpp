//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Some numerical algorithms taken from "Numerical Recipes in C":
// @article{teukolsky1992numerical,
//   title={Numerical recipes in C},
//   author={Teukolsky, Saul A and Flannery, Brian P and Press, WH and Vetterling, WT},
//   journal={SMR},
//   volume={693},
//   number={1},
//   pages={59--70},
//   year={1992}
// }

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "config.h"
#include "blas.h"
#include "conv.h"
#include "nrc.h"
#include "trace.h"

using namespace std;

double diffNorm (int n, double *a, double *b);
std::string striter (double const* v, int nr);

static void eval (int n, double x, double h, Fdfcn fcn, double* y);

double
NRC::
root(double x0, double x1, double tol, Rootfcn fcn) {
// Ridders' method for finding a root of a function known to
// lie between x0 and x1; from Ref. 1 section 9.2
	Trace trc(1,"NRC::root");
	auto takesign = [](double a, double b) { return (b >= 0.0 ? abs(a) : -abs(a)); };
	double fl = fcn(x0);
	double fh = fcn(x1);
	if (fl == 0.0)
		return x0;
	if (fh == 0.0)
		return x1;
	if (fl*fh > 0.0)
		return -1.0;
		//!! throw runtime_error(vastr("there is no root between ",x0," and ",x1));
	double xl = x0;
	double xh = x1;
	double rval = std::numeric_limits<double>::max();
	int maxiter{50};
	for (int j=0; j<maxiter; j++) {
		double xm = 0.5*(xl + xh);
		double fm = fcn(xm);
		double s = sqrt(fm*fm - fl*fh);
		if (s == 0.0) return rval;
		double xnew = xm + (xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
		if (abs(xnew-rval) <= tol) {
			trc.dprint("returning ",rval,", ",j+1," iterations, f=",fm);
			return rval;
		}
		rval = xnew;
		double fnew = fcn(rval);
		if (fnew == 0.0)
			return rval;
		if (takesign(fm, fnew) != fm) {
			xl = xm;
			fl = fm;
			xh = rval;
			fh = fnew;
		} else if (takesign(fl,fnew) != fl) {
			xh = rval;
			fh = fnew;
		} else if (takesign(fh,fnew) != fh) {
			xl = rval;
			fl = fnew;
		} else
			throw runtime_error("internal error in NRC::root");
		if (abs(xh-xl) < tol) {
			trc.dprint("returning ",rval,", ",j+1," iterations");
			return rval;
		}
	}
	throw runtime_error("NRC::root: no convergence");
}

bool
NRC::
fdapprox (int n, double x, double xmin, double xmax,
		double* y, Fdfcn fcn, double& err) {
//------------------------------------------------------------------
// Estimate the derivative of a real, vector-valued function
//			y = f(x)
// with respect to a scalar (x) using Ridders' implementation
// of Richardson "deferred approach to the limit" and Neville's
// algorithm for constructing the interpolating polynomial.
//
// Adapted from "Numerical Recipes in C" sect. 5.7 pg. 188
//
// Input:
//		n			number of elements in vector y
//		x			(double) scalar value at which the derivative is to
//					be estimated
//		fcn		function object which evaluates y = f(x) at x.
//		         fcn must have
//					   void operator()(int n, double x, double *y)
//					where
//						n		number of elements in y
//						x		the (scalar) value of x where y
//								is to be evaluated
//						y		double vector to receive the output y = f(x)
//		y			(double*) pointer to vector of n doubles
//
// Output:
//		y			(double*) on successful output y contains the
//		         estimate of the derivative.
//		err      estimate of the error in the derivative
//
// Convergence is when the largest change between the current
// estimate and the previous estimate does not add significantly to
// the largest of MAXY and the abs of the input value of y.
//------------------------------------------------------------------
	Trace trc(2,"NRC::fdapprox");
	bool rval{true};
	double big{std::numeric_limits<double>::max()};
	int nred = 0;
	int iter{0};
	int maxiter{20};
	// double safe{2.0};
	double safe{100.0};
	double* G[maxiter][maxiter];
	double oldrelerr{0.0};

	trc.dprint("n ",n,", x ",x,", [",xmin,':',xmax,"]");

	err = big;

	double macheps = std::numeric_limits<double>::epsilon();
	double eps = 10.0*(macheps);
	double reltol = eps;

	// Initialize the G array
	for (int i=0; i<maxiter; i++)
		for (int j=0; j<maxiter; j++)
			G[i][j] = nullptr;

	std::vector<double> ans(n, 0.0);
	G[0][0] = new double[n];

	// Get an initial h which is large enough to cause significant
	// change in the function. Start with very small h, double it until
	// we get reasonable change in the function. If we don't get a
	// reasonable difference return zero for the derivative.
	// Don't exceed 2^10 increase
	double h = std::sqrt(eps)*(1.0 + std::abs(x));
	std::vector<double> fxplus(n, 0.0);
	std::vector<double> fxminus(n, 0.0);
	fcn(n, x, &fxminus[0]);
	double fx1 = blas_snrm2(n, &fxminus[0], 1);
	double hfac{64.0};  // increase h by hfac each step
	for (int i=0; i<10; i++) {
		fcn(n, x+h, &fxplus[0]);
		double fx2 = blas_snrm2(n, &fxplus[0], 1);
		double fxnorm = std::max(fx1, fx2);
		double diff = diffNorm (n, &fxminus[0], &fxplus[0]);
		trc.dprint("h ",h,", fxnorm ",fxnorm,", f(x+h)-f(x-h) ",diff);
		if (diff > 0.01*std::max(fxnorm, 0.1))
			break;
		h *= hfac;
	}
	// Following the advice on pg 189 of ref.:
	// h = std::max(0.1*std::abs(x), sqrt(eps));
	double hmin = std::max(eps*std::abs(x), macheps);
	trc.dprint("using initial h = ",h,", hmin = ",hmin);


	// try the algorithm with smaller initial h if the error
	// from the iter loop is too large
	while(h > hmin) {
		double x0 = x;
		if (x-h < xmin) {
			trc.dprint("bumped into lower limit (",xmin,")");
			x0 = xmin + h;
		}
		if (x+h > xmax) {
			trc.dprint("bumped into upper limit(",xmax,")");
			x0 = xmax - h;
		}
		eval (n, x0, h, fcn, G[0][0]);

		double con{1.4};		// stepsize decreased by con at each step
		for (iter=1; iter<maxiter; iter++) {
			h /= con;
			x0 = x;
			if (x-h < xmin)
				x0 = xmin + h;
			if (x+h > xmax)
				x0 = xmax - h;
			trc.dprint("------------- iteration ",iter," h = ",h," ----------");
			if (G[0][iter] == nullptr)
				G[0][iter] = new double[n];
			eval (n, x0, h, fcn, G[0][iter]);
			trc.dprint("G[0][",iter,"]\n",striter(G[0][iter],n));
			double fac = con*con;
			double errt{0.0};
			// fill in this row of the tableau
			for (int j=1; j<=iter; j++) {
				if (G[j][iter] == nullptr)
					G[j][iter] = new double[n];
				blas_copy (n, G[j-1][iter], 1, G[j][iter], 1);
				blas_scal (n, fac, G[j][iter], 1);
				assert(G[j-1][iter-1]);
				blas_axpy (n, -1.0, G[j-1][iter-1], 1, G[j][iter], 1);
				blas_scal (n, 1.0/(fac-1.0), G[j][iter], 1);
				fac *= con*con;
				trc.dprint("G[",j,"][",iter,"]\n",striter(G[j][iter],n));
				double erra = diffNorm (n, G[j][iter], G[j-1][iter]);
				double errb = diffNorm (n, G[j][iter], G[j-1][iter-1]);
				trc.dprint("diff between g[",j,"][",iter, "] and g[",j-1,"][",iter,"] = ",erra);
				trc.dprint("diff between g[",j,"][",iter, "] and g[",j-1,"][",iter-1,"] = ",errb);
				errt = std::max(erra, errb);
				if (errt <= err) {
					trc.dprint("max diff (",errt,") is < ",err,": copy g[",iter, "][",j,"] to ans");
					err = errt;
					blas_copy (n, G[j][iter], 1, &ans[0], 1);
				}
			}
			errt = diffNorm (n, G[iter][iter], G[iter-1][iter-1]);
			if (errt >= safe*err) {
				trc.dprint("finished: rel diff between g[",iter,',',iter,"] and g[", iter-1,',',iter-1,"] = ",errt," > ",safe,'*',err);
				break;
			}
		}

		double anorm = blas_snrm2(n, &ans[0], 1);
		double relerr = err/std::max(eps, anorm);
		trc.dprint("rel error ",relerr,", ",nred," reductions");
		// converged?
		if (relerr < reltol || (std::abs(relerr-oldrelerr) < 0.001*relerr)) {
			break;
		} else {
			oldrelerr = relerr;
			h /= 2.0;
			nred++;
			trc.dprint("reduce (",nred,") start interval to ",h);
		}
	}

	if (h <= hmin) {
		rval = false;
		trc.dprint("returning false: h (",h,") less than hmin (",hmin,")");
	}

	// copy the last estimate to y
	if (rval)
		blas_copy (n, &ans[0], 1, y, 1);

	for (int i=0; i<maxiter; i++)
		for (int j=0; j<maxiter; j++)
			if (G[i][j] != nullptr)
				delete[]  G[i][j];

	trc.dprint("returning ",rval?"true":"false", ", iter ",iter+1);
	return rval;
}

static void
eval (int n, double x, double h, Fdfcn fcn, double* y) {
	std::vector<double> fxm(n, 0.0);
	fcn(n, x-h, &fxm[0]);
	fcn(n, x+h, y);
	blas_axpy (n, -1.0, &fxm[0], 1, y, 1);
	double t = 1.0/(2.0*h);
	blas_scal (n, t, y, 1);
}

double
diffNorm (int n, double *a, double *b) {
// Returns the infinity norm of the difference between two n-vectors
	double rval{0.0};

	for (int i=0; i<n; i++) {
		double diff = std::abs(a[i] - b[i]);
		if (diff > rval)
			rval = diff;
	}
	return rval;
}


string
striter (double const* v, int nr) {
	ostringstream os;

	for (int i=0; i<nr; i++) {
		os << v[i] << endl;
	}
	return os.str();
}


#ifdef MAIN

#include "conv.h"

int
test_fdapprox() {
	double omega{1.0};
	double pi{flaps::pi};
	int n = 1;

	double y = 0.0;
	double tmin = -2.0*pi/omega;
	double tmax = 2.0*pi/omega;
	double ta = 0.0;
	double tb = pi/omega;
	int nstep = 11;
	double deltat = (tb - ta)/(nstep-1);
	double maxerr{0.0};
	double terr{0.0};
	double rerr{0.0};  // error estimate from ridders

	for (int i=0; i<nstep; i++) {
		double t = ta + i*deltat;
		double exact = omega*cos(omega*t);

		// using a lambda for the function
		NRC::fdapprox(n, t, tmin, tmax, &y,
			[&omega](int n, double t, double *ft) {
				*ft = sin (omega*t);
				cerr << "sin(" << t << ") = " << *ft << endl;
			}, rerr);

		double error = std::abs(exact-y);
		if (error > maxerr) {
			maxerr = error;
			terr = t;
		}
		cout << "deriv of sin(" << omega*t << ") = " << y
			<< ", ridders error est " << rerr << ", actual error " << error << endl;
	}
	cout << "max error = " << maxerr << " at t = " << terr << endl;
	return 0;
}

int
test_root() {
	double x0{0.0};
	double x1 = flaps::pi;
	double r{numeric_limits<double>::min()};
	double tol{0.000001};
	try {
		r = NRC::root(x0, x1, tol, [](double x) { return cos(x); });
	} catch (runtime_error& s) {
		cerr << "caught: " << s.what() << endl;
		return 1;
	}
	cout << "root returned " << r << endl;
	return 0;
}

static void
usage(char const* prog) {
	cerr << "Usage: " << prog << " option\n";
	cerr << "  options:\n";
	cerr << "    -f:   test fdapprox\n";
	cerr << "    -r:   test root\n";
	exit(1);
}

int
main (int argc, char **argv) {

	if (argc  < 2 || argv[1][0] != '-')
		usage(argv[0]);
	char opt = argv[1][1];
	if (opt == 'f')
		return test_fdapprox();
	else if (opt == 'r')
		return test_root();
	else
		usage(argv[0]);

	return 0;
}
#endif
