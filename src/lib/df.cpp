//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "config.h"

#include <mutex>
#include <exception>
#ifdef HAVE_FFTW3_H
#include <fftw3.h>
#endif

#include "conv.h"
#include "df.h"
#include "interp.h"
#include "pset.h"
#include "trace.h"

using namespace std;

Interp* df_build(const string& dfid, const vector<double>& params,
		Qfcn qs_fcn, Ffcn fq_fcn, Interpfcn interp_fcn);
Interp* df_fft(const string& dfid, const vector<double>& params,
		Qfcn qs_fcn, Ffcn fq_fcn, Interpfcn interp_fcn);

// 1) built in functions

// bilinear
static
Ad
bldf (Ad gamma, Ad r) {
// following the notation in meyer2016continuation (and nl2.fp):
// in nl2.fp I define gamma = |q_j|/\delta
	T_(Trace trc(2,"bldf");)
	Ad rval{1.0};
	if (gamma.value() > 1.0) {
		Ad pi{flaps::pi};
		Ad xi = asin(1.0/gamma);
		Ad one(1.0);
		Ad two(2.0);
		Ad gamsq = gamma*gamma;
		rval = r + (two/pi)*(one - r)*(xi + sqrt(gamsq-one)/gamsq);
		if (isnan(rval.value())) {
			throw runtime_error(vastr("bldf nan at gamma ",gamma.value()));
		}
	}
	T_(trc.dprint("gamma ",gamma,", ratio ",r,", returning ",rval);)
	return rval;
}

Ad
df::
blfcn(Ad const& gamma, Ad const& ratio, double smooth) {
// gamma:   absgc/x0
	T_(Trace trc(1,"blfcn");)
	Ad rval{0.0};
	Ad pi(flaps::pi);

	T_(trc.dprint("gamma=",gamma,", ratio=",ratio,", smooth=",smooth);)

	double g{gamma.value()};
	// smoothing between [gm:gp]
	double gp = 1.0 + smooth;
	double gm = 1.0 - smooth;

	// smoothing if gm <= g <= gp
	// let t(g) = (g - gm)/(2*smooth)
	//     t0 = 0, t1 = 1
	// then f'(t) = df/dt = (df/dg)*(dg/dt) = (df/dg)/(dt/dg) = 2*smooth*df/dg
	// Fit a cubic polynomial in terms of t
	// so that for t0 <= t <= t1
	//   f(t) = c3*t^3 + c2*t^2 + c1*t + c0
	// 4 equations:
	//   f(0) = 1 = c0
	//   f'(0) = 0 = c1
	//   f(1) = c3 + c2 + 1 = ft1
	//   f'(1) = 2c2 + 3c3 = ftp1
	// solving for c2 and c3:
	//   c2 = 3*f(1) - f'(1) - 3
	//   c3 = f'(1) - 2*f(1) + 2

	if (g <= gm) {
		rval = Ad(1.0);
	} else if ((gm < g) && (g < gp)) {
		vector<double> c(4);
		Ad gam1(gp);
		gam1.der(0,1.0); // make der(0) df/dgamma
		Ad ft1 = bldf (gam1, ratio);
		c[0] = 1.0;
		c[1] = 0.0;
		c[2] = 3.0*ft1.value() - 2.0*smooth*ft1.der(0) - 3.0;
		c[3] = 2.0*smooth*ft1.der(0) - 2.0*ft1.value() + 2.0;
		T_(trc.dprint("c0=",c[0],", c1=",c[1],", c2=",c[2],", c3=",c[3]);)
		Ad tad = (gamma - gm)/(2.0*smooth);
		rval = c[3]*tad*tad*tad + c[2]*tad*tad
						+ c[1]*tad + c[0];

		T_(trc.dprint("returning smoothed f{",tad,"} = ",rval);)
	} else {
		rval = bldf(gamma, ratio);
	}
	T_(trc.dprint("returning ",rval.value());)
	return rval;
}

Ad
df::
bilinear(pset& plt, Ad const& x0, Ad const& ratio, int gcno) {
// Describing function for a bilinear stiffness - a number between 1 and "ratio"
// which multiplies the nominal stiffness to give an equivalent stiffness.
// The transition between k1 and k2 may be smoothed by creating a parameter
// "blsmooth" in param; for example
//    param {..., blsmooth=0.05", ...}
// will produce a cubic between x0*0.95 and x0*1.05
// Input:
//    x0    (Ad) transition between k1 and k2
//    ratio ratio of the regions: k2/k1
//    gcno  (int) g.c. number (1b) which has the bilinear stiffness
// If x0 is zero:       return 1 regardless of abs(gcno)
// if abs(gcno) <= x0:  return 1
//    abs(gcno) >> x0:  approaches "ratio" asymtotically
	T_(Trace trc(2,"df::bilinear");)

	Ad absgc = plt.absgcv(gcno);
	Ad gamma = absgc/x0;
	Ad rval{0.0};

	double rx0 =  x0.value();

	// smooth the transition between x0-alpha and x0+alpha
	// where alpha = blsmooth*x0. Default: blsmooth=0.1
	double blsmooth{0.1};
	double alpha{0.0};
	const Par* blsmoothp = plt.findp("blsmooth");
	if (blsmoothp != nullptr) {
		blsmooth = blsmoothp->value();
		if (blsmooth < 0.0 || blsmooth > 1.0) {
			throw runtime_error(vastr("blsmooth (", blsmooth,
						") is out of the range 0-1"));
		}
		alpha = blsmooth*rx0;
		static int visit{0};
		if (visit == 0)
			flaps::info("got blsmooth = ",blsmooth,": smoothing between ",
				rx0-alpha, " and ",rx0+alpha," (gamma[",1.0-blsmooth,":",
				1.0+blsmooth, "]) width ",2.0*alpha);
		visit++;
	}
	rval = blfcn(gamma, ratio, blsmooth);

	T_(trc.dprint("returning ",rval);)
	return rval;
}

// freeplay
Ad
df::
freeplay(Ad const& gap, Ad const& absgc, double smooth) {
// Describing function for freeplay - a number between 0 and 1 which
// multiplies the nominal (gapless) stiffness to give an equivalent
// stiffness. The transition between the gap and engagement may
// be smoothed by passing a value for "smooth" between 0 and 1.
// for example smooth=0.05 will produce a cubic between gap*0.95
// and gap*1.05
// Input:
//   gap   (Ad) gap size
//   absgc (Ad) mag of gc which has the gap
//   smooth (double)  smoothing factor between 0-1
// If gap is zero then return 1 regardless of absgc
// if absgc <= gap return 0
//    absgc >> gap -> 1 asymtotically
	T_(Trace trc(2,"freeplay");)
	double x = absgc.value();
	double x0 =  gap.value();

	if (x0 == 0.0) {
		T_(trc.dprint("quick return: zero gap");)
		return Ad(1.0);
	}

	// smooth the transition between x0-alpha and x0+alpha
	// where alpha = smooth*x0
	double alpha{smooth*x0};

	Ad pi(flaps::pi);

	T_(trc.dprint("alpha=",alpha,", gap=",gap,", abs gc=",absgc);)

	double x0p = x0 + alpha;
	double x0m = x0 - alpha;

	if (x <= x0m) {
		T_(trc.dprint("quick return: gc<",x0m);)
		return Ad(0.0);
	}

	// smoothing if x0m <= x <= x0p
	// let t = (x - x0m)/2*alpha,
	//     t0 = 0, t1 = 1
	// then f'(t) = df/dt = df/dx*dx/dt = 2*alpha*df/dx
	// Fit a cubic polynomial in terms of t
	// so that for t0 <= t <= t1
	//   f(t) = c3*t^3 + c2*t^2 + c1*t + c0
	// 4 equations:
	//   f(0) = 0 = c0
	//   f'(0) = 0 = c1
	//   f(1) = c2 + c3 = ft1
	//   f'(1) = 2c2 + 3c3 = ftp1
	// solving for c2 and c3:
	//   c2 = 3*f(1) - f'(1)
	//   c3 = f'(1) - 2*f(1)
	if ((x0m < x) && (x < x0p)) {
// #define CORNER 1
#ifdef CORNER // use Corner
		constexpr double pi{flaps::pi};
		double u = gap.value()/x0p;
		double du = -u/x0p;
		double xi = 2.0*asin(u);
		double dxi = (2.0/sqrt(1.0 - u*u))*du;
		double fx0p = 1.0 - (xi + sin(xi))/pi;
		double dfx0p = -(dxi + cos(xi)*dxi)/pi;
		corner c(x0m, 0.0, 0.0, x0p, fx0p, dfx0p);
		Ad rval = c.eval(absgc);
		T_(trc.dprint("f(x0p) = ",fx0p,", returning smoothed ",rval);)
		return rval;
#else
		// it is more efficient to compute c0,c1,c2,c3 directly
		vector<double> c(4);
		Ad x0pad(x0p);
		x0pad.der(0,1.0); // make der(0) df/dx
		Ad gapad{gap.value()};  // zero derivs in case gap has derivs
		Ad adxi = 2.0*asin(gapad/x0pad);
		Ad ft1 = 1.0 - (adxi + sin(adxi))/pi;
		//!! Ad ft1 = 1.0 - 2.0*asin(gap/x0pad)
		//!! 	- sin(2.0*asin(gap/x0pad)))/pi;
		c[0] = 0.0;
		c[1] = 0.0;
		c[2] = 3.0*ft1.value() - 2.0*alpha*ft1.der(0);
		c[3] = 2.0*alpha*ft1.der(0) - 2.0*ft1.value();
		T_(trc.dprint("c0=",c[0],", c1=",c[1],", c2=",c[2],", c3=",c[3]);)
		Ad tad = (absgc - gap + alpha)/(2.0*alpha);
		Ad rval = c[3]*tad*tad*tad + c[2]*tad*tad
						+ c[1]*tad + c[0];

		T_(trc.dprint("returning smoothed f(",tad,") = ",rval);)
		return rval;
#endif // CORNER : use Corner
	}

  //!! Ad rval = (pi - Ad(2.0)*asin(gap/absgc)
	//!!		- sin(Ad(2.0)*asin(gap/absgc)))/pi;
	Ad xi = 2.0*asin(gap/absgc);
	Ad rval = 1.0 - (xi + sin(xi))/pi;

	T_(trc.dprint("returning ",rval);)
	return rval;
}

Ad
df::
freeplay(pset& plt, Ad const& gap, int gcno) {
// improved gap fcn: takes gcno instead of absgc
// gcno - 1b complex gc number
// evaluate absgc = gcnorm*ev[gcno], and fpsmooth,
// then call freeplay(...Ad absgc)
// The transition between the gap and engagement may
// be smoothed by creating a parameter "fpsmooth" in
// param; for example
// param {..., fpsmooth=0.05", ...}
// will produce a cubic between gap*0.95 and gap*1.05
// or 0.95 <= gamma <= 1.05
	Ad absgc = plt.absgcv(gcno);
	if (gap.value() <= 0.0)
		return 1.0;
	// check for a parameter called "fpsmooth"
	const Par* fpsmoothp = plt.findp("fpsmooth");
	double fpsmooth{0.0};
	if (fpsmoothp != nullptr) {
		fpsmooth = fpsmoothp->value();
		if (fpsmooth < 0.0 || fpsmooth > 1.0) {
			throw runtime_error(vastr("fpsmooth (", fpsmooth,
					") is out of the range 0-1"));
		}
		static int visit{0};
		if (visit == 0) {
			double x0{gap.value()};
			double alpha{x0*fpsmooth};
			flaps::info("got fpsmooth = ",fpsmooth,": smoothing between ",
				x0-alpha, " and ",x0+alpha," (gamma[",1.0-fpsmooth,":",
				1.0+fpsmooth, "]) width ",2.0*alpha);
		}
		visit++;
	}

	Ad rval = df::freeplay(gap, absgc, fpsmooth);

	//!! plt.setpar(rval, vastr("fpfac",gcno));

	return rval;
} // freeplay


//------------------------------------------------------------------
// 2) time-domain description method
//------------------------------------------------------------------
string
df_iid(const string& dfid, const vector<double>& params) {
// given a dfid which identifies a block of user-written code describing
// a describing function, and a set of parameter/value pair, create a
// string which can be used to identify a df interpolation
	T_(Trace trc(1,"df_iid");)
	string rval{dfid};
	for (auto& pi : params)
		rval += vastr('@',pi);
	T_(trc.dprint("returning ",rval);)
	return rval;
}

bool
userdf (const string& dfid, Qfcn& qs_fcn, Ffcn& fq_fcn, Interpfcn& interp_fcn) {
// placeholder to avoid undefined references while linking programs with
// this library. When flaps run the first thing done is to create userdf.so
// containing the user's code. This is done in src/bin/flaps/load.cpp
	return false;
}

Ad
dfval(pset& plt, const string& dfid, int gcno, const vector<double>& params) {
// Evaluate a user DF at the abs value of gc 'gcno' and parameters
// 'params' and 'plt'. A list of parameter values may be specified with
// an initializer-list, like this:
//   Ad fac = dfval(plt, dfid, gcno, {0.001, 2.0});
	T_(Trace trc(2,"dfval");)
	Ad absgc = plt.absgcv(gcno);
	// get the interpolant id
	string iid = df_iid(dfid, params);
	T_(trc.dprint("dfid ",dfid,", iid ",iid,", gc ",gcno,", absgc ",absgc);)
	// DF interpolations that have been seen are kept in a static vector
	// so care must be taken when modifying the vector to make it thread-safe
	static vector<Interp*> interps;
	Interp* theinterp{nullptr};

	// first check to see if they are in the vector of Interps...
	for (auto ip : interps) {
		if (ip->iid() == iid) {
			theinterp = ip;
			break;
		}
	}
	// ... then try fetching Interp "iid", add it to interps,
	//     guarding the static vector "interps"
	if (theinterp == nullptr) {
		theinterp = Interp::fetch(iid);
		if (theinterp != nullptr) {
			std::mutex mtx;
			std::lock_guard<std::mutex> guard(mtx);
			interps.push_back(theinterp);
		}
	}
	// ... finally, it doesn't exist so build it
	// guard the list of Interps while adding to it, and store the coefficients
	if (theinterp == nullptr) {
		T_(trc.dprint("describing function '",iid,"' has not been built yet");)
		// get pointers to the 3 user functions
		//   1) vector<double> qs_fcn(vector<double>&)
		//   2) double fq_fcn(double y,vector<double> params)
		//   3) Interp* interp_fcn(vector<double> qs, vector<double> fs)
		Qfcn qs_fcn;
		Ffcn fq_fcn;
		Interpfcn interp_fcn;
		if (!userdf(dfid, qs_fcn, fq_fcn, interp_fcn))
			throw runtime_error(vastr("DF ",dfid," has not been loaded"));
		theinterp = df_fft(dfid, params, qs_fcn, fq_fcn, interp_fcn);
		{
			std::mutex mtx;
			std::lock_guard<std::mutex> guard(mtx);
			interps.push_back(theinterp);
			theinterp->store();
		}
	}
	// check that absgc is in range: if not use limits
	double lo, hi;
	theinterp->limits(lo,hi);
	if (absgc.value() < lo)
		absgc.value(lo);
	if (absgc.value() > hi)
		absgc.value(hi);
	// evaluate
	Ad rval = theinterp->eval(absgc);
	T_(trc.dprint("returning ",rval);)
	return rval;
}

// myfcn is for creating function objects which are passed to the
// integrator which then calls it with a single argument (wt)
// XXX replace with a lambda
class myfcn {
	double q_;	// abs value of the gc: amplitude of the time-domain oscillation
	vector<double> params_;	// parameter values used in df_fq
	double (*fcn_)(double y, const vector<double>& params);
public:
	myfcn(double q, const vector<double>& params,
			double (*fq)(double y, const vector<double>& params)) :
				q_(q), params_(params), fcn_(fq) {}
	// operator() is called from the integrator with the single arg
	// double wt [0 : 2pi]
	double operator()(double wt) {
		double sinwt = sin(wt);
		double y = q_*sinwt;
		return fcn_(y,params_)*sinwt;
	}
};

#ifndef HAVE_LIBFFTW3
void
four1(double* data, size_t nn, int isign) {
// Replaces 2*nn vector data by its discrete Fourier transform, if isign is
// input as 1, or replaces data by nn times its inverse discrete
// Fourier transform if isign is -1. data is a complex array of length
// nn or, equivalently, a real array of length 2*nn. nn MUST be
// an integer power of 2
	size_t n = nn<<1;
	size_t j{1};

	for (size_t i=1; i<n; i+=2) {
		if (j > i) {
			swap(data[j], data[i]);
			swap(data[j+1], data[i+1]);
		}
		size_t m = n >> 1;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	// Danielson-Lanczos section
	size_t mmax{2};
	while (n > mmax) {
		size_t istep = mmax << 1;
		double theta = isign*(2.0*flaps::pi/mmax);
		double wtemp = sin(0.5*theta);
		double wpr = -2.0*wtemp*wtemp;
		double wpi = sin(theta);
		double wr = 1.0;
		double wi = 0.0;
		for (size_t m=1; m < mmax; m += 2) {
			for (size_t i=m; i<=n; i += istep) {
				j = i + mmax;
				double tempr = wr*data[j] - wi*data[j+1];
				double tempi = wr*data[j+i] + wi*data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp = wr)*wpr - wi*wpi + wr;
			wi = wi*wpr + wtemp*wpi + wi;
		}
		mmax = istep;
	}
}
#endif // !HAVE_LIBFFTW3

Interp*
df_fft(const string& dfid, const vector<double>& params,
		Qfcn qs_fcn, Ffcn fq_fcn, Interpfcn interp_fcn) {
// Given a DF id, and pointers to the 3 user-written functions which
// define a DF, and create a set of interpolation coefficients that are
// a function of the abs value of a gc. The DF is also a function of
// fixed parameters 'params'.
	T_(Trace trc(1,"df_fft");)

	// 1) get the user's list of displacement magnitudes
	vector<double> qs = qs_fcn(params);

	// 2) evaluate f(q) at each one: 2nd fft coef
	size_t ncoef{64};
	std::vector<double> fs;
	double pi{flaps::pi};
	vector<double> in(ncoef, 0.0);
	vector<double> df;
	for (auto qi : qs) {
		T_(trc.dprint("fft df_fq(",qi,"*sin(wt)");)
		// compute f(q) from 0 -> 2pi
		double dt = 2.0*pi/ncoef;
		for (size_t i=0; i<ncoef; i++) {
			double wt = i*dt;
			double sinwt = sin(wt);
			double y = qi*sinwt;
			in[i] = fq_fcn(y,params);
		}
		// fft the in
#ifdef HAVE_LIBFFTW3
		int nc = ncoef/2 + 1;
		// use std::vector instead of fftw_malloc, cast to fftw_complex*
		vector<complex<double>> out(nc);
		fftw_plan p = fftw_plan_dft_r2c_1d(ncoef, &in[0],
				reinterpret_cast<fftw_complex*>(out.data()), FFTW_ESTIMATE);
		fftw_execute(p);
		double t = 2.0*abs(out[1].imag())/(qi*ncoef);
#else
		size_t nn = ncoef/2;
		four1(&in[0], nn, 1);
		double t = 2.0*abs(in[3])/(qi*ncoef);
#endif

		T_(trc.dprint("df(",qi,") = ",t);)
		df.push_back(t);
	}
	T_(trc.dprintvv(qs, df, "q   df");)

	// 3) interpolate df as a function of qs
	string iid{df_iid(dfid, params)};
	T_(trc.dprint("built df \'",dfid,"\'");)
	return interp_fcn(iid, params, qs, df);
}

void
trapzd(std::function<double(double)> fcn, double a, double b, int n, double& s) {
// Compute the nth stage of refinement of an extended trapezoidal rule.
	// When n=1 return the crudest estimate
	if (n == 1) {
		s = 0.5*(b-a)*(fcn(a) + fcn(b));
		return;
	}
	int it{1};
	for (int j=1; j<n-1; j++) it <<= 1;
	double tnm = it;
	double del = (b-a)/tnm; // spacing of points to be added
	double x = a + 0.5*del;
	double sum{0.0};
	for (int j=1; j <= it; j++) {
		sum += fcn(x);
		x += del;
	}
	 s = 0.5*(s + (b-a)*sum/tnm);
	 return;
}

double
trapezoidal(std::function<double(double)> fcn, double a, double b) {
// Compute the integral of the function fcn from a to b using the trapezoid rule.
// From "Numerical Recipes in C" pg 137
	double olds{std::numeric_limits<double>::max()};
	int maxstep{20};
	double s{0.0};
	double eps{1024*std::numeric_limits<double>::epsilon()};
	for (int j=1; j<=maxstep; j++) {
		trapzd(fcn, a, b, j, s);
		if (abs(s-olds) < eps*abs(olds)) return s;
		if (s == 0.0 && olds == 0.0 && j > 6) return s;
		olds = s;
	}
	throw runtime_error("trapezoidal did not converge");
	return 0.0;
}

Interp*
df_build(const string& dfid, const vector<double>& params,
		Qfcn qs_fcn, Ffcn fq_fcn, Interpfcn interp_fcn) {
// Given a DF id, get pointers to the 3 user-written functions which
// define a DF, and create a set of interpolation coefficients that are
// a function of the abs value of a gc. The DF is also a function of
// fixed parameters 'params'.
	T_(Trace trc(1,"df_build");)

	// 1) create a list of displacement magnitudes
	vector<double> qs = qs_fcn(params);

	// 2) evaluate f(q) at each one: integrate 0-2pi
	std::vector<double> fs;
	double a{0.0};
	double pi{flaps::pi};
	double b = pi*2.0;
	for (auto qi : qs) {
		T_(trc.dprint("integrating df_fq(",qi,"*sin(wt)");)
		// integrate f*sin(wt) from 0-2*pi
		// double t = boost::math::quadrature::trapezoidal(myfcn(qi,params,fq_fcn), a, b);
		double t = trapezoidal([&](double wt) {
						double sinwt = sin(wt);
						double y = qi*sinwt;
						return fq_fcn(y,params)*sinwt;}, a, b);
		// scale by 1/pi*qi
		if (qi > 0.0)
			t *= 1.0/(pi*qi);
		T_(trc.dprint("df(",qi,") = ",t);)
		fs.push_back(t);
	}
	T_(trc.dprintvv(qs, fs, "q   f");)

	// 3) interpolate fs as a function of qs
	string iid{df_iid(dfid, params)};
	T_(trc.dprint("built df \'",dfid,"\'");)
	return interp_fcn(iid, params, qs, fs);
}

void
dfplot(const string& dfid) {
	// fetch Interp "dfid"
	Interp* ip = Interp::fetch(dfid);
	if (ip == nullptr) {
		throw runtime_error(vastr("describing function '",dfid,"' is not available"));
	}
	int nstep{1001};
	ip->plot(dfid, nstep);
}


#ifdef MAIN
#define MAIN 1

#include "Ad.h"
#include "exim.h"
#include "nrc.h"
#include "tyme.h"

using namespace std;

namespace gap {
vector<double> df_qs(const vector<double>& params);
double df_fq(double y, const vector<double>& params);
Interp* df_interp(const string& iid, const vector<double>& params,
	const vector<double>& qs, const vector<double>& fs);
}

namespace bilinear {
vector<double> df_qs(const vector<double>& params);
double df_fq(double y, const vector<double>& params);
Interp* df_interp(const string& iid, const vector<double>& params,
	const vector<double>& qs, const vector<double>& fs);
}

vector<double>
gap::
df_qs(const vector<double>& params) {
   vector<double> rval;
   int nstep{1001};  // use many steps to take significant time
	assert(params.size() == 1);
	double gap = params[0];
   // put a bunch of points between gap and 100*gap...
   double qmax{100.0*gap};
   double del{(qmax-gap)/(double)(nstep-1)};
   for (int i=0; i<nstep; i++)
      rval.push_back(gap+i*del);
   // ... and a single point at 3
   // rval.push_back(3.0);
   return rval;
}

double
gap::
df_fq(double y, const vector<double>& params) {
// returns the factor which multiplies K to get force at displacement 'y'
   T_(Trace trc(1,"df_fs");)
	double gap = params[0];
   double rval{0.0};
   if (abs(y) > gap) {
      if (y > 0.0)
         rval = y-gap;
      else
         rval = y+gap;
   }
   T_(trc.dprint("returning ",rval);)
   return rval;
}

Interp*
gap::
df_interp(const string& iid, const vector<double>& params,
	const vector<double>& qs, const vector<double>& fs) {

	double gap = params[0];
   // fit a spline to qs/fs unsmoothed
   spline* sp = new spline(qs, fs);
   // f(q[0:gap]) = 0: fit a piecewise linear
   vector<double> q0{{0.0}, {gap}};
   vector<double> f0{{0.0}, {0.0}};
   plinear* pl = new plinear(q0,f0);
   // smoothed transition from the piecewise linear to the spline
   double width{0.1};
   double x0{gap*(1.0-width)};
   double y0{0.0};
   double yp0{0.0};
   double x1{gap*(1.0+width)};
   double y1 = sp->eval(x1);
   double yp1 = sp->der(x1);
   corner* cp = new corner(x0, y0, yp0, x1, y1, yp1);
   // first interpolant must be the corner
	Interp* coef = new Interp(iid, "q", "f");
   coef->add(cp);
   coef->add(pl);
   coef->add(sp);
   //!! coef->plot(iid,2001);
	return coef;
}


vector<double>
bilinear::
df_qs(const vector<double>& params) {
   vector<double> rval;
   int nstep{101};
   double delta = params[0];
   // put a bunch of points between gap and 100*gap
   // double qmax{1.0};
   double qmax{100.0*delta};
   double del{(qmax-delta)/(double)(nstep-1)};
   for (int i=0; i<nstep; i++)
      rval.push_back(delta+i*del);
   // ... and a single point at 3
   // rval.push_back(3.0);
   return rval;
}

double
bilinear::
df_fq(double y, const vector<double>& params) {
// returns the factor which multiplies K to get force
// at dimensionless time "wt"
   T_(Trace trc(1,"df_fq");)
   double rval = y;
   double delta = params[0];
   double ratio = params[1];
   if (abs(y) > delta) {
      if (y > 0.0)
         rval = delta + ratio*(y-delta);
      else
         rval = -delta + ratio*(y+delta);
   }
   T_(trc.dprint("returning ",rval);)
   return rval;
}

Interp*
bilinear::
df_interp(const string& iid, const vector<double>& params,
      const vector<double>& qs, const vector<double>& fs) {
   // fit a spline to qs/fs unsmoothed
   double delta = params[0];
   spline* sp = new spline(qs, fs);
   // f(q[0:delta]) = 1: fit a piecewise linear
   vector<double> q0{{0.0}, {delta}};
   vector<double> f0{{1.0}, {1.0}};
   plinear* pl = new plinear(q0,f0);
   // smoothed transition from the piecewise linear to the spline
   double width{0.1};
   double x0{delta*(1.0-width)};
   double y0{1.0};
   double yp0{0.0};
   double x1{delta*(1.0+width)};
   double y1 = sp->eval(x1);
   double yp1 = sp->der(x1);
   corner* cp = new corner(x0, y0, yp0, x1, y1, yp1);
   // first interpolant must be the corner
   Interp* coef = new Interp(iid, "q", "f");
   coef->add(cp);
   coef->add(pl);
   coef->add(sp);
   coef->plot(iid,2001);
   return coef;
}

void
plot_force (const string& plotfile, Ffcn fq_fcn,
		vector<double> params, double qmax) {
// plot over one cycle (0 - 2pi):
// - displacement qmax*sin(wt)
// - force coefficent: K(y)y/K_0
// - df approximation: c_1 * sin(wt)
	// double qmax{10};
	int ncoef{128};
	double pi{flaps::pi};
	double dt = 2.0*pi/ncoef;
	vector<double> wts;
	vector<double> fs;
	vector<double> ys;
	double fmax{0.0};
	for (int i=0; i<ncoef; i++) {
		double wt = i*dt;
		wts.push_back(wt);
		double y = qmax*sin(wt);
		ys.push_back(y);
		double f = fq_fcn(y, params)/qmax;
		fs.push_back(f);
		fmax = std::max(abs(f), fmax);
	}
	// now fft the fs
	// use std::vector instead of fftw_malloc, but it requires a
	// C-style cast!
	//!! fftw_complex* out = (fftw_complex*)fftw_malloc(nc*sizeof(fftw_complex));
	int nc = ncoef/2 + 1;
	vector<complex<double>> out(nc);
	fftw_plan p = fftw_plan_dft_r2c_1d(ncoef, &fs[0],
			(fftw_complex*)&out[0], FFTW_ESTIMATE);
	fftw_execute(p);
	//!! double t = 2.0*abs(out[1].imag())/(qmax*ncoef);
	double t = 2.0*abs(out[1].imag())/ncoef;

	// last point
	wts.push_back(2.0*pi);
	ys.push_back(ys[0]);
	fs.push_back(fs[0]);
	// plot displacement vs wt
	// scale the y's
	normalize_inf(ys.size(), &ys[0], 1, true, fmax/2.0);
	//
	string cid{"disp"};
	vector<vector<double>> xs;
	xs.push_back(wts);
	xs.push_back(ys);
	vector<string> names{"wt", "y"};
	bool append{false};
	Apf::plot(plotfile, cid, xs, names, append);
	// the force
	append = true;
	cid = "force";
	xs.clear();
	xs.push_back(wts);
	xs.push_back(fs);
	names[1] = "f";
	Apf::plot(plotfile, cid, xs, names, append);

	// and plot t*sin wt
	fs.clear();
	for (auto wt : wts) {
		fs.push_back(t*sin(wt));
	}
	xs.clear();
	xs.push_back(wts);
	xs.push_back(fs);
	names[1] = "df";
	append = true;
	cid = "approx";
	Apf::plot(plotfile, cid, xs, names, append);
}

void
test_fftgap() {
	string dfid{"fftgap"};
	vector<double> params{1.0};  // delta
	double qmax{2.0};
	// plot the force and its df approx
	plot_force ("gap-force.apf", gap::df_fq, params, qmax);
	Interp* interp;
	{
		Tyme tyme("df_fft");
		interp = df_fft(dfid, params, gap::df_qs, gap::df_fq, gap::df_interp);
	}
	string plotfile = dfid + ".apf";
	int nstep{101};
	interp->plot(plotfile, nstep);
	// repeat with df_build
	{
		// Tyme tyme("df_build");
		// df_build(dfid, params, gap::df_qs, gap::df_fq, gap::df_interp);
	}
}

void
test_fftbilinear() {
	string dfid{"fftbilinear"};
	vector<double> params{1.0, 8.0};  // delta, ratio
	double qmax{1.5};  // assuming delta==1
	// plot the force and its df approx
	plot_force ("bilinear-force.apf", bilinear::df_fq, params, qmax);

	{
		Tyme tyme("df_fft");
		df_fft(dfid, params, bilinear::df_qs, bilinear::df_fq, bilinear::df_interp);
	}
}

void
test_freeplay(double delta) {
// plot value from the builtin function "freeplay"
	double qmax = 10.0*delta;
	int nstep = 501;
	// smoothing?
	double fpsmooth{0.5};
	// create some new parameters for testing the derivatives
	// of freeplay wrt gap, absgc
	Par gap("gap",vastr(delta));
	Par* gap_par = gpset::get().add(&gap);
	Par absgc("absgc","0");
	Par* absgc_par = gpset::get().add(&absgc);
	vector<string> names{"gap", "absgc"};
	Ad::initialize(names);
	gpset::get().realloc();
	// create a set of amplitudes to plot; if there is no smoothing,
	// start at the gap since we know that below that the value is zero
	double start{delta};
	if (fpsmooth > 0.0)
		start = 0.0;
	double del = (qmax - start)/(double)(nstep-1);
	vector<double> qs;
	for (int i=0; i<nstep; i++)
		qs.push_back(start + i*del);

	Ad adgap = gap_par->advalue();
	Ad absgc_ad;

	vector<double> fgap;
	vector<double> df_dgap_fd;
	vector<double> df_dgap_ad;
	vector<double> df_dq_fd;
	vector<double> df_dq_ad;
	for (auto q : qs) {
		absgc_par->value(q);
		absgc_ad = absgc_par->advalue();
		Ad t = df::freeplay(adgap, absgc_ad, fpsmooth);
		fgap.push_back(t.value());
		// compute df/dgap by finite diff to compare with AD...
		double der{0.0};
		//!! double err{0.0};
		//!! NRC::fdapprox(1,delta,0.0,1.0,&der,gap_deriv,err);
		df_dgap_fd.push_back(der);
		df_dgap_ad.push_back(t.der(0));
		// ... and also df/dq
		//!! NRC::fdapprox(1,q,0.0,2.0*qmax,&der,q_deriv,err);
		df_dq_fd.push_back(der);
		df_dq_ad.push_back(t.der(1));
	}

	string filename{vastr("freeplay",delta,".apf")};
	string cid{"freeplay"};
	vector<vector<double>> xs;
	xs.push_back(qs);
	xs.push_back(fgap);
	//!! xs.push_back(df_dgap_fd);
	xs.push_back(df_dgap_ad);
	//!! xs.push_back(df_dq_fd);
	xs.push_back(df_dq_ad);
	//!! vector<string> plot_names{"q", "f", "dfdgap_fd", "dfdgap_ad", "dfdq_fd", "dfdq_ad"};
	vector<string> plot_names{"q", "f", "dfdgap_ad", "dfdq_ad"};
	bool append{false};
	Apf::plot(filename, cid, xs, plot_names, append);
}

void
test_bl(double ratio) {
// plot value from the builtin function "blfcn"
	// gamma = absgc/delta [0:gmax]
	double gamma{0.0};
	double gmax = 10.0;
	int nstep = 501;
	// smoothing?
	double blsmooth{0.5};
	// create some new parameters for testing the derivatives
	// of bilinear wrt gamma, ratio
	Par gamma_p("gamma",vastr(gamma));
	Par* gamma_par = gpset::get().add(&gamma_p);
	Par ratio_p("ratio",vastr(ratio));
	Par* ratio_par = gpset::get().add(&ratio_p);
	vector<string> names{"gamma", "ratio"};
	Ad::set_adpar_names(names);
	gpset::get().realloc();
	// create a set of gammas to plot: gamma = absgc/delta
	double start{0.0};
	double del = (gmax - start)/(double)(nstep-1);
	vector<double> gs;
	for (int i=0; i<nstep; i++)
		gs.push_back(start + i*del);

	Ad gamma_ad = gamma_par->advalue();
	Ad ratio_ad;

	vector<double> f;
	vector<double> df_dgap_fd;
	vector<double> df_dgamma_ad;
	vector<double> df_dq_fd;
	vector<double> df_dratio_ad;
	for (auto g : gs) {
		gamma_par->value(g);
		gamma_ad = gamma_par->advalue();
		Ad t = blfcn(gamma_ad, ratio_ad, blsmooth);
		f.push_back(t.value());
		// compute df/dgap by finite diff to compare with AD...
		double der{0.0};
		//!! double err{0.0};
		//!! NRC::fdapprox(1,delta,0.0,1.0,&der,gap_deriv,err);
		df_dgap_fd.push_back(der);
		df_dgamma_ad.push_back(t.der(0));
		// ... and also df/dq
		//!! NRC::fdapprox(1,q,0.0,2.0*qmax,&der,q_deriv,err);
		df_dq_fd.push_back(der);
		df_dratio_ad.push_back(t.der(1));
	}

	string filename{vastr("bilinear",ratio,".apf")};
	string cid{"bilinear"};
	vector<vector<double>> xs;
	xs.push_back(gs);
	xs.push_back(f);
	//!! xs.push_back(df_dgap_fd);
	xs.push_back(df_dgamma_ad);
	//!! xs.push_back(df_dq_fd);
	xs.push_back(df_dratio_ad);
	//!! vector<string> plot_names{"q", "f", "dfdgap_fd", "dfdgap_ad", "dfdq_fd", "dfdq_ad"};
	vector<string> plot_names{"gamma", "f", "dfdgamma_ad", "dfdratio_ad"};
	bool append{false};
	Apf::plot(filename, cid, xs, plot_names, append);
}

int
main (int argc, char** argv) {
	T_(Trace trc(1,"df");)
	double delta{0.001};
	double ratio{2.0};
	bool bl{false};
	bool fp{false};
	bool df{false};
	string usage{"df [-b ratio][-f gap][-d]"};

	if (argc > 1) {
		string arg(argv[1]);
		if (arg == "-b") {
			bl = true;
			str2double(arg, ratio);
		} else if (arg == "-f") {
			fp = true;
			str2double(arg, delta);
		} else if (arg == "-d") {
			df = true;
		} else {
			cerr << usage << endl;
			return 1;
		}
	} else {
		cerr << usage << endl;
		return 1;
	}

	// test the fft gap & bilinear dfs
	if (df) {
		try {
			test_fftgap();
		} catch (runtime_error& s) {
			cerr << "test_fftgap failed: " << s << endl;
		}
		try {
			test_fftbilinear();
		} catch (runtime_error& s) {
			cerr << "test_fftbilinear failed: " << s << endl;
		}
	}
	// test the builtin freeplay df
	if (fp) {
		try {
			test_freeplay(delta);
		} catch (runtime_error& s) {
			cerr << "test_freeplay failed: " << s << endl;
		}
	}

	// 1 argument (ratio): test the builtin bilinear df
	if (bl) {
		try {
			test_bl(ratio);
		} catch (runtime_error& s) {
			cerr << "test_bl failed: " << s << endl;
		}
	}
}
#endif // MAIN
