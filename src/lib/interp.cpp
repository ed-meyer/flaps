//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


#include <iomanip>
#include <iostream>
#include <limits>

#include "config.h"
#include "conv.h"
#include "blas.h"
#include "exim.h"
#include "fio.h"
#include "interp.h"
#include "lapack.h"
#include "trace.h"

using namespace std;

int
cubgcv (double const* x, double* f, int n, double rho, double* c);

// class interpolant implementation

void
interpolant::
put(Sender& s) const {}
//!! moved to interp.h interpolant::~interpolant() {}

// class corner implementation
// register class for Fio::get/put serialization
bool corner::regd = Fio::register_class(corner::id(), corner::get);

corner::
corner (double xa, double ya, double ypa, double xb, double yb, double ypb) {
	T_(Trace trc(1,"corner constructor");)
	int n{4};

	T_(trc.dprint("between (x,y,yp) = (",xa,',',ya,',', ypa,") and (",xb,',',yb,',',ypb,')');)

	x0 = xa;
	x1 = xb;
	c = vector<double>(n, 0.0);

	// It is ok if the width is zero (x0 == x1) => sharp corner
	// but if y0 != y1 we have a discontinuity!
	if (x0 == x1) {
		T_(trc.dprint("returning sharp corner");)
		return;
	}

	vector<double> a(n*n, 0.0);
	vector<double> af(n*n, 0.0);
	vector<double> b(n, 0.0);

	double t1 = x1 - x0;

	// compute the interp coeff array in terms of t = x - x0
	// so the interpolant is
	//   y(t) = c0 + t*c1 + t^2*c2 + t^3*c3
	//   y'(t) = c1 + t*c2 + t^2*c3
	//   Ac = b
	//   b = {y(0) y'(0) y(t1) y'(t1)}^T
	//   c = {c0 c1 c2 c3}^T
	// A row 1: y(0) = c0 = ya
	a[0] = 1.0;
	b[0] = ya;

	// row 2: y'(0) = c1 = ypa
	a[1+n] = 1.0;
	b[1] = ypa;

	// row 3: y(t1) = c0 + t1*c1 + t1^2*c2 + t1^3*c3 = yb
	a[2] = 1.0;
	a[2+n] = t1;
	a[2+(2*n)] = t1*t1;
	a[2+(3*n)] = t1*t1*t1;
	b[2] = yb;

	// row 4: y'(t1) = c1 + t1*c2 + t1^2*c3 = ypb
	a[3+n] = 1.0;
	a[3+(2*n)] = 2.0*t1;
	a[3+(3*n)] = 3.0*t1*t1;
	b[3] = ypb;

	T_(trc.dprintm(4,4,4,a,"A = ");)

	vector<int> ipiv(n, 0);
	vector<double> rowsf(n);
	vector<double> colsf(n);
	double rcond, ferr, berr;
	int info;
	char equed[2];

	vector<double> work(4*n);
	vector<int> iwork(n);
	info = lapack::dgesvx("n", "n", n, 1, &a[0], n, &af[0], n, &ipiv[0],
			equed, &rowsf[0], &colsf[0], &b[0], n, &c[0], n, &rcond, &ferr, &berr);
	if (info != 0) {
		ostringstream os;
		os << "cannot interpolate between (" <<
			xa << ',' << ya << ',' << ypa << ") and (" <<
			xb << ',' << yb << ',' << ypb << ") sgesvx returned rcond " << rcond
			<< ", info " << info << endl;
		os << strArray(4,4,4,&a[0]) << endl;
		T_(trc.dprint("throwing exception: ",os.str());)
		throw runtime_error(os.str());
	}
	T_(trc.dprint("rcond = ",rcond," c = ",c[0],", ",c[1],", ",c[2],", ",c[3]);)
}

corner::
corner(Receiver& s) {
	s.serialize(x0);
	s.serialize(x1);
	s.serialize(c);
}

void
corner::
put(Sender& s) const {
	s.serialize(x0);
	s.serialize(x1);
	s.serialize(c);
}

bool
corner::
contains (double t) const {
	if (x0 <= t && t <= x1)
		return true;
	return false;
}

void
corner::
limits(double& min, double& max) const {
	min = x0;
	max = x1;
}

double
corner::
eval(double x) const {
	T_(Trace trc(1,"corner::eval ",x);)
	if (!this->contains(x))
		throw runtime_error("corner::eval called with x out-of-range");

	double t = x - x0;
	double rval = ((c[3]*t + c[2])*t + c[1])*t + c[0];
	T_(trc.dprint("returning ",rval);)
	return rval;
}

Ad
corner::
eval(Ad const& x) const {

	if (!this->contains(x.value()))
		throw runtime_error("corner::eval called with x out-of-range");

	Ad t = x - x0;
	return ((c[3]*t + c[2])*t + c[1])*t + c[0];
}

double
corner::
deriv(double x) const {
	double t = x - x0;
	return (c[3]*t + c[2])*t + c[1];
}


// plinear implementation
// register class for Fio::get/put serialization
bool plinear::regd = Fio::register_class(plinear::id(), plinear::get);


plinear::
plinear(Receiver& s) {
	s.serialize(x);
	s.serialize(y);
}

void
plinear::
put(Sender& s) const {
	s.serialize(x);
	s.serialize(y);
}

bool
plinear::
contains (double t) const {
	for (size_t i=0; i<x.size()-1; i++) {
		if (x[i] <= t && t <= x[i+1])
			return true;
	}
	return false;
}

void
plinear::
limits(double& min, double& max) const {
	min = x[0];
	max = x[x.size()-1];
}

double
plinear::
eval(double xa) const {
	T_(Trace trc(1,"plinear::eval ",xa);)
	size_t nm1{x.size()-1};
	// check for extrapolation
	if (xa < x[0]) {
		double ya = y[0] + (xa-x[0])*(y[1]-y[0])/(x[1]-x[0]);
		T_(trc.dprint("returning extrapolation y(",xa,") = ",ya);)
		return ya;
	}
	if (xa > x[nm1]) {
		double ya = y[nm1] + (xa-x[nm1])*(y[nm1]-y[nm1-1])/(x[nm1]-x[nm1-1]);
		T_(trc.dprint("returning extrapolation y(",xa,") = ",ya);)
		return ya;
	}
	// check each segment
	for (size_t i=0; i<nm1; i++) {
		if (x[i] <= xa && xa <= x[i+1]) {
			double t = xa - x[i];
			double yap = (y[i+1] - y[i])/(x[i+1] - x[i]);
			double ya = y[i] + t*yap;
			T_(trc.dprint("returning ",ya);)
			return ya;
		}
	}
	throw runtime_error("plinear::eval arg out of range");
}

Ad
plinear::
eval(Ad const& xa) const {
	T_(Trace trc(2,"plinear::eval ",xa);)
	size_t nm1{x.size()-1};
	double xv{xa.value()};
	// check for extrapolation
	if (xv < x[0]) {
		Ad ya = y[0] + (xa-x[0])*(y[1]-y[0])/(x[1]-x[0]);
		T_(trc.dprint("returning extrapolation y(",xa,") = ",ya);)
		return ya;
	}
	if (xv > x[nm1]) {
		Ad ya = y[nm1] + (xa-x[nm1])*(y[nm1]-y[nm1-1])/(x[nm1]-x[nm1-1]);
		T_(trc.dprint("returning extrapolation y(",xa,") = ",ya);)
		return ya;
	}
	// check each segment
	for (size_t i=0; i<nm1; i++) {
		if (x[i] <= xv && xv <= x[i+1]) {
			Ad t = xa - x[i];
			Ad yap{(y[i+1] - y[i])/(x[i+1] - x[i])};
			Ad ya = y[i] + t*yap;
			return ya;
		}
	}
	throw runtime_error(vastr("plinear::eval arg (",xa,") out of range"));
}

double
plinear::
deriv(double xa) const {
	T_(Trace trc(1,"plinear::deriv ",xa);)
	size_t nm1{x.size()-1};
	for (size_t i=0; i<nm1; i++) {
		// should never be called with xa a breakpoint
		if (x[i] <= xa && xa <= x[i+1]) {
			if (xa == x[i] || xa == x[i+1])
				throw runtime_error(vastr("plinear::deriv(",xa,") called with breakpoint"));
			double yap = (y[i+1] - y[i])/(x[i+1] - x[i]);
			T_(trc.dprint("returning ",yap);)
			return yap;
		}
	}
	throw runtime_error("plinear::eval arg out of range");
}

ostream&
operator<<(ostream& s, const plinear& t) {
	s << "piecwise-linear interpolation with " << t.x.size() << " breakpoints:\n";
	for (size_t i=0; i<t.x.size(); i++)
		s << "   " << t.x[i] << "   " << t.y[i] << endl;
	return s;
}


// class spline implementation
// register class for Fio::get/put serialization
bool spline::regd = Fio::register_class(spline::id(), spline::get);

#define IJ(i,j,ld1) ((i) + (ld1)*(j))

spline::
spline (const int nsa, const double *x, const double *y, double rho) {
// scalar cubic spline with optional smoothing.
// rho  controls the smoothing:
//      0:   no smoothing
//      >0   smoothing; values between 1-10 yield the best smoothing.
//           values >10 approach a least-squares linear fit
	T_(Trace trc(1,"spline(rho) constructor");)
	int i;
	ostringstream os;

	T_(trc.dprint("rho ", rho);)

	int nint = nsa - 1;

	vector<double> c(4*nint, 0.0);
	// save x in xbp, y in f and ybp
	xbp = vector<double>(nsa, 0.0);
	vector<double> f(nsa);
	ybp = vector<double>(nsa, 0.0);

	for (i=0; i<nsa; i++) {
		xbp[i] = x[i];
		f[i] = y[i];
		ybp[i] = y[i];
	}

	int ier = cubgcv (&x[0], &f[0], nsa, rho, &c[0]);

	if (ier == 131) {
		os << "attempt to interpolate the following data where "
			"the x-values are not strictly increasing:\n";
		for (i=0; i<nsa; i++)
			os << x[i] << "   " << y[i] << endl;
		T_(trc.dprint("throwing exception \"",os.str(),"\"");)
		throw runtime_error(os.str());
	}

	if (ier != 0) {
		throw runtime_error(vastr("interpolation failed: ier = ", ier));
	}

	T_(trc.dprintm(nint,4,nint,c,"coef");)

	// Note: c0 may not be exactly the same as ybp if rho > 0
	c0 = vector<double>(nsa);
	c1 = vector<double>(nint);
	c2 = vector<double>(nint);
	c3 = vector<double>(nint);

	for (i=0; i<nint; i++) {
		c0[i] = c[IJ(i,3,nint)];
		c1[i] = c[IJ(i,2,nint)];
		c2[i] = c[IJ(i,1,nint)];
		c3[i] = c[IJ(i,0,nint)];
	}
	c0[nint] = y[nint];
}

spline::
spline(const vector<double>& x, const vector<double>& y, double rho) {
// spline(vector) constructor: just call pointer version
	T_(Trace trc(1,"spline(vector,vector) constructor");)
	T_(trc.dprintvv(x,y,"input x,y");)
	if (x.size() != y.size())
		throw runtime_error(vastr("spline constructor called with ",
				x.size()," x, but ", y.size(), "y"));
	*this = spline(x.size(), &x[0], &y[0], rho);
}

spline::
spline(Receiver& s) {
	s.serialize(xbp);
	s.serialize(ybp);
	s.serialize(c0);
	s.serialize(c1);
	s.serialize(c2);
	s.serialize(c3);
}

void
spline::
put(Sender& s) const {
	s.serialize(xbp);
	s.serialize(ybp);
	s.serialize(c0);
	s.serialize(c1);
	s.serialize(c2);
	s.serialize(c3);
}

bool
spline::
contains (double t) const {
	for (size_t i=0; i<xbp.size()-1; i++) {
		if (xbp[i] <= t && t <= xbp[i+1])
			return true;
	}
	return false;
}

void
spline::
limits(double& min, double& max) const {
	min = xbp[0];
	max = xbp[xbp.size()-1];
}

double
spline::
eval (double t) const {
	T_(Trace trc(1,"spline::eval ",t);)
	double d{0.0};
	double rval{0.0};

	int nint = xbp.size() - 1;

	int i = interval(t);
	if (i < 0 || i >= nint) {
		string exc{vastr("interval ", i, " is out of the range 0-", nint-1,
			", t = ", t, ", x[", nint, "] = ", xbp[nint])};
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}

	d = t - xbp[i];
	rval = ((c3[i]*d + c2[i])*d + c1[i])*d + c0[i];
	T_(trc.dprint("returning ",rval);)
	return rval;
}

Ad
spline::
eval (Ad const& t) const {
	T_(Trace trc(3,"eval<Ad>");)
	Ad rval;

	T_(trc.dprint("t = ",t);)

	int nint = xbp.size() - 1;
	int i = interval(t.value());
	if (i < 0 || i >= nint) {
		string exc = vastr("interval ", i, " is out of the range 0-", nint-1,
					", t = ", t, ", x[", nint, "] = ", xbp[nint]);
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}

	Ad d = t - xbp[i];
	rval = ((c3[i]*d + c2[i])*d + c1[i])*d + c0[i];
	T_(trc.dprint("returning ",rval);)
	return rval;
}

double
spline::
deriv (double t) const {
	T_(Trace trc(1,"spline::deriv ",t);)
	// If t is out of range extrapolate using the fact that
	// since this is a natural spline the second deriv at the endpoint is zero...
	// Take advantage of the fact that we saved all nint+1 values
	// of the original y in c0.
	int nint = xbp.size() - 1;
	if (t < xbp[0])
		t = xbp[0];
	if (t > xbp[nint])
		t = xbp[nint];

	int i = interval(t);
	double d = t - xbp[i];
	double rval = (3.0*c3[i]*d + 2.0*c2[i])*d + c1[i];
	T_(trc.dprint("returning ",rval);)
	return rval;
}

int
spline::
interval(double t) const {
// Returns the (zero-base) interval that "t" is in
// If t is below the range, returns 0.
// If t is above the range, returns nint
	int nint = xbp.size() - 1;
	for (int i=0; i<nint; i++) {
		if (t <= xbp[i+1]) {
			return i;
		}
	}
	return nint;
}


// vspline implementation
// cubic interpolation of vectors with optional smoothing

// register class for Fio::get/put serialization
bool vspline::regd = Fio::register_class(vspline::id(), vspline::get);

vspline::
vspline (vector<double> const& x, vector<vector<double> > const& y,
		int ext, double rh) {
// x    (nx) vector of breakpoints
// y    nx vectors of length ny: array of data
// Coefficients are such that
//   f(t) = c0[i] + d*(c1[i] + d*(c2[i] + d*c3[i]))
// where d = t-bp[i] and bp[i] <= t <= bp[i+1]
// Note the coefficients are stored reversed from cubgcv
	T_(Trace trc(1,"vspline constructor");)
	size_t i, j;

	size_t nx = x.size();
	// the number of intervals is the number of breakpoints - 1
	size_t nint = nx - 1;
	size_t ny = y[0].size();
	assert(y.size() == nx);
	extrap = ext;
	rho = rh;

	T_(trc.dprint("nx ",nx,", ny ",ny,", rho ",rho,", extrap ",ext);)

	// (nint,4) arrary of coeff
	vector<double> c(nint*4, 0.0);
	vector<double> yj(nx);

	// copy orig data values to yorig
	yorig = y;
	// ... and the breakpoints to bp
	bp = x;
	// initialize c1, c2, c3 (nx-1,ny)
	vector<double> cj(ny, 0.0);
	for (i=0; i<nint; i++) {
		c0.push_back(cj);
		c1.push_back(cj);
		c2.push_back(cj);
		c3.push_back(cj);
	}
	// if there are only 2 data points (nint == 1) do a linear interp
	if (nint == 1) {
		c0[0] = y[0];
		//!! c0[1] = y[1];
		double dt = bp[1] - bp[0];
		for (j=0; j<ny; j++)
			c1[0][j] = (y[1][j] - y[0][j])/dt;
		T_(trc.dprintm(nint,3,nint,c,"coefficients");)
		return;
	}

	// interpolate each element of the y vectors separately
	for (j=0; j<ny; j++) {
		for (i=0; i<nx; i++) {
			yj[i] = y[i][j];
		}
		T_(trc.dprintvv(x, yj, "x     y",j);)

		int ier = cubgcv (&x[0], &yj[0], nx, rho, &c[0]);
		T_(trc.dprint("cubgcv returned ",ier);)

		if (ier != 0) {
			string exc = vastr("vspline interpolation failed: ier = ", ier);
			T_(trc.dprint("throwing exception: ",exc);)
			throw runtime_error(exc);
		}
		T_(trc.dprintm(nint,3,nint,c,"coefficients");)

		// copy the coeff reversed from cubgcv
		for (i=0; i<nint; i++) {
			c0[i][j] = c[IJ(i,3,nint)];
			c1[i][j] = c[IJ(i,2,nint)];
			c2[i][j] = c[IJ(i,1,nint)];
			c3[i][j] = c[IJ(i,0,nint)];
		}
	}
}

vspline::
vspline (Receiver& s) {
	T_(Trace trc(2,"vspline Receiver constructor");)
	s.serialize(bp);
	s.serialize(yorig);
	s.serialize(c0);
	s.serialize(c1);
	s.serialize(c2);
	s.serialize(c3);
	s.serialize(extrap);
	s.serialize(rho);
}

void
vspline::
put(Sender& s) const {
	T_(Trace trc(2,"vspline::put");)
	s.serialize(bp);
	s.serialize(yorig);
	s.serialize(c0);
	s.serialize(c1);
	s.serialize(c2);
	s.serialize(c3);
	s.serialize(extrap);
	s.serialize(rho);
}

bool
vspline::
contains (double t) const {
	for (size_t i=0; i<bp.size()-1; i++) {
		if (bp[i] <= t && t <= bp[i+1])
			return true;
	}
	return false;
}

void
vspline::
limits(double& min, double& max) const {
	min = bp[0];
	max = bp[bp.size()-1];
}

// vspline has a different eval from the other interpolants: CAdvector arg
bool
vspline::
eval (Ad const& t, vector<complex<Ad>>& result) const {
	T_(Trace trc(2,"vspline::eval", " at t = ",t);)
	bool rval{true};
	size_t nx = bp.size();
	size_t nint = nx - 1;  // # of intervals

	// check sizes: orig matrices must have been cast to complex if they
	// are real, so the orig data (real) is twice the (complex) result
	// XXX allow for real data (no casting)?
	if (yorig[0].size() != 2*result.size())
		throw runtime_error(vastr("size mismatch in vspline::eval: ",result.size(),
			" result, ",yorig[0].size()," in coeff"));

	// if t is > max+eps (bp[nint]) return the last breakpoint
	// if it is closer to min or max extrapolate (with Ad derivatives)
	// so that the Jacobian has nearly correct derivatives and the
	// tangent angle won't be too large, causing a retry of the step
	double tv = t.value();
	double meps = numeric_limits<double>::epsilon();
	double eps = sqrt(meps)*(1.0 + tv);
	//!! if (is_equal(tv, bp[nint], 8) || tv > bp[nint])
	if (tv > bp[nint]+eps) {
		for (size_t i=0; i<result.size(); i++)
			result[i] = complex<Ad>(yorig[nint][2*i], yorig[nint][2*i+1]);
		T_(trc.dprint("t is at its max: returning yorig[nint]");)
		return rval;
	}
	// if t is < min-eps (bp[0]) return the first breakpoint
	//!! if (is_equal(tv, bp[0], 8) || tv < bp[0])
	if (tv < bp[0]-eps) {
		for (size_t i=0; i<result.size(); i++)
			result[i] = complex<Ad>(yorig[0][2*i], yorig[0][2*i+1]);
		T_(trc.dprint("t is at its min: returning yorig[0]");)
		return rval;
	}

	int iseg = segment(t);

	// if t is out of range check if extrapolation is ok
	if (iseg < 0 || iseg >= (int)nint) {
		if (extrap) {
			if (iseg < 0)
				iseg = 0;
			if (iseg >= (int)nint)
				iseg = nint-1;
			T_(trc.dprint("extrapolating: t ",t," bp ",bp[iseg]);)
		} else {
			string exc{vastr("interval ",iseg," is out of the range 0-",
					nint-1,", t = ", t, ", x[", nint, "] = ", bp[nint])};
			throw runtime_error(exc);
		}
	}

	// result = c0 + d*c1 + d^2*c2 + d^3*c3
	Ad d(t);
	double dmbp = d.value() - bp[iseg];
	d.value(dmbp);
	vector<Ad> ds(4);
	ds[0].value(1.0);
	ds[1] = d;
	ds[2] = d*d;
	ds[3] = ds[2]*d;
	T_(trc.dprint("evaluating vspline at t-bp = ",t," - ",bp[iseg]," = ",d);)
	//!! Ad d2 = d*d;
	//!! Ad d3 = d2*d;
	vector<double> cr(4);
	vector<double> ci(4);
	for (size_t i=0; i<result.size(); i++) {
		int k = 2*i;
		cr[0] = c0[iseg][k];
		cr[1] = c1[iseg][k];
		cr[2] = c2[iseg][k];
		cr[3] = c3[iseg][k];
		ci[0] = c0[iseg][k+1];
		ci[1] = c1[iseg][k+1];
		ci[2] = c2[iseg][k+1];
		ci[3] = c3[iseg][k+1];
		blas::dot(4, &cr[0], 1, &ds[0], 1, Ad::real(result[i])); 
		blas::dot(4, &ci[0], 1, &ds[0], 1, Ad::imag(result[i])); 
	}
	return rval;
}

int
vspline::
segment(double t) const {
// Returns the (zero-based) interval that "t" is in
// If t is below the range, returns 0.
// If t is above the range, returns nint
	T_(Trace trc(2,"vspline::segment");)
	int nx = bp.size();
	int nint = nx - 1;

	if (t < bp[0]) {
		T_(trc.dprint("t (",t,") < bp[0] (",bp[0],") returning -1");)
		return -1;
	}

	for (int i=0; i<nint; i++) {
		if (t <= bp[i+1]) {
			T_(trc.dprint("returning ",i," bp ",bp[i+1]);)
			return i;
		}
	}
	T_(trc.dprint("t (",t,") is > bp[",nint,"]: returning ",nint);)
	return nint;
}

int
vspline::
segment(Ad const& t) const {
	return vspline::segment(t.value());
}

ostream&
operator<<(ostream& s, const vspline& t) {
	s << "B-spline: " << t.getnx() << " breakpoints, "
		" vector size = " << t.getny() << endl;
	return s;
}

//------------------------------------------------------------------
// vspline2 implementation
//------------------------------------------------------------------

// register class for Fio::get/put serialization
bool vspline2::regd = Fio::register_class(vspline2::id(), vspline2::get);


vspline2::
vspline2(const vector<double>& x, const vector<double>& y,
		const vector<vector<double> >& data) :
		x_basis(x), y_basis(y) {
	size_t nx = x.size();
	size_t ny = y.size();
	size_t nd = data[0].size();
	size_t nxy = data.size();
	if (nx*ny != nxy) {
		string exc{vastr("data size ",nxy," != nx (",nx,") * ny (",ny,")")};
		throw runtime_error(exc);
	}
	assert(nx*ny == nxy);
	vector<double> z(nxy);
	for (size_t j=0; j<nd; j++) {
		for (size_t i=0; i<nxy; i++) {
			z[i] = data[i][j];
		}
		splines.push_back(NUBspline_2d(z,x_basis,y_basis));
	}
}

vspline2::
vspline2 (Receiver& s) {
	T_(Trace trc(2,"vspline2 Receiver constructor");)
	// grids
	s.serialize(x_basis.grid);
	s.serialize(x_basis.xVals);
	s.serialize(x_basis.dxInv);

	s.serialize(y_basis.grid);
	s.serialize(y_basis.xVals);
	s.serialize(y_basis.dxInv);

	size_t nspline{0};
	s.serialize(nspline);
	splines = vector<NUBspline_2d>(nspline);
	for (auto& spline : splines) {
		s.serialize(spline.coefs);
		s.serialize(spline.x_stride);
	}
}

void
vspline2::
put(Sender& s) const {
	T_(Trace trc(2,"vspline2::put");)
	// grid
	s.serialize(x_basis.grid);
	s.serialize(x_basis.xVals);
	s.serialize(x_basis.dxInv);

	s.serialize(y_basis.grid);
	s.serialize(y_basis.xVals);
	s.serialize(y_basis.dxInv);

	s.serialize(splines.size());
	for (auto& spline : splines) {
		s.serialize(spline.coefs);
		s.serialize(spline.x_stride);
	}
}

bool
vspline2::
contains (double t) const {
// not implemented: is t x or y?
	return false;
}

void
vspline2::
limits(double& min, double& max) const {
// not implemented: need x & y limits
}

// vspline2 has a different eval from the other interpolants: vector<complex<Ad>>& arg
bool
vspline2::
eval (vector<Ad> const& t, vector<complex<Ad>>& result) {
	T_(Trace trc(2,"vspline2::eval", " at t = ",t);)
	bool rval{true};
	//!! size_t nr = result.rsize();
	//!! size_t nc = result.csize();
	// sanity check
	if (splines.size() != 2*result.size())
		throw runtime_error(vastr("size mismatch in vspline2: ",splines.size(),
			" splines, ",result.size()," complex result"));
	Ad x{t[0]};
	Ad y{t[1]};
	Ad a[4];
	Ad b[4];
	Ad cb[4];
	int ix = get_NUBasis_funcs_d(x_basis, x, a);
	int iy = get_NUBasis_funcs_d(y_basis, y, b);
	//!! size_t is{0};
	//!! for (size_t j=1; j<=nc; j++) {    // 1b
		//!! for (size_t i=1; i<=nr; i++) 
			//!! result(i,j) = complex<Ad>(splines[is].eval(a,b,ix,iy),
				//!! 	splines[is+1].eval(a,b,ix,iy));
	for (size_t i=0; i<result.size(); i++)
		result[i] = complex<Ad>(splines[2*i].eval(a,b,cb,ix,iy),
				splines[2*i+1].eval(a,b,cb,ix,iy));
	return rval;
}

ostream&
operator<<(ostream& s, const vspline2& t) {
	s << "2D B-spline with grid: " << t.x_basis.grid.size()
		<< " by " << t.y_basis.grid.size() <<
		" data size = " << t.ndata << endl;
	return s;
}

//------------------------------------------------------------------
// vspline3 implementation
//------------------------------------------------------------------

// register class for Fio::get/put serialization
bool vspline3::regd = Fio::register_class(vspline3::id(), vspline3::get);


vspline3::
vspline3(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& z,
		const vector<vector<double> >& data) :
		x_basis(x), y_basis(y), z_basis(z) {
	size_t nx = x.size();
	size_t ny = y.size();
	size_t nz = z.size();
	size_t nd = data[0].size();
	size_t nxyz = data.size();
	if (nx*ny*nz != nxyz) {
		string exc{vastr("data size ",nxyz," != nx (",nx,") * ny (",
				ny,") * nz (",nz,")")};
		throw runtime_error(exc);
	}
	vector<double> w(nxyz);
	for (size_t j=0; j<nd; j++) {
		for (size_t i=0; i<nxyz; i++) {
			w[i] = data[i][j];
		}
		splines.push_back(NUBspline_3d(w,x_basis,y_basis,z_basis));
	}
}

vspline3::
vspline3 (Receiver& s) {
	T_(Trace trc(2,"vspline3 Receiver constructor");)
	// grids
	s.serialize(x_basis.grid);
	s.serialize(x_basis.xVals);
	s.serialize(x_basis.dxInv);

	s.serialize(y_basis.grid);
	s.serialize(y_basis.xVals);
	s.serialize(y_basis.dxInv);

	s.serialize(z_basis.grid);
	s.serialize(z_basis.xVals);
	s.serialize(z_basis.dxInv);

	size_t nspline{0};
	s.serialize(nspline);
	splines = vector<NUBspline_3d>(nspline);
	for (auto& spline : splines) {
		s.serialize(spline.coefs);
		s.serialize(spline.x_stride);
		s.serialize(spline.y_stride);
	}
}

void
vspline3::
put(Sender& s) const {
	T_(Trace trc(2,"vspline3::put");)
	// grid
	s.serialize(x_basis.grid);
	s.serialize(x_basis.xVals);
	s.serialize(x_basis.dxInv);

	s.serialize(y_basis.grid);
	s.serialize(y_basis.xVals);
	s.serialize(y_basis.dxInv);

	s.serialize(z_basis.grid);
	s.serialize(z_basis.xVals);
	s.serialize(z_basis.dxInv);

	s.serialize(splines.size());
	for (auto& spline : splines) {
		s.serialize(spline.coefs);
		s.serialize(spline.x_stride);
		s.serialize(spline.y_stride);
	}
}

bool
vspline3::
contains (double t) const {
// not implemented: is t x or y?
	return false;
}

void
vspline3::
limits(double& min, double& max) const {
// not implemented: need x & y limits
}

// vspline3 has a different eval from the other interpolants: vector<complex<Ad>>& arg
bool
vspline3::
eval (vector<Ad> const& t, vector<complex<Ad>>& result) {
	T_(Trace trc(2,"vspline3::eval", " at t = ",t);)
	bool rval{true};
	// sanity check
	if (splines.size() != 2*result.size())
		throw runtime_error(vastr("size mismatch in vspline3: ",splines.size(),
			" splines, ",result.size()," complex result"));
	Ad x{t[0]};
	Ad y{t[1]};
	Ad z{t[2]};
	Ad a[4];
	Ad b[4];
	Ad c[4];
	int ix = get_NUBasis_funcs_d(x_basis, x, a);
	int iy = get_NUBasis_funcs_d(y_basis, y, b);
	int iz = get_NUBasis_funcs_d(z_basis, z, c);
	for (size_t i=0; i<result.size(); i++)
		result[i] = complex<Ad>(splines[2*i].eval(a,b,c,ix,iy,iz),
				splines[2*i+1].eval(a,b,c,ix,iy,iz));
	return rval;
}

ostream&
operator<<(ostream& s, const vspline3& t) {
	s << "2D B-spline with grid: " << t.x_basis.grid.size()
		<< " by " << t.y_basis.grid.size() <<
		" data size = " << t.ndata << endl;
	return s;
}

// class Interp implementation
// a collection of different types of interpolation

// register class for Fio::get/put serialization
bool Interp::regd = Fio::register_class(Interp::id(), Interp::get);

Interp::
Interp(Receiver& s) {
// Interp deserialize constructor
	s.serialize(iid_);
	s.serialize(parnames_);
	s.serialize(dataname_);
	// read the vector of interpolants
	vector<Fio*> tmp;
	getVFio(s, tmp);
	for (auto tp : tmp) {
		interpolant* ip = dynamic_cast<interpolant*>(tp);
		if (ip == nullptr) {
			throw runtime_error(vastr("expecting an interpolant, got a <",tp->vid(),">"));
		} else {
			segments.push_back(ip);
		}
	}
}

Interp::
Interp(const std::string& id) {
// Interp fetch constructor
	vector<string> ids = fio::catalog(id);
	if (ids.empty()) {
		string exc = vastr ("Interp \"",id,"\" is not available");
		throw runtime_error(exc);
	}

	Receiver file(id);
	Fio* op = Fio::get(file);
	Interp* rval = dynamic_cast<Interp*>(op);
	// check it is the right type
	if (rval == nullptr) {
		string exc{vastr("fetching ",id,", expecting a Interp, got a ",op->vid())};
		throw runtime_error(exc);
	}
	// need move
	*this = *rval;
}

Interp*
Interp::
fetch(const string& iid) {
// fetch an Interp named "iid", return nullptr if not available
	T_(Trace trc(1,"Interp::fetch ",iid);)
	Interp* rval{nullptr};
	try {
		rval = new Interp(iid);
	} catch (runtime_error& s) {
		T_(trc.dprint("caught exception: ",s.what());)
		return nullptr;
	}
	return rval;
}

void
Interp::
put(Sender& s) const {
	s.serialize(iid_);
	s.serialize(parnames_);
	s.serialize(dataname_);
	vector<Fio*> tmp;
	for (auto ip : segments)
		tmp.push_back(static_cast<interpolant*>(ip));
	putVFio(s, tmp);
}

void
Interp::
store() {
	Sender file(iid_);
	Fio::put(file, *this);
}

void
Interp::
add (vector<interpolant*> ts) {
	for (auto tp : ts) {
		segments.push_back(tp);
	}
}

void
Interp::
limits(double& lo, double& hi) const {
	T_(Trace trc(1,"Interp::limits");)
	double big{std::numeric_limits<double>::max()};
	lo = big;
	hi = -big;
	double mini;
	double maxi;
	for (auto ip : segments) {
		ip->limits(mini, maxi);
		T_(trc.dprint(ip->vid()," segment limits: ",mini,":",maxi);)
		lo = std::min(lo, mini);
		hi = std::max(hi, maxi);
	}
}

double
Interp::
eval(double t) const {
	double rval{0.0};
	// evaluate the first segment containing t
	for (auto sp : segments) {
		if (sp->contains(t))
			return sp->eval(t);
	}
	return rval;
}

Ad
Interp::
eval(Ad const& t) const {
	T_(Trace trc(2,"Interp::eval ",segments.size()," segments");)
	Ad rval{0.0};
	// evaluate the first segment containing t
	for (auto sp : segments) {
		if (sp->contains(t.value()))
			return sp->eval(t);
	}
	double lo;
	double hi;
	this->limits(lo,hi);
	string exc = vastr(this->iid()," called with ",t,
			" which is out of the range ",lo," to ",hi);
	T_(trc.dprint("throwing exception: ",exc);)
	throw std::range_error(exc);
	return rval;
}

bool
Interp::
eval (vector<Ad> const& t, vector<complex<Ad>>& result) {
// special eval function for vspline2/3
	for (auto sp : segments) {
		vspline* vp = dynamic_cast<vspline*>(sp);
		if (vp != nullptr) {
			return vp->eval(t[0], result);
		}
		vspline2* vp2 = dynamic_cast<vspline2*>(sp);
		vspline3* vp3 = dynamic_cast<vspline3*>(sp);
		if (vp2 != nullptr) {
			return vp2->eval(t, result);
		} else if (vp3 != nullptr) {
			return vp3->eval(t, result);
		}
	}
	// no segment in limits - throw an exception
	double lo;
	double hi;
	this->limits(lo,hi);
	string exc = vastr(this->iid()," called with ",t,
			" which is out of the range ",lo," to ",hi);
	throw std::range_error(exc);
	//  return false;
}


void
Interp::
plot(const std::string& plotfile,  int nstep) {
// plot each segment in an Interp, one curve per segment
	T_(Trace trc(1,"Interp::plot");)

	// limits for each segment to ensure adjacent segment end
	// points are included in a segment's curve
	double t0{0.0};
	double t1{3.0};
	vector<double> lim;
	for (auto sp : segments) {
		sp->limits(t0, t1);
		lim.push_back(t0);
		lim.push_back(t1);
	}

	vector<string> names;
	vector<string> dep{dependson()};
	names.push_back(dep[0]);
	names.push_back(dataname());
	
	// each segment is a separate curve
	int i;
	int nc{0}; // corner number
	bool append{false};
	for (auto sp : segments) {
		// curve id: the type of interpolation
		string cid = sp->vid();
		// if multiple segments of same type add an int to the name
		if (cid == "corner")
			cid = vastr("corner",++nc);
		// if this is a piecewise linear segment only plot the breakpoints
		vector<double> x;
		vector<double> y;
		if (cid == "plinear") {
			plinear* pl = dynamic_cast<plinear*>(sp);
			x = pl->x;
			y = pl->y;
		} else {
			// limits for this segment
			sp->limits(t0, t1);
			// make up a list of points to plot: nstep pts + end pts of adjacent
			// segments (only points within this segment's limits are plotted)
			vector<double> ts = lim;
			double del{(t1-t0)/(nstep-1)};
			for (i=0; i<nstep; i++)
				ts.push_back(t0 + i*del);
			std::sort(ts.begin(), ts.end());
			for (auto t : ts) {
				if (sp->contains(t)) {
					x.push_back(t);
					y.push_back(sp->eval(t));
				}
			}
		}
		vector<vector<double>> xs;
		xs.push_back(x);
		xs.push_back(y);
		// Apf::plot writes one curve
		Apf::plot(plotfile, cid, xs, names, append);
		append = true;
	}
}

ostream&
operator<<(ostream& s, Interp& t) {
	s << "Interp " << t.iid() << ": ";
	string sep("");
	for (auto pi : t.dependson()) {
		s << sep << pi;
		sep = ", ";
	}
	s << " vs " << t.dataname_;
	if (t.segments.size() > 1)
		s << " with " << t.segments.size() << " segments";
	return s;
}

// convenience functions
// plsmooth(x,y,width):   piecewise-linear interpolation with smooothed corners
Interp*
plsmooth(const vector<double>& x, const vector<double>& y,
		double width, const string& iid, const string& xname, const string& yname) {
	Interp* rval{nullptr};
	plinear* pl = new plinear(x,y);
	vector<corner*> corners = plinear_smooth(*pl, width);
	string id = iid;
	if (id.empty())
		id = "smoothed piecewise-linear";
	string xnm = xname;
	if (xnm.empty())
		xnm = "x";
	string ynm = yname;
	if (ynm.empty())
		ynm = "y";
	vector<string> pnames{xnm};
	rval = new Interp(id, pnames, ynm);
	// add the corners first so they take precedence over plinear
	for (auto ip : corners)
		rval->add(ip);

	rval->add(pl);
	return rval;
}

vector<corner*>
plinear_smooth(const plinear& pl, double width) {
// create a corner at each internal junction between piecewise-linear
// segments, with width "width"
	vector<corner*> rval;

	int ncorn = pl.x.size() - 2;
	if (ncorn < 1)
		return rval;

	for (int i=0; i<ncorn; i++) {
		// skip if a previous corner overlaps this point
		bool skip{false};
		for (auto cp : rval) {
			if (cp->contains(pl.x[i+1])) {
				skip = true;
				break;
			}
		}
		if (skip)
			continue;
		double eps{1048.0*std::numeric_limits<double>::epsilon()};
		double x0 = std::max(pl.x[i+1] - width, pl.x[i]+eps);
		double y0 = pl.eval(x0);
		// double yp0 = (pl.y[i+1]-pl.y[i])/(pl.x[i+1]-pl.x[i]);
		double yp0 = pl.deriv(x0);
		double x1 = std::min(pl.x[i+1] + width, pl.x[i+2]-eps);
		if (pl.x[i+1] == pl.x[i+2])
			x1 = pl.x[i+1]+width;
		double y1 = pl.eval(x1);
		// double yp1 = (pl.y[i+2]-pl.y[i+1])/(pl.x[i+2]-pl.x[i+1]);
		double yp1 = pl.deriv(x1);
		rval.push_back(new corner(x0,y0,yp0,x1,y1,yp1));
	}
	return rval;
}


static int CUBGCV_F77 (double const* x, double const* f, double* df, int const* n,
		double* y, double* c, int const* ic, double* var, double* rho,
		int const* job, double* se, double* wk, int* ier);

int
cubgcv (double const* x, double* f, int n, double rho, double* c) {
/*
 * 1D interpolation with optional smoothing of f(x) at
 * n breakpoints. If rho=0 => interpolation, if rho = infinity
 * the result is a least-squares line
 * Output:
 *    c   double (n-1, 4) array of interpolation coefficients
 *        such that
 *           f(t) = ((c(i,0)*d + c(i,1))*d + c(i,2))*d + c(i,3)
 *              d = t - x_i     x_i <= t <= x_{i+1}
 *    f       new values of the breakpoints, moved to smooth
 *            the curve
 */
	double stddev = 1.0;
	double var = -1.0;
	int nm = n - 1;

	vector<double> df(n, stddev);
	vector<double> coef(nm*3, 0.0);
	vector<double> y(n, 0.0);
	vector<double> se(n, 0.0);
	vector<double> wk(7*(n+2), 0.0);

	int f_n = n;
	int f_ic = nm;
	int f_job = 0;
	int ier = 0;

	if (rho == 0.0)
		var = 0.0;

	CUBGCV_F77(x, f, &df[0], &f_n, &y[0], &coef[0], &f_ic, &var, &rho,
			&f_job, &se[0], &wk[0], &ier);

	// put y and coef into the return coefficient array c
	// so that
	//    a(t) = ((c(i,0)*d + c(i,1))*d + c(i,2))*d + c(i,3)
	//
	// compute largest difference between breakpoints before
	// and after smoothing relative to the largest breakpoint
	/*
	 * coef    a (n-1, 4, nix) array
	 * cubgcv stores the Y array backwards from the ABCD splines.
	 * cubgcv:
	 *    s(t) = ((c(i,2)*d + c(i,1))*d + c(i,0))*d + y(i)   d = t - x(i)
	 * ABCD:
	 *    a(t) = ((c(i,0)*d + c(i,1))*d + c(i,2))*d + c(i,3)
	 */
	double dmag = blas::normalize_inf(nm, &y[0], 1);
	double eps = sqrt(std::numeric_limits<double>::epsilon());
	vector<double> absdiff(nm, 0.0);
	vector<double> reldiff(nm, 0.0);
	double maxreldiff = 0.0;

	for (int i=0; i<nm; i++) {
		absdiff[i] = std::abs(f[i] - y[i]);
		if (dmag > eps)
			reldiff[i] = absdiff[i]/dmag;
		else
			reldiff[i] = 0.0;
		if (reldiff[i] > maxreldiff) {
			maxreldiff = reldiff[i];
		}
		c[IJ(i,3,nm)] = y[i];
		for (int j=0; j<3; j++) {
			c[IJ(i,j,nm)] = coef[IJ(i,2-j,nm)];
		}
	}
	// replace f with the smoothed breakpoints (y)
	blas::copy(n, &y[0], 1, f, 1);

	return ier;
}

/* cubgcv.f -- translated by f2c (version 20021022).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

// #include "f2c.h"

static int
CUBGCV_F77 (double const* x, double const* f, double* df, int const* n,
		double* y, double* c__, int const* ic, double* var, double* rho,
		int const* job, double* se, double* wk, int* ier) {

/*     ALGORITHM 642 COLLECTED ALGORITHMS FROM ACM. */
/*     ALGORITHM APPEARED IN ACM-TRANS. MATH. SOFTWARE, VOL.12, NO. 2, */
/*     JUN., 1986, P. 150. */
/*   SUBROUTINE NAME     - CUBGCV */

/* -------------------------------------------------------------------------- */

/*   COMPUTER            - VAX/DOUBLE */

/*   AUTHOR              - M.F.HUTCHINSON */
/*                         CSIRO DIVISION OF MATHEMATICS AND STATISTICS */
/*                         P.O. BOX 1965 */
/*                         CANBERRA, ACT 2601 */
/*                         AUSTRALIA */

/*   LATEST REVISION     - 15 AUGUST 1985 */

/*   PURPOSE             - CUBIC SPLINE DATA SMOOTHER */

/*   USAGE               - CALL CUBGCV (X,F,DF,N,Y,C,IC,VAR,JOB,SE,WK,IER) */

/*   ARGUMENTS    X      - VECTOR OF LENGTH N CONTAINING THE */
/*                           ABSCISSAE OF THE N DATA POINTS */
/*                           (X(I),F(I)) I=1..N. (INPUT) X */
/*                           MUST BE ORDERED SO THAT */
/*                           X(I) .LT. X(I+1). */
/*                F      - VECTOR OF LENGTH N CONTAINING THE */
/*                           ORDINATES (OR FUNCTION VALUES) */
/*                           OF THE N DATA POINTS (INPUT). */
/*                DF     - VECTOR OF LENGTH N. (INPUT/OUTPUT) */
/*                           DF(I) IS THE RELATIVE STANDARD DEVIATION */
/*                           OF THE ERROR ASSOCIATED WITH DATA POINT I. */
/*                           EACH DF(I) MUST BE POSITIVE.  THE VALUES IN */
/*                           DF ARE SCALED BY THE SUBROUTINE SO THAT */
/*                           THEIR MEAN SQUARE VALUE IS 1, AND UNSCALED */
/*                           AGAIN ON NORMAL EXIT. */
/*                           THE MEAN SQUARE VALUE OF THE DF(I) IS RETURNED */
/*                           IN WK(7) ON NORMAL EXIT. */
/*                           IF THE ABSOLUTE STANDARD DEVIATIONS ARE KNOWN, */
/*                           THESE SHOULD BE PROVIDED IN DF AND THE ERROR */
/*                           VARIANCE PARAMETER VAR (SEE BELOW) SHOULD THEN */
/*                           BE SET TO 1. */
/*                           IF THE RELATIVE STANDARD DEVIATIONS ARE UNKNOWN, */
/*                           SET EACH DF(I)=1. */
/*                N      - NUMBER OF DATA POINTS (INPUT). */
/*                           N MUST BE .GE. 3. */
/*                Y,C    - SPLINE COEFFICIENTS. (OUTPUT) Y */
/*                           IS A VECTOR OF LENGTH N. C IS */
/*                           AN N-1 BY 3 MATRIX. THE VALUE */
/*                           OF THE SPLINE APPROXIMATION AT T IS */
/*                           S(T)=((C(I,3)*D+C(I,2))*D+C(I,1))*D+Y(I) */
/*                           WHERE X(I).LE.T.LT.X(I+1) AND */
/*                           D = T-X(I). */
/*                IC     - ROW DIMENSION OF MATRIX C EXACTLY */
/*                           AS SPECIFIED IN THE DIMENSION */
/*                           STATEMENT IN THE CALLING PROGRAM. (INPUT) */
/*                VAR    - ERROR VARIANCE. (INPUT/OUTPUT) */
/*                           IF VAR IS NEGATIVE (I.E. UNKNOWN) THEN */
/*                           THE SMOOTHING PARAMETER IS DETERMINED */
/*                           BY MINIMIZING THE GENERALIZED CROSS VALIDATION */
/*                           AND AN ESTIMATE OF THE ERROR VARIANCE IS */
/*                           RETURNED IN VAR. */
/*                           IF VAR IS NON-NEGATIVE (I.E. KNOWN) THEN THE */
/*                           SMOOTHING PARAMETER IS DETERMINED TO MINIMIZE */
/*                           AN ESTIMATE, WHICH DEPENDS ON VAR, OF THE TRUE */
/*                           MEAN SQUARE ERROR, AND VAR IS UNCHANGED. */
/*                           IN PARTICULAR, IF VAR IS ZERO, THEN AN */
/*                           INTERPOLATING NATURAL CUBIC SPLINE IS CALCULATED. */
/*                           VAR SHOULD BE SET TO 1 IF ABSOLUTE STANDARD */
/*                           DEVIATIONS HAVE BEEN PROVIDED IN DF (SEE ABOVE). */
/*                JOB    - JOB SELECTION PARAMETER. (INPUT) */
/*                         JOB = 0 SHOULD BE SELECTED IF POINT STANDARD ERROR */
/*                           ESTIMATES ARE NOT REQUIRED IN SE. */
/*                         JOB = 1 SHOULD BE SELECTED IF POINT STANDARD ERROR */
/*                           ESTIMATES ARE REQUIRED IN SE. */
/*                SE     - VECTOR OF LENGTH N CONTAINING BAYESIAN STANDARD */
/*                           ERROR ESTIMATES OF THE FITTED SPLINE VALUES IN Y. */
/*                           SE IS NOT REFERENCED IF JOB=0. (OUTPUT) */
/*                WK     - WORK VECTOR OF LENGTH 7*(N + 2). ON NORMAL EXIT THE */
/*                           FIRST 7 VALUES OF WK ARE ASSIGNED AS FOLLOWS:- */

/*                           WK(1) = SMOOTHING PARAMETER (= RHO/(RHO + 1)) */
/*                           WK(2) = ESTIMATE OF THE NUMBER OF DEGREES OF */
/*                                   FREEDOM OF THE RESIDUAL SUM OF SQUARES */
/*                           WK(3) = GENERALIZED CROSS VALIDATION */
/*                           WK(4) = MEAN SQUARE RESIDUAL */
/*                           WK(5) = ESTIMATE OF THE TRUE MEAN SQUARE ERROR */
/*                                   AT THE DATA POINTS */
/*                           WK(6) = ESTIMATE OF THE ERROR VARIANCE */
/*                           WK(7) = MEAN SQUARE VALUE OF THE DF(I) */

/*                           IF WK(1)=0 (RHO=0) AN INTERPOLATING NATURAL CUBIC */
/*                           SPLINE HAS BEEN CALCULATED. */
/*                           IF WK(1)=1 (RHO=INFINITE) A LEAST SQUARES */
/*                           REGRESSION LINE HAS BEEN CALCULATED. */
/*                           WK(2) IS AN ESTIMATE OF THE NUMBER OF DEGREES OF */
/*                           FREEDOM OF THE RESIDUAL WHICH REDUCES TO THE */
/*                           USUAL VALUE OF N-2 WHEN A LEAST SQUARES REGRESSION */
/*                           LINE IS CALCULATED. */
/*                           WK(3),WK(4),WK(5) ARE CALCULATED WITH THE DF(I) */
/*                           SCALED TO HAVE MEAN SQUARE VALUE 1.  THE */
/*                           UNSCALED VALUES OF WK(3),WK(4),WK(5) MAY BE */
/*                           CALCULATED BY DIVIDING BY WK(7). */
/*                           WK(6) COINCIDES WITH THE OUTPUT VALUE OF VAR IF */
/*                           VAR IS NEGATIVE ON INPUT.  IT IS CALCULATED WITH */
/*                           THE UNSCALED VALUES OF THE DF(I) TO FACILITATE */
/*                           COMPARISONS WITH A PRIORI VARIANCE ESTIMATES. */

/*                IER    - ERROR PARAMETER. (OUTPUT) */
/*                         TERMINAL ERROR */
/*                           IER = 129, IC IS LESS THAN N-1. */
/*                           IER = 130, N IS LESS THAN 3. */
/*                           IER = 131, INPUT ABSCISSAE ARE NOT */
/*                             ORDERED SO THAT X(I).LT.X(I+1). */
/*                           IER = 132, DF(I) IS NOT POSITIVE FOR SOME I. */
/*                           IER = 133, JOB IS NOT 0 OR 1. */

/*   PRECISION/HARDWARE  - DOUBLE */

/*   REQUIRED ROUTINES   - SPINT1,SPFIT1,SPCOF1,SPERR1 */

/*   REMARKS      THE NUMBER OF ARITHMETIC OPERATIONS REQUIRED BY THE */
/*                SUBROUTINE IS PROPORTIONAL TO N.  THE SUBROUTINE */
/*                USES AN ALGORITHM DEVELOPED BY M.F. HUTCHINSON AND */
/*                F.R. DE HOOG, 'SMOOTHING NOISY DATA WITH SPLINE */
/*                FUNCTIONS', NUMER. MATH. (IN PRESS) */

/* ----------------------------------------------------------------------- */
    /* Initialized data */

    double ratio = 2.;
    double tau = 1.618033989;
    double zero = 0.;
    double one = 1.;

    /* System generated locals */
    int c_dim1, c_offset, wk_dim1, wk_offset, i__1;

    /* Local variables */
    int i__;
    double p, q, r1, r2, r3, r4, gf1, gf2, gf3, gf4, avh, err, 
	    avdf, avar, stat[6], delta;
    extern int spcof1_(double const*, double *, double const*, double *,
			 int const*, double *, double *,
	     double *, double *, int const*, double *, double *
	    ), spfit1_(double const*, double *, double *, int const*, 
	    double *, double *, double *, double *, 
	    double *, double *, double *, double *, int const*,
	     double *, double *, double *, double *),
		 sperr1_(double const*, double *, double *, int const*, double *,
	     double *, double *, double *),
		 spint1_(double const*, double *, double const*, double *, double *, int const*,
	     double *, double *, int const* , double *, double *, int*);


/* ---SPECIFICATIONS FOR ARGUMENTS--- */

/* ---SPECIFICATIONS FOR LOCAL VARIABLES--- */

    /* Parameter adjustments */
    wk_dim1 = *n + 1 - 0 + 1;
    wk_offset = 0 + wk_dim1;
    wk -= wk_offset;
    --se;
    --y;
    --df;
    --f;
    --x;
    c_dim1 = *ic;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */

/* ---INITIALIZE--- */
    *ier = 133;
    if (*job < 0 || *job > 1) {
	goto L140;
    }
    spint1_(&x[1], &avh, &f[1], &df[1], &avdf, n, &y[1], &c__[c_offset], ic, &
	    wk[wk_offset], &wk[wk_dim1 * 4], ier);
    if (*ier != 0) {
	goto L140;
    }
    avar = *var;
    if (*var > zero) {
	avar = *var * avdf * avdf;
    }

	 // Modification: new arg (rho)
	if (*rho != zero) {
		r1 = *rho;
		goto L90;
	}
/* ---CHECK FOR ZERO VARIANCE--- */
	if (*var == zero) {
		r1 = zero;
		goto L90;
	}

/* ---FIND LOCAL MINIMUM OF GCV OR THE EXPECTED MEAN SQUARE ERROR--- */
    r1 = one;
    r2 = ratio * r1;
    spfit1_(&x[1], &avh, &df[1], n, &r2, &p, &q, &gf2, &avar, stat, &y[1], &
	    c__[c_offset], ic, &wk[wk_offset], &wk[wk_dim1 * 4], &wk[wk_dim1 *
	     6], &wk[wk_dim1 * 7]);
L20:
    spfit1_(&x[1], &avh, &df[1], n, &r1, &p, &q, &gf1, &avar, stat, &y[1], &
	    c__[c_offset], ic, &wk[wk_offset], &wk[wk_dim1 * 4], &wk[wk_dim1 *
	     6], &wk[wk_dim1 * 7]);
    if (gf1 > gf2) {
	goto L30;
    }

/* ---EXIT IF P ZERO--- */
    if (p <= zero) {
	goto L100;
    }
    r2 = r1;
    gf2 = gf1;
    r1 /= ratio;
    goto L20;
L30:
    r3 = ratio * r2;
L40:
    spfit1_(&x[1], &avh, &df[1], n, &r3, &p, &q, &gf3, &avar, stat, &y[1], &
	    c__[c_offset], ic, &wk[wk_offset], &wk[wk_dim1 * 4], &wk[wk_dim1 *
	     6], &wk[wk_dim1 * 7]);
    if (gf3 > gf2) {
	goto L50;
    }

/* ---EXIT IF Q ZERO--- */
    if (q <= zero) {
	goto L100;
    }
    r2 = r3;
    gf2 = gf3;
    r3 = ratio * r3;
    goto L40;
L50:
    r2 = r3;
    gf2 = gf3;
    delta = (r2 - r1) / tau;
    r4 = r1 + delta;
    r3 = r2 - delta;
    spfit1_(&x[1], &avh, &df[1], n, &r3, &p, &q, &gf3, &avar, stat, &y[1], &
	    c__[c_offset], ic, &wk[wk_offset], &wk[wk_dim1 * 4], &wk[wk_dim1 *
	     6], &wk[wk_dim1 * 7]);
    spfit1_(&x[1], &avh, &df[1], n, &r4, &p, &q, &gf4, &avar, stat, &y[1], &
	    c__[c_offset], ic, &wk[wk_offset], &wk[wk_dim1 * 4], &wk[wk_dim1 *
	     6], &wk[wk_dim1 * 7]);

/* ---GOLDEN SECTION SEARCH FOR LOCAL MINIMUM--- */
L60:
    if (gf3 > gf4) {
	goto L70;
    }
    r2 = r4;
    gf2 = gf4;
    r4 = r3;
    gf4 = gf3;
    delta /= tau;
    r3 = r2 - delta;
	 if (r3 == r2) {
		 *ier = 1;
		 return 0;
	}
    spfit1_(&x[1], &avh, &df[1], n, &r3, &p, &q, &gf3, &avar, stat, &y[1], &
	    c__[c_offset], ic, &wk[wk_offset], &wk[wk_dim1 * 4], &wk[wk_dim1 *
	     6], &wk[wk_dim1 * 7]);
    goto L80;
L70:
    r1 = r3;
    gf1 = gf3;
    r3 = r4;
    gf3 = gf4;
    delta /= tau;
    r4 = r1 + delta;
    spfit1_(&x[1], &avh, &df[1], n, &r4, &p, &q, &gf4, &avar, stat, &y[1], &
	    c__[c_offset], ic, &wk[wk_offset], &wk[wk_dim1 * 4], &wk[wk_dim1 *
	     6], &wk[wk_dim1 * 7]);
L80:
    err = (r2 - r1) / (r1 + r2);
    if (err * err + one > one && err > 1e-6) {
	goto L60;
    }
    r1 = (r1 + r2) * .5;

/* ---CALCULATE SPLINE COEFFICIENTS--- */
L90:
    spfit1_(&x[1], &avh, &df[1], n, &r1, &p, &q, &gf1, &avar, stat, &y[1], &
	    c__[c_offset], ic, &wk[wk_offset], &wk[wk_dim1 * 4], &wk[wk_dim1 *
	     6], &wk[wk_dim1 * 7]);
L100:
    spcof1_(&x[1], &avh, &f[1], &df[1], n, &p, &q, &y[1], &c__[c_offset], ic, 
	    &wk[wk_dim1 * 6], &wk[wk_dim1 * 7]);

/* ---OPTIONALLY CALCULATE STANDARD ERROR ESTIMATES--- */
    if (*var >= zero) {
	goto L110;
    }
    avar = stat[5];
    *var = avar / (avdf * avdf);
L110:
    if (*job == 1) {
	sperr1_(&x[1], &avh, &df[1], n, &wk[wk_offset], &p, &avar, &se[1]);
    }

/* ---UNSCALE DF--- */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	df[i__] *= avdf;
/* L120: */
    }

/* --PUT STATISTICS IN WK--- */
    for (i__ = 0; i__ <= 5; ++i__) {
	wk[i__ + wk_dim1] = stat[i__];
/* L130: */
    }
    wk[wk_dim1 + 5] = stat[5] / (avdf * avdf);
    wk[wk_dim1 + 6] = avdf * avdf;
    goto L150;

/* ---CHECK FOR ERROR CONDITION--- */
L140:
/*     IF (IER.NE.0) CONTINUE */
L150:
    return 0;
} /* cubgcv_ */

int
spint1_(double const* x, double *avh, double const* y, 
	double *dy, double *avdy, int const* n, double *a, 
	double *c__, int const* ic, double *r__, double *t, int* ier)
{
    /* Initialized data */

    double zero = 0.;

    /* System generated locals */
    int c_dim1, c_offset, r_dim1, r_offset, t_dim1, t_offset, i__1;

    /* Builtin functions */

    /* Local variables */
    double e, f, g, h__;
    int i__;


/* INITIALIZES THE ARRAYS C, R AND T FOR ONE DIMENSIONAL CUBIC */
/* SMOOTHING SPLINE FITTING BY SUBROUTINE SPFIT1.  THE VALUES */
/* DF(I) ARE SCALED SO THAT THE SUM OF THEIR SQUARES IS N */
/* AND THE AVERAGE OF THE DIFFERENCES X(I+1) - X(I) IS CALCULATED */
/* IN AVH IN ORDER TO AVOID UNDERFLOW AND OVERFLOW PROBLEMS IN */
/* SPFIT1. */

/* SUBROUTINE SETS IER IF ELEMENTS OF X ARE NON-INCREASING, */
/* IF N IS LESS THAN 3, IF IC IS LESS THAN N-1 OR IF DY(I) IS */
/* NOT POSITIVE FOR SOME I. */

/* ---SPECIFICATIONS FOR ARGUMENTS--- */

/* ---SPECIFICATIONS FOR LOCAL VARIABLES--- */
    /* Parameter adjustments */
    t_dim1 = *n + 1 - 0 + 1;
    t_offset = 0 + t_dim1;
    t -= t_offset;
    r_dim1 = *n + 1 - 0 + 1;
    r_offset = 0 + r_dim1;
    r__ -= r_offset;
    --a;
    --dy;
    --y;
    --x;
    c_dim1 = *ic;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */

/* ---INITIALIZATION AND INPUT CHECKING--- */
    *ier = 0;
    if (*n < 3) {
	goto L60;
    }
    if (*ic < *n - 1) {
	goto L70;
    }

/* ---GET AVERAGE X SPACING IN AVH--- */
    g = zero;
    i__1 = *n - 1;
	for (i__ = 1; i__ <= i__1; ++i__) {
		h__ = x[i__ + 1] - x[i__];
		if (h__ <= zero) {
	    goto L80;
		}
		g += h__;
	}
	*avh = g / (*n - 1);

/* ---SCALE RELATIVE WEIGHTS--- */
    g = zero;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (dy[i__] <= zero) {
	    goto L90;
	}
	g += dy[i__] * dy[i__];
/* L20: */
    }
    *avdy = sqrt(g / *n);

    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	dy[i__] /= *avdy;
/* L30: */
    }

/* ---INITIALIZE H,F--- */
    h__ = (x[2] - x[1]) / *avh;
    f = (y[2] - y[1]) / h__;

/* ---CALCULATE A,T,R--- */
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	g = h__;
	h__ = (x[i__ + 1] - x[i__]) / *avh;
	e = f;
	f = (y[i__ + 1] - y[i__]) / h__;
	a[i__] = f - e;
	t[i__ + t_dim1] = (g + h__) * 2. / 3.;
	t[i__ + (t_dim1 << 1)] = h__ / 3.;
	r__[i__ + r_dim1 * 3] = dy[i__ - 1] / g;
	r__[i__ + r_dim1] = dy[i__ + 1] / h__;
	r__[i__ + (r_dim1 << 1)] = -dy[i__] / g - dy[i__] / h__;
/* L40: */
    }

/* ---CALCULATE C = R'*R--- */
    r__[*n + (r_dim1 << 1)] = zero;
    r__[*n + r_dim1 * 3] = zero;
    r__[*n + 1 + r_dim1 * 3] = zero;
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	c__[i__ + c_dim1] = r__[i__ + r_dim1] * r__[i__ + r_dim1] + r__[i__ + 
		(r_dim1 << 1)] * r__[i__ + (r_dim1 << 1)] + r__[i__ + r_dim1 *
		 3] * r__[i__ + r_dim1 * 3];
	c__[i__ + (c_dim1 << 1)] = r__[i__ + r_dim1] * r__[i__ + 1 + (r_dim1 
		<< 1)] + r__[i__ + (r_dim1 << 1)] * r__[i__ + 1 + r_dim1 * 3];
	c__[i__ + c_dim1 * 3] = r__[i__ + r_dim1] * r__[i__ + 2 + r_dim1 * 3];
/* L50: */
    }
    return 0;

/* ---ERROR CONDITIONS--- */
L60:
    *ier = 130;
    return 0;
L70:
    *ier = 129;
    return 0;
L80:
    *ier = 131;
    return 0;
L90:
    *ier = 132;
    return 0;
} /* spint1_ */

int
spfit1_(double const* x, double *avh, double *dy, 
	int const* n, double *rho, double *p, double *q, double* fun,
	double *var, double *stat, double *a, double *c__,
	int const* ic, double *r__, double *t, double *u, double *v) {
    /* Initialized data */

    double zero = 0.;
    double one = 1.;
    double two = 2.;

    /* System generated locals */
    int c_dim1, c_offset, r_dim1, r_offset, t_dim1, t_offset, i__1;
    double d__1;

    /* Local variables */
    double e, f, g, h__;
    int i__;
    double rho1;


/* FITS A CUBIC SMOOTHING SPLINE TO DATA WITH RELATIVE */
/* WEIGHTING DY FOR A GIVEN VALUE OF THE SMOOTHING PARAMETER */
/* RHO USING AN ALGORITHM BASED ON THAT OF C.H. REINSCH (1967), */
/* NUMER. MATH. 10, 177-183. */

/* THE TRACE OF THE INFLUENCE MATRIX IS CALCULATED USING AN */
/* ALGORITHM DEVELOPED BY M.F.HUTCHINSON AND F.R.DE HOOG (NUMER. */
/* MATH., IN PRESS), ENABLING THE GENERALIZED CROSS VALIDATION */
/* AND RELATED STATISTICS TO BE CALCULATED IN ORDER N OPERATIONS. */

/* THE ARRAYS A, C, R AND T ARE ASSUMED TO HAVE BEEN INITIALIZED */
/* BY THE SUBROUTINE SPINT1.  OVERFLOW AND UNDERFLOW PROBLEMS ARE */
/* AVOIDED BY USING P=RHO/(1 + RHO) AND Q=1/(1 + RHO) INSTEAD OF */
/* RHO AND BY SCALING THE DIFFERENCES X(I+1) - X(I) BY AVH. */

/* THE VALUES IN DF ARE ASSUMED TO HAVE BEEN SCALED SO THAT THE */
/* SUM OF THEIR SQUARED VALUES IS N.  THE VALUE IN VAR, WHEN IT IS */
/* NON-NEGATIVE, IS ASSUMED TO HAVE BEEN SCALED TO COMPENSATE FOR */
/* THE SCALING OF THE VALUES IN DF. */

/* THE VALUE RETURNED IN FUN IS AN ESTIMATE OF THE TRUE MEAN SQUARE */
/* WHEN VAR IS NON-NEGATIVE, AND IS THE GENERALIZED CROSS VALIDATION */
/* WHEN VAR IS NEGATIVE. */

/* ---SPECIFICATIONS FOR ARGUMENTS--- */

/* ---LOCAL VARIABLES--- */
    /* Parameter adjustments */
    t_dim1 = *n + 1 - 0 + 1;
    t_offset = 0 + t_dim1;
    t -= t_offset;
    r_dim1 = *n + 1 - 0 + 1;
    r_offset = 0 + r_dim1;
    r__ -= r_offset;
    --a;
    --dy;
    --x;
    --stat;
    c_dim1 = *ic;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */

/* ---USE P AND Q INSTEAD OF RHO TO PREVENT OVERFLOW OR UNDERFLOW--- */
    rho1 = one + *rho;
    *p = *rho / rho1;
    *q = one / rho1;
    if (rho1 == one) {
	*p = zero;
    }
    if (rho1 == *rho) {
	*q = zero;
    }

/* ---RATIONAL CHOLESKY DECOMPOSITION OF P*C + Q*T--- */
    f = zero;
    g = zero;
    h__ = zero;
    for (i__ = 0; i__ <= 1; ++i__) {
	r__[i__ + r_dim1] = zero;
/* L10: */
    }
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	r__[i__ - 2 + r_dim1 * 3] = g * r__[i__ - 2 + r_dim1];
	r__[i__ - 1 + (r_dim1 << 1)] = f * r__[i__ - 1 + r_dim1];
	r__[i__ + r_dim1] = one / (*p * c__[i__ + c_dim1] + *q * t[i__ + 
		t_dim1] - f * r__[i__ - 1 + (r_dim1 << 1)] - g * r__[i__ - 2 
		+ r_dim1 * 3]);
	f = *p * c__[i__ + (c_dim1 << 1)] + *q * t[i__ + (t_dim1 << 1)] - h__ 
		* r__[i__ - 1 + (r_dim1 << 1)];
	g = h__;
	h__ = *p * c__[i__ + c_dim1 * 3];
/* L20: */
    }

/* ---SOLVE FOR U--- */
    u[0] = zero;
    u[1] = zero;
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	u[i__] = a[i__] - r__[i__ - 1 + (r_dim1 << 1)] * u[i__ - 1] - r__[i__ 
		- 2 + r_dim1 * 3] * u[i__ - 2];
/* L30: */
    }
    u[*n] = zero;
    u[*n + 1] = zero;
    for (i__ = *n - 1; i__ >= 2; --i__) {
	u[i__] = r__[i__ + r_dim1] * u[i__] - r__[i__ + (r_dim1 << 1)] * u[
		i__ + 1] - r__[i__ + r_dim1 * 3] * u[i__ + 2];
/* L40: */
    }

/* ---CALCULATE RESIDUAL VECTOR V--- */
    e = zero;
    h__ = zero;
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	g = h__;
	h__ = (u[i__ + 1] - u[i__]) / ((x[i__ + 1] - x[i__]) / *avh);
	v[i__] = dy[i__] * (h__ - g);
	e += v[i__] * v[i__];
/* L50: */
    }
    v[*n] = dy[*n] * (-h__);
    e += v[*n] * v[*n];

/* ---CALCULATE UPPER THREE BANDS OF INVERSE MATRIX--- */
    r__[*n + r_dim1] = zero;
    r__[*n + (r_dim1 << 1)] = zero;
    r__[*n + 1 + r_dim1] = zero;
    for (i__ = *n - 1; i__ >= 2; --i__) {
	g = r__[i__ + (r_dim1 << 1)];
	h__ = r__[i__ + r_dim1 * 3];
	r__[i__ + (r_dim1 << 1)] = -g * r__[i__ + 1 + r_dim1] - h__ * r__[i__ 
		+ 1 + (r_dim1 << 1)];
	r__[i__ + r_dim1 * 3] = -g * r__[i__ + 1 + (r_dim1 << 1)] - h__ * r__[
		i__ + 2 + r_dim1];
	r__[i__ + r_dim1] = r__[i__ + r_dim1] - g * r__[i__ + (r_dim1 << 1)] 
		- h__ * r__[i__ + r_dim1 * 3];
/* L60: */
    }

/* ---CALCULATE TRACE--- */
    f = zero;
    g = zero;
    h__ = zero;
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	f += r__[i__ + r_dim1] * c__[i__ + c_dim1];
	g += r__[i__ + (r_dim1 << 1)] * c__[i__ + (c_dim1 << 1)];
	h__ += r__[i__ + r_dim1 * 3] * c__[i__ + c_dim1 * 3];
/* L70: */
    }
    f += two * (g + h__);

/* ---CALCULATE STATISTICS--- */
    stat[1] = *p;
    stat[2] = f * *p;
    stat[3] = *n * e / (f * f);
    stat[4] = e * *p * *p / *n;
    stat[6] = e * *p / f;
    if (*var >= zero) {
	goto L80;
    }
    stat[5] = stat[6] - stat[4];
    *fun = stat[3];
    goto L90;
L80:
/* Computing MAX */
    d__1 = stat[4] - two * *var * stat[2] / *n + *var;
    stat[5] = std::max(d__1,zero);
    *fun = stat[5];
L90:
    return 0;
} /* spfit1_ */

int
sperr1_(double const* x, double *avh, double *dy, 
	int const* n, double *r__, double *p, double *var, 
	double *se) {
    /* Initialized data */

    double zero = 0.;
    double one = 1.;

    /* System generated locals */
    int r_dim1, r_offset, i__1;
    double d__1;

    /* Builtin functions */

    /* Local variables */
    double f, g, h__;
    int i__;
    double f1, g1, h1;


/* CALCULATES BAYESIAN ESTIMATES OF THE STANDARD ERRORS OF THE FITTED */
/* VALUES OF A CUBIC SMOOTHING SPLINE BY CALCULATING THE DIAGONAL ELEMENTS */
/* OF THE INFLUENCE MATRIX. */

/* ---SPECIFICATIONS FOR ARGUMENTS--- */

/* ---SPECIFICATIONS FOR LOCAL VARIABLES--- */
    /* Parameter adjustments */
    --se;
    r_dim1 = *n + 1 - 0 + 1;
    r_offset = 0 + r_dim1;
    r__ -= r_offset;
    --dy;
    --x;

    /* Function Body */

/* ---INITIALIZE--- */
    h__ = *avh / (x[2] - x[1]);
    se[1] = one - *p * dy[1] * dy[1] * h__ * h__ * r__[r_dim1 + 2];
    r__[r_dim1 + 1] = zero;
    r__[(r_dim1 << 1) + 1] = zero;
    r__[r_dim1 * 3 + 1] = zero;

/* ---CALCULATE DIAGONAL ELEMENTS--- */
    i__1 = *n - 1;
    for (i__ = 2; i__ <= i__1; ++i__) {
	f = h__;
	h__ = *avh / (x[i__ + 1] - x[i__]);
	g = -f - h__;
	f1 = f * r__[i__ - 1 + r_dim1] + g * r__[i__ - 1 + (r_dim1 << 1)] + 
		h__ * r__[i__ - 1 + r_dim1 * 3];
	g1 = f * r__[i__ - 1 + (r_dim1 << 1)] + g * r__[i__ + r_dim1] + h__ * 
		r__[i__ + (r_dim1 << 1)];
	h1 = f * r__[i__ - 1 + r_dim1 * 3] + g * r__[i__ + (r_dim1 << 1)] + 
		h__ * r__[i__ + 1 + r_dim1];
	se[i__] = one - *p * dy[i__] * dy[i__] * (f * f1 + g * g1 + h__ * h1);
/* L10: */
    }
    se[*n] = one - *p * dy[*n] * dy[*n] * h__ * h__ * r__[*n - 1 + r_dim1];

/* ---CALCULATE STANDARD ERROR ESTIMATES--- */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing MAX */
	d__1 = se[i__] * *var;
	se[i__] = sqrt((std::max(d__1,zero))) * dy[i__];
/* L20: */
    }
    return 0;
} /* sperr1_ */

int
spcof1_(double const* x, double *avh, double const* y, double *dy, int const* n,
		double *p, double *q, double *a, double *c__, int const* ic, double *u, double *v) {
    /* System generated locals */
    int c_dim1, c_offset, i__1;

    /* Local variables */
    double h__;
    int i__;
    double qh;


/* CALCULATES COEFFICIENTS OF A CUBIC SMOOTHING SPLINE FROM */
/* PARAMETERS CALCULATED BY SUBROUTINE SPFIT1. */

/* ---SPECIFICATIONS FOR ARGUMENTS--- */

/* ---SPECIFICATIONS FOR LOCAL VARIABLES--- */

/* ---CALCULATE A--- */
    /* Parameter adjustments */
    --a;
    --dy;
    --y;
    --x;
    c_dim1 = *ic;
    c_offset = 1 + c_dim1;
    c__ -= c_offset;

    /* Function Body */
    qh = *q / (*avh * *avh);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	a[i__] = y[i__] - *p * dy[i__] * v[i__];
	u[i__] = qh * u[i__];
/* L10: */
    }

/* ---CALCULATE C--- */
    i__1 = *n - 1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	h__ = x[i__ + 1] - x[i__];
	c__[i__ + c_dim1 * 3] = (u[i__ + 1] - u[i__]) / (h__ * 3.);
	c__[i__ + c_dim1] = (a[i__ + 1] - a[i__]) / h__ - (h__ * c__[i__ + 
		c_dim1 * 3] + u[i__]) * h__;
	c__[i__ + (c_dim1 << 1)] = u[i__];
/* L20: */
    }
    return 0;
} /* spcof1_ */



#ifdef MAIN

#include "conv.h"

void
test_sin(int argc, char** argv) {
	size_t i;
	size_t nx = 21;
	// size_t ny = 1;
	// smoothing?
	double rho = 4.0;
	Interp myinterp("sin", "x", "f");

	// Usage:  interp [rho]
	if (argc > 1)
		str2double(argv[1], rho);

	// interpolate a sin
	vector<double> x(nx);
	vector<double> y;
	double pi{flaps::pi};
	double xmin = 0.0;
	double xmax = 2.0*pi;
	double del = (xmax - xmin)/(double(nx-1));
	for (i=0; i<nx; i++) {
		x[i] = xmin + del*i;
		double yi = sin(x[i]);
		y.push_back(yi);
	}
	try {
		myinterp.add(new spline(x,y,rho));
		myinterp.plot("spline.apf", 101);
	} catch (runtime_error &s) {
		cerr << "interpolation failed: " << s << endl;
	}
}

void
test_discont() {
// smoothed piecewise linear discontinuous fcn
	vector<double> x;
	vector<double> y;
	Interp myinterp("discontinuous", "x", "f");

	x.push_back(0.2);
	y.push_back(0.0);
	x.push_back(1.0);
	y.push_back(0.0);
	x.push_back(1.0);
	y.push_back(1.0);
	x.push_back(2.0);
	y.push_back(0.0);
	x.push_back(2.8);
	y.push_back(0.0);
	plinear* pl = new plinear(x,y);
	double width{0.1};
	vector<corner*> cl = plinear_smooth(*pl, width);
	myinterp.clear();
	for (auto cp : cl)
		myinterp.add(cp);
	myinterp.add(pl);
	myinterp.plot("plinear.apf", 101);
}

int
main (int argc, char **argv) {
	T_(Trace trc(1,"interp");)

	test_sin(argc, argv);
	test_discont();

}
#endif /* MAIN */
