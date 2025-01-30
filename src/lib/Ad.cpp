//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// classes for automatic differentiation (Ad)
//   Ad::initialize(const vector<string>& names)
//      this is the *only* way to set the names of parameters to be
//      used as autodiff derivatives and it must be called only once.
//      It will also reallocate the Ad members of each parameter in the
//      gpset, and set parameter derivatives to 1 in the Ad member.
//   set_adpar_derivs(pset& pl)
//      for each Adparname, find the parameter in pl and set the
//      corresponding Ad derivative to 1.0
// namespace Ad comprises functions for maintaining a list
// of parameter names for the automatic-differentiation class Ad

#include <cassert>
#include <cstddef>

#include "Ad.h"
#include "atmos.h"
#include "blas.h"
#include "conv.h"
#include "fptype.h"
#include "matrix.h"
#include "message.h"
#include "trace.h"

using namespace std;

// *the* list of Ad parameters, accessible only in this file
// with the Ad::* functions, public with adnames()
static vector<string> Adparameters;

bool Ad::initialize_called{false};

void
Ad::
initialize(const vector<string>& names) {
// Set the list of Ad parameter names; the list may be accessed with adnames().
// This function must only be called once in a process so that all
// Ad instances have the same derivative parameters.
// After calling this function is may be necessary to call Ad::realloc()
// for any Ad that were created prior to calling this function, e.g.
// the global pset: gpset::get().realloc()
	T_(Trace trc(1,"Ad::initialize");)

	if (initialize_called) {
		string exc = vastr("Ad::initialize has already been called: ",
				Adparameters.size()," parameters set");
		throw runtime_error(exc);
	}
	if (!Adparameters.empty()) {
		string exc = vastr("Ad::initialize: Adparameters is not empty "
				"it has ",Adparameters.size()," derivative parameters: ",
				Ad::toString());
		throw runtime_error(exc);
	}
	if (names.empty())
		throw runtime_error("Ad::initialize called with no parameter names");

	Adparameters = names;
	T_(trc.dprint("autodiff names:",Adparameters);)
	initialize_called = true;
}

void
Ad::
realloc() {
// reallocate data_ to reflect the number of ad parameters; should
// only be neceesary once: right after calling Ad::initialize
	double v{data_[0]};	// save only the value, not derivatives
	delete[] data_;
	data_ = new double[Ad::ndata()];
	std::fill_n(data_, Ad::ndata(), 0.0);
	data_[0] = v;	// restore value
}

int
Ad::
nder() { return Adparameters.size(); }

int
Ad::
ndata() {
// Returns the number of double elements in a real variable, i.e.
// the value and all derivative values
	return Adparameters.size() + 1; }

std::vector<std::string> const&
Ad::
adnames() {
// returns a const reference to the vector of Ad parameter names
	return Adparameters;
}

std::string
Ad::
toString() {
// Returns a string containing the Ad parameter names (comma separated)
// in one line; you may want to run it through roff:
//        roff(rval, indent, linelen);
	std::ostringstream os;

	string sep("");
	for (auto& name : Adparameters) {
		os << sep << name;
		sep = ", ";
	}
	return os.str();
}

int
Ad::
find(std::string const& name) {
// Returns the 0b index into data() of parameter "name";
// if name is empty or == "0" => value
// If "name" is not found in Adparameters, return -1
	if (name.empty() || name == "0")
		return 0;
	for (size_t i=0; i<Adparameters.size(); i++) {
		if (Adparameters[i] == name)
			return (int)(i+1);
	}
	return -1;
}

std::string
Ad::
live() {
	// return vastr(Adarray::live()," Adarray (",Adarray::peak()," peak), ",
	// 		Ad::live()," Ad (",Ad::peak()," peak), ",
	// 		vector<complex<Ad>>::live()," vector<complex<Ad>> (",vector<complex<Ad>>::peak()," peak), ",
	// 		cAd::live()," cAd (",cAd::peak()," peak)");
	return "";
}

int
Ad::
nlive() {
	int rval{0};
	// rval = Adarray::live() + Ad::live();
	// rval += vector<complex<Ad>>::live() + cAd::live();
	return rval;
}


double
Ad::
der(int d) const {
// return the value of derivative "d" (0b) 0-(nder()-1)
	if (d < 0 || d >= (int)Adparameters.size())
		throw runtime_error(vastr("Ad::deriv illegal arg: ",d));
	return data_[d+1];
}

void
Ad::
der(int d, double v) {
// set the value of derivative "d" (0b) 0-(nder()-1) to "v"
	if (d < 0 || d >= (int)Adparameters.size())
		throw runtime_error(vastr("Ad::deriv illegal arg: ",d));
	data_[d+1] = v;
}

std::ostream&
operator<< (std::ostream& s, Ad const& a) {
	s << a.data_[0];
	for (int i=1; i<=Ad::nder(); i++) {
		s << "[" << a.data_[i] << ']';
	}
	return s;
}

bool
is_equal(Ad const& a, Ad const& b, int sigfig) {
	double const* ad = a.data();
	double const* bd = b.data();
	for (int i=0; i<Ad::ndata(); i++) {
		if (!is_equal(ad[i], bd[i], sigfig))
			return false;
	}
	return true;
}

// math functions

Ad
abs (const Ad& x) {
	Ad rval{x};
	for (int i=0; i<Ad::ndata(); i++)
		rval.data_[i] = std::abs(rval.data_[i]);
	return rval;
}

Ad
acos (const Ad& x) {
	Ad rval;
	double xv(x.value());
	rval.value(acos(xv));
	double d = sqrt(1.0 - xv*xv);
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i,-x.der(i)/d);
	return rval;
}

Ad
asin (const Ad& x) {
	Ad rval;
	double xv(x.value());
	rval.value(asin(xv));
	// d/du asin u = 1/sqrt(1-u^2)
	// watch out for 1-u^2 = 0
	double d = sqrt(double(1.0) - xv*xv);
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i,x.der(i)/d);
	return rval;
}

Ad
atan (Ad const& x) {
	Ad rval;
	double xv(x.value());
	rval.value(atan(xv));
	double d = 1.0 + xv*xv;
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i, x.der(i)/d);
	return rval;
}

Ad
atan2 (Ad const& a, Ad const& b) {
	double av{a.value()};
	double bv{b.value()};
	Ad rval;
	rval.value(atan2(av,bv));
	double d(av*av + bv*bv);
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i, (a.der(i)*bv - av*b.der(i))/d);
	return rval;
}

Ad
cos (Ad const& x) {
	Ad rval;
	double xv(x.value());
	rval.value(cos(xv));
	double sx{-sin(xv)};
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i, x.der(i)*sx);
	return rval;
}

Ad
exp (Ad const& x) {
	Ad rval;
	rval.value(exp(x.value()));
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i, x.der(i)*rval.value());
	return rval;
}

Ad
log (Ad const& x) {
	Ad rval;
	rval.value(log(x.value()));
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i, x.der(i)/x.value());
	return rval;
}

Ad
log10 (const Ad& x) {
	Ad rval;
	rval.value(log10(x.value()));
	double le = log10(exp(1.0));
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i,le*x.der(i)/x.value());
	return rval;
}

Ad
min (Ad const& a, Ad const& b) {
	if (a.value() < b.value())
		return a;
	return b;
}

Ad
max (Ad const& a, Ad const& b) {
	if (a.value() > b.value())
		return a;
	return b;
}

Ad
pow (Ad const& x, Ad const& y) {
	Ad rval;
	double u = x.value();
	double v = y.value();
	// u^v	
	double t = pow(u, v);
	// derivatives: v*u^(v-1)(du/dp) + ln(u)*u^v*(dv/dp)
	rval.value(t);
	double dt = v*pow(u, v-1.0);
	double lnu = 0.0;
	if (u > 0.0)
		lnu = log(u);
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i, t*y.der(i)*lnu + dt*x.der(i));
	return rval;
}

Ad
sin (Ad const& x) {
	Ad rval;
	double xv(x.value());
	rval.value(sin(xv));
	double cx{cos(xv)};
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i, x.der(i)*cx);
	return rval;
}

Ad
sqrt (const Ad& x) {
	Ad rval;
	double xv = x.value();
	if (xv < 0.0)
		throw runtime_error(vastr("attempt to take sqrt of ",xv));
	double v(sqrt(xv));
	rval.value(v);
	if (v == 0.0) {
		for (int i=0; i<Ad::nder(); i++)
			rval.der(i, 0.0);
	} else {
		for (int i=0; i<Ad::nder(); i++)
			rval.der(i, x.der(i)/(2.0*v));
	}
	return rval;
}

Ad
tan (Ad const& x) {
	Ad rval;
	double xv(x.value());
	rval.value(tan(xv));
	double sec = 1.0/cos(xv);
	double t(sec*sec);
	for (int i=0; i<Ad::nder(); i++)
		rval.der(i, x.der(i)*t);
	return rval;
}

double
mag(complex<Ad> const& x) {
	double xr = Ad::real(x).data()[0];
	double xi = Ad::imag(x).data()[0];
	return sqrt(xr*xr + xi*xi);
}

//------------------------------------- end of Ad math fcns

// 4 versions of extract:
//    Ad, name -> double
//    Ad, idx -> double
//    complex Ad, name -> complex
//    complex Ad, idx -> complex

bool
extract (const vector<Ad>& from, std::string const& name, double* val) {
// extract from "from" the derivative "name" or the value
// if "name" is empty or "0" (zero), put the result in "val"
	T_(Trace trc(2,"extract(Ad,string)");)
	int idx{0};	// default: value

	T_(trc.dprint("name: ",name);)

	if (!name.empty() && name != "0") {
		idx = Ad::find(name);	// 0b index in data
		if (idx < 0)
			throw runtime_error(vastr(name," is not an automatic-derivative parameter"));
	}

	return extract (from, idx, val);
}

bool
extract (const vector<Ad>& from, int idx, double* val) {
// extract from "from" data element "idx" (0b): the value if idx==0,
// the 1b derivative number (1-Ad::nder()) otherwise.
	T_(Trace trc(2,"extract(Ad,int)");)

	T_(trc.dprint("index ",idx);)

	size_t n = from.size();

	if (idx < 0 || idx >= (int)Ad::ndata())
		throw runtime_error(vastr("attempt to extract index ", idx));

	for (size_t i=0; i<n; i++)
		val[i] = from[i].data()[idx];

	return true;
}

bool
extract (const vector<complex<Ad>>& from, std::string const& name, complex<double>* val) {
// extract from "from" the derivative "name" or the value
// if "name" is empty or "0" (zero), put the result in "val"
	T_(Trace trc(2,"extract(complex<Ad>,string)");)
	int idx{0};	// default: the value

	T_(trc.dprint("name: ",name);)

	if (!name.empty() && name != "0") {
		idx = Ad::find(name);
		if (idx < 0)
			throw runtime_error(vastr(name,
				" is not an automatic-derivative parameter"));
	}

	return extract (from, idx, val);
}

bool
extract (const vector<complex<Ad>>& from, int idx, complex<double>* val) {
// extract from "from" data element "idx" (0b): the value if idx==0,
// the 1b derivative number (1-Ad::nder()) otherwise.
	T_(Trace trc(2,"extract(complex<Ad>,int)");)

	T_(trc.dprint("index ",idx);)

	if (idx < 0 || idx >= (int)Ad::ndata())
		throw runtime_error(vastr("attempt to extract idx(0b) ", idx));

	complex<double>* vp = val;
	for (auto& fi : from)
		*vp++ = complex<double>(fi.real().data()[idx], fi.imag().data()[idx]);

	return true;
}


template<>
template<>
complex<Ad>::
complex(const complex<Ad>& z) {
// copy constructor with fewer Ad ctors
   this->real(Ad::real(z));
   this->imag(Ad::imag(z));
}

void
Ad::
multad (const Ad& a, const Ad& b, Ad& c) {
// multiply c = a*b  works for a = a*b
// This is a helper function for various overrides and blas specializations
   double av = a.data_[0];
   double bv = b.data_[0];
   c.data_[0] = av*bv;
   for (int i=1; i<Ad::ndata(); i++)
      c.data_[i] = a.data_[i]*bv + b.data_[i]*av;
}

void
Ad::
multcad(const complex<Ad>& a, const complex<Ad>& b, complex<Ad>& c, Ad& wk) {
// multiply 2 complex Ad's: c = a*b, without creating new Ad's
// This is a helper function for various overrides and blas specializations
   const Ad& ar = Ad::real(a);
   const Ad& ai = Ad::imag(a);
   const Ad& br = Ad::real(b);
   const Ad& bi = Ad::imag(b);
   Ad& cr = Ad::real(c);
   Ad& ci = Ad::imag(c);
   Ad::multad(ar, br, cr);
   Ad::multad(ai, bi, wk);
   cr -= wk;
   Ad::multad(ar, bi, ci);
   Ad::multad(ai, br, wk);
   ci += wk;
}

template<>
template<>
complex<Ad>&
complex<Ad>::
operator*=(const complex<Ad>& rhs) {
// override std::complex<T>::operator*= to avoid copy ctors
   Ad& ar = Ad::real(*this);
   Ad& ai = Ad::imag(*this);
   const Ad& br = Ad::real(rhs);
   const Ad& bi = Ad::imag(rhs);
   Ad wk(0.0);
   Ad r;          // XXX how to avoid?
   Ad::multad(ar, br, r);
   Ad::multad(ai, bi, wk);
   r -= wk;
   Ad::multad(ai, br, ai);
   Ad::multad(ar, bi, wk);
   ai += wk;
   ar = r;
   return *this;
}

template<>
template<>
complex<Ad>&
complex<Ad>::
operator+=(const complex<Ad>& rhs) {
// override std::complex<T>::operator+= to avoid copy constructors
   Ad& ar = Ad::real(*this);
   Ad& ai = Ad::imag(*this);
   const Ad& br = Ad::real(rhs);
   const Ad& bi = Ad::imag(rhs);
   ar += br;
   ai += bi;
   return *this;
}


adfcn1p
adfcn1(const std::string& fcn) {
// Return a pointer to single-Ad-arg function "fcn"
	T_(Trace trc(3,"adfcn1");)
	adfcn1p rval{nullptr};
	// XXX how can we specify *exactly* which function?
	if (fcn == "abs")
		rval = abs;
	else if (fcn == "acos")
		rval = acos;
	else if (fcn == "asin")
		rval = asin;
	else if (fcn == "atan")
		rval = atan;
	else if (fcn == "cos")
		rval = cos;
	else if (fcn == "exp")
		rval = exp;
	else if (fcn == "log")
		rval = log;
	else if (fcn == "log10")
		rval = log10;
	else if (fcn == "sin")
		rval = sin;
	else if (fcn == "sqrt")
		rval = sqrt;
	else if (fcn == "tan")
		rval = tan;
	else if (fcn == "atmos::rho")
		rval = atmos::rho;
	else if (fcn == "atmos::press")
		rval = atmos::press;
	else if (fcn == "atmos::temp")
		rval = atmos::temp;
	else if (fcn == "atmos::vsonic")
		rval = atmos::vsonic;
	T_(trc.dprint("returning ",rval);)
	return rval;
}

adfcn2p
adfcn2(const std::string& fcn) {
// Return a pointer to double-Ad-arg function "fcn"
	T_(Trace trc(3,"adfcn2");)
	adfcn2p rval{nullptr};
	if (fcn == "atan2")
		rval = atan2;
	else if (fcn == "max")
		rval = max;
	else if (fcn == "min")
		rval = min;
	else if (fcn == "pow")
		rval = pow;
	else if (fcn == "atmos::vcas")
		rval = atmos::vcas;
	else if (fcn == "atmos::cdpress")
		rval = atmos::cdpress;
	T_(trc.dprint("returning ",rval);)
	return rval;
}


#ifdef MAIN
#undef MAIN

#include <complex>
#include <cstdlib>  // for rand()
#include <iostream>
#include <iomanip>

#include "nrc.h"

double
relerrAd(const Ad& a, const Ad& b) {
	double t = abs(a.value() - b.value());
	double am = abs(a.value());
	double bm = abs(b.value());
	double ei = t/max(max(am, bm), 0.1);
	double rval = ei;

	for (int i=0; i<Ad::nder(); i++) {
		double t = abs(a.deriv(i) - b.deriv(i));
		double am = abs(a.deriv(i));
		double bm = abs(b.deriv(i));
		double ei = t/max(max(am, bm), 0.1);
		rval = max(rval, ei);
	}
	return rval;
}

class ad_app {
	Ad (*fcn)(const Ad& a);
public:
	ad_app(Ad (*f)(const Ad&)) : fcn(f) {}
	void operator()(int n, double t, double* ft) {
		Ad a(t);
		Ad at = fcn(a);
		ft[0] = at.value();
	}
};

double
ad_fcn (const string& title, Ad (*fcn)(const Ad& a),
		Ad (*inverse)(const Ad&)) {
// test Ad single-arg functions: sin, cos, etc
	int nstep = 8;
	double pi{flaps::pi};
	double amin = 0.0;
	double amax = pi/2.0;
	double del = (amax - amin)/(double)(nstep-1);
	double norm, maxnorm = 0.0;
	int nd = Ad::nder();
	vector<double> der(nd, 1.0);

	cerr << "------------------- check " << title << endl;

	nstep--;  // avoid 0, pi/2
	for (int i=1; i<nstep; i++) {
		double ai = amin + i*del;
		Ad a(ai);
		a.deriv(0, 1.0);
		Ad b = fcn(a);
		Ad c = inverse(b);
		Ad err = c - a;
		cerr << "      fcn(a) = " << b << endl;
		cerr << "      inverse(fcn(a)) = " << c << endl;
		cerr << title << " err @ " << a << " = " << err << endl;
		// T norm = inner_product(err.begin(), err.end(), err.begin(), T(0.0));
		norm = relerrAd(a,c);
		maxnorm = max(norm, maxnorm);
	}

	// check derivative
	Ad a((amax+amin)/2.0);
	a.deriv(0,1.0);
	Ad exact = fcn(a);
	double dexact = exact.deriv(0);
	double y{0.0};
	double ridders_err;
	NRC::fdapprox(1, a.value(), amin, amax, &y, ad_app(fcn), ridders_err);
	double err = std::abs(y-dexact);
	double relerr = err/std::max(0.1, std::abs(y));
	maxnorm = std::max(maxnorm, relerr);
	cerr << "exact & approx deriv @ " << a.value() << ": " << dexact << ", " << y
		<< ", rel err " << relerr << ", ridders err " << ridders_err << endl;

	// watch for large error
	double eps = sqrt(numeric_limits<double>::epsilon());
	if (maxnorm > eps)
		cerr << "ERROR: ";
	cerr << title << " largest err norm = " << maxnorm << endl;
	return maxnorm;
}


double
ad_fcn2(const string& title, Ad (*fcn)(const Ad& a, const Ad& b),
		Ad (*inverse)(const Ad& a, const Ad& b)) {
// test Ad two-argument functions: operator+, etc
	int nstep = 4;
	double pi{flaps::pi};
	double amin = 0.0;
	double amax = pi/2.0;
	double norm, maxnorm = 0.0;
	double del = (amax - amin)/(double)(nstep-1);
	int nd = Ad::nder();
	vector<double> der(nd, double(0.5));

	cerr << "------------------- check " << title << endl;

	nstep--;  // avoid 0, pi/2
	for (int i=1; i<nstep; i++) {
		double ai = amin + i*del;
		Ad a(ai);
		Ad b(10.0*ai);
		b.deriv(0, 0.5);
		Ad c = fcn(a,b);

		cerr << endl << title << " at " << a << ", " << b << endl;
		cerr << "     fcn = " << c << endl;
		Ad d = inverse(c, b);
		cerr << "     inverse = " << d << endl;
		Ad err{d - a};
		norm = relerrAd(a,d);
		maxnorm = max(norm, maxnorm);
		cerr << " err = " << err << ", rel err = " << norm << endl;
	}

	// watch for large error
	double eps = sqrt(numeric_limits<double>::epsilon());
	if (maxnorm > eps)
		cerr << "ERROR: ";
	cerr << title << " largest err norm = " << maxnorm << endl;
	return maxnorm;
}


double
xrand() {
// Returns a random number between -1 and 1
	int k = rand();
	double rval = (double)k/(double)RAND_MAX;
	rval = 2.0*rval - 1.0;
	return rval;
}

Ad
random_Ad() {
// return an Ad with random value and derivatives
	int nd = Ad::nder();
	Ad rval{xrand()};
	for (int i=0; i<nd; i++)
		rval.deriv(i, xrand());
	return rval;
}

double
ad_mixed_mode() {
// test mixed Ad-double arithmetic
	double rval{0.0};
	Ad a = random_Ad();
	Ad b = random_Ad();

	rval = std::max(rval, relerrAd((a+b)-b, a));
	rval = std::max(rval, relerrAd((a*b)/b, a));
	// Ad and double
	double c = xrand();
	rval = std::max(rval, relerrAd((a+c)-c, a));
	rval = std::max(rval, relerrAd((a*c)/c, a));

	cerr << "returning rel err = ";
	cerr << std::setw(20) << std::setprecision(14) << rval << endl;
	return rval;
}
	
double
test_arith() {
	T_(Trace trc(1,"test_arith");)
	double rval{0.0};
	cerr << "--------------------------------------------------------\n";
	cerr << " Testing Ad arithmetic operations\n";
	cerr << "sizeof(Ad) = " << sizeof(Ad) << endl;

	rval = ad_mixed_mode();
	cerr << "Ad mix-mode err = " << rval << endl;
	return rval;
}

Ad
tan2 (const Ad& c, const Ad& b) {
// the inverse function for atan2: tan(c)*b
	return tan(c)*b;
}

Ad
unary_minus (const Ad& a) {
// for testing unary minus: -x
	return -a;
}

Ad
square(const Ad& a) {
// the inverse function for sqrt: a*a
	return a*a;
}

double
sgemv_test() {
	T_(Trace trc(1,"sgemv_test");)
	int n = 3;
	vector<Ad> a(n*n);
	vector<Ad> x(n);
	Ad p(2.0);
	p.deriv(0, 0.5);
	for (int i=0; i<n; i++) {
		a[i*(n+1)] = 1.0;
		x[i] = p;
	}
	vector<Ad> y(n);

	T_(trc.dprint("test_sgemv a:\n",a);)
	T_(trc.dprint("test_sgemv x:\n",x);)
	blas_sgemv("n", n, n, Ad(1.0), &a[0], n, &x[0], 1, Ad(0.0), &y[0], 1);
	T_(trc.dprint("test_sgemv y:\n",y);)
	Ad t{0.0};
	for (int i=0; i<n; i++) {
		Ad ti = x[i] - y[i];
		t += ti*ti;
	}
	double rval = sqrt(t).value();
	return rval;
}

double
cgemv_test() {
	T_(Trace trc(1,"cgemv_test");)
	int n = 3;
	vector<complex<Ad>> a(n*n);
	vector<complex<Ad>> x(n);
	// note: we cannot say
	//    p.real().deriv(0, 0.5);
	// must do this:
	Ad pr(2.0);
	pr.deriv(0,0.5);
	complex<Ad> p(pr);

	for (int i=0; i<n; i++) {
		a[i*(n+1)] = 1.0;
		x[i] = p;
	}
	vector<complex<Ad>> y(n);

	T_(trc.dprint("cgemv_test a:\n",a);)
	T_(trc.dprint("cgemv_test x:\n",x);)
	blas_cgemv("n", n, n, complex<Ad>(1.0), &a[0], n, &x[0], 1,
		complex<Ad>(0.0), &y[0], 1);
	T_(trc.dprint("cgemv_test y:\n",y);)
	// return err
	complex<Ad> t{0.0};
	for (int i=0; i<n; i++) {
		complex<Ad> ti = x[i] - y[i];
		t += ti*ti;
	}
	double rval = std::abs(sqrt(t)).value();
	return rval;
}

bool
expression_test() {
	double dpr{1.0};
	double dpi{0.0};
	complex<double> dp(dpr,dpi);
	Ad pr{dpr};
	pr.deriv(1,1.0);	// set deriv 1 to 1: p is par2
	Ad pi{dpi};
	complex<Ad> p(pr,pi);
	double R0{2.0};
	double R1{4.0};
	double R2{6.0};
	complex<Ad> result;

	result = R0 + p*R1 + p*p*R2;
	cerr << "--------------------------------------------------------\n";
	cout << "test expression: R0 + p*R1 + p*p*R2\n";
	cout << "  with p = " << p << ", R0=" << R0 << ", R1 = " << R1 << ", R2 = " << R2 << endl;
	cout << "result: " << result << endl;
	cout << "should be " << R0 + dp*R1 + dp*dp*R2 << endl;
	cout << " derivative should be R1 + 2*R2 = " << R1+2.0*R2 << endl;
	return true;
}

int
main(int argc, char** argv) {
	T_(Trace trc(1,"Ad testing");)
	int nder = 2;  // test with 2 derivative parameters
	ostringstream os;
	bool check_arith = false;
	bool check_fcns = false;
	bool check_cgemv = false;
	bool check_sgemv = false;
	bool check_expression = true;
	double maxnorm{0.0};
	Ad Adsc(1.0);		// before initialize().

	// initial Ad
	cerr << "before initialize(): " << Adsc << endl;

	// set nder parameter names as Ad parameters, then
	// Ad constructors will use the number of parameters
	// to allocate space for derivatives
	vector<string> mypar;
	for (int i=0; i<nder; i++) {
		mypar.push_back(vastr("par", i+1));
	}
	Ad::initialize(mypar);
	cerr << Ad::nder() << " automatic differentiation parameters: "
		<< Ad::toString() << endl;

	// test the assignment op after initialize
	Ad newAdsc = Adsc;
	cerr << "assignment: " << newAdsc << endl;

	try {
		if (check_expression)
			expression_test();

		if (check_arith)
			maxnorm = std::max(maxnorm, test_arith());

		if (check_fcns) {
			// test functions like sin, cos, etc
			// fcn/inverse_fcn and cfcn/cinverse_fcn take one arg
			// fcn2/inverse_fcn2 and cfcn2/cinverse_fcn2 take 2 args
			Ad (*fcn)(const Ad&);
			Ad (*inverse_fcn)(const Ad&);
			Ad (*fcn2)(const Ad&, const Ad&);
			Ad (*inverse_fcn2)(const Ad&, const Ad&);

			fcn = sin;
			inverse_fcn = asin;
			maxnorm = std::max(maxnorm, ad_fcn(string("sin, asin"), fcn, inverse_fcn));

			fcn = cos;
			inverse_fcn = acos;
			maxnorm = std::max(maxnorm, ad_fcn(string("cos, acos"), fcn, inverse_fcn));

			fcn = unary_minus;
			inverse_fcn = unary_minus;
			maxnorm = std::max(maxnorm, ad_fcn(string("unary minus, unary minus"),
					fcn, inverse_fcn));

			fcn = tan;
			inverse_fcn = atan;
			maxnorm = std::max(maxnorm, ad_fcn(string("tan, atan"), fcn, inverse_fcn));

			fcn2 = atan2;
			inverse_fcn2 = tan2;
			maxnorm = std::max(maxnorm, ad_fcn2(string("tan2, atan2"),
					fcn2, inverse_fcn2));

			fcn = exp;
			inverse_fcn = log;
			maxnorm = std::max(maxnorm, ad_fcn(string("exp, log"), fcn, inverse_fcn));

			fcn = sqrt;
			inverse_fcn = square;
			maxnorm = std::max(maxnorm, ad_fcn(string("sqrt, square"), fcn, inverse_fcn));
		}

		if (check_sgemv)
			maxnorm = std::max(maxnorm, sgemv_test());
		if (check_cgemv)
			maxnorm = std::max(maxnorm, cgemv_test());

		cerr << "Largest error: " << maxnorm << "\n";
	} catch (const std::exception& err) {
		cerr << "caught exception: " << err.what() << endl;
	}
}
#endif // MAIN
