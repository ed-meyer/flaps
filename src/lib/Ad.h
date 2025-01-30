//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// classes for automatic differentiation (AD)
//   Ad::initialize(const vector<string>& names)
//      this is the *only* way to set the names of parameters to be
//      used as autodiff derivatives and it may be called only once
//      It will also reallocate the Ad members of each parameter in the
//      gpset, and set parameter derivatives to 1 in the Ad member.
//   set_adpar_derivs(pset& pl)
//      for each adparname, find the parameter in pl and set the
//      corresponding AD derivative to 1.0
//   std::vector<std::string> const& adnames()
//      returns a reference to a vector of the Ad derivative parameter names
// namespace ad comprises functions for maintaining a list
// of parameter names for the automatic-differentiation class Ad



#ifndef Ad_h
#define Ad_h

#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <exception>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

class Ad {
public:
	// std::unique_ptr<double> data_{nullptr};
	double* data_{nullptr};

	// constructors
	// default
	Ad() {
		data_ = new double[Ad::ndata()];
		std::fill_n(data_, Ad::ndata(), 0.0);
	}
	// a double
	Ad(double a) {
		data_ = new double[Ad::ndata()];
		std::fill_n(data_, Ad::ndata(), 0.0);
		data_[0] = a;
	}

	// copy constructor
	Ad (Ad const& from) {
		data_ = new double[Ad::ndata()];
		std::copy_n (from.data_, Ad::ndata(), data_);
	}
	// move copy constructor
	Ad (Ad&& from) noexcept {
		data_ = from.data_;
		from.data_ = nullptr;
	}
	// assignment operator
	Ad& operator=(Ad const& rhs) {
		std::copy_n (rhs.data_, Ad::ndata(), data_);
		return *this;
	}
	// move assignment operator
	Ad& operator=(Ad&& rhs) noexcept {
		std::swap (data_, rhs.data_);
		return *this;
	}

	// destructor
	~Ad () {
		if (data_ != nullptr)
			delete[] data_;
		data_=nullptr;
	}

	// member access
	// value: return or set just the value, not derivatives
	double value() const { return this->data_[0]; }
	void value(double v) { data_[0] = v; }

	// reference to the real or imag parts - not avail in std::complex
	// cplusplus.github.io/LWG/issue 387: complex over-encapsulated
	// XXX do this also for complex<double>?
	static Ad& real(std::complex<Ad>& x) {
		return reinterpret_cast<Ad(&)[2]>(x)[0];
	}
	static const Ad& real(std::complex<Ad> const& x) {
		return reinterpret_cast<const Ad(&)[2]>(x)[0];
	}
	static Ad& imag(std::complex<Ad>& x) {
		return reinterpret_cast<Ad(&)[2]>(x)[1];
	}
	static const Ad& imag(std::complex<Ad> const& x) {
		return reinterpret_cast<const Ad(&)[2]>(x)[1];
	}
	// return the "magnitude" of a complex<Ad>: a double
	static double mag(std::complex<Ad> const& x);

	// zero values & derivatives of this Ad
	void zero() { std::fill_n(data_, Ad::ndata(), 0.0); }

	// get/set derivative "der" 0b
	double der(int d) const; // { return data_[der+1]; }
	void der(int d, double v); // { data_[der+1] = v; }

	// get a pointer to an Ad data array: 1 + nder doubles
	double* data() { return data_; }
	const double* data() const { return data_; }

	// reallocate data_ to reflect the number of derivative parameters
	void realloc();

	// Assignment
	Ad& operator=(double rhs) {
		for (int i=0; i<Ad::ndata(); i++)
			data_[i] = 0.0;
		data_[0] = rhs;
		return *this;
	}

	//------------------------------------------------------------------
	// arithmetic assignment operators
	// addition
	inline Ad& operator+= (const Ad& rhs) {
		for (int i=0; i<Ad::ndata(); i++)
			data_[i] += rhs.data_[i];
		return *this;
	}
	inline Ad& operator+= (double rhs) {
		data_[0] += rhs;
		return *this;
	}

	// subtraction
	inline Ad& operator-= (const Ad& rhs) {
		for (int i=0; i<Ad::ndata(); i++)
			data_[i] -= rhs.data_[i];
		return *this;
	}
	inline Ad& operator-= (double rhs) {
		data_[0] -= rhs;
		return *this;
	}

	// multiplication a = a*b, a' = a'*b + a*b'
	inline Ad& operator*= (const Ad& rhs) {
		double a{data_[0]};
		double b{rhs.data_[0]};
		for (int i=1; i<Ad::ndata(); i++)
			data_[i] *= b;
		for (int i=1; i<Ad::ndata(); i++)
			data_[i] += a*rhs.data_[i];
		data_[0] = a*b;
		return *this;
	}
	inline Ad& operator*= (double rhs) {
		for (int i=0; i<Ad::ndata(); i++)
			data_[i] *= rhs;
		return *this;
	}

	// division: c = a/b
	//           c' = a'/b - (a/b^2)b' = (a'b - ab')/b^2
	inline Ad& operator/=(const Ad& a) {
		double a2 = a.data_[0];
		double r = this->data_[0]/a2;
		this->data_[0] = r;
		for (int i=1; i<Ad::ndata(); i++)
			this->data_[i] = (this->data_[i] - r*a.data_[i])/a2;
		return *this;
	}

	inline Ad& operator/=(double a) {
		for (int i=0; i<Ad::ndata(); i++)
			this->data_[i] /= a;
		return *this;
	}

	// end of arithmetic assignment operators
	//------------------------------------------------------------------
	// Relational operators
	// XXX should these also check derivatives? or should there be
	// a special function to do that?
	// provide is_equal to test with a tolerance?

	bool
	operator== (const Ad& b) const { return this->value() == b.value(); }

	bool
	operator!= (const Ad& b) const { return !this->operator==(b); }

	bool
	operator> (const Ad& b) const {
		return this->value() > b.value();
	}

	bool
	operator< (const Ad& b) const {
		return this->value() < b.value();
	}

	bool
	operator<= (const Ad& b) const {
		return this->value() <= b.value();
	}

	bool
	operator>= (const Ad& b) const {
		return this->value() >= b.value();
	}
	//------------------------------------------------------------------
	// Unary operators

	Ad
	operator+() // non-const
		{ return *this; }

	Ad
	operator-() const {
		Ad rval(*this);
		double* rp = rval.data();
		for (int i=0; i<Ad::ndata(); i++)
			rp[i] = -rp[i];
		return rval;
	}
	//------------------------------------------------------------------
	// this is the *only* way to set the names of parameters to be
	// used as autodiff derivatives and it may be called only once
	// It will also reallocate the Ad members of each parameter in the
	// stdpl, and set parameter derivatives to 1 in the Ad member.
	static bool initialize_called;
	static void initialize(const std::vector<std::string>& names);

	static int nder();

	static int ndata();

	// returns a const reference to the vector of AD parameter names
	static std::vector<std::string> const& adnames();

	static std::string toString();

	// return the 0b index into data_ of parameter "name" or -1
	static int find(std::string const& name);

	static std::string live();
	static int nlive();

	// multiplication c = a*b for Ad and complex<Ad>
   static void multad (const Ad& a, const Ad& b, Ad& c);
   static void multcad(const std::complex<Ad>& a, const std::complex<Ad>& b,
         std::complex<Ad>& c, Ad& wk);

}; // class Ad ------------------------------------------------------------------

// Specializations for complex<Ad> to get around the stdlib real() and
// imag() member fcns which return an Ad value instead of a reference
// implementations in Ad.cpp
template<>
template<>
std::complex<Ad>::
complex<Ad>(const std::complex<Ad>& z);

template<>
template<>
std::complex<Ad>&
std::complex<Ad>::operator*=(const std::complex<Ad>& rhs);

template<>
template<>
std::complex<Ad>&
std::complex<Ad>::operator+=(const std::complex<Ad>& rhs);


// Non-member functions

bool
is_equal(Ad const& a, Ad const& b, int sigfig);

// adfcn1p and adfcn2p are pointers to functions taking 1
// or 2 arguments of type Ad; these are mainly for adeqn
using adfcn1p = Ad (*)(const Ad&);
adfcn1p
adfcn1(const std::string& fcn);

using adfcn2p = Ad (*)(const Ad&, const Ad&);
adfcn2p
adfcn2(const std::string& fcn);

std::ostream& operator<<(std::ostream& s, const Ad& a);

// arithmetic operators: create an Ad, call arithmetic assignment op
inline
Ad
operator+(const Ad& a, const Ad& b) {
	Ad c{a};
	return c.operator+=(b);
}

inline
Ad
operator+(const Ad& a, double b) {
	Ad c{a};
	return c.operator+=(b);
}

inline
Ad
operator+(double a, const Ad& b) {
// commutative
	return operator+(b,a);
}

// XXX complex arith ops are necessary - get 'no match for operator+'
// need to do the same for all arith ops
inline
std::complex<Ad>
operator+(const std::complex<Ad>& a, double b) {
	Ad cr = a.real()+b;
	Ad ci = a.imag();
	return std::complex<Ad>(cr,ci);
}
inline
std::complex<Ad>
operator+(double a, const std::complex<Ad>& b) {
// commutative
	return operator+(b,a);
}

inline
Ad
operator-(const Ad& a, const Ad& b) {
	Ad c{a};
	return c.operator-=(b);
}

inline
Ad
operator-(const Ad& a, double b) {
	Ad c{a};
	return c.operator-=(b);
}

inline
Ad
operator-(double a, const Ad& b) {
// non-commutative
	Ad c{a};
	return c.operator-=(b);
}

inline
std::complex<Ad>
operator-(const std::complex<Ad>& a, double b) {
	Ad cr = a.real()-b;
	Ad ci = a.imag();
	return std::complex<Ad>(cr,ci);
}
inline
std::complex<Ad>
operator-(double a, const std::complex<Ad>& b) {
// non-commutative
	Ad cr(a);
	std::complex<Ad> c(cr);
	return operator-(c,b);
}

inline
Ad
operator*(const Ad& a, const Ad& b) {
	Ad c{a};
	return c.operator*=(b);
}

inline
Ad
operator*(const Ad& a, double b) {
// commutative
	return operator*(b,a);
}

inline
Ad
operator*(double a, const Ad& b) {
// commutative
	Ad c{b};
	return c.operator*=(a);
}

inline
std::complex<Ad>
operator*(const std::complex<Ad>& a, double b) {
	Ad ar = a.real();
	Ad ai = a.imag();
	ar.operator*=(b);
	ai.operator*=(b);
	return std::complex<Ad>(ar,ai);
	// std::complex<Ad> c{a};
	// c.real() *= b;
	// c.imag() *= b;
	// return c;
}

inline
std::complex<Ad>
operator*(double a, const std::complex<Ad>& b) {
// commutative
	return operator*(b,a);
}

// division: non-commutative
inline
Ad
operator/(const Ad& a, const Ad& b) {
	Ad c{a};
	return c.operator/=(b);
}

inline
Ad
operator/(const Ad& a, double b) {
	Ad c{a};
	return c.operator/=(b);
}

inline
Ad
operator/(double a, const Ad& b) {
// non-commutative
	Ad c(a);
	return c.operator/=(b);
}


// math function declarations

// single ad argument fcns
Ad
abs (Ad const& x);

Ad
acos (Ad const& x);

Ad
asin (Ad const& x);

Ad
atan (Ad const& x);

Ad
cos (Ad const& x);

Ad
exp (Ad const& x);

Ad
log (Ad const& x);

Ad
log10 (Ad const& x);

Ad
sin (Ad const& x);

Ad
sqrt (Ad const& x);

Ad
tan (Ad const& x);

// 2-arg Ad functions
Ad
atan2 (Ad const& x, const Ad& y);

Ad
max (Ad const& x, const Ad& y);

Ad
min (Ad const& x, const Ad& y);

Ad
pow (Ad const& x, const Ad& y);

// extract a named derivative from arrays of Ad or complex<Ad>
bool extract (const std::vector<Ad>& from, std::string const& name, double* val);
bool extract (const std::vector<Ad>& from, int der, double* val);
bool extract (const std::vector<std::complex<Ad>>& from,
		std::string const& name, std::complex<double>* val);
bool extract (const std::vector<std::complex<Ad>>& from,
		int der, std::complex<double>* val);

#endif // Ad_h
