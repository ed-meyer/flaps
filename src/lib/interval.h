//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//------------------------------------------------------------------
// References:
//  [1] Moore, R.E., Methods and Applications of Interval Analysis,
//      SIAM, 1979
//  [2] Alfeld, G. and Herzberger, J., Introduction to Interval
//      Analysis, Academic Press, 1983
//  [3] Neumaier, A., Interval Methods for Systems of Equations,
//      Cambridge Univ. Press, 1990
//  [4] Hammer, R., Hocks, M., Kulisch, U., and Ratz, D.,
//      C++ Toolbox for Verified Computing I, Springer, 1995
//  [5] Kearfott, R. Baker, Rigorous Global Search: Continuous
//      Problems, Kluwer Academic Publishers, 1996
//  [6] Aberth, O., Precise Numerical Methods Using C++,
//      Academic Press, 1998.
//  [7] Hamzo, Chadi and Kreinovich, Vladik, "On Average Bit
//      Complexity of Interval Arithmetic", Bulletin of the European
//      Association for Theoretical Computer Science (EATCS), 1999,
//      Vol. 68, pp. 153-156.
//      (http://www.cs.utep.edu/vladik/1999/tr99-20.ps.gz)
//      Show how at most 3 multiplications are necessary for
//      an interval multiplication.
//------------------------------------------------------------------

#ifndef Interval_h
#define Interval_h
 
#include <algorithm>
#include <iostream>
#include <cmath>
#include <string>

#include "message.h"

void round_up();
void round_down();
void round_near();

class Interval {
public:
	// static const String& desc();
	// static const int type();

		// constructors:
		//  - let copy constructor default to bitwise
		//  - must include a single-argument constructor because
		//    the default sup in this case is inf, not 0.0
	Interval() : l(0.0), u(0.0) { }		// no initializer
	Interval(double m) : l(m), u(m) { }	// single double
	Interval(double inf, double sup) {
		l = std::min(inf,sup);
		u = std::max(inf,sup);
	}
	Interval(const Interval& x) { l = x.l; u = x.u; }
	~Interval() {}
	// default memberwise assignment is ok
	Interval& operator=(double x) { l = u = x; return *this; }

	bool innerIncludes (const Interval& a) const {
		if (a.lower() == a.upper())
			return true;
		if (l == a.lower() && u > a.upper())
			return true;
		if (l < a.lower() && u == a.upper())
			return true;
		return l < a.lower() && u > a.upper(); }
	bool innerIncludes (double a) const { return l < a && u > a; }
	bool contains (const Interval& a) const { return l <= a.lower() && u >= a.upper(); }
	bool contains (double a) const { return l <= a && u >= a; }
	bool intersects (const Interval& a) const {
		return contains(a.lower()) || contains(a.upper()); }
	bool disjoint (const Interval& a) const;

	Interval intersection (const Interval& a);
	Interval Union (const Interval& a);
	double upper() const { return u; }
	double lower() const { return l; }
	double midpoint() const { return ((u + l)/2.0); }
	double relDiam() const;
	double absMin() const;
	double width() const { return std::abs(u-l); }
	double radius() const { return (std::abs(u - l)/2.0); }
	double mig() const;
	double mag() const;
	// NOTE: abs() is usually ([1],[2], and [3]) defined
	// as max(abs(u),abs(l)), but [4] defines it as an interval!
	// I would like to return to the more common definition
	// However - the STL defines template<FLT> abs(complex<FLT>& a)
	// ie it returns the datatype of its argument; to be consistent
	// with the STL lets define abs to return Interval
	Interval abs() const { return Interval(mig(), mag()); }

	Interval square() const;
	Interval sqrt() const;

	// These are so we can write, e.g. midpoint(a) instead of a.midpoint()
	friend double midpoint(const Interval&);
	friend double width(const Interval&);
	friend double radius(const Interval&);
	friend double upper(const Interval&);
	friend double lower(const Interval&);
	friend double mig(const Interval&);
	friend double mag(const Interval&);

	friend Interval safeDivide(const Interval& a, const Interval& b);

	Interval& operator += (const Interval&);
	Interval& operator += (double b) { return operator+=(Interval(b)); }
	Interval& operator -= (const Interval&);
	Interval& operator -= (double b) { return operator-=(Interval(b)); }
	Interval& operator *= (const Interval&);
	Interval& operator *= (double b) { return operator*=(Interval(b)); }
	Interval& operator /= (const Interval&);
	Interval& operator /= (double b) { return operator/=(Interval(b)); }

	// Unary minus XXX was wrong: modified *this
 	Interval operator-() { return Interval(-u,-l); }

	bool operator==(const Interval& b) const {
		return lower() == b.lower() && upper() == b.upper(); }
	bool operator==(double b) const { return lower() == b && upper() == b; }

	bool operator<=(const Interval& b) const {
		return lower() <= b.lower() && upper() <= b.upper(); }
	bool operator<=(double b) const { return lower() <= b && upper() <= b; }

	bool operator>=(const Interval& b) const {
		return lower() >= b.lower() && upper() >= b.upper(); }
	bool operator>=(double b) const { return lower() >= b && upper() >= b; }

	bool operator<(const Interval& b) const { return (upper() < b.lower()); }
	bool operator<(double b) const { return (upper() < b); }

	bool operator>(const Interval& b) const { return (lower() > b.upper()); }
	bool operator>(double b) const { return (lower() > b); }


#ifdef NEVER
        friend  double  norm(Interval);
        friend  Interval cosh(Interval);
        friend  Interval log(Interval);
        friend  complex pow(complex, int);
        friend  complex pow(complex, double);
        friend  complex pow(complex, complex);
        friend  complex polar(double, double = 0);
        friend  complex sin(complex);
        friend  complex sinh(complex);
        friend  complex sqrt(complex);
#endif
	
private:
	double l, u;
	// static String* description;
			
};			// end of class Interval

//------------------------------------------------------------------
// Operators

inline Interval operator+(const Interval& a, const Interval& b) {
	Interval r = a;
	r += b;
	return r;
}
inline Interval operator+(double a, const Interval& b) {
	Interval r(a);
	r += b;
	return r;
}
inline Interval operator+(const Interval& a, double b) {
	Interval r = a;
	r += Interval(b);
	return r;
}

inline Interval operator-(const Interval& a, const Interval& b) {
	Interval r = a;
	r -= b;
	return r;
}
inline Interval operator-(double a, const Interval& b) {
	Interval r(a);
	r -= b;
	return r;
}
inline Interval operator-(const Interval& a, double b) {
	Interval r = a;
	r -= Interval(b);
	return r;
}

inline Interval operator*(const Interval& a, const Interval& b) {
	Interval r = a;
	r *= b;
	return r;
}
inline Interval operator*(double a, const Interval& b) {
	Interval r(a);
	r *= b;
	return r;
}
inline Interval operator*(const Interval& a, double b) {
	Interval r = a;
	r *= Interval(b);
	return r;
}

inline Interval operator/(const Interval& a, const Interval& b) {
	Interval r = a;
	r /= b;
	return r;
}
inline Interval operator/(double a, const Interval& b) {
	Interval r(a);
	r /= b;
	return r;
}
inline Interval operator/(const Interval& a, double b) {
	Interval r = a;
	r /= Interval(b);
	return r;
}

#ifdef NEVER
int     operator==(const Interval&, const Interval&);
int     operator>(const Interval&, const Interval&);
int     operator<(const Interval&, const Interval&);
int     operator!=(const Interval&, const Interval&);
#endif // NEVER
// double distance (const Interval& a, const Interval& b);

//------------------------------------------------------------------
Interval cos(const Interval& a);
Interval acos(const Interval& a);
Interval sin (const Interval& a);
Interval asin (const Interval& a);
Interval tan (const Interval& a);
Interval atan (const Interval& a);
inline Interval atan2 (const Interval& x, const Interval& y) { return Interval(atan2(x.lower(),y.upper()), atan2(x.upper(),y.lower())); }
Interval exp(const Interval&);
Interval sqrt(const Interval&);
Interval pow (const Interval&, double);
Interval pow (const Interval&, const Interval&);

Interval hypot (const Interval&, const Interval&);

// Hyperbolic functions
//  XXXX These need some arg checking!!!

inline Interval sinh (const Interval& x) { return Interval(sinh(x.lower()), sinh(x.upper())); }
inline Interval cosh (const Interval& x) {
	if (x.upper() < 0.0) {
		return Interval(cosh(x.upper()), cosh(x.lower()));
	} else if (x.lower() > 0.0) {
		return Interval(cosh(x.lower()), cosh(x.upper()));
	} else {
		if (-x.lower() > x.upper())
			return Interval(0.0, cosh(x.lower()));
		else
			return Interval(0.0, cosh(x.upper()));
	}
}
inline Interval tanh (const Interval& x) {return Interval(tanh(x.lower()), tanh(x.upper()));}

inline Interval log (const Interval& x) { return Interval(log(x.lower()), log(x.upper())); }



//------------------------------------------------------------------
// Inline functions

inline double midpoint(const Interval& a)  { return (a.upper() + a.lower())*double(0.5); }
inline double width(const Interval& a) { return fabs(a.u - a.l); }
inline double radius(const Interval& a) { round_up();return fabs(a.u - a.l)*0.5; }
inline double upper(const Interval& a) { return a.u; }
inline double lower(const Interval& a) { return a.l; }

inline double mig(const Interval& x) { return x.mig(); }

namespace ax {
inline double mag(const Interval& x) { return x.mag(); }
}

//------------------------------------------------------------------
// mag always returns a double, whereas abs returns the same
// datatype (AR in this case)
// See comment above in the class defn regarding abs()
//------------------------------------------------------------------

inline Interval
abs(const Interval& a) { return Interval(a.mig(), a.mag()); }

Interval safeDivide(const Interval& a, const Interval& b);

std::ostream& operator<<(std::ostream& s, const Interval& a);

//------------------------------------------------------------------
// XInterval: extended interval
//------------------------------------------------------------------

class XInterval {
	Interval low, hi;
	bool minusInf;
	bool plusInf;
	bool isSplit;
public:

	XInterval() : low(0.0), hi(0.0), minusInf(false), plusInf(false), isSplit(false) {}
	XInterval(const Interval& a);
	XInterval(const Interval& a, const Interval& b);

	bool split() const { return isSplit; }
	bool hasNegInf() const { return minusInf; }
	bool hasPosInf() const { return plusInf; }
	XInterval intersection (const Interval& a) const;
	bool disjoint(const Interval& a) const {
		bool rval = low.disjoint(a);
		if (isSplit)
			rval = rval && hi.disjoint(a);
		return rval;
	}
	Interval lowerInterval() const { return low; }
	Interval upperInterval() const { return hi; }
	static double infinity() { return std::numeric_limits<double>::max(); }
	friend XInterval xdivide(const Interval& a, const Interval& b);
	friend XInterval xsubtract(double a, const XInterval& b);
};

std::ostream& operator<<(std::ostream& s, const XInterval& a);

XInterval xdivide (const Interval& x, const Interval& y);

#endif

