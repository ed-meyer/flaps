//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Implementation of class Interval

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <string>
#include <exception>
#include <stdexcept>

#include "config.h"
#include "blas.h"
#include "conv.h"
#include "interval.h"
#include "text.h"
#include "trace.h"

using namespace std;

// should use std::numeric_limits<double>::round_style
void round_up() {}
void round_down() {}
void round_near() {}

constexpr double Pi = flaps::pi;
constexpr double Eps = std::numeric_limits<double>::epsilon();

Interval
Interval::
intersection (const Interval& a) {
	if (disjoint(a))
		throw runtime_error(vastr("no intersection between ",*this, " and ",a));
	return Interval(std::max(a.lower(),lower()), std::min(a.upper(),upper()));
}


bool
Interval::
disjoint (const Interval& a) const {
	double tol = (1.0 + width() + a.width())*Eps;
	return (((l - a.upper()) > tol) || ((a.lower() - u) > tol));
}

Interval
Interval::
Union (const Interval& a) {
	if (disjoint(a))
		throw runtime_error(vastr("no intersection between ",*this, " and ",a));
	return Interval(std::min(a.lower(),lower()), std::max(a.upper(),upper()));
}


double
Interval::
absMin() const {
// XXXX This function should probably be replaced by the more
//      standard mig()
	if (contains(0.0))
		return 0.0;
	else if (lower() > 0.0)
		return lower();
	else
		return -upper();
}

double
Interval::
relDiam() const {
	double d = mag();

		// XXXX note toolbox uses absMin

	if (d > 0.0)
		return width()/d;
	else
		return width();
}

Interval
Interval::
square() const {
	if (contains(0.0)) {
		return Interval(0.0, u*u);
	} else {
		if (u < 0.0)
			return Interval(u*u, l*l);
		else
			return Interval(l*l, u*u);
	}
}

Interval
Interval::
sqrt () const {
// sqrt: throw an exception if interval is < or contains zero
	if (l < 0.0)
		throw runtime_error(vastr("sqrt(",*this,"): invalid interval"));
	else
		return Interval(::sqrt(l), ::sqrt(u));
}

ostream& operator<<(ostream& s, const Interval& a) {
	if (a.lower() == a.upper())
		return s << "[" << a.lower() << "]";
	else if (a.width() < std::max(Eps, Eps*a.mag()))
		return s << '[' << a.midpoint() << 'R' << a.radius() << ']';
	return s << '[' << a.lower() << ':' << a.upper() << ']';
}

// const String&
// Interval::desc() { return *description; }

// const int
// Interval::type() { return Datatype::Interval; }

Interval&
Interval::operator += (const Interval& a) {
	round_down();
	l += a.l;
	round_up();
	u += a.u;
	round_near();
	return *this;
}

Interval&
Interval::operator-=(const Interval& a) {
	round_down();
	l -= a.u;
	round_up();
	u -= a.l;
	round_near();
	return *this;
}

Interval&
Interval::operator*=(const Interval& a) {
// 9 different cases to consider; see R.E. Moore, Methods and
// Applications of Interval Analysis, p. 12, SIAM, 1979
	double tmp;
   if (a.l >= 0.0) {
      if (l >= 0.0) {                // case (1) both positive
			round_down();
         l = a.l*l;
			round_up();
         u = a.u*u;
      } else {
         if (u > 0.0) {             // case (4)
				round_down();
            l = a.u*l;
				round_up();
            u = a.u*u;
         } else {                    // case (6)
				round_down();
            l = a.u*l;
				round_up();
            u = a.l*u;
         }
      }
   } else {
      if (a.u > 0.0) {
         if (l >= 0.0) {          // case (2)
				round_down();
            l = a.l*u;
				round_up();
            u = a.u*u;
         } else {
            if (u <= 0.0) {        // case 7)
					round_down();
               tmp = a.u*l;
					round_up();
               u = a.l*l;
					l = tmp;
						// The case (9) where both contain zero requires 4 mults
            } else {
					round_down();
               double lu = a.l*u;
               double ul = a.u*l;
					round_up();
               double ll = a.l*l;
               double uu = a.u*u;
               l = lu < ul ? lu : ul;
               u = ll > uu ? ll : uu;
            }
         }
      } else {
         if (l >= 0.0) {               // case (3)
				round_down();
            tmp = a.l*u;
				round_up();
            u = a.u*l;
				l = tmp;
         } else {
            if (u > 0.0) {            // case (5)
					round_down();
               tmp = a.l*u;
					round_up();
               u = a.l*l;
					l = tmp;
            } else {                   // case (8)
					round_down();
               tmp = a.u*u;
					round_up();
               u = a.l*l;
					l = tmp;
            }
         }
      }
   }
   return *this;
}

Interval&
Interval::operator/=(const Interval& y) {
// Ref: R. Hammer, M. Hocks, U. Kulisch, D. Ratz,
// C++ Toolbox for Verified Computing, Springer-Verlag, 1995
// p. 33
// 6 cases to consider:
//       x/y              0 < y.l             y.u < 0
//  ---------------------------------------------------------
//     0 <= x.l    |  (x.l/y.u, x.u/y.l)   (x.u/y.u, x.l/y.l)
//  x.l < 0 < x.u  |  (x.l/y.l, x.u/y.l)   (x.u/y.u, x.l/y.u)
//     x.u <= 0    |  (x.l/y.l, x.u/y.u)   (x.u/y.l, x.l/y.u)
//                         1                     4
//                         2                     5
//                         3                     6
	if (y.l > 0.0) {
		if (l >= 0.0) {
			round_down();
			l /= y.u;
			round_up();
			u /= y.l;
		} else if (u > 0.0 && l < 0.0) {
			round_down();
			l /= y.l;
			round_up();
			u /= y.l;
		} else {
			round_down();
			l /= y.l;
			round_up();
			u /= y.u;
		}
	} else if (y.u < 0.0) {
		if (l >= 0.0) {
			double lo;
			round_down();
			lo = u/y.u;
			round_up();
			u = l/y.l;
			l = lo;
		} else if (u > 0.0 && l < 0.0) {
			round_up();
			double up = l/y.u;
			round_down();
			l = u/y.u;
			u = up;
		} else {
			round_up();
			double up = l/y.u;
			round_down();
			l = u/y.l;
			u = up;
		}
	} else
		throw runtime_error(vastr("attempt to divide by ",y));
	round_near();
	return *this;
}


// Equality tests
// ARE THESE NECESSARY since some are also member functions?
// Also - they return int, not bool!

#ifdef NEVER
int operator==(const Interval& a, const Interval& b) {
   return (a.lower() == b.lower() && a.upper() == b.upper());
}

int operator!=(const Interval& a, const Interval& b) {
   return !(a == b);
}

int
operator<(const Interval& a, const Interval& b) {
	return (a.upper() < b.lower());
}

int
operator>(const Interval& a, const Interval& b) {
	return (a.lower() > b.upper());
}
#endif // NEVER

//double
//distance (const Interval& a, const Interval& b) {
//	return Max(fabs(a.lower() - b.lower()), fabs(a.upper() - b.upper()));
//}

// Math functions

double
Interval::
mag() const {
	return std::max(fabs(lower()), fabs(upper()));
}

double
Interval::
mig() const {
	if (lower() < 0.0 && upper() > 0.0) {
		return 0.0;
	} else {
		return std::min(fabs(lower()), fabs(upper()));
	}
}

Interval
cos (const Interval& a) {
	double cl = cos(a.lower());
	double cu = cos(a.upper());
	double sl = sin(a.lower());
	double su = sin(a.upper());
	double l=0.0, u=0.0;

	if (a.width() > 2.0*Pi) {
		return Interval(-1.0, 1.0);
	}
	if (a.width() > Pi) {
		if (sl >= 0.0) {
			l = -1.0;
			u = std::max(cl,cu);
		} else {
			u = 1.0;
			l = std::min(cl, cu);
		}
		return Interval(l,u);
	}
	if (sl*su < 0.0) {
		if (sl >= 0.0) {
			l = -1.0;
			u = std::max(cl,cu);
		} else {
			u = 1.0;
			l = std::min(cl,cu);
		}
		return Interval(l,u);
	} else {
		l = std::min(cl,cu);
		u = std::max(cl,cu);
		return Interval(l,u);
	}
	return Interval(l,u);
}

Interval
acos (const Interval& a) {
	if (a.lower() < -1.0 || a.upper() > 1.0)
		throw runtime_error(vastr("acos(",a,"): invalid interval"));
	return atan(sqrt(1.0 - a.square())/a);
}

Interval
sin (const Interval& a) {
	Interval b(a.lower()-Pi/2.0, a.upper()-Pi/2.0);
	return cos(b);
}

Interval
asin (const Interval& a) {
	Interval b(a.lower()-Pi/2.0, a.upper()-Pi/2.0);
	return cos(b);
}

Interval
tan (const Interval& a) {
	return sin(a)/cos(a);
}

Interval
atan (const Interval& a) {
	return Interval(atan(a.lower()), atan(a.upper()));
}

Interval
exp (const Interval& a) {
	return Interval(exp(a.lower()), exp(a.upper()));
}

Interval
sqrt (const Interval& a) {
// sqrt: throw an exception if "a" is < or contains zero
	if (a.lower() < 0.0)
		throw runtime_error(vastr("sqrt(",a,"): invalid interval"));
	else
		return Interval(::sqrt(a.lower()), ::sqrt(a.upper()));
}

Interval
pow (const Interval& x, double y) {
   double a = ::pow(double(lower(x)), double(y));
   double b = ::pow(double(upper(x)), double(y));
		// Let the constructor sort out min & max...
   return Interval (a, b);
}

Interval
pow (const Interval& x, const Interval& y) {
	double ll, lu, ul, uu, lo, hi;

	ll = ::pow (double(lower(x)), double(lower(y)));
	lu = ::pow (double(lower(x)), double(upper(y)));
	ul = ::pow (double(upper(x)), double(lower(y)));
	uu = ::pow (double(upper(x)), double(upper(y)));
	lo = std::min(ll, std::min(lu, std::min(ul, uu)));
	hi = std::max(ll, std::max(lu, std::max(ul, uu)));
   return Interval (lo, hi);
}

Interval
hypot (const Interval& a, const Interval& b) {
	return sqrt(a.square() + b.square());
}

Interval
safeDivide(const Interval& a, const Interval& b) {
/*------------------------------------------------------------------
 * Divide two Intervals "safely"
 *------------------------------------------------------------------*/
	double den = b.mig();
	double num = a.mig();
	static double tol = 0.0;
	static double big = 0.0;
	Interval rval(0.0);

	if (tol == 0.0) {
		tol = sqrt(sqrt(Eps));
		big = 1.0/tol;
	}
	if (abs(den) > abs(num*tol)) {
		return a/b;
	} else if (den >= 0.0) {
		if (num > 0.0)
			rval = Interval(big);
		else if (num == 0.0)
			rval = Interval(0.0);
		else
			rval = Interval(-big);
	} else {
		if (num > 0.0)
			rval = Interval(-big);
		else if (num == 0.0)
			rval = Interval(0.0);
		else
			rval = Interval(big);
	}
	return rval;
}

//------------------------------------------------------------------
// Complex Interval stuff
//------------------------------------------------------------------
// ostream& operator<<(ostream& s, const CInterval& a) {
// 	return s << "(" << a.real() << ", " << a.imag() << ")";
// }

//------------------------------------------------------------------
// Extended Interval stuff
//------------------------------------------------------------------

XInterval::
XInterval (const Interval& a) {
	low = a;
	isSplit = false;
	minusInf = false;
	plusInf = false;
}

XInterval::
XInterval (const Interval& a, const Interval& b) {
	low = a;
	hi = b;
	isSplit = true;
	minusInf = false;
	plusInf = false;
}

XInterval
XInterval::
intersection (const Interval& a) const {
	bool emptyLower, emptyUpper;
	double lb, ub;
	Interval lowerInt, upperInt;

	if (a.lower() < low.upper()) {
		lb = std::max(a.lower(), low.lower());
		ub = std::min(a.upper(), low.upper());
		lowerInt = Interval(lb,ub);
		emptyLower = false;
	} else {
		emptyLower = true;
	}

	// check upper interval if there is one...

	if (isSplit) {
		if (a.upper() > hi.lower()) {
			ub = std::min(a.upper(), hi.upper());
			lb = std::max(a.lower(), hi.lower());
			upperInt = Interval(lb, ub);
			emptyUpper = false;
		} else {
			emptyUpper = true;
		}
	} else {
		emptyUpper = true;
	}

	if (emptyLower && emptyUpper)
		throw runtime_error(vastr("no intersection between ",*this," and ",a));

	if (!emptyLower && !emptyUpper)
		return XInterval(lowerInt, upperInt);
	
	if (!emptyLower)
		return XInterval(lowerInt);
	else
		return XInterval(upperInt);
}

XInterval
xdivide (const Interval& x, const Interval& y) {
// Extended interval division: x/y where y.contains(0)
// Ref: C++ Toolbox for Verified Computing, pg 40
	XInterval rval;
	double inf = XInterval::infinity();

	// If divisor y does not contain zero return vanilla interval division
	if (!y.contains(0.0)) {
		rval.low = x/y;
		rval.isSplit = false;
		rval.minusInf = rval.plusInf = false;
		return rval;
	}
	// If dividend (numerator) contains zero or the divisor (denominator) y
	// is a thin zero the result is +- infinity
	if (x.contains(0.0) || (y.lower() == 0.0 && y.upper() == 0.0)) {
		rval.low = Interval(-inf, inf);
		rval.minusInf = rval.plusInf = true;
		return rval;
	}
	// If a limit of y is zero the result is a single interval with infinity...
	if (y.upper() == 0.0) {
		if (x.upper() <= 0.0) {
			rval.plusInf = true;
			rval.low = Interval(x.upper()/y.lower(), inf);
			rval.isSplit = false;
			return rval;
		} else if (x.lower() >= 0.0) {
			rval.minusInf = true;
			rval.low = Interval(-inf, x.lower()/y.lower());
			rval.isSplit = false;
			return rval;
		}
	}

	if (y.lower() == 0.0) {
		if (x.upper() <= 0.0) {
			rval.minusInf = true;
			rval.low = Interval(-inf, x.upper()/y.upper());
			return rval;
		} else if (x.lower() >= 0.0) {
			rval.plusInf = true;
			rval.low = Interval(x.lower()/y.upper(), inf);
			return rval;
		}
	}

	// dividend (numerator) x does not contain zero: two
	// intervals result, "hi" and "low"
	rval.isSplit = true;
	rval.minusInf = rval.plusInf = true;
	if (x.upper() <= 0.0) {
		rval.low = Interval(-inf, x.upper()/y.upper());
		rval.hi = Interval(x.upper()/y.lower(), inf);
	} else {
		rval.low = Interval(-inf, x.lower()/y.lower());
		rval.hi = Interval(x.lower()/y.upper(), inf);
	}
	return rval;
}

XInterval
xsubtract (double x, const XInterval& y) {
// Extended interval subtraction: x-y
// Ref: C++ Toolbox for Verified Computing, pg 40
	XInterval rval;
	double oo = XInterval::infinity();

	rval.minusInf = true;
	rval.plusInf = true;

	if (y.isSplit) {
		rval.hi = Interval(x-y.low.upper(), oo);
		rval.low = Interval(-oo, x - y.hi.lower());
		rval.isSplit = true;
	} else {
		if (y.minusInf) {
			if (y.plusInf) {
				rval.low = Interval(-oo, oo);
			} else {
				rval.low = Interval(x - y.low.upper(), oo);
				rval.minusInf = false;
			}
		} else {
			if (y.plusInf) {
				rval.low = Interval(-oo, x - y.low.lower());
				rval.plusInf = false;
			} else {
				rval.low = Interval(x) - y.low;
				rval.minusInf = false;
				rval.plusInf = false;
			}
		}
	}
	return rval;
}

ostream& operator<<(ostream& s, const XInterval& a) {
	if (a.split()) {
		if (a.hasNegInf() && a.hasPosInf()) {
			return s << "[-inf, " << a.lowerInterval().upper() <<
				"]U[" << a.upperInterval().lower() << ", inf]";
		} else if (!a.hasNegInf() && a.hasPosInf()) {
			return s << "[" << a.lowerInterval().lower() << ", " << a.lowerInterval().upper() <<
				"]U[" << a.upperInterval().lower() << ", inf]";
		} else if (a.hasNegInf() && !a.hasPosInf()) {
			return s << "[-inf, " << a.lowerInterval().upper() <<
				"]U[" << a.upperInterval().lower() << ", " << a.upperInterval().upper() << "]";
		} else if (!a.hasNegInf() && !a.hasPosInf()) {
			return s << "[" << a.lowerInterval().lower() << ", " << a.lowerInterval().upper() <<
				"]U[" << a.upperInterval().lower() << ", " << a.upperInterval().upper() << "]";
		}
	} else {
		if (a.hasNegInf() && a.hasPosInf()) {
			return s << "[-inf, inf]";
		} else if (!a.hasNegInf() && a.hasPosInf()) {
			return s << "[" << a.lowerInterval().lower() << ", inf]";
		} else if (!a.hasNegInf() && a.hasPosInf()) {
			return s << "[-inf, " << a.lowerInterval().upper() << "]";
		}
	}
	return s;	// hush -Wall
}
			

#ifdef MAIN
#undef MAIN

using std::endl;
using std::cout;

int
main (int argc, char* argv[]) {
	const int n = 50000;
	Interval a[n];
	Interval b = 1.0, d(2.0,3.0);
	double t;
	int i;

	cout << "\n\nTesting Interval " << sizeof(Interval) <<
		" bytes" << endl;
	cout << "test of operator<<: a[0] = " << a[0] << ", b = " << b
		<< ", d = " << d << endl;
	cout << "length of " << n << " Intervals = " << sizeof(a) << endl;
	cout << "b*d = " << b*d << endl;
	cout << "b/d = " << b/d << endl;
	t = 1.0/sqrt((double)n);
	for (i = 0; i<n; i++)
		a[i] = Interval(t);
	b = Interval(0.0, 0.0);
	for (i = 0; i<n; i++)
		b += a[i]*a[i];
	cout << "dot-product of 2 " << n << "-vectors of " <<
		a[0] << " = " << b << endl;

	b = Interval(2.0);
	d = Interval(3.0-0.01, 3.0+0.01);
	cout << b << " ^ " << d << " = " << pow(b,d) << endl;
	if (b < d)
		cout << b << " < " << d << endl;
	else
		cout << b << " > " << d << endl;

	b = Interval(-1.0, 2.0);
	cout << "abs (" << b << ") = " << abs(b) << endl;

	cout << "exp(-1,2) = " << exp(b) << endl;

	{
		Interval a(33.17,33.43);
		Interval b(-1.0,-0.99);
		Interval c(0.0,0.0);

		c = a*b;
		cout << "c = " << a << " * " << b << " = " << c << endl;
	}

	{
		Interval x(4.0, 5.0), y(-1.0, 2.0);
		XInterval xi = xdivide(x, y);
		cout << "XInterval divide: " << x << "/" << y << " = " << xi << endl;
	}
	{
		Interval a[6];
		Interval tmp(0.0);
		int i;
		for (i=0; i<6; i++)
			a[i] = Interval(-1.0, 1.0);
		a[0] = Interval(0.0, 1.0);
		a[1] = Interval(0.0);
		for (i=0; i<6; i++)
			tmp += a[i].square();
		cout << "a = " << a[2] << ", a*a = " << a[2].square() << endl;
		cout << "a*a = " << tmp << endl;
	}
}
#endif	// MAIN
