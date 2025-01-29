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
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "fptype.h"
#include "trace.h"

using namespace std;

std::string
toString(double const& a) {
// Returns a string representation of a double with near max precision
	std::ostringstream os;
	os << setprecision(18) << a;
	return os.str();
}

double
flaps::
xrand() {
// Returns a random number between -1 and 1
	int k = rand();
	double rval;

	rval = (double)k/(double)RAND_MAX;
	rval = 2.0*rval - 1.0;
	return rval;
}

bool
flaps::
is_valid(double const& t) {
// is "t" a valid double, i.e. check for NAN, infinity
	if (std::numeric_limits<double>::has_infinity && t == std::numeric_limits<double>::infinity())
		return false;
	// XXX the IEEE std says that comparing nan's will always fail!!
	// if (std::numeric_limits<double>::has_quiet_NaN && t == std::numeric_limits<double>::quiet_NaN())
	if (isnan(t))
		return false;
	// if (std::numeric_limits<double>::has_signaling_NaN && t == std::numeric_limits<double>::signaling_NaN())
	// 	return false;
	return true;
}

double
safe_divide (double num, double den) {
/*------------------------------------------------------------------
 * Divide two doubles "safely"
 *------------------------------------------------------------------*/
	static double eps = std::numeric_limits<double>::epsilon();
	static double big = 1.0/eps;
	double rval = 0.0;
	// test to make sure den > eps
	if (std::abs(den) > std::abs(num*eps) && std::abs(den) > eps) {
		rval = num/den;
	} else {
		if (num == 0.0) {
			rval = 0.0;
		} else {
			if (den < 0.0) {
				if (num > 0.0)
					rval = -big;
				else
					rval = big;
			} else {
				if (num > 0.0)
					rval = big;
				else
					rval = -big;
			}
		}
	}
	return rval;
}


double
ulp_frac (double a, double b) {
/*
 * Return the fraction of binary digits in the mantissa which
 * agree between two doubles. Thus if 1.0 is returned they agree
 * exactly, and if 0 is returned they do not agree at all.
 */
	double t = 2.0;		// 2^ulp
	double epsilon = std::numeric_limits<double>::epsilon();
	double eps = epsilon;
	int ulp, maxulp;
	double mantissa = std::numeric_limits<double>::digits; // digits in the mantissa

	maxulp = (int)mantissa;
	for (ulp=1; ulp<maxulp; ulp++) {
		if (fcmp(a, b, eps) == 0) {
			break;
		}
		t *= 2.0;
		eps = (t - 1.0)*epsilon;
	}
	if (ulp == maxulp)
		return 0.0;
	else
		return (double)(maxulp-ulp)/mantissa;
}

int
ulp_double (double a, double b) {
/*
 * Return the number of units in the last place difference between two
 * doubles; that is, the number of binary digits in the mantissa in the
 * difference between two doubles. If two doubles differ by only 1-3 ulp
 * they can be considered equal
 */
	double t = 2.0;		// 2^ulp
	double epsilon = std::numeric_limits<double>::epsilon();
	double eps = epsilon;
	int ulp;
	int maxulp = 40; // no point going further than this

	for (ulp=1; ulp<maxulp; ulp++) {
		if (fcmp(a, b, eps) == 0) {
			return ulp;
		}
		t *= 2.0;
		eps = (t - 1.0)*epsilon;
	}
	return ulp;
}

bool
is_equal (double const& a, double const& b, int ndigits) {
/*
 * Are a and b equal to within "ndigits" significant figures?
 * Note: this function is different from fcmp in a very important
 * way. fcmp always does a comparison similar to
 *       a - b < eps*max(abs(a), abs(b))
 * so that if a and b are tiny the comparison could fail. For example
 *       fcmp(1.e-19, 0.0, 1.e-8) = 1
 * so 1.e-19 is not considered to be zero within a 1.e-8 interval!
 *
 * In this function I take a different approach. First I compare the
 * two numbers and if abs(a-b) is less than epsilon = 10^-ndigits I declare
 * them equal; otherwise I call fcmp(a,b,epsilon) and if it returns
 * zero I call them equal.
 */
	double eps = std::numeric_limits<double>::epsilon();
	double epsilon = pow(10.0, -(double)ndigits);
	epsilon = std::max(epsilon, eps);

	if (std::abs(a) < eps && std::abs(b) < eps)
		return true;
	if (std::abs(a-b) < epsilon)
		return true;

	if (ndigits <= 0)
		throw runtime_error("is_equal called with zero significant digits");


	return (fcmp(a,b,epsilon) == 0);
}

bool
is_lessthan (double a, double b, int ndigits) {
	return (a < b && !is_equal(a, b, ndigits));
}

bool
is_greaterthan (double a, double b, int ndigits) {
	return (a > b && !is_equal(a, b, ndigits));
}

double
closest_10(double x) {
// Given a number "x", returns a number that is
// the closest power of 10 to that number
	Trace trc(2,"closest_10");
	double rval = 1.0;

	if (x <= 0.0) {
		trc.dprint("returning 1.0: x<=0");
		return rval;
	}

	if (x > rval) {
		x /= 2.0;
		while (10.0*rval < x) {
			// rval *= 2.0;
			rval *= 10.0;
		}
	} else {
		x *= 2.0;
		while (rval/10.0 > x) {
			// rval /= 2.0;
			rval /= 10.0;
		}
	}
	trc.dprint("returning ",rval);
	return rval;
}


/*------------------------------------------------------------------
fcmp
Copyright (c) 1998-2000 Theodore C. Belding
University of Michigan Center for the Study of Complex Systems
<mailto:Ted.Belding@umich.edu>
<http://www-personal.umich.edu/~streak/>                

This file is part of the fcmp distribution. fcmp is free software; you
can redistribute and modify it under the terms of the GNU Library
General Public License (LGPL), version 2 or later.  This software
comes with absolutely no warranty. See the file COPYING for details
and terms of copying.

File: README 

Description: README documentation file
************************************************************************

FCMP: SAFER COMPARISON OF FLOATING-POINT NUMBERS

It is generally not wise to compare two floating-point values for
exact equality, for example using the C == operator.  The fcmp package
implements Knuth's [1] suggestions for safer floating-point comparison
operators as a C function.

FCMP HOMEPAGE AND CURRENT RELEASES

The fcmp homepage is at <http://fcmp.sourceforge.net/>. There you can
always find the current release of fcmp, as well as the fcmp CVS
source code repository, which is available for anonymous read-only
access.

To receive announcements of new fcmp releases, you can subscribe to
the fcmp-announce email list by visiting the web page
<http://lists.sourceforge.net/mailman/listinfo/fcmp-announce> or by
sending an email to <fcmp-announce-request@lists.sourceforge.net> with
the body "subscribe". fcmp-announce is a very low volume list for fcmp
announcements only.

REPORTING BUGS 

To report a bug in fcmp, visit the fcmp bug tracking system on the web
at <http://sourceforge.net/bugs/?group_id=1799> or send an email to
<fcmp-bugs@lists.sourceforge.net> (the web page is the preferred
method). Please give full details about the bug, including what fcmp
version, platform, OS, and compiler you used.

USING FCMP

The fcmp function prototype is:

int fcmp(double x1, double x2, double epsilon);

where x1 and x2 are two floating-point numbers to be compared, and
epsilon determines the tolerance of the comparison (see the section
BACKGROUND for details). 

fcmp returns -1 if x1 is less than x2, 0 if x1 is equal to x2, and 1
if x1 is greater than x2 (relative to the tolerance).

For example, the following program should print "the result is: 1"

#include "fcmp.h"
#include <float.h>
#include <stdio.h>

int main() {
  int result;
  result = fcmp(3.0, 2.0, DBL_EPSILON);
  printf("the result is: %d\n", result);
  return 0;
}

BACKGROUND 

Floating point numbers are inherently inexact, and comparisons between
them are also inexact.  Equality or `==' comparisons are particularly
unreliable.  In the following code snippet, x is unlikely to be
exactly equal to y, even if both variable nominally have the same
value. This is due to sources of error such as truncation error and
rounding error.

float x;
float y;


if (x == y) {  unlikely to work!
  code to be executed if x == y
}

For example, the following program (from Priest [3]) may produce the
output "Equal" or "Not Equal", depending on the platform, compiler,
and compile-time options that are used:

#include <stdio.h>

int main() {
    double q;

    q = 3.0/7.0;
    if (q == 3.0/7.0) printf("Equal\n");
    else printf("Not Equal\n");
    return 0;
}

The same problem holds for `<', `>', `<=', `>=', and `!=' comparisons.
Since equality cannot be exactly determined with floating point
numbers, inequality cannot either.  

What is needed is a comparison operator that takes into account a
certain amount of uncertainty:

if (fabs(x - y) <= epsilon) {
  code to be executed if x == y
}

if (x - y > epsilon) {
  code to be executed if x > y
}

if (x - y < -epsilon) {
  code to be executed if x < y
}

In the above code, a neighborhood is defined that extends a distance
epsilon to either side of y on the real number line.  If x falls
within epsilon of y, x is declared to be equal to y (the first case,
above).  If x is greater than y by an amount that is greater than
epsilon, x is declared to be greater than y (the second case, above).
If x is less than y by an amount that is greater than epsilon, x is
declared to be less than y (the third case, above).

The problem then becomes to determine an appropriate value of epsilon.
A fixed value of epsilon would not work for all x and y; epsilon
should be scaled larger or smaller depending on the magnitudes of the
numbers to be compared.

A floating point number is represented by two numbers, the significand
(also called the fraction or mantissa) and the exponent, and a sign,
where

0 <= significand < 1 

and 

number = sign * significand * pow(2, exponent).

Knuth's suggestion is to scale epsilon by the exponent of the larger of the
two floating point numbers to be compared:

delta = epsilon * maxExponent,

where maxExponent is the exponent of max(x, y).  Delta can then be
substituted for epsilon in the code snippets above.

The routine fcmp() in this package implements Knuth's comparison
operators.  Given a value for epsilon, and two float or double numbers
x and y, it returns 1 if x is determined to be greater than y, 0 if x
equals y, and -1 if x is less than y.  (The routine automatically
converts floats to doubles, since there is no single-precision
equivalent for the standard C routines frexp() and ldexp(). However,
it works for both floats and doubles, without modification.)

DETERMINING EPSILON

Now that we have found a way to scale epsilon to work with a wide
range of x and y, we still need to choose an appropriate epsilon,
before scaling.  

If the number of binary digits of error, e, is known, then epsilon
can be calculated as follows:

epsilon = (pow(2, e) - 1) * FLT_EPSILON         (for floats)
epsilon = (pow(2, e) - 1) * DBL_EPSILON         (for doubles)

FLT_EPSILON and DBL_EPSILON are equivalent to 1 ulp for single- and
double-precision numbers, respectively; they are defined in the
standard C header file <float.h>. (An ulp is one unit in the last
place of the significand, or fraction part, of a floating point
number; see Knuth[1] for more details.)

TIMING RESULTS

To measure the performance penalty caused by using fcmp() instead
of normal floating-point comparisons such as <, I ran 4 different
benchmarks. In each benchmark except the first, which simply 
measured the overhead of the non-comparison portions of the program, 
60000000 comparisons were executed. The results are:

comparison type         total user seconds      seconds/comparison

no comparisons          31.55                   N/A
< (normal)              35.32                   6.28e-8
fcmp                    66.36                   5.80e-7
fcmp2                   66.62                   5.85e-7

The platform was a 266 MHz Pentium II, running Redhat 5.1 and egcs 1.1
(installed for i386-redhat-linux). All programs were compiled with
-O3. "fcmp2" is a variant of fcmp, where the formula

epsilon * max(fabs(x1), fabs(x2)) 

is used instead of 

epsilon * b^e(max(fabs(x1), fabs(x2)))

In summary, fcmp() is 9.23 times slower than normal floating point
comparisons (<).

This may seem like a huge performance hit, but it will not actually
have a major impact on program speed unless fcmp() is used in an
inner loop or another bottleneck within the program. As always, you 
should profile your program to see where the bottlenecks are before
attempting to optimize the program by hand for speed.

At some point in the future, I hope to make it possible to use fcmp as
a C macro, to avoid the overhead of calling it as a function.

BIBLIOGRAPHY

[2] Goldberg, David. (1991). What every computer scientist should know
about floating-point arithmetic. ACM Computing Surveys 23(1): 5-48.
Online at <http://www.validgh.com/goldberg/paper.ps>. The online
version includes an addendum by Doug Priest [3].

Goldberg, David. (1996). Computer arithmetic. Appendix A in Hennessy,
John L., and David A. Patterson, Computer Architecture: A Quantitative
Approach. Second edition. pp. A1-A77. San Francisco: Morgan
Kaufmann. ISBN 1-55860-329-9.

[1] Knuth, Donald E. (1998). The Art of Computer Programming.  Volume
2: Seminumerical Algorithms. Third edition. Section 4.2.2,
p. 233. Reading, MA: Addison-Wesley.  ISBN 0-201-89684-2.

Patterson, David A., and John L. Hennessy. (1998). Arithmetic for
computers. Chapter 4 in Computer Organization and Design: The
Hardware/Software Interface. Second edition. pp. 208-335. San
Francisco: Morgan Kaufmann. ISBN 1-55860-428-6. A more basic
introduction than Goldberg or Knuth.

[3] Priest, Doug. (1997). Differences among IEEE 754 implementations.
Addendum to Goldberg's [2] paper. Online at
<http://www.validgh.com/goldberg/addendum.html> and
<http://www.validgh.com/goldberg/paper.ps>; the latter version
includes Goldberg's [2] paper.

Summit, Steve. (1996). C Programming FAQs: Frequently Asked
Questions. Question 14.5, pp. 250-251. Reading, MA:
Addison-Wesley. ISBN 0-201-84519-9. Another algorithm for comparing
floating-point numbers.
 *------------------------------------------------------------------*/

int
fcmp(double x1, double x2, double epsilon) {
  int exponent;
  double delta;
  double difference;
  
  /* Get exponent(max(fabs(x1), fabs(x2))) and store it in exponent. */

  /* If neither x1 nor x2 is 0, */
  /* this is equivalent to max(exponent(x1), exponent(x2)). */

  /* If either x1 or x2 is 0, its exponent returned by frexp would be 0, */
  /* which is much larger than the exponents of numbers close to 0 in */
  /* magnitude. But the exponent of 0 should be less than any number */
  /* whose magnitude is greater than 0. */
  
  /* So we only want to set exponent to 0 if both x1 and */
  /* x2 are 0. Hence, the following works for all x1 and x2. */

  frexp(std::abs(x1) > std::abs(x2) ? x1 : x2, &exponent);

  /* Do the comparison. */

  /* delta = epsilon * pow(2, exponent) */

  /* Form a neighborhood around x2 of size delta in either direction. */
  /* If x1 is within this delta neighborhood of x2, x1 == x2. */
  /* Otherwise x1 > x2 or x1 < x2, depending on which side of */
  /* the neighborhood x1 is on. */
  
  delta = ldexp(epsilon, exponent); 
  
  difference = x1 - x2;

  if (difference > delta)
    return 1; /* x1 > x2 */
  else if (difference < -delta) 
    return -1;  /* x1 < x2 */
  else /* -delta <= difference <= delta */
    return 0;  /* x1 == x2 */
}


#ifdef MAIN

int
main(int argc, char** argv) {
// Usage:
//   fptype a b nsig  (test is_equal)
//   fptype           (test fcmp)
//   fptype str       (test str2double)
	double r;
	double num, den;
	double t;
	int ulp, maxulp = 40;

	// convert argv to strings
	vector<string> args = cstr2tok(argc, argv);

	// 3 arguments: test is_equal
	if (argc == 4) {
		double a;
		str2double(args[1], a);
		double b;
		str2double(args[2], b);
		int ndig = stoi(args[3]);
		cout << a << " and " << b << " are ";
		if (!is_equal(a,b,ndig)) {
			cout << "not";
		}
		cout << " equal\n";
		exit(0);
	}
	// 1 argument: test str2double
	if (argc == 2) {
		double x;
		str2double(args[1], x);
		cout << "double(" << args[1] << ") = " << x << endl;
	}

	// no arguments: test safe_divide
	if (argc == 1) {
		double smin = std::numeric_limits<double>::min();
		cout << "min: " << smin << endl;
		num = 1.0;
		den = -smin/2.0;
		cout << num << '/' << den << " = " << safe_divide(num,den) << endl;
		num = 0.0;
		den = -smin/2.0;
		cout << num << '/' << den << " = " << safe_divide(num,den) << endl;
		num = -smin/2.0;
		den = -smin/2.0;
		cout << num << '/' << den << " = " << safe_divide(num,den) << endl;
		den = 0.0;
		cout << num << '/' << den << " = " << safe_divide(num,den) << endl;

		// test fcmp
		double a1 = 1234.002;
		double a2 = 1234.0;
		double c = 1000.0;
		r = (c*(a1 - a2)+1.0)/7.0;
		if (r != 3.0/7.0) {
			cout << "r = " << c << "*(" << a1 << " - " << a2 << ")+1.0)/7.0"
				<< " = " << r << " != 3/7 = " << 3.0/7.0 << "\n";
		} else {
			cout << "r == 3/7\n";
		}

		t = 2.0;		// 2^ulp
		double eps = std::numeric_limits<double>::epsilon(); // (2^ulp - 1)*DBL_EPSILON
		double epsilon(eps);

		for (ulp=1; ulp<maxulp; ulp++) {
			if (fcmp(r, 3.0/7.0, epsilon) == 0) {
				cout << "r == 3/7 with " << ulp << " ulp error (" << epsilon << ")\n";
				break;
			}
			t *= 2.0;
			epsilon = (t - 1.0)*eps;
		}
		if (ulp == maxulp)
			cout << "r (" << r << ") != 3/7 (" <<  3.0/7.0 << ") with " << ulp
				<< " ulp error (" << epsilon << ")!\n";

		cout << "fcmp(1.e-19, 0, 1.e-8) = " << fcmp(1.e-19, 0.0, 1.e-8) << endl;

		// Test ulp_double
		ulp = ulp_double(r, 3.0/7.0);
		cout << "diff between r (" << r << " and 3/7 is " << (r-3.0/7.0) << ", " << ulp
			<< " ulp\n";
		int mantissa = std::numeric_limits<double>::digits; // digits in the mantissa
		cout << "mach eps = " << eps << ", should be 2^(-" << mantissa << ')'
			<< " = " << pow(2.0, -(double)mantissa) << endl;
		cout << "diff between r and 3/7 should be 2^-(" << mantissa
			<< '-' << ulp << ") = " << pow(2.0, (double)(ulp-mantissa)) << endl;
		cout << "r and 3/7 agreee in " << 100.0*ulp_frac(r, 3.0/7.0)
			<< " percent of their digits\n";
		// test is_valid
		double a{2.2};
		double b{2.2};
		double t = 1.0/(a-b);
		if (!flaps::is_valid(t))
			cout << t << " is not valid\n";
		t = std::numeric_limits<double>::quiet_NaN();
		if (!flaps::is_valid(t))
			cout << t << " is not valid\n";
		else
			cout << t << " is valid\n";
	}

}

#endif
