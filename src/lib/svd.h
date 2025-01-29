//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef SVD_H
#define SVD_H

// minimum 2-norm solutions of underdetermined linear systems
// as described in Algorithm 5.6.2 in
//    Golub, G.H.; Van Loan, C.F. Matrix Computations, 4th ed; The
//    Johns Hopkins University Press: Baltimore, MD, USA, 2013

// factorize an (m,n) double, real matrix A = U \sigma V^t
// assuming m <= n:
//    u: (m,m)
//    s: (m,1)
//    vt: (n,n)

#include <string>
#include <vector>

#include "interval.h"
#include "message.h"

// exception classes thrown by SVD methods:
class SVDErr : public std::runtime_error {
public:
	SVDErr (const std::string& m) : std::runtime_error(m) {}
};

class SVDSingular : public SVDErr {
public:
	SVDSingular (const std::string& m) : SVDErr(m) {}
};

class SVD {
public:
	const double* a;
	int nr, nc; // dimensions of the input matrix A
	double rcond;
	int rankdef;
	std::vector<double> u; // (nr,nr)
	std::vector<double> vt;  // (nc,nc)
	std::vector<double> s;   // (nr,1) assuming nr <= nc

	// constructor(s)
	SVD() : nr(0), nc(0), rcond(0.0), rankdef(0) {}
	// factorize A = U*S*V(t)
	SVD (int nr, int nc, double const* A);
	~SVD () {}
	// minimum 2-norm solution to Ax = b
	void solve (std::vector<double> const& b, std::vector<double>& x);
	std::vector<Interval> solve (std::vector<Interval> const& b);

	std::vector<double> nullProj(std::vector<double> const& proj);

	int rank_def();

	int perfIndex (double const* x, double const* b);
	int perfIndex (double* x, double* dir);
};


/*
 * Some misc matrix  utilities...
 */
std::ostream&
operator<<(std::ostream& s, const SVD& t);
int performanceIndex (int nr, int nc, double const* a, double const* x, double const* b);
void analyze (FILE* stream, char const* title, int nr, int nc, double const* a);

#endif	/* SVD_H */
