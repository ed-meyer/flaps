//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Describing functions. 2 methods for computing DF's:
// 1) built in functions:
//	   Ad df::freeplay(pset& plt, Ad const& gap, int gcno);
//	   Ad df::bilinear(pset& plt, Ad const& x0,
//			Ad const& k1k2, int gcno);
//
// 2) Describing Function creation from a time-domain description of a
//    nonlinear "force" or a factor which multiplies a matrix element
//    to produce the time-domain force.
// XXX replace with df_q(), etc
//    The describing function is specified in a user-written class which has
//    the following form:
// class F {
//   // private data, e.g. gap size
//   double q;   // amplitude of oscillation: the displacement is
//               // y(t) = q sin wt
// public:
//   double operator()(double wt) {
//     // returns the factor at \tau = wt
//     double y = q*sin(wt);
//         ...
//     // we are integrating f*sin(wt) so mult by sin(wt):
//     return f*sin(wt);
//   }
//
//   Interp* interpolate(vector<double>& qs, const vector<double>& f) {
//     // create an Interp object: see interp.h
//   }
// };
//
// This class is known as a "function object" because it has operator() which
// is needed for the integration function

#ifndef DF_h
#define DF_h 1

//!! #include <boost/math/quadrature/trapezoidal.hpp>
#include <vector>

#include "interp.h"
#include "pset.h"

//	1) builtin functions
namespace df {
Ad freeplay(pset& plt, Ad const& gap, int gcno);
Ad freeplay(Ad const& gap, Ad const& absgc, double smooth);
Ad bilinear(pset& plt, Ad const& x0, Ad const& k1k2, int gcno);
Ad blfcn(Ad const& gamma, Ad const& ratio, double smooth);
}


// 2) time-domain descriptions
// Public interface:
Ad
dfval(pset& plt, const std::string& dfid, int gcno, const std::vector<double>& params);
void
dfplot(const std::string& dfid);

using Qfcn = std::vector<double> (*)(const std::vector<double>&);
using Ffcn = double (*)(double,const std::vector<double>&);
using Interpfcn = Interp* (*)(const std::string& iid,
		const std::vector<double>& params,
		const std::vector<double>& qs, const std::vector<double>& fs);

struct Df_fcns {
	std::string dfid;
	Qfcn qs_fcn;
	Ffcn fq_fcn;
	Interpfcn interp_fcn;

	Df_fcns(const std::string& d, Qfcn q, Ffcn f, Interpfcn i) : 
		dfid{d}, qs_fcn{q}, fq_fcn{f}, interp_fcn{i} {}
};

// userdf is created from the prototype lib/userdf.c by adding pointers
// to the 3 functions defining a DF. It is compiled into $FTMP/userdf.so
// so that it may be called from any program to get the function pointers
// and build a set of interpolation coefficients (df_build).
// df_build is called any time dfval is called (usually from user-written
// functions) and the interpolation coefficients have not been built yet.
bool
userdf (const std::string& dfid, Qfcn& qs_fcn, Ffcn& fq_fcn, Interpfcn& interp_fcn);

#endif // DF_h
