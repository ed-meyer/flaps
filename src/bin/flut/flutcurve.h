//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef FLUTCURVE_H
#define FLUTCURVE_H 1

#include <iostream>
#include <mutex>
#include <string>
#include <vector>

class Flutcurve;
class Fstack;

#include "Curve.h"
#include "pac.h"
#include "specs.h"

// 3 different sets of units:
// 1) "continuation" - (CU) used in the continuation method so that the
//                     variables have approx the same size
// 2) "equation"     - (EU) so that all parameters have a consistent set
//                     of units in the equations; by default inches, lbf
// 3) "presentation" - (PU) units that make sense to the user, e.g. knots or m/s
// Conversions:
// EU = CU*cu2eu where cu2eu = Flutcurve::ivar_.second
// PU = EU*Par::Conv

//------------------------------------------------------------------
// Class for creating function objects to add constraints to a Pac
//------------------------------------------------------------------
class Cfo {
public:
	// Constraint function object to be used as a Pac::Constraintfcn
	Target* target;
	Flutcurve* fcv{nullptr};
	// constructors
	Cfo() = default;
	Cfo(Target* t, Flutcurve* f) : target(t), fcv(f) {}
	// assignment op
	// Cfo& operator=(Cfo&) = default;
	//!! Cfo& operator=(Cfo&&) = default;
#ifdef NEVER // mv to .cpp
	double operator()(std::vector<double> const& x, std::vector<double>& c) {
		return fcv->constraintfcn(target, x, c);
	}
#else // NEVER // mv to .cpp
	double operator()(std::vector<double> const& x, std::vector<double>& c);
#endif // NEVER // mv to .cpp
};

class Ivar {
public:
	Par* par;
	double cu2eu_;
	Ivar(Par* p, double c) : par(p), cu2eu_(c) {}
	double cu2eu() const { return cu2eu_; }
	// get/set the parameter value in EU (not AD derivatives)
	double value() const { return par->value(); }
	double value(double v) { return par->value(v); }
	double cuvalue() const { return par->value()/cu2eu_; }
	double cuvalue(double cuv) { return par->value(cuv*cu2eu_)/cu2eu_; }
};

std::ostream&
operator<<(std::ostream& s, const Ivar& t);

class Discontinuity {
public:
	int coord;
	std::vector<double> angles;
	Discontinuity(double c, double a) : coord(c) { angles.push_back(a); }
};

class Flutcurve : public Curve {
	int phase_index_{0};			// phase-normalization element in an Ivar 0b
	int evstart_{0};				// start of eigenvector in ivar_ 0b
	std::vector<Ivar> ivar;		// indep variables: EU=CU*ivar.cu2eu
	double time{0.0};				// time to trace
	Cfo cfo;							// constraint
public:
	//!! Uncert* uncert{nullptr};
	Pac origin;
	std::string atlimit;
	Specs specs;		// so this curve can have different specs?
	int thread{0};					// tracked on thread #
	int ncev{0};					// number of complex eigenvector components
	// workspace for pac_fjac
	std::vector<std::complex<Ad>> adev;	// eigenvector as complex AD
	std::vector<std::complex<Ad>> dynmat;	// dynamic matrix
	std::vector<std::complex<Ad>> work;	// workspace for dmatrix
	std::vector<std::complex<Ad>> Dy;	// dynmat*adev
	std::vector<Pac> bif_origins;
	// map between Ad derivatives and ivar index
	std::vector<int> ad2ivar;
	Discontinuity* discon{nullptr};

	// constructors
	Flutcurve(std::string const& analysis_id, std::string const& cid, std::string const& vzid);
	Flutcurve(std::string const& analysis_id, std::string const& cid,
			std::string const& vzid, Pac& origin);
	// delete the copy constructor, assignment op
	Flutcurve(Flutcurve&) = delete;
	void operator=(const Flutcurve&) = delete;

	// clear the "solns" vector
	void clear_solns() { this->params.clear_solns(); }

	// nx: the number of independent variables
	int get_nx() { return ivar.size(); }

	// update my params
	void update(const std::vector<double>& x);
	void update(const Pac& v);
	void e9n_update(Pac const&);

	// trace the curve starting at "origin"
	bool pac_track (double newcoord); 

	// compute the dynamic matrix (static for dependson)
	static void dmatrix(pset& plt, std::vector<std::complex<Ad>>& result,
			std::vector<std::complex<Ad>>& work);
	// compute the jacobian and function
	int pac_fjac (std::vector<double> const& x, std::vector<double>& f,
		std::vector<double>& jac);
	// get the names of all parameters used
	static std::vector<std::string> dependson(pset& plt);

	// access elements of the Ivar vector
	int index(const std::string& parname) const;
	int phase_index() { return phase_index_; }
	int evstart() { return evstart_; }

	// compute a scale factor to give an indep variable a reasonable
	// range for the continuation method; the scaled variables are in
	// "continuation units" or CU, so EU = cu2eu*CU
	double cu2eu(Par* pp);

	// create a constraint...
	Cfo pac_constraint(Pac& to, std::string& atlim);
	// ... and evaluate it
	double constraintfcn(Target* t, const std::vector<double>& x, std::vector<double>& c);

	// process the results from Pac::go: call from a lambda
	int processfcn(Pac& from, Pac& to, bool& looped, Issue*);

	// check determinant signs between 2 points: a change means
	// a possible bifurcation
	int checkdsc(Pac& a, Pac& b);

	// get a summary of the current (toprint) values
	std::string current_values();
	// ... and a summary of all values
	std::string print_solns(const std::vector<std::string>& leaders, bool firstlast);
	// a header for either of the above
	std::string make_header() const;
	// summarize an nx-vector of indep variables (Ivar)
	std::string summarize(std::vector<double> x);

	// to plot all stored flutcurves with "aid" or in a .apf file with fcurves:
	static void plot(std::string const& aid, std::string const& plotfile,
			const std::vector<Flutcurve*>& fcurves);
	void normalize (int n, std::complex<double>* cv, int& indexReal);

	// compute lcostab
	void lcostability(Pac& p);

}; // Flutcurve

// an Fstack is a vector of Flutcurve pointers with one member function:
// pop() returns a pointer to the next Flutcurve* to be traced, guarded
// for thread safety, unlike std::stack
class Fstack : public std::vector<Flutcurve*> {
	size_t next{0};
	std::mutex mtx;
public:
	Fstack() : next{0} {}
	Fstack(std::vector<Flutcurve*>& t);
	Flutcurve* pop();
};

// Trace a set of curves by getting the next available from curves->pop()
void trace_curves(Fstack* curves, int threadno);
// correct the parameter values in Flutcurve by homotopy
//!! void homotopy (std::vector<Flutcurve*>& curves, const std::string& constpar);


#endif // FLUTCURVE_H
