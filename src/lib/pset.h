//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef PSET_H
#define PSET_H

class pset;

#include <map>

#include "Ad.h"
#include "fio.h"
#include "Par.h"
#include "Stack.h"

//------------------------------------------------------------------
// pset is a container (map) for parameters (Par*); it inherits
// from Fio so it has put/get
//
// Eigenvectors
// ------------
// Eigenvector components are kept as parameters: the real
// and imag parts of each component are kept in a separate parameter
// with a name like "ev1.real" or "ev22.imag". The values
// are changed like any other parameter, with Par::advalue(x).
// updateEv(std::vector<double> const& ev) is
// provided to allow the ev parameters to be updated from a vector.
// The eigenvector is always normalized to 2-norm 1.0.
//
// Generalized coordinates are the product of the parameter "gcnorm"
// times the eigenvector. In a linear flutter analysis gcnorm is
// 1.0 so the generalized-coordinate vector and the eigenvector
// are the same.
//
// The eigenvector components are accessed with:
// 	getadev(vector<complex<Ad>>& rval):
// 	   fill "rval" with eigenvector as complex AD 
// 	vector<double> getev():
// 	   current eigenvector as doubles
//    vector<complex> get_eigenvector(size_t idx):
//       the idx-th Par::soln of each eigenvector component
//    vector<const Par*> getnlev() const
//    vector<Par*> getnlev():
//       returns pointers to eigenvector Pars which others are dependent on
//    vector<vector> get_evmatrix() const:
//       all "Par::soln" values of the eigenvector Pars collected into a
//       vector of vectors
//    void updateEv(vector<double> const& ev)
//    void updateEv(size_t n, double const* ev)
//       update the eigenvector components from a vector of doubles
//------------------------------------------------------------------

class pset : public Fio {
   std::string desc_;
	bool eqns_made{false};
	std::deque<std::vector<std::string>> monitors;
	std::map<std::string,Par*> pmap_;
	std::vector<Eigpar*> eigv_;
public:

	// mostly use default (map) constructors...
	pset() : eqns_made{false} { }
	pset(std::string const& ident);
	pset(const pset& from);  // copy constructor
	pset& operator=(const pset& rhs); // assignment operator
	// XXX move ctor & assignment?

	// pset destructor deletes all contained Par
	~pset();


	// to access the maps:
	std::map<std::string,Par*>& pmap() { return pmap_; }
	const std::map<std::string,Par*>& pmap() const { return pmap_; }
	std::vector<Eigpar*>& eigv() { return eigv_; }
	const std::vector<Eigpar*>& eigv() const { return eigv_; }
	size_t size() const { return pmap_.size()+eigv_.size(); }
	bool empty() const { return size() == 0; }
	void clear() { pmap_.clear(); eigv_.clear(); }
	// get all parameters of various states:
	std::vector<Par*> all();	// all except nostate? even Eigpar?
	std::vector<Par*> get_nostate();
	std::vector<Par*> get_indep();
	std::vector<Par*> get_derived();
	std::vector<Par*> get_fixed();
	std::vector<Par*> get_mvf();
	std::vector<Eigpar*> get_eigv();
	std::vector<Par*> get_aux();

	// set up equations for each parameter by trying each
	// equation in "candidates"
	void equations();
	int choose_eqn (Par& par); // give "par" an evaluator
	// has equations() been called?
	bool has_eqns() const { return eqns_made; }
	bool has_eqns(bool state) { eqns_made = state; return true; }

	// eigenvector handling stuff: eigenvector names (evname) in Par.h
	// put pointers to all eigenvector components in member "eigenvector"
	// get the eigenvector as complex autodiffs...
	void getadev(std::vector<std::complex<Ad>>&);
	// ... or just the values...
	std::vector<double> getev();
	// ... or the complex idx (0b) solns...
	std::vector<std::complex<double>> get_eigenvector(size_t idx);
	// ... or as a matrix of the gc's XXX return a vector<complex>?
	std::vector<std::vector<double> > get_evmatrix();
	// update the eigenvector from a double* or a vector<double>
	void updateEv(size_t n, double const* ev);
	void updateEv(std::vector<double> const& ev);
	// create the real and imag parameters for components 1-nev...
	void make_eigv(int nev, double init=1.0);  // in pset.cpp
	// ... or just component k
	//!! std::vector<Par*> make_eigvk(int k, double init=1.0);
	std::vector<std::string> make_eigvk(int k, double init=1.0);

	// push current values from list "from" onto the "solns"
	// array of "this" either in the front or the back
	void back_solns(pset& from);
	void front_solns(pset& from);

	// set my parameter values to those in "pl"; 2 versions:
	// 1) update changes the advalue and throws an exception on
	//    attempt to change Fixed Par
	// 2) setValues  sets only the value part of the advalue of
	//    each Par, and will change Fixed parameters
	bool update(pset const& pl);
	bool setValues (const pset& pl);
	void freshen();		// set all Par fresh to true

	// return a vector of pointers to all parameters in this pset
	std::vector<Par*> allpar();
	// get a pointer to a particular parameter
	const Par* findp (const std::string& name) const;
	Par* findp (const std::string& name);
	Par* findg (const std::string& name);
	// CAUTION when using findRe: you could match, e.g.
	// "cdpress" when you are looking for "dpress" (it happened).
	// Probably better to anchor the re: "^dpress$"
	std::vector<Par*> findRe(const std::string& pattern);
	// find parameters matching a wildcard patter
	std::vector<Par*> findwildcard(const std::string& pattern);
	Par* findic(const std::string& name);

	// interpolate each Par's solns:
	//   t = solns[start] + t*(solns[start+1]-solns[start])
	// set current value to t
	void interp(std::pair<size_t, double>);

	// add a new parameter or members to an existing parameter,
	// return a pointer to the one on the list which may be a copy
	Par* add (const Par* pp);
	// ...similar, but replace the entire Par if existing
	Par* replace (const Par* pp);

	// check all parameters that have limits for in-range, return
	// an out-of-range (or empty) message
#ifdef NEVER // use inrange
	std::string inRange();	// XXX deprecate: use inrange
#else // NEVER // use inrange
	std::vector<std::string> inrange();
#endif // NEVER // use inrange

	// printing stuff
	//   desc():          short string, e.g. "mode 2,target=growth"
	//   desc(string):    set the description
	//   summary()                desc() + number of parameters
	//   summary(vector<string>)  one-line current values of "toprint" parameters
	//   summary_solns    one-line for each element of "solns" array
	//   display          summary of all parameters, one line per
	//   operator<<       long summary of each parameter
	std::string desc() const { return desc_; }
	void desc(std::string const& des) { desc_ = des; }
	std::string summary();
	// one-line summary of current values
	std::string summary(std::vector<std::string> const& toprint, int maxchar=160);
	// one-line summary of each set in the "solns" arrays
	// Note: summary_solns  modifies "this" so it cannot be const
	std::string summary_solns(std::vector<std::string> const& toprint,
		const std::string& header, std::vector<std::string> const& leader,
		bool firstlast);
	std::string display (std::string const& path, std::string const& title) const;
	std::string summaryTable() const; // one line (name=val) for each param

	// plot specified parameters
	void plot (std::string const& path, std::string const& title,
			std::string const& runid,
			std::vector<std::string> const& toPlot, bool append=false);
	// simple plot: all parameters, use desc() as path, title & curve id
	void plot();

	// fetching, storing
	static pset* fetch(const std::string& mid);
	void store(std::string const& mid);
	//! static pset pop(std::string const& mid);

	// set the appropriate spot in the advalue member of each parameter
	// which is an AD parameter to 1
	void set_adpar_derivs();

	void touch();		// set all dependent params out-of-date
	void eval();      // evaluate all parameters
	//!! void evalf();      // force evaluation of all parameters

	size_t nsolns() const;  // throws exception if not all same
	virtual void clear_solns();
	// Set the values of all (except Fixed) parameters in this pset from element
	// "index" of the "solns" arrays of my Par (from == nullptr)
	// or the Par in "from".
	void update_from_solns(size_t index);

	bool verify (const pset& from);

	// Fio/serialization
	static std::string id() { return std::string("pset"); }
	std::string vid() const { return id(); }
	void put (Sender&) const;
	static Fio* get (Receiver&);
	static bool regd;  // has been registered?
	Fio* clone() const { return new pset(*this); }

	// parval(name, inrange):  return the AD value of Par "name", inrange
	// setpar(x,name):  change the AD value of Par "name"
	Ad parval (std::string const& name, bool inrange=false);
	Ad absgcv(int gcno);  // abs gc gcno (1b)
	Ad absgcp(int gc1, int gc2);  // abs(gc1-gc2) (1b)
	void setpar (Ad const& x, std::string const& name);
	// x is either a parameter name or a double
	Ad par_or_double(const std::string& x);

	// monitoring calls to findp: save names in dependents, return
	// when monitoring is turned off
	// monitor():   turn on monitoring
	// monitor_add(string nm):  add nm if monitoring and nm is not already in
	// rotinom():  stop the most recent monitor, return the dependent parameter names
	void monitor();
	void monitor_add(const std::string& nm);
	std::vector<std::string> rotinom();
	bool monitoring() const { return !monitors.empty(); }

};  // pset

/*------------------------------------------------------------------
* gpset: the global parameter set - available to all programs,
* stored in-between programs.
* The only way to access the only instance of gpset is by calling
* gpset::get() which return a reference to the gpset.
* The first time gpset::get() is called in a program, either
* an existing set is fetched or a new one is created and filled
* with make_gpset().
*------------------------------------------------------------------*/

class gpset : public pset {
protected:
	gpset(); // one and only constructor: private

public:
	~gpset(); // stores the gpset on exit

	static std::string persistentName();    // specialialized in pset.c

	// reallocate all parameter advalues with new AD parameter names,
	// intended to be called immediatly after calling Ad::initialize()
	void realloc();

	// access the list like this:
	//    gpset& gps = gpset::get();
	// the first time get() is called the static instance
	// (thegpset) is constructed
	static gpset& get() {
		static gpset thegpset;
		return thegpset;
	}
	static Par* find(const std::string& name); // in pset.c

	// special store: no solns (in pset.c) XXX necessary now that solns/altval diff?
	static void store();

	// print a table of all Par split into states
	std::ostream&  table(std::ostream& s, const std::string& title="");

}; // gpset


// for finding an element of getIndep(), etc functions above
// 2 versions: for const and non-const lists
const Par*
find(const std::vector<const Par*>& list, const std::string& nm);

Par*
find(std::vector<Par*>& list, const std::string& nm);

std::ostream&
operator<<(std::ostream& s, const pset& pl);

//!! std::ostream&
//!! operator<<(std::ostream& s, gpset& pl);

#endif // PSET_H
