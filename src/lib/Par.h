//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef Par_h
#define Par_h

// 3 different sets of units:
// 1) "continuation" - (CU) used in the continuation method so that the
//                     variables have approx the same size
// 2) "equation"  -    (EU) so that all parameters have a consistent set
//                     of units in the equations; SI by default:
//                     meters, Newtons, seconds
// 3) "presentation" - (PU) units that make sense to the user, e.g. knots or m/s
// Conversions:
// EU = CU*Flutcurve::ivar[i].cu2eu()
// PU = EU*Par::Conv

// A Par has the general definition
//	   name(desc)[ min:max ]<conv eunits/iunits> = value units {eqn} state
//	where
//	  name    string with no blanks
//	  desc    arbitrary-length string, blanks ok
//	  min     limits
//	  max
//	  conv    (double) conversion factor: multiply value in EU to get PU
//	  eunits  (string) PU description (e.g. ft/s)
//	  iunits  (string) EU description (e.g. in/s)
//	  value   (double) current value
//	  units   (string) either the equation or presentation units description.
//	  eqn     (string) equation definition, parsed in adeqn
//	  state   current state of the parameter:
//	          nostate  undefined
//	          indep    independent
//	          derived  derived from other parameters
//	          fixed    declared constant
//	          aux      declared "auxilliary": from a side-calculation,
//	                   e.g. generalized-coordinates or stepsize parameters
//------------------------------------------------------------------
//
// The Par(string,string) constructor is the basic way
// to parse a user's parameter definition.
// Then you usually want to add it to the standard parameter list,
// while at the same time converting input values and units to
// equation units:
//    Par* newpp = gpset::get().add(&pp)
// which returns a pointer to the parameter in the gpset (a copy of pp).
//
// Evaluation of a parameter:
//    Par::eval()    evaluate the parameter using its Evaluator, putting
//                   the result in advalue_ member.
//------------------------------------------------------------------

#include <memory>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <vector>

void AttemptToUpdate(std::string);

// forward declarations
class pset;
class Conv;
class Par;
class Eigpar;


#include "Ad.h"
#include "message.h"
#include "fio.h"
#include "fptype.h"
#include "Regex.h"

class Conv {

	void parse (const std::string& des, double& factor,
			std::string& eu, std::string& iu);
	void splitUnits(std::string const& str, std::string& eu, std::string& iu);
public:
	double factor{1};			// multiply a quantity in "equation units" (EU)
									// by "factor" to get "presentation units" (PU)
									// XXX rename factor to eu2pu?
	std::string eunits;     // description of equation units
	std::string punits; 		// description of presentation units

	// constructors
	Conv () noexcept = default;
	Conv (const Conv&) noexcept = default;
	Conv (Conv&&) noexcept = default;
	// construct from a string like "factor eunits/iunits"
	Conv (const std::string& str);
	Conv (double fac, char const* pu, char const* eu) :
			factor(fac), eunits(eu), punits(pu) {}
	// default copy constructor, assignment op ok
	Conv& operator=(Conv const& rhs) = default; // assignment operators
	Conv& operator=(Conv&& rhs) = default;

	// serialization
	Conv(Receiver& s);
	void put(Sender& s) const;
	
};  // class Conv

// class State
enum class State { nostate, indep, derived, fixed, aux};

// a Par is a serializable parameter with automatic-differentiation
//
class Par : public Fio {
protected:
	std::string desc_;	// e.g. "Velocity (ktas)"
	Ad advalue_{0.0};		// current value as an AD in EU
public:
	std::string name;		// e.g. "vtas", "dpress"
	int pref{0};			// was the value set as a preference in this process?
	State state;				// indep|derived|fixed|aux
	bool constant{false};	// all my dependents are constant?
	bool fresh{false};		// up-to-date wrt my dependents?
	//!! bool active{false};		// active in an analysis?	XXX deprecate
	Conv conv;
	std::optional<double> minp;   // in equation units (EU)
	std::optional<double> maxp;   // in equation units (EU)
	// the parameter is evaluated in eval() using equation
	std::vector<std::string> candidates;  // candidate eqns
	std::string equation;
	// "soln" values are kept in a deque so we can push_front or push_back
	std::deque<double> solns; // in EU
	// "altval" are for multiple-fixed-valued parameters
	std::vector<double> altval; // in EU
	int work{0};

	// Constructors
	Par() noexcept = default;
	Par (const Par& ap) = default; // copy constructor: use clone()
	Par (Par&& ap) = default;      // move constructor
	// initializer-list constructor: defn contains lhs & srhs
	Par(std::string const& defn, std::initializer_list<std::string> candidate_eqns);
	// parsing constructors (with single or multiple rhs);
	// rhs may both be in lhs. A single rhs may be an equation
	Par (std::string const& lhs, std::string const& srhs);
	Par (std::string const& lhs, std::vector<std::string> const& svec);
	// copy constructors are private (above): use clone
	Par& operator=(Par const& rhs) = default; // assignment operators
	Par& operator=(Par&& rhs) = default;

	~Par() = default;		// dtor

	// virtual constructor or factory method
	virtual Par* clone() const { return new Par(*this); }

	// parse() parses parameter definitions in lhs=srhs
	static void parse(std::string const& lhs, std::string const& srhs,
			std::string& name, std::string& desc, std::string& minstr,
			std::string& maxstr, Conv& conv, std::string& eqn);
	static bool parseLimits (std::string const& limits,
			std::string& minstr, std::string& maxstr);

	// Accessing parameter value: an autodiff scalar (advalue_)
	// value()      returns a double: the current value (no derivs),
	//              throws if derived && stale
	// value(x)     sets the value to x, leaves derivatives untouched
	// valuef()     returns the current value even if stale
	// valuef(x)    sets the value to x, even if fixed
	// deriv(j)     returns a double: the value of the jth 0b AD derivative
	//              i.e. the partial of this par wrt AD par j
	// deriv(j,x)   sets the value of AD derivative j (0b) to x
	// deriv(pl,w)  returns a double: the partial wrt parameter "w" in "pl"
	// advalue      returns an Ad: the value+derivatives
	// advalue(x)   set the value+derivatives to x
	// iradvalue()  returns an Ad that is within the min/max
	double value() const;
	double value(double x);
	double valuef() const { return advalue_.value(); }
	void valuef(double x) { advalue_.value(x); }
	double deriv(size_t j) const;
	void deriv(size_t j, double x);
	Ad advalue() const;
	void advalue(const Ad& x);
	Ad iradvalue();
#ifdef NEVER // deprecate addata
	const double* addata() const { return advalue_.data(); }
#endif // NEVER // deprecate addata
	// min()        returns the minimum value; Caution: returns -std::max() if none,
	//              so check with has_min() first
	// max()        returns the max value, numeric_limits<double>::max() if none,
	//              so check with has_max() first
	// min(x)       sets the minimum value to x
	// max(x)       sets the maximum value to x
	bool has_min() const { if (minp) return true; else return false; }

	double min() const;
	bool has_max() const { if (maxp) return true; else return false; }
	double max() const;
	void min(double x) { minp = x; }
	void max(double x) { maxp = x; }
	void minmax (double& min, double& max) const;

	// set or get the value in presentation units:
	double prevalue() const;
	double prevalue(double val);

	// get the limits in presentation units
	double premin() const;
	double premax() const;
	double predefault() const;

	// AD elements must be re-allocated if the number of derivatives changes
	void adrealloc();

	// stuff for printing...
	std::string parSum(int full) const;
	std::string stringPreValue(int nsig=0) const;
	std::string stringPreDefault(int nsig=0) const;
	// limits in presentation units, e.g. [-1:]
	std::string stringLimits() const;
	std::string shortsummary() const;
	std::string summary() const;
	std::string longsummary() const;
	std::string adsummary() const;
	// 2 versions of desc(): 1st returns a const string, 2nd has an argument
	std::string const desc() const {
		if (desc_.empty()) return name;
		else return desc_; }
	void desc(std::string const& desc) { desc_ = desc; }

	bool has_conv() { return !(conv.factor == 1.0 && conv.punits.empty()); }
	std::string presUnits() const;
	std::string compUnits() const;
	std::string table () const;
	// a description of the current State:
	std::string state_desc() const;

	// Factor to convert from EU to PU
	double convFactor() const { return conv.factor; }	// XXX deprecate
	double conv_factor() const { return conv.factor; }
	void addConversion (const Conv& cp);

	void update(const Par& from);

	// change all my members depending on member "pref"
	void upgrade(const Par* from);

	// functions for handling the "solns" array
	void update_from_solns(size_t index, const Par* from=nullptr);
	double min_solns() const;
	double max_solns() const;
	size_t nsolns() const { return solns.size(); }
	void clear_solns() { solns.clear(); }
	std::vector<std::pair<size_t, double> > find_solns(double val) const;
	std::vector<std::pair<size_t, double> > find_closest(double val) const;
	// interpolate solns[start] + t*(solns[start+1]-solns[start])
   //        = (1-t)*solns[start] + t*solns[start+1]
	// where pair(start,t)
	void interp(std::pair<size_t, double>);

	bool operator==(const Par& ap) const { return name == ap.name; }
	bool operator!=(const Par& ap) const { return !operator==(ap); }

	// State functions
	bool is_nostate() const { return state == State::nostate; }
	bool is_indep() const { return state == State::indep; }
	bool is_derived() const { return state == State::derived; }
	bool is_fixed() const { return state == State::fixed; }
	bool is_eigv() const; // special function: test for Eigpar
	bool is_aux() const { return state == State::aux; }
	// a Par's constness is determined in make_equations
	bool is_constant() const { return (constant || is_fixed()); }

	// State handling (set/get) functions
	// XXX what's wrong with setting Aux to Indep/Derived/Fixed?
	void set_nostate() {
		if (state == State::aux) throw std::runtime_error("attempt to set Aux to nostate");
		state = State::nostate; }
	void set_indep() {
		//!! if (state == State::aux)
			//!! throw std::runtime_error("attempt to set Aux to Indep");
		state = State::indep; }
	void set_derived() {
		if (state == State::aux)
			throw std::runtime_error("attempt to set Aux to Derived");
		state = State::derived; }
	void set_fixed() {
		state = State::fixed; }
	void set_aux() {
		state = State::aux; }
	State get_state() const { return state; }
	void set_state(State const& s) { state = s; }

	// fio stuff
	// serialize() must be static so an instance is not necessary
	// static void serialize (Xdr*, Par<Type>**);
	static std::string id() { return "Par"; }
	std::string vid() const { return id(); }
	// Par (Receiver& s);
	static Fio* get (Receiver& s);
	void put (Sender& s) const;
	static bool regd;

	// Evaluation functions
	Ad& eval(pset& plt);
	bool inRange(int sigfig=8) const;
	double inrange(int sigfig=8) const;

	// return a vector of parameter names this parameter depends on
	std::vector<std::string> dependson(pset&);

}; // Par

// an Eigpar is a component (real or imaginary) of an eigenvector
// there is minimal difference between Eigpar and Par, just a way
// to make the distinction clear
class Eigpar : public Par {
public:
	Eigpar() : Par() {}
	// initializer-list constructor: defn contains lhs & srhs
	Eigpar(std::string const& defn,
			std::initializer_list<std::string> candidate_eqns) :
		 Par(defn, candidate_eqns) {}
	// parsing constructors (with single or multiple rhs);
	// rhs may both be in lhs. A single rhs may be an equation
	Eigpar (std::string const& lhs, std::string const& srhs) : Par (lhs, srhs) {}
	Eigpar (std::string const& lhs, std::vector<std::string> const& svec) :
		Par (lhs, svec) {}
	Eigpar(const Par& pp) : Par(pp) {}

	// virtual constructor or factory method
	virtual Eigpar* clone() const { return new Eigpar(*this); }
};


// ostream output
std::ostream& operator<<(std::ostream& s, const State&);
std::ostream& operator<<(std::ostream& s, const Conv&);
std::ostream& operator<<(std::ostream& s, const Par&);

// eigenvector and generalized-coordinate names
// evparse: returns the 0-based component # in the equivalent real
//          vector; e.g. ev4.imag = 7
//          Returns -1 if the name is not of an eigenvector component
// gcparse: returns the 0-based component # in the equivalent real
//          vector; e.g. gc4.imag = 7
//          Returns -1 if the name is not of a gc component
int
evparse(const std::string& name);
int
gcparse(const std::string& name);
// 2 ways to get the name of an element of a complex eigenvector:
//   evname:  given the 1-based component # in the complex vector and real/imag
//   revname: given the 0-based component # in the equivalent real vector
// rgcno  - 0-based component number in the equivalent real vector
std::string
evname(int k1b, bool imag);
std::string
revname (size_t rgcno);   // name of complex element given 0b element #
// return true/false if name is/is not an gc component name.
// if it is, rgcno is the 0b element number in the real
// representation of the complex eigenvector, so for example
//    ev4.imag -> rgcno = 7
bool is_gcname(std::string const& name, int& rgcno);  // parse a gccomp name
bool is_absgcname(std::string const& name, int& cgcno);  // parse a absgc comp name
// to convert a 0b rgcno to a 1b cgcno - i.e. convert a 0b component of
// the real representation to the 1b component of the complex representation
int rgc2cgc(int rgcno, bool& is_imag);

#endif // Par_h
