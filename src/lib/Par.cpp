//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#include <map>

#include "config.h"
#include "conv.h"
#include "adeqn.h"
#include "matrix.h"
#include "Par.h"
#include "trace.h"

using namespace std;


// register class Par for fio serialization
bool Par::regd = Fio::register_class(Par::id(), Par::get);

string
delimitedSubstr (string const& str, char open, char close,
		string::size_type& end);


//------------------------------------------------------------------
// class State
//------------------------------------------------------------------
std::string
Par::
state_desc() const {
	switch (this->state) {
		case State::nostate:
			return "no state";
		case State::indep:
			return "Independent";
		case State::derived:
			return "Derived";
		case State::fixed: {
			string rval{"Fixed"};
			if (!altval.empty())
				rval += vastr(altval.size()," values");
			return rval;
		}
		case State::aux:
			return "Auxilliary";
	}
	return "illegal state";
}

/*------------------------------------------------------------------
 * class Par implementation
 *------------------------------------------------------------------*/

Par::
Par(const string& defn, initializer_list<string> eqns) {
// Par parser constructor: parse a standard Par definition:
//	   name(desc)[ min:max ]<conv eunits/iunits> = value units {eqn} state
//	and a set of candidate equations in a std::initializer_list, e.g.
//	   {"mach/vsound", "freq/rf"}
//	see pset.c:make_biparams() for examples
	T_(Trace trc(1,"Par initializer_list");)

	// candidate equations
	for (auto& s : eqns)
		candidates.push_back(s);

	T_(trc.dprint(candidates.size()," candidates");)

	std::string minstr;
	std::string maxstr;
	std::string srhs;
	std::string eqn;

	// get name, desc, minstr, etc as strings
	Par::parse(defn, srhs, name, desc_, minstr, maxstr, conv, eqn);

	double factor = conv.factor;
	// if eqn is just a floating-point number, set advalue_
	// and set the state to fixed; otherwise set equation=eqn
	double x;
	if (str2double(eqn, x)) {
		advalue_ = x;
		state = State::fixed;
	} else {
		equation = eqn;
	}
	// parse the limits string, convert from PU to EU
	if (!minstr.empty()) {
		double min;
		if (str2double(minstr, min))
			minp = min/factor;
	}
	if (!maxstr.empty()) {
		double max;
		if (str2double(maxstr, max))
			maxp = max/factor;
	}

	pref = 0;
	state = State::nostate;
	fresh = false;
	T_(trc.dprint("returning ",this->longsummary());)
}


Par::
Par (std::string const& lhs, std::vector<std::string> const& rhsvec) {
// "Parser" constructor (multiple rhs).
	T_(Trace trc(1,"Par parser constructor ",rhsvec.size()," rhs");)
	std::string minstr;
	std::string maxstr;
	std::string srhs;
	std::string eqn;

	// only consider the first rhs in parsing, deal with rhsvec below
	if (!rhsvec.empty())
		srhs = rhsvec[0];

	// get name, desc, minstr, etc
	// XXX eqn is just a copy of srhs?
	Par::parse(lhs, srhs, name, desc_, minstr, maxstr, conv, eqn);

	double factor = conv.factor;

	// parse the limits string, convert from PU to EU
	if (!minstr.empty()) {
		double min;
		if (str2double(minstr, min))
			minp = min/factor;
	}
	if (!maxstr.empty()) {
		double max;
		if (str2double(maxstr, max))
			maxp = max/factor;
	}

	// set some defaults
	pref = 0;
	state = State::nostate;
	fresh = false;

	// check if the rhsvec are floats
	double eqnval{};
	if (!rhsvec.empty()) {
		string dblchar("+-.0123456789dDeE");
		bool isfloat = (rhsvec[0].find_first_not_of(dblchar) == string::npos);
		if (isfloat && str2double(rhsvec[0], eqnval)) {
			eqnval /= factor;
			advalue_ = Ad(eqnval);
			// put the rhs values into "altval" if there are more than one
			if (rhsvec.size() > 1) {
				for (auto& s : rhsvec) {
					double x;
					if (!str2double(s, x))
						throw runtime_error(vastr("in the definition of \"",lhs,
							"\", \"",s," is not a valid floating-point number"));
					x /= factor;
					altval.push_back(x);
				}
				T_(trc.dprint("now have ",altval.size()," altvalues");)
			}
		} else {
			// not floats - take them as equation candidates
			for (auto& s : rhsvec)
				candidates.push_back(s);
			// if there was only one rhs set it to the equation and
			// make the parameter Derived
			if (rhsvec.size() == 1) {
				equation = rhsvec[0];
				state = State::derived;
			}
		}
	}

	T_(trc.dprint("returning ",this->summary());)
}

Par::
Par (const std::string& lhs, const std::string& srhs) {
// "Parser" constructor (single rhs). Just put the rhs in
// a vector and call the multiple rhs version.
	T_(Trace trc(2,"Par parser(",lhs,") = ",srhs);)

	std::vector<std::string> svec;
	if (!srhs.empty())
		svec.push_back(srhs);
	*this = Par(lhs, svec);
}

void
Par::
adrealloc() {
// reallocate the advalue_ member of this parameter in order that
// it have the current number of AD derivatives. If there
// is currently an advalue_ transfer it's value to the new one.
	T_(Trace trc(2,"adrealloc ",this->name," ",Ad::nder()," derivs");)
	advalue_.realloc();
}

// functions for accessing & setting a parameter's value

double
Par::
value() const {
// Return the current value (no derivatives)
#ifdef NEVER // disable for Eduardo
	if (is_derived() && !fresh) {
		cerr << backtrace() << endl;
		throw runtime_error(vastr("attempt to return a stale value for ",name));
	}
#endif // NEVER // disable for Eduardo
	return advalue_.value();
}

double
Par::
value( double x) {
// Set the current value (no derivatives), return the current value
	double rval = advalue_.value();
	if (this->is_fixed())
		throw runtime_error(vastr("attempt to set fixed parameter ",name));
	advalue_.value(x);
	return rval;
}

double
Par::
deriv(size_t j) const {
// Return my derivative wrt the 0b jth Ad derivative
#ifdef NEVER // disable for Eduardo
	if (is_derived() && !fresh) {
		cerr << backtrace() << endl;
		throw runtime_error(vastr("attempt to return a stale derivative for ",name));
	}
#endif // NEVER // disable for Eduardo
	return advalue_.der(j);
}

void
Par::
deriv(size_t j, double x) {
// set my Ad derivative "j" (0-Ad::nder()-1) to "x"
	advalue_.der(j,x);
}

Ad
Par::
advalue() const {
// Return the current value as an Ad
	if (is_derived() && !fresh) {
		cerr << backtrace() << endl;
		throw runtime_error("attempt to return a stale advalue");
	}
	return advalue_;
}

void
Par::
advalue(const Ad& x) {
// Set my value as an Ad
	// it is an error to call this function for a fixed or indep parameter
	if (is_fixed() || is_indep()) {
		if (!is_equal(x.value(), value(), 8)) {
			cerr << backtrace() << endl;
			throw runtime_error(vastr("attempt to change ",summary(),
				" - use valuef() to change a ",state_desc()," value"));
		}
	}
	advalue_ = x;
}

Ad
Par::
iradvalue() {
// returns an Ad that is "inrange": no greater than max() or less than min();
// if it's value is out of range the Ad that is returns has zero derivatives
	double v = value();
	if (has_max() && v > max())
		return max();
	if (has_min() && v < min())
		return min();
	return advalue();
}

double
Par::
min() const {
	if (!minp)
		return -std::numeric_limits<double>::max();
	else
		return *minp;
}

double
Par::
max() const {
	if (!maxp)
		return std::numeric_limits<double>::max();
	else
		return *maxp;
}

Ad&
Par::
eval(pset& plt) {
// Evaluate this parameter by passing its equation to adeqn::eval
	T_(Trace trc(2,"Par::eval ",name);)

	T_(trc.dprint("state ",state_desc(),", uptodate ",fresh);)
	// quick return if not derived, constant, or has nostate
	// the State may not be set yet so check for an equation
	//!! if (equation.empty() || is_constant()) XXX cannot rely on constant
	if (equation.empty()) {
		T_(trc.dprint("quick return (",advalue_,"): derived? ",is_derived(),", constant? ",is_constant(),", state? ",state);)
		this->fresh = true;
	} else if (!equation.empty()) {
		bool noexc{false};	// throw exception on error
		if (plt.monitoring()) noexc = true;
		try {
			advalue_ = adeqn::eval(plt, equation, noexc);
		} catch (runtime_error& s) {
			throw runtime_error(vastr("evaluating ",this->name,": ",s.what()));
		}
		this->fresh = true;
	} else {
		if (!plt.monitoring())
			throw runtime_error(vastr(name," is Derived but has no equation"));
	}
	T_(trc.dprint("returning ",advalue_);)
	return advalue_;
}

bool
Par::
inRange(int sigfig) const {
// returns true if this parameter is in range
	T_(Trace trc(2,"Par::inrange");)
	double v{value()};
	double mn{min()};
	double mx{max()};
	T_(trc.dprint(name," value ",v,", min ",mn,", max ",mx);)
	if (has_min() && is_lessthan(v, mn, sigfig)) {
		T_(trc.dprint("returning false: ",name," is < ",mn," to ",sigfig," places");)
		return false;
	}
	if (has_max() && is_greaterthan(v, mx, sigfig)) {
		T_(trc.dprint("returning false: ",name," is > ",mx," to ",sigfig," places");)
		return false;
	}
	T_(trc.dprint("returning true");)
	return true;
}

double
Par::
inrange(int sigfig) const {
// Returns the difference between this parameter's min or max
// if it's value is out of range
	double val = value();
	double vmin = min() - val;
	double vmax = max() - val;
	double tol = pow(10.0, -sigfig)*val;
	if (vmin > tol)
		return vmin;
	if (vmax < -tol)
		return vmax;
	return 0.0;
}

void
Par::
parse(string const& lhs, string const& rhs, string& name, string& desc,
		string& minstr, string& maxstr, Conv& conv, string& eqn) {
/*------------------------------------------------------------------
 * parse an string definition of a parameter of the form
 *   name (desc) [ min : max ] <conv eunits/iunits> = value units {eqn} state
 *      ------ - --- - --- - -------------------- ------- -----
 * where
 *   name     name (e.g. "vtas", "rsf", etc)
 *   desc     description of the parameter (e.g. "Velocity ktas")
 *   min      (string) minimum value for the parameter in presentation units
 *   max      (string) maximum value for the parameter in presentation units
 *   conv     conversion factor: multiply quantities in EU to get PU
 *   eqn      (string) parameter value or equation in presentation units
 *
 * Input:
 *   lhs      string containing either the entire defn
 *            (if rhs is empty) or everything left of the '='
 *   rhs      string containing everything right of the '='
 *            or empty if lhs contains the whole mess
 * Output:
 *   no numerical values for min, max or rhs are returned, only
 *   strings because we do not know what datatype to expect.
 *
 * Throws runtime_error exception if parsing fails
 *------------------------------------------------------------------*/
	T_(Trace trc(2,"Par::parse(string rhs)");)

	if (lhs.empty())
		throw runtime_error("no left-hand side");

	string lhsbuf{lhs};
	stripquotes(lhsbuf);
	string rhsbuf{rhs};
	
	// Copy the lhs and (possibly) rhs into buffers (lhsbuf, rhsbuf)
	// The rhs may be contained in the lhs; if so split
	// it into lhs and rhs. Watch out for quoted lhs
	if (rhs.empty()) {
		if (lhsbuf.find("=") != string::npos) {
			throw runtime_error(vastr("lhs has an equal-sign: ",lhsbuf));
		}
	}
	T_(trc.dprint("lhs<",lhsbuf,"> = rhs<",rhsbuf,">");)

	string::size_type idx;
	string::size_type current = 0;
	string::size_type end;

	// Left-hand side: first look for an open parenthesis which
	// will end the name and start the description
	idx = lhsbuf.find('(', current);
	if (idx != string::npos) {
		name = lhsbuf.substr(current, idx-current);
		desc = delimitedString(lhsbuf, '(', ')', idx, end);
		if (desc.empty()) {
			string exc = vastr("no closing ) in ",lhsbuf.substr(idx));
			T_(trc.dprint("throwing exception: ",exc);)
			throw runtime_error(exc);
		}
		current = end + 1;
	}

	// next look for limits
	idx = lhsbuf.find('[', current);
	if (idx != string::npos) {
		string limitsStr = delimitedString(lhsbuf, '[', ']', idx, end);
		if (limitsStr.empty()) {
			string exc = vastr("no closing ] in ",lhsbuf.substr(idx));
			T_(trc.dprint("throwing exception: ",exc);)
			throw runtime_error(exc);
		}
		// if parsing the limits fails this is part of the name
		// if there was no description
		if (parseLimits(limitsStr, minstr, maxstr)) {
			if (name.empty())
				name = lhsbuf.substr(0, idx);
		} else {
			if (current != 0) {
				string exc = vastr("illegal limits: ",limitsStr);
				T_(trc.dprint("throwing exception: ",exc);)
				throw runtime_error(exc);
			}
		}
		current = end + 1;
	}

	// next look for a conversion
	idx = lhsbuf.find('<', current);
	if (idx != string::npos) {
		if (name.empty())
			name = lhsbuf.substr(0, idx);
		string convStr = delimitedString(lhsbuf, '<', '>', idx, end);
		if (convStr.empty()) {
			string exc = vastr("no closing > in ",lhsbuf.substr(idx));
			T_(trc.dprint("throwing exception: ",exc);)
			throw runtime_error(exc);
		}
		conv = Conv(convStr);
		current = end + 1;
	}

	// all done with lhs - take the whole lhs as name if empty
	if (name.empty())
		name = lhsbuf;

	name = stripwhitespace(name);

	// check the name for illegal characters - maybe this was
	// not intended to be a parameter
	if (name.find_first_of(" ={}()") != string::npos) {
		string exc = vastr("illegal parameter name: \"",name,"\"");
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}

	// Rhs: just return it as a string an let the caller
	// parse it because we don't know the datatype
	eqn = rhsbuf;

	T_(trc.dprint("returning name<",name,"> desc<",desc,"> minstr<", minstr,"> maxstr<",maxstr,">  eqn<",eqn,">");)

	return;
}

bool
Par::
parseLimits (string const& limits, string& minstr, string& maxstr) {
	T_(Trace trc(1,"Par::parseLimits");)
	ostringstream os;

	T_(trc.dprint("parsing \"",limits,"\"");)

	if (limits.empty()) {
		T_(trc.dprint("returning false: no limits");)
		return false;
	}

	std::string::size_type idx = 0;
	while (limits[idx] == ' ') idx++;

	if (limits[idx] == '[') idx++;

	std::string::size_type last = limits.find(':', idx);
	if (last == std::string::npos) {
		T_(trc.dprint("returning false: \"",limits, "\" is not a limits string: missing colon");)
		return false;
	}
	if (last != idx)
		minstr = limits.substr(idx, last-idx);

	idx = last+1;
	last = limits.find(']', idx);
	if (last == std::string::npos) {
		maxstr = limits.substr(idx);
	} else {
		maxstr = limits.substr(idx, last - idx);
	}

	T_(trc.dprint("returning true: min(",minstr,") max(",maxstr,")");)
	return true;
}


void
Par::
addConversion (const Conv& cp) {
// Give a parameter a conversion factor and convert its
// value, altval and limits from their current units (assumed
// to be presentation) to the new equation units.
// Throws an exception if this Par already has a Conv
	T_(Trace trc(1,"Par::addConversion ",cp);)
	std::ostringstream os;
	size_t i;

	// error if already have a conv
	if (!conv.eunits.empty())
		throw runtime_error(vastr("attempt to add conversion to ", *this));

	conv = cp;

	// convert values and limits from presentation (PU) to
	// equation units (EU)
	advalue_ /= conv.factor;
	if (has_min())
		minp = *minp/conv.factor;
	if (has_max())
		maxp = *maxp/conv.factor;
	if (!altval.empty())
		for (i=0; i<altval.size(); i++)
			altval[i] /= conv.factor;
	T_(trc.dprint("returning converted parameter: ",*this);)
}

void
Par::
update_from_solns (size_t index, const Par* from) {
// Set the current value of a parameter to the index-th (0b)
// element of its (or "from"s) solns array. Do this even if the
// parameter is Fixed
// "from" will be a nullptr if we are called with 1 arg
	T_(Trace trc(2,"Par::update_from_solns ",index);)
	if (index >= solns.size()) {
		string exc = vastr("attempt to access value ",
				index, "(0b) of ", name, " - only ",
				solns.size(), " are available");
		throw runtime_error(exc);
	}

	double x = this->solns[index];

	// if "from" is not a nullptr take the value from it instead of "this"
	if (from != nullptr) {
		if (name != from->name) {
			throw runtime_error(vastr("attempt to update ", name,
					" from ", from->name));
		}

		if (index >= from->solns.size()) {
			throw runtime_error(vastr("attempt to access value ",
					index, "(0b) of ", from->name, " - only ",
					from->solns.size(), " are available"));
		}
		x = from->solns[index];
	}

	if (is_fixed()) {
		if (!is_equal(x, this->value(), 8)) {
			T_(trc.dprint("changing fixed parameter ",name," from ",this->value()," to ",x);)
			valuef(x);
		}
	} else {
		value(x);
	}
	fresh = true;
}

double
Par::
min_solns() const {
	if (solns.empty())
		return 0.0;
	return *std::min_element(solns.begin(), solns.end());
}

double
Par::
max_solns() const {
	if (solns.empty())
		return 0.0;
	return *std::max_element(solns.begin(), solns.end());
}

void
Par::
update(const Par& from) {
// Set my advalue to from.advalue
// Throws exception if my state is Fixed
// and the values are not the same

	if (name != from.name) {
		throw runtime_error(vastr("attempt to update ",name," from ",from.name));
	}

	// Do not allow fixed parameters to change...
	if (this->is_fixed()) {
		double x = from.value();
		double y = this->value();
		if (!is_equal(x, y, 5)) {
			throw runtime_error(vastr("attempt to change ",summary(),
					" from ", x, " to ", y));
		}
	} else {
		advalue_ = from.advalue();
	}
}

void
Par::
upgrade(const Par* from) {
// Copy certain members from "from" if missing, and if the defnloc
// of "from->pref" is true copy certain values.
	T_(Trace trc(2,"Par::upgrade ", from->desc());)
	double factor = 1.0;

	// Description
	if (desc().empty() && !from->desc().empty()) {
		desc_ = from->desc();
		T_(trc.dprint("taking desc ",desc_);)
	}
	// Conversion factor: error if there already is one and it
	// is different...
	if (!from->conv.eunits.empty()) {
		if (conv.eunits.empty()) {
			if (has_min() || has_max()) {
				flaps::warning("giving \"", name, "\" a set of units and conversion (",
						from->conv, "). Currently it is: ", desc());
			}
			addConversion(from->conv);
		} else {
			conv = from->conv;
		}
	// Input par has no conv: assume it is in PU and convert it to EU
	// before taking values
	} else {
		if (!conv.eunits.empty())
			factor = conv.factor;
	}
	// If I have no limits but the input does, take them...
	if (!has_min() && from->has_min())
		min(from->min()/factor);
	if (!has_max() && from->has_max())
		max(from->max()/factor);

	// ... some members are replaced only if the precedence
	// of the input parameter is equal to or higher than the existing,
	// e.g. value, altval, min, max.
	if (from->pref > 0) {
		T_(trc.dprint("input has higher precedence");)
		this->valuef(from->valuef()/factor);
		T_(trc.dprint("took input value ", advalue_);)

		if (from->has_min())
			min(from->min()/factor);
		if (from->has_max())
			max(from->max()/factor);

		// input has altval...
		if (!from->altval.empty()) {
			altval.clear();
			for (size_t i=0; i<from->altval.size(); i++)
				altval.push_back(from->altval[i]/factor);
			T_(trc.dprint("taking ",from->altval.size()," alt values");)
		}
		// state...
		set_state(from->get_state());
		// input has equation: change state to derived
		if (!from->equation.empty()) {
			altval.clear();
			equation = from->equation;
			set_derived();
			candidates = from->candidates;
			if (candidates.empty())
				candidates.push_back(equation);
		}
		pref = 1;
	}
	T_(trc.dprint("returning ", *this);)
}  // upgrade()

vector<pair<size_t, double> >
Par::
find_solns(double val) const {
// Search my "solns" array for 2 values that enclose
// "val", return a vector of pair<size_t, double>, one for each
// enclosure found, such that for each pair<start,t>:
//    val = solns[start] + t*(solns[start+1]-solns[start])
//        = (1-t)*solns[start] + t*solns[start+1]
// The interpolated value for any other Par can be set by a call
// to Par::interp(pair<start,t>)
	T_(Trace trc(2,"Par::find_solns");)
	std::vector<pair<size_t, double> > rval;
	double eps{sqrt(numeric_limits<double>::epsilon())};

	T_(trc.dprint("searching for ",name," = ",val," in ",solns.size(),"-value array");)
	
	if (solns.empty()) {
		T_(trc.dprint("returning empty: solns is empty");)
		return rval;
	}
	if (solns.size() == 1) {
		if (is_equal(val, solns[0], 8)) {
			rval.push_back(std::pair<size_t,double>(0,0.0));
		}
		T_(trc.dprint("returning ",rval.size()," locations");)
		return rval;
	}

	size_t nm = solns.size() - 1;
	for (size_t i=0; i<nm; i++) {
		double a = val - solns[i];
		double b = solns[i+1] - val;
		if (std::abs(a) < eps) {
			rval.push_back(make_pair(i,0.0));
			T_(trc.dprint("found one at solns ",i,": ",solns[i]);)
			// continue until solns[i] - val > eps
			while(i<nm && abs(val-solns[i]) < eps) i++;
			continue;
		}
		if (std::abs(b) < eps) {
			rval.push_back(make_pair(i+1,0.0));
			T_(trc.dprint("found one at solns ",i+1,": ",solns[i+1]);)
			// continue until solns[i] - val > eps
			while (i<nm && abs(val-solns[i+1]) < eps) i++;
			continue;
		}

		if (a*b >= 0.0) {
			double t = 0.5;
			double den = solns[i+1] - solns[i];
			if (std::abs(den) > eps)
				t = (val - solns[i])/den;
			t = std::max(t, 0.0);
			t = std::min(t, 1.0);
			// ignore if this was found at upper end of previous interval
			if (i == 0 || t > 0.0) {
				rval.push_back(make_pair(i,t));
				T_(trc.dprint("found one between solns ",i," and ",i+1,", t = ",t);)
			}
		}
	}
	// check final point if it was not gotten above
	if (rval.empty() || rval[rval.size()-1].first < nm-1) {
		if(is_equal( val, solns[nm],8)) {
			rval.push_back(make_pair(nm,0));
			T_(trc.dprint("found one at ",nm+1);)
		}
	}
	T_(trc.dprint("returning ",rval.size()," locations");)
	return rval;
}

vector<pair<size_t, double> >
Par::
find_closest(double val) const {
// find all places in my "solns" array where the (approx) derivative
// wrt milepost is zero, i.e. where the slope changes sign:
//     (v[j] - v[j-1])*(v[j+1]-v[j]) < 0
// and return a vector of pair<idx,t> such that solns[idx] is approx
// where the slope changes (t is zero).
	T_(Trace trc(1,"Par::find_closest");)
	std::vector<pair<size_t,double> > rval;

	T_(trc.dprint("search for closest to ",val," in ",solns.size(),"-value array");)
	
	if (solns.size() <= 3) {
		T_(trc.dprint("returning empty: only ",solns.size()," solns");)
		return rval;
	}

	size_t nm1 = solns.size() - 1;
	for (size_t i=1; i<nm1; i++) {
		if ((solns[i] - solns[i-1])*(solns[i+1]-solns[i]) < 0.0) {
			rval.push_back({i,solns[i]});
		}
	}

	// sort by increasing abs(value) of second with a lambda criterion
	std::sort(rval.begin(), rval.end(),
		[](const pair<size_t,double> a, const pair<size_t,double> b) {
			return (std::abs(a.second) < std::abs(b.second)); });
			
	// set all seconds to zero
	for (auto& vp : rval)
		vp.second = 0.0;

	return rval;
}

void
Par::
interp(pair<size_t,double> loc) {
// set the current value (advalue_.value, no derivatives) to the value
// interpolated from the "solns" array at "loc":
//   value = solns[loc.first] + loc.second*(solns[loc.first+1]-solns[loc.first])
	T_(Trace trc(2,"Par::interp ",name);)
	size_t start = loc.first;
	double t = loc.second;
	double v1{solns[start]};
	double v{v1};
	if (start+1 < solns.size()) {
		double v2{solns[start+1]};
		v = v1 + t*(v2 - v1);
		T_(trc.dprint(name," = ",v1," + ",t,"*",v2-v1," = ",v);)
	}
	advalue_.value(v);
	this->fresh = true;
}

double
Par::
prevalue() const {
// Returns a parameters value in presentation units (PU)
// The parameter is *NOT* evaluated!
	return this->valuef()*conv.factor;
}

double
Par::
prevalue(double x) {
// Set a parameters value to "x", first converting it to equation units (EU)
	double rval{x/convFactor()};
	value(rval);
	return rval;
}

double
Par::
premin() const { return min()*conv.factor; }

double
Par::
premax() const { return max()*conv.factor; }

string
Par::
adsummary() const {
	return vastr(name," = ",advalue());
}

std::string
Par::
shortsummary() const {
	return parSum(0);
}

std::string
Par::
summary() const {
	return parSum(1);
}

std::string
Par::
longsummary() const { return parSum(3); }


std::string
Par::
presUnits () const {
	static std::string none("dimensionless");
	static std::string illegal("illegal");
	std::string rval;
	
	if (conv.punits.empty())
		return none;

	if (!conv.punits.empty())
		rval = conv.punits;
	else
		rval = illegal;
	return rval;
}

std::string
Par::
compUnits () const {
	static std::string none("dimensionless");
	static std::string illegal("illegal");
	std::string rval;
	
	if (conv.eunits.empty())
		return none;

	if (!conv.eunits.empty())
		rval = conv.eunits;
	else
		rval = illegal;
	return rval;
}

std::string
Par::
stringPreValue(int nsig) const {
	std::ostringstream os;
	if (nsig == 0)
		nsig = 6;

	os << std::showpoint << std::left << std::setprecision(nsig) << prevalue();

	return os.str();
}

std::string
Par::
parSum (int full) const {
/*------------------------------------------------------------------
 * Write a text description of a parameter to a string;
 * this string should be in the same form expected by
 * the Param(lhs,srhs) constructor.
 * If "full" is
 *   0  - name = value
 *   1  - name(desc) = value
 *   2  - name(desc)[limits] = value
 *   3  - name(desc)<conv> = value {eqn} state
 * the string returned when full=1 should result in the exact same
 * parameter if passed to the Param(lhs,srhs) constructor.
 *------------------------------------------------------------------*/
	T_(Trace trc(3,"Par::parSum");)
	std::string minstr, maxstr;
	std::ostringstream os;

	T_(trc.dprint("full = ",full);)

	os << name;

	if (full > 0 && !desc_.empty()) {
		os << " (" << desc_ << ')';
	}

	if (full > 1 && (has_min() || has_max()))
		os << stringLimits();

	if (full > 2 && !conv.punits.empty())
		os << " <" << conv << '>';

	os << " = " << this->advalue_*this->convFactor() << " " << this->presUnits();

	// candidate equations
	if (!candidates.empty())
		os << " (" << candidates.size() << " candidates)";

	//  equation
	if (full > 2) {
		if (!equation.empty())
			os << " {" << equation << "}";
		else
			os << " {no equation}";
	}

	// if (full > 2)
	os << " " << state_desc();

	if (full > 1 && solns.size() > 1) {
		os << " {" << solns.size() << " solns}";
	}

	T_(trc.dprint("returning ",os.str());)
	return os.str();
}

ostream&
operator<<(ostream& s, const State& t) {
// writes a string representation of a parameter's state
// enum class State { nostate, indep, derived, fixed, aux};
	switch (t) {
		case State::nostate:
			s << "no state";
			break;
		case State::indep:
			s << "independent";
			break;
		case State::derived:
			s << "derived";
			break;
		case State::fixed:
			s << "fixed";
			break;
		case State::aux:
			s << "auxilliary";
			break;
		default:
			s << "unrecognized state";
	}
	return s;
}

ostream&
operator<<(ostream& s, const Par& t) {
	s << t.parSum(3);
	return s;
}

std::string
Par::
table () const {
	std::string des;

	if (desc_.empty()) {
		if (name.empty())
			des = "<no name>";
		else
			des = name;
	} else {
		des = desc_;
	}
	std::ostringstream os;
	os << des;
	os << "    " << prevalue();

	if (!conv.punits.empty())
		os << ' ' << presUnits();
	return os.str();
}


std::string
Par::
stringLimits() const {
/*------------------------------------------------------------------
 * Write parameter limits (in presentation units) to a string
 * in the standard format, e.g.:
 *    [ -0.1 : 0.3 ]
 *------------------------------------------------------------------*/
	double fac = convFactor();
	std::ostringstream os;

	os << '[';
	if (has_min()) {
		double x = fac*min();
		os << std::setprecision(4) << x;
	}
	os << ':';
	if (has_max()) {
		double x = fac*max();
		os << std::setprecision(4) << x;
	}
	os << ']';
	return os.str();
}

// serialization functions

// get method: Receiver constructor
Fio*
Par::
get(Receiver& s) {
	std::ostringstream os;
	T_(Trace trc(2,"Par::get");)

	Par* rval = new Par();
	s.serialize(rval->name);
	T_(trc.dprint("got name \"", rval->name, "\"");)
	s.serialize(rval->desc_);
	T_(trc.dprint("got desc \"", rval->desc_, "\"");)
	s.serialize(rval->pref);
	T_(trc.dprint("got pref \"", rval->pref, "\"");)
	int ist;
	s.serialize(ist);
	rval->state = static_cast<State>(ist);
	T_(trc.dprint("got state \"", rval->state, "\"");)

	// conv
	rval->conv = Conv(s);

	rval->fresh = false;

	// only serialize the value (0th derivative)
	// XXX Ad::initialize() may not have been called yet
	double v;
	s.serialize(v);
	rval->valuef(v);

	// min, max, and default are optional<double>
	s.serialize(rval->minp);
	s.serialize(rval->maxp);

	// candidates
	s.serialize(rval->candidates);
	T_(trc.dprint("got ",rval->candidates.size()," candidates");)

	s.serialize(rval->equation);

	// solns
	s.serialize(rval->solns);
	T_(trc.dprint("got ", rval->solns.size(), " solns");)

	// do not serialize the altval

	return  rval;
}

// put() method
void
Par::
put(Sender& s) const {
	T_(Trace trc(2,"Par::put ", *this);)
	s.serialize(name);
	s.serialize(desc_);
	s.serialize(pref);
	int ist = static_cast<int>(state);
	s.serialize(ist);

	// conv
	conv.put(s);

	// only serialize the value (0th derivative)
	s.serialize(valuef());

	// min, max, and default are optional<double>
	s.serialize(minp);
	s.serialize(maxp);

	// candidates & solns: send using serialize(vector<T>)
	T_(trc.dprint(name, " has ",candidates.size(), " candidates");)
	s.serialize(candidates);
	// the chosen candidate
	s.serialize(equation);
	// put solns but not altval
	s.serialize(solns);
}

vector<string>
Par::
dependson(pset& plt) {
// returns a vector of the names of all parameters this parameter
// depends on including dependents-of-dependents
// Note: this only works properly after pset::equations()
	vector<string> rval;
	// no equation? no dependents
	if (this->equation.empty())
		return rval;

	// get direct dependents from the equation
	bool is_complex;
	vector<string> dep = adeqn::dependson(plt, this->equation, is_complex);
	// deep-dive recursively to get *all* dependents
	for (auto& di : dep) {
		Par* pp = plt.findp(di);
		if (pp == nullptr)
			continue;
		//!! if (pp->is_indep() || pp->is_eigv()) {
			if (std::find(rval.begin(), rval.end(), di) == rval.end())
				rval.push_back(di);
		//!! } else {
			vector<string> depi = pp->dependson(plt);		// recurse
			for (auto& dj : depi) {
				if (std::find(rval.begin(), rval.end(), dj) == rval.end())
					rval.push_back(dj);
			}
		//!! }
	}
	return rval;
}

// end of Par implementation
//------------------------------------------------------------------
/*--------------------------------------------------------------------
 * Conv implementation
 *------------------------------------------------------------------*/

Conv::
Conv (const string& des) {
/*
 * Parse a string like "<factor pu eu>"
 * or <factor pu/eu>.
 * If the second form is used pu and eu must be enclosed in parentheses
 * if they contain a slash (/); e.g.  0.049 knots/(in/s)
 * Also legal:
 *    <factor>
 *    <pu>
 *    <pu/eu>
 *    <pu eu>
 * If factor is not included it defaults to 1.0
 * If only pu is included eu defaults to pu
 */
	T_(Trace trc(2,"Conv parser constructor");)
	Conv::parse (des, factor, punits, eunits);
}

void
Conv::
splitUnits(string const& s, string& pu, string& eu) {
// given a string with one of the forms
//   pu
//   (pu)
//   pu/eu
//   (pu)/eu
//   (pu)/(eu)
//   pu/(eu)
// where eu are equation units and pu are presentation units.
// Returns pu and eu
	string str = stripwhitespace(s);
	std::string::size_type idx = str.find('/');
	std::string::size_type end;
	ostringstream os;

	if (idx == string::npos) {
		if (str[0] == '(') {
			pu = delimitedSubstr (str, '(', ')', end);
		} else {
			pu = str;
		}
		return;
	}
	if (str[0] == '(') {
		pu = delimitedSubstr (str, '(', ')', end);
		if (end != string::npos) {
			while(str[end] == ' ') end++;
			if (str[end] != '/') {
				throw runtime_error(vastr("missing / between presentation and "
						"equation units in ", str));
			}
			eu = str.substr(end+1);
			if (eu[0] == '(')
				eu = delimitedSubstr(eu, '(', ')', end);
		}
	} else {
		pu = str.substr(0,idx);
		eu = stripwhitespace(str.substr(idx+1));
		if (eu[0] == '(')
			eu = delimitedSubstr(eu, '(', ')', end);
	}
}

void
Conv::
parse (const string& des, double& factor, string& pu, string& eu) {
// Given a string description of a conversion factor like
//    <2.54 (furlongs/fortnight)/(m/s)>
// split it into
//   double factor (2.54)
//   pu: presentation units (furlongs/fortnight)
//   eu: equation units (m/s)
	T_(Trace trc(2,"Conv::parse ",des);)
	std::string::size_type end;
	string thestring(delimitedSubstr(des, '<', '>', end));
	ostringstream os;

	T_(trc.dprint("thestring ",thestring);)
	// check for some standard conversions, string versions of
	// the builtin conversion factors
	if (thestring == "m2ft") {
		factor = flaps::m2ft;
		eu = "m";
		pu = "ft";
		return;
	} else if (thestring == "kg2slug") {
		factor = flaps::kg2slug;
		eu = "kg";
		pu = "slug";
		return;
	} else if (thestring == "kg2lbm") {
		factor = flaps::kg2lbm;
		eu = "kg";
		pu = "lbm";
		return;
	} else if (thestring == "radps2Hz") {
		factor = flaps::radps2Hz;
		eu = "rad/s";
		pu = "Hz";
		T_(trc.dprint("returning eu<",eu,"> pu<",pu,"> fac ",factor);)
		return;
	} else if (thestring == "radps2rpm") {
		factor = flaps::radps2rpm;
		eu = "rad/s";
		pu = "rpm";
		return;
	} else if (thestring == "rad2deg") {
		factor = flaps::rad2deg;
		eu = "rad";
		pu = "deg";
		return;
	} else if (thestring == "mps2knot") {
		factor = flaps::mps2knot;
		eu = "m/s";
		pu = "knot";
		return;
	} else if (thestring == "Pa2psf") {
		factor = flaps::Pa2psf;
		eu = "Pa";
		pu = "psf";
		return;
	}

	// First token is either the factor, or a set of units; if it is
	// just a set of units they are the equation and presentation units
	// and the factor is 1.0
	vector<string> quotes;
	quotes.push_back("()");
	vector<string> tok = string2tok(thestring, " ", quotes);
	T_(trc.dprint("split into tokens: ",tok);)
	double x;
	if (tok.size() == 1) {
		if (str2double(thestring, x)) {
			factor = x;
		} else {
			factor = 1.0;
			eu = thestring;
			pu = thestring;
		}
		return;
	} else if (tok.size() == 2) {
		if (str2double(tok[0], x)) {
			factor = x;
			splitUnits(tok[1], pu, eu);
			if (eu.empty())
				eu = pu;
		} else {
			factor = 1.0;
			pu = tok[0];
			eu = tok[1];
		}
	} else if (tok.size() == 3) {
		if (!str2double(tok[0], x)) {
			string exc{"conversion factors must have the form "
				"<factor presentation_units equation_units>"};
			T_(trc.dprint("throwing exception: ",exc);)
			throw runtime_error(exc);
		}
		factor = str2double(tok[0], factor);
		pu = tok[1];
		string::size_type end;
		if (pu[0] == '(')
			pu = delimitedSubstr (pu, '(', ')', end);
		eu = tok[2];
		if (eu[0] == '(')
			eu = delimitedSubstr (eu, '(', ')', end);
	} else if (tok.size() == 4) {
		if (tok[2] == "/") {
			os << tok[0] << " " << tok[1] << tok[2] << tok[3];
			Conv::parse (os.str(), factor, pu, eu);
			return;
		} else {
			string exc{vastr("illegal conversion factor: \"", des, '\"')};
			T_(trc.dprint("throwing exception: ",exc);)
			throw runtime_error(exc);
		}
	} else {
		string exc{vastr("illegal conversion factor: \"", des, '\"')};
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}
	T_(trc.dprint("returning eu<",eu,"> pu<",pu,"> fac ",factor);)
	return;	
}

string
delimitedSubstr (string const& str, char open, char close,
		string::size_type& end) {
/*------------------------------------------------------------------
 * Extract a delimited string out of a string. str points to
 * either the first char of the delimited string or to the open
 * delimiter
 * The string may contain embedded open/close pair as long as
 * they are paired (AND open and close are different characters!)
 * Returns:  a copy of the delimited string without the delimiters
 * Returns:
 *   end    index of the character just past the "close" char or
 *          string::npos if "close" is the last char in the string
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"delimitedSubstr");)
	bool needclose = false;
	string rval;
	int level = 0;

	T_(trc.dprint("str ",str,", open ",open,", close ",close);)
	
	std::string::size_type first = 0;
	if (str[0] == open) {
		needclose = true;
		first = 1;
	}

	T_(trc.dprint("first<",first,">");)
	/*
	 * go through the string a char at a time, increasing
	 * "level" if an open delimiter is found, decreasing
	 * "level" if a close is found until a close is found
	 * and level is zero
	 */
	std::string::size_type last = first;
	while (last != str.size()) {
		if (str[last] == open && open != close)
			level++;
		if (str[last] == close) {
			if (level == 0)
				break;
			level--;
		}
		last++;
	}

	if (needclose && str[last] != close) {
		string exc = vastr("no matching ", close, " in ", str);
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	} else {
		rval = str.substr(first, last-first);
	}

	if (last == str.size()) {
		end = string::npos;
	} else {
		end = last+1;
	}
	T_(trc.dprint("returning ",rval);)
	return rval;
}

Conv::
Conv (Receiver& s) {
	T_(Trace trc(2,"Conv receiver constructor");)
	s.serialize(factor);
	s.serialize(eunits);
	s.serialize(punits);
	T_(trc.dprint("returning ",*this);)
}

void
Conv::
put(Sender& s) const {
	s.serialize(factor);
	s.serialize(eunits);
	s.serialize(punits);
}

ostream&
operator<<(ostream& s, const Conv& t) {
// writes a string representation of a conversion factor like
// 12.0 knots/(in/s). Note that if units contain a slash they
// are enclosed in parentheses for readability.
// First write the data to a local string, then write it to
// s in case s has format flags which should apply to the
// whole string, not just each piece of the Conv.

	// if this is a no-op conversion, e.g. "1 ft/ft" just write "ft"
	if (t.factor == 1.0 && t.punits == t.eunits) {
		s << t.punits;
		return s;
	}

	ostringstream os;
	os << t.factor;

	if (!t.punits.empty()) {
		if (t.punits.find('/') == string::npos) {
			os << " " << t.punits;
		} else {
			os << " (" << t.punits << ")";
		}
	}
		
	if (!t.eunits.empty()) {
		if (!t.punits.empty())
			os << "/";

		if (t.eunits.find('/') == string::npos) {
			os << t.eunits;
		} else {
			os << "(" << t.eunits << ")";
		}
	}
	s << os.str();
	return s;
}
//------------------------------------------------------------------
// end of Conv implementation
//------------------------------------------------------------------
//------------------------------------------------------------------
// Each eigenvector component has a name like "ev55.real"
// These functions are where the names are defined and tested for
//------------------------------------------------------------------
bool
Par::
is_eigv() const { return dynamic_cast<const Eigpar*>(this) != nullptr; }

std::string
evname(int k1b, bool imag) { return vastr("ev",k1b, imag?".imag":".real");}

std::string
revname(size_t rgcno) {
// create the name of a complex eigenvector component parameter,
// given the 0b component number of the real vector equivalent.
// For example revname(1) = ev[1].imag and
// revname(4) = ev[3].real
	int gcno = rgcno/2 + 1;
	bool isimag = false;
	if (rgcno%2 != 0)
		isimag = true;
	return evname(gcno, isimag);
}

int
gcparse (std::string const& name) {
// returns the 0b element number in the real representation
// of the complex gc vector, or -1 if it is not a gc component.
// so for example  gc4.imag -> 7
	// just replace "gc" with "ev" and call evparse
	if (name.substr(0,2) != "gc")
		return -1;
	string nm{name};
	nm[0] = 'e';
	nm[1] = 'v';
	return evparse(nm);
}

bool
is_absgcname (std::string const& name, int& cgcno) {
// return true/false if name is/is not the abs value of the abs value of
// a gc component. if it is, cgcno is the 1b element number of the
// complex gc vector, so for example
//    absgc41 -> cgcno = 41
	// just replace "gc" with "ev" and call is_evname
	if (name.substr(0,5) != "absgc")
		return false;
	string::size_type pos = name.find_first_of("0123456789");
	cgcno = stod(name.substr(pos));
	return true;
}


int
rgc2cgc (int rgcno, bool& is_imag) {
// given a 0-based component number of the real representation
// of a complex vector, return the 1-based component number of
// the complex representation (cgcno). So, for example
//   rgcno         cgcno
//     0            1
//     1            1
//     2            2
//     3            2
	int rval{0};
	if (rgcno%2 == 1) {		// imaginary part
		rval = (rgcno-1)/2 + 1;
		is_imag = true;
	} else {
		rval = rgcno/2 + 1;
		is_imag = false;
	}
	return rval;
}

int
evparse (const string& name) {
// does "name" refer to an eigenvector component? if so
// return the 0-based component number in the real equivalent
// vector, otherwise return -1
	int rval{-1};
	static regex* rx{nullptr};
	static regex* ix{nullptr};
	if (rx == nullptr)
		rx = new regex(R"(ev([1-9][0-9]*).real)");
	if (ix == nullptr)
		ix = new regex(R"(ev([1-9][0-9]*).imag)");
	smatch mch;
	if (regex_match(name, mch, *rx) || regex_match(name, mch, *ix)) {
		rval = 2*(stoi(mch[1])-1);
		if (rsubstr(name,4) == "imag")
			rval++;
		return rval;
	}
	return rval;
}

#ifdef MAIN
#undef MAIN

Par
testmove(string& lhs, string& rhs) {
	Par* rval = new Par(lhs, rhs);
	return move(*rval);
}

int
main (int argc, char* argv[]) {
	if (argc < 2) {
		cout << "Usage: Par par_defn\n";
		exit(1);
	}
	Ftmpdir ftmp;
	string lhs{argv[1]};
	string rhs;
	if (argc > 2)
		rhs = argv[2];
	// Par par(lhs,rhs);
	Par par{testmove(lhs,rhs)};
	cout << "got \"" << par << "\"\n";
	// test put/get
	string file{"Partest"};
	Sender s(file);
	par.put(s);
	Receiver t(file);
	Fio* fio = Par::get(t);
	Par* pp = dynamic_cast<Par*>(fio);
	if (pp == nullptr) {
		cerr << "Par::get failed\n";
		exit(1);
	}
	cout << "put/get returned " << *pp << endl;
}
#endif // MAIN
