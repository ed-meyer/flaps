//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <sys/types.h> // for getpid()
#include <unistd.h> // for getpid()

#include "config.h"
#include "conv.h"
#include "adeqn.h"
#include "atmos.h"
#include "Curve.h"
#include "exim.h"
#include "matrix.h"
#include "./pset.h"
#include "settings.h"
#include "Stack.h"

using namespace std;
using namespace flaps;

//------------------------------------------------------------------
// pset implementation
//------------------------------------------------------------------
pset::
~pset() {
// pset destructor deletes all contained Par
	Trace trc(1,"pset destructor ",desc());
	for (auto pp : this->pmap_)
		delete pp.second;
	for (auto pp : this->eigv_)
		delete pp;
}

// register class pset for fio serialization
bool pset::regd = Fio::register_class(pset::id(), pset::get);

Fio*
pset::
get(Receiver& s) {
// read a pset, each Par and insert into our map
	size_t n; // the number of parameters
	pset* rval = new pset();
	s.serialize(rval->desc_);
	s.serialize(n);
	for (size_t i=0; i<n; i++) {
		Fio* np = Fio::get(s);
		Par* newpar = dynamic_cast<Par*>(np);
		if (newpar == nullptr)
			throw runtime_error(vastr("expecting Par, got ", np->vid()));
		// all parameters are "fresh" until equations are made...
		newpar->fresh = true;
		// ...and they have not been changed by a preference
		newpar->pref = 0;
		// insert in pmap or eigv
		if (evparse(newpar->name) >=0)
			rval->eigv().push_back(new Eigpar(*newpar));
		else
			rval->pmap_.insert(make_pair(newpar->name, newpar));
	}
	return rval;
}

void
pset::
put(Sender& s) const {
	s.serialize(desc_);
	s.serialize(this->size());
	for (auto& pp : this->pmap_)
		Fio::put(s, *pp.second);
	for (auto& pp : this->eigv_)
		Fio::put(s, *pp);
}

// constructors
pset::
pset(std::string const& ident) {
// no parameters are added, only the copy constructor does that
	Trace trc(1,"pset constructor ",ident);
	desc_ = ident;
	eqns_made = false;
}

pset::
pset(const pset& from) {
// copy constructor (deep: creates new Par)
	Trace trc(1,"pset copy constructor ",from.desc());
	for (const auto& i : from.pmap())
		this->pmap_.insert(make_pair(i.first, i.second->clone()));
	for (const auto& i : from.eigv()) {
		this->eigv_.push_back(i->clone());
	}
	desc_ = from.desc();
}

pset&
pset::
operator=(const pset& rhs) {
// deep copy assignment operator
	Trace trc(1,"pset::operator= ",rhs.desc());
	
	if (this == &rhs) {
		return *this;
	}
	// delete all contained Par
	for (auto pp : this->pmap_)
		delete pp.second;
	for (auto pp : this->eigv_)
		delete pp;
	this->clear();

	for (const auto& i : rhs.pmap())
		this->pmap_.insert(make_pair(i.first, i.second->clone()));
	for (const auto& i : rhs.eigv())
		this->eigv_.push_back(i->clone());

	this->desc(rhs.desc());
	return *this;
}

bool
pset::
update (const pset& pl) {
// Set my param advalues to the advalues in corresponding (same name)
// parameters in "pl", but only if the pl parameter was not fixed by
// a preference.
// Returns true if all my parameters were set; false if there were
// parameters in pl that are not in mine.
	Trace trc(2,"pset::update(pset)", desc());
	bool rval{true};
	
	// pmap
	for (const auto& i : pl.pmap()) {
		auto j = this->pmap_.find(i.first);
		if (j != this->pmap_.end()) {
			Par* pp = j->second;
			if (pp->is_fixed()) {
				// do not change fixed Par if set by user
				if (pp->pref == 1)
					continue;
				else
					pp->valuef(i.second->value());
			} else
				pp->update(*i.second);
		} else {
			rval = false;
		}
	}
	// eigv
	//!! assert(pl.eigv().size() == this->eigv().size());
	if (pl.eigv().size() != this->eigv().size())
		throw runtime_error(vastr("attempt to update ",eigv().size(),
			"-eigenvector with ",pl.eigv().size(),"-eigenvector"));
	for (size_t i=0; i<eigv_.size(); i++)
		eigv_[i]->update(*pl.eigv()[i]);

	return rval;
}

bool
pset::
setValues (const pset& pl) {
// Set my param values (not advalues) to the values in corresponding
// (same name) parameters in "pl", even if the parameter is Fixed.
// The diff between this and update() is 1) checking for Fixed
// and 2) update() changes advalue, this function only changes value()
	Trace trc(2,"pset::setValues ",this->desc()," from ",pl.desc());
	bool rval = true;

	for (const auto& par : pl.pmap()) {
		const Par* pi = par.second;
		trc.dprintn("working on \"",pi->summary(),"\"...");
		Par* pp = this->findp(pi->name);
		if (pp) {
			trc.dprint("found in list: set value to ",pi->stringPreValue());
			pp->valuef(pi->valuef());
		} else {
			trc.dprint("not in this param list");
			rval = false;
		}
	}

	// eigv
	assert(pl.eigv().size() == this->eigv().size());
	for (size_t i=0; i<eigv_.size(); i++)
		eigv_[i]->valuef(pl.eigv()[i]->value());

	return rval;
}

void
pset::
freshen() {
// set all Par->fresh to true to avoid potential errors caused
// by stale derived Par. USE WITH CAUTION
	for (auto& par : this->pmap_) {
		Par* pi = par.second;
		pi->fresh = true;
	}
}


// monitoring calls to findp: save names in dependents, return
// when monitoring is turned off
// monitor():   turn on monitoring by adding a vector<string> to monitors
// monitor_add(string nm):  add nm to all vectors on the stack if not already in
// rotinom():  turn off and return the dependent parameter names
void
pset::
monitor() {
// put a new vector<string> on the stack
	Trace trc(2,"monitor");
	vector<string> empty;
	monitors.push_back(empty);
	trc.dprint("now have ",monitors.size()," monitors");
}

void
pset::
monitor_add(const string& nm) {
// add nm to all vectors in the monitors stack if not already there
	Trace trc(2,"monitor_add ",nm);
	for (auto& j : monitors) {
		if (std::find(j.begin(), j.end(), nm) == j.end()) {
			j.push_back(nm);
			trc.dprint("added ",nm);
		}
	}
}

vector<string>
pset::
rotinom() {
// pop the top vector off the stack, return its list
	Trace trc(2,"rotinom");
	vector<string> rval(monitors.back());
	monitors.pop_back();
	trc.dprint("now have ",monitors.size()," monitors, returning ",rval);
	return rval;
}

vector<Par*>
pset::
allpar() {
// XXX same as all()?
// return a vector of pointers to all parameters in this pset
	vector<Par*> rval;
	for (auto par : pmap_)
		rval.push_back(par.second);
	return rval;
}

// various functions for finding a parameter
Par*
pset::
findp(const std::string& nm) {
	// cast this to const, call the const version
	Par* rval = const_cast<Par*>(static_cast<const pset&>(*this).findp(nm));

	// add it to the monitored names
	if (rval != nullptr && this->monitoring())
		this->monitor_add(nm);
	return rval;
}

const Par*
pset::
findp(const std::string& nm) const {
// Return a pointer to parameter named "nm" or nullptr
	Trace trc(3,"findp ",nm);

	const Par* rval{nullptr};

	// search the maps
	auto p = pmap_.find(nm);
	if (p == pmap_.end()) {
		for (auto i : eigv_) {
			if (i->name == nm) {
				rval = i;
				break;
			}
		}
	} else
		rval = p->second;

	trc.dprint("returning ",rval==nullptr?"nullptr":rval->name);
	return rval;
}

Par*
pset::
findg(const std::string& nm) {
// find or throw exception
	Par* rval = findp(nm);
	if (rval == nullptr)
		throw runtime_error(vastr(nm, " is not a parameter"));
	else
		return rval;
}

std::vector<Par*>
pset::
findRe(const std::string& pattern) {
// find all Par with names matching an ECMA regular expression,
// ignoring case
	regex re(pattern, regex_constants::icase);
	std::vector<Par*> rval;
	for (auto& par : this->pmap_) {
		Par* pp = par.second;
		if (regex_match(pp->name, re)) {
			// names only get added to monitor list if
			// they already exist
			this->monitor_add(pp->name);
			rval.push_back(pp);
		}
	}
	// maybe an eigenvector name
	for (auto& pp : this->eigv_) {
		if (regex_match(pp->name, re)) {
			// names only get added to monitor list if
			// they already exist
			this->monitor_add(pp->name);
			rval.push_back(pp);
		}
	}
	return rval;
}

Par*
pset::
findic(const std::string& name) {
// find a Par with "name" ignoring case, ignore eigenvector
	regex re(name, regex_constants::icase);
	Par* rval{nullptr};
	for (auto& par : this->pmap_) {
		Par* pp = par.second;
		if (regex_match(pp->name, re)) {
			// names only get added to monitor list if
			// they already exist, unlike findp() and findg()
			this->monitor_add(pp->name);
			return pp;
		}
	}
	return rval;
}

void
pset::
interp(std::pair<size_t, double> loc) {
	Trace trc(2,"pset interp");
// interpolate each Par's solns:
//   t = solns[start] + t*(solns[start+1]-solns[start])
// set current value to t
	for (auto p : this->pmap_) {
		Par* pp = p.second;
		pp->interp(loc);
	}
	for (auto e : this->eigv_) {
		e->interp(loc);
	}
}

Par*
pset::
add(const Par* newpar) {
/*------------------------------------------------------------------
 * Add a *copy* of a Par to the parameter list if it does not
 * already exist. Return a pointer to which ever ends up on the
 * parameter list.
 *
 * If the parameter already exists, certain members of the
 * Param will be copied to the existing parameter if they
 * do not exist there also:
 *		desc
 *		altval
 *		min
 *		max
 * value, min, and max are assumed to be in presentation units.
 *------------------------------------------------------------------*/
	Trace trc(2,"pset::add");
	Par* rval{nullptr};
	std::string newconv;
	std::string oldconv;
	std::string exc;

	assert(newpar != nullptr);
	trc.dprint("adding \"",*newpar,"\" to pset ",this->desc());

	// Search for the parameter in existing list
	// If found, take stuff from the new parameter
	Par* existing{nullptr};
	existing = findp(newpar->name);

	if (existing == nullptr) {
		// insert a copy of the input parameter
		rval = newpar->clone();
		if (rval->is_eigv())
			this->eigv_.push_back(dynamic_cast<Eigpar*>(rval));
		else
			this->pmap_.insert(make_pair(rval->name, rval));
		trc.dprint("added ",rval->summary(),", now have ",this->size());
	} else if (existing == newpar) {
		trc.dprint("input and existing ptrs identical: no data transfer");
	} else {
		existing->upgrade(newpar);
	}
	// return a pointer to the Par in "this"
	rval = this->findp(newpar->name);
	return rval;
}

Par*
pset::
replace(const Par* newpar) {
/*------------------------------------------------------------------
 * Add a *copy* of a Par to the parameter list if it does not
 * already exist, replace the existing parameter if it does.
 * Return a pointer to which ever ends up on the parameter list.
 *------------------------------------------------------------------*/
	Trace trc(2,"pset::replace");
	Par* rval{nullptr};

	assert(newpar != nullptr);
	trc.dprint("adding \"",*newpar,"\" to pset ",this->desc());

	// Search for the parameter in existing list
	// If found, take stuff from the new parameter
	Par* oldpar = this->findp(newpar->name);
	if (oldpar != nullptr) {
		*oldpar = *newpar;
		rval = oldpar;
	} else {
		// insert a copy of the input parameter
		rval = newpar->clone();
		if (rval->is_eigv())
			this->eigv_.push_back(dynamic_cast<Eigpar*>(rval));
		else
			this->pmap_.insert(make_pair(rval->name, rval));
		trc.dprint("added ",rval->summary(),", now have ",this->size());
	}

	return rval;
}

#ifdef NEVER // use inrange
string
pset::
inRange() {
// check all parameters except "rf" for in-range,
// return a description of the first out-of-range
	Trace trc(2,"pset::inRange");
	string rval;
	int sigfig{4};
	for (auto& i : this->pmap_) {
		Par* pp = i.second;
		if (pp->name == "rf" && !pp->is_indep())
			continue;
		if (!pp->inRange(sigfig)) {
			rval = vastr(pp->name, " is out of the range ",
				pp->stringLimits(), ": ", pp->prevalue());
			break;
		}
	}
	trc.dprint("returning \"",rval,"\"");
	return rval;
}
#else // NEVER // use inrange

vector<string>
pset::
inrange() {
// Check all parameter limits, return the names of
// all parameters that are out-of-range
	Trace trc(2,"pset::inrange");
	vector<string> rval;
	int sigfig{4};
	for (auto& i : this->pmap_) {
		Par* pp = i.second;
		if (!pp->inRange(sigfig))
			rval.push_back(pp->name);
	}
	trc.dprint("returning ",rval.size()," out of range: ",rval);
	return rval;
}
#endif // NEVER // use inrange

std::string
pset::
summary () {
// description + number of parameters
	string rval = vastr(this->desc(), "  ",this->size()," parameters");
	return rval;
}

std::string
pset::
summary (std::vector<std::string> const& toprint, int maxchar) {
// One-line summary of a parameter list in "presentation" units.
// Only include parameters in "toprint"; if toprint is empty
// just return the desc and number of Par.
// Each parameter gets a min width of 14: 2 spaces plus a "setw(12)"
// which sets the min (right-adjusted) width. Note this does not limit
// the width of the gc print - only pagewidth does that
// The width of the line is determined by page_width() which
// can be set with env variable PAGEWIDTH or settings{pagewidth=x}
//
// If a parameter is out-of-range print an asterisk at it's right edge
	Trace trc(3,"pset::summary");
	std::ostringstream os;
	int pw = page_width();
	bool printev{false};

	trc.dprint(this->size()," par, ",toprint.size()," printed, pagewidth ",pw);

	if (toprint.empty()) {
		os << desc() << " " << this->size() << " parameters";
		return os.str();
	}

	for (const auto& nm : toprint) {
		trc.dprint("toprint: ", nm);
		if (nm == "gc" || nm == "ev" || nm == "eigenvector") {
			printev = true;
			continue;
		}
		const Par* pp = findp(nm);
		if (pp) {
			os << "  " << std::setw(12) << pp->stringPreValue();
			if (!pp->inRange(8))
				os << '*';
		}
	}

	// print the eigenvector if requested and if there is enough room
	if (printev) {
		int wid = maxchar - os.str().size();
		static bool visit{false};
		if (wid < 8 && !visit) {
			flaps::warning("the max char (",maxchar,") is insufficient to print",
					" a summary of the eigenvector");
			visit = true;
		} else {
			// normalize it first so that the largest element is 1
			vector<double> ev = this->getev();  // make a copy
			blas::normalize_inf(ev.size()/2, (complex<double>*)&ev[0], 1, true, 1.0);
			os << "   " << flaps::summarize(ev, wid, true);
		}
	}
	trc.dprint("returning ",os.str());
	return os.str();
}

string
pset::
summary_solns (vector<string> const& toprint, const string& header,
		vector<string> const& leaders, bool firstlast) {
// Summary of a parameter list which has the same number of solns for
// each parameter. Returns a one-line summary for each value in the "solns"
// array (if the "solns" array is empty the current values are used).
// Only include parameters in "toprint". If "leader" is not empty, start
// each line with it. If "firstlast" is true only print the first & last solns
// with an ellipsis in between. If "header" is not empty, print it.
// Note: this function modifies "this" so it cannot be const
	Trace trc(2,"pset::summary_solns");
	std::ostringstream os;
	size_t nval = nsolns();
	size_t i;

	trc.dprint(this->size()," par, ",nval," solns, ",toprint.size()," to print");

	// determine the widest leader so we know how to offset the
	// header column labels to line up with the values, and
	// how wide to make the leader field so all columns line up...
	// If "leader" is empty use pset::desc()
	size_t lwidth = 0;
	if (leaders.empty()) {
		if (nval == 0) {
			lwidth = desc().size()+3;
		} else {
			for (i=0; i<nval; i++) {
				if (firstlast && i != 0 && i != nval-1)
					continue;
				os.str("");
				os << std::setw(7) << i+1 << ")  ";
				lwidth = std::max(lwidth, os.str().size());
			}
		}
	} else {
		for (auto& s : leaders)
			lwidth = std::max(lwidth, s.size());
	}

	if (!header.empty()) {
		os.str("");
		os << string(lwidth, ' ') << header << endl;
	}

	// Print the summaries for either the first and last,
	// or for all solns (firstlast=false)
	// If nval is zero only print the current values
	if (nval == 0) {
		if (leaders.empty()) {
			os << std::setw(lwidth) << desc() << ":  "
				<< summary(toprint) << std::endl;
		} else {
			os << leaders[0] << summary(toprint) << std::endl;
		}
	} else {
		for (i=0; i<nval; i++) {
			if (firstlast) {
				if (i == 1 && nval > 2)
					os << "     ...\n";
				if (i != 0 && i != nval-1)
					continue;
			}
			update_from_solns(i);
			if (leaders.size() < i+1) {
				os << std::setw(7) << i+1 << ")  " << summary(toprint) << std::endl;
			} else {
				os << leaders[i] << summary(toprint) << std::endl;
			}
		}
	}

	// If the return string contains any astrisks print a footnote
	// explaining what this means
	if (os.str().find('*') != std::string::npos)
		os << "  *out of range\n";

	return os.str();
} // summary_solns

std::string
pset::
display(std::string const& path, std::string const& title) const {
	Trace trc(3,"pset::display");
	std::string rval;
	std::ofstream file;
	std::ostream* out = &std::cout;
	if (!path.empty()) {
		if (path == "cout")
			out = &std::cout;
		else if (path == "cerr" || path == "debug")
			out = &std::cerr;
		else {
			file.open(path.c_str());
			out = &file;
			rval = path;
		}
	}
	(*out) << *this << std::endl;
	return rval;
}

std::string
pset::
summaryTable () const {
	std::ostringstream os;
	for (const auto& i : this->pmap_)
		os << i.first << "  " << i.second->stringPreValue() << endl;

	return os.str();
}

void
pset::
set_adpar_derivs() {
// set the appropriate AD derivative to 1.0 for each AD parameter
	Trace trc(2,"set_adpar_derivs ",this->desc());
	vector<string> const& adparnames = Ad::adnames();

	// for each ad par in adparnames find it in this, set deriv i to 1
	for (size_t i=0; i<adparnames.size(); i++) {
		trc.dprint("setting der ",i," for ",adparnames[i]);
		Par* pp = this->findp(adparnames[i]);
		assert(pp != nullptr);
		pp->deriv(i, 1.0);
		// check that all other deriv are zero
		for (size_t j=0; j<adparnames.size(); j++) {
			if (j != i && pp->deriv(j) != 0.0)
				throw runtime_error(vastr("deriv ",j," of ",pp->name," is ",pp->deriv(j)));
		}
	}
}

void
pset::
touch() {
// Mark all dependent (derived) parameters out-of-date (!fresh) so that
// they will be evaluated the next time eval() is called.
// Also, Constant parameters are changed to Derived, evaluated,
// and checked for Constant-ness.
// It is important to call this function after changing Fixed parameter values.
	Trace trc(3,"pset::touch");
	for (auto& i : this->pmap_)
		i.second->fresh = false;
}

void
pset::
eval() {
// evaluate all derived parameters
	Trace trc(2,"pset::eval ", this->desc());

	// put a 1 in the appropriate AD spots
	this->set_adpar_derivs();

	// evaluate all derived parameters
	vector<Par*> derived = this->get_derived();
	for (auto i : derived)
		i->eval(*this);
}

//------------------------------------------------------------------
// Functions for handling "solns" member
//------------------------------------------------------------------
size_t
pset::
nsolns() const {
// Returns the number of elements in the "solns" member of
// the contained parameter.
// Throws a runtime_error exception if not all parameters have the same
// number of parameters
	Trace trc(2,"pset::nsolns()");
	size_t rval{};
	int nvi{-1};
	std::string first_name;

	for (const auto& k : this->pmap_) {
		if (nvi == -1) {
			rval = k.second->nsolns();
			first_name = k.first;
		}
		nvi = k.second->nsolns();

		if (nvi != (int)rval)
			throw runtime_error(vastr("nsolns() expecting a parameter list with the ",
				"same number of values, but ", first_name,
				" has ", rval, " and ", k.first, " has ", nvi));
	}

	// eigv
	for (const auto& k : this->eigv_) {
		nvi = k->nsolns();
		if (nvi != (int)rval)
			throw runtime_error(vastr("nsolns() expecting a parameter list with the ",
				"same number of values, but ", first_name,
				" has ", rval, " and ", k->name, " has ", nvi));
	}

	trc.dprint("returning ", rval);
	return rval;
}

void
pset::
clear_solns() {
// empty the "solns" member of each parameter
	for (auto& i : this->pmap_)
		i.second->clear_solns();
	for (auto& i : this->eigv_)
		i->clear_solns();
}

void
pset::
update_from_solns(size_t index) {
// Set the values of all (even Fixed) parameters in this pset
// from element "index" of the "solns" arrays of my Par
	for (auto& i : this->pmap_)
		i.second->update_from_solns(index);
	for (auto& i : this->eigv_)
		i->update_from_solns(index);
}

bool
pset::
verify (const pset& pl) {
// check to see that "pl" has the same parameters and in
// the same order as "this"
	if (pl.size() == this->size()) {
		return true;
	} else
		throw runtime_error(vastr("pset ",pl.desc()," has ",pl.size(),
			" parameters, but ",this->desc()," has ",this->size()));
}

//------------------------------------------------------------------
// functions for getting lists of parameters of various states
//------------------------------------------------------------------

vector<Par*>
pset::
all () {
// Returns pointers to all parameters except nostate
	Trace trc(3,"pset::all ", this->desc());
	vector<Par*> rval;

	for (auto& i : this->pmap_) {
		if (!i.second->is_nostate())
			rval.push_back(i.second);
	}
	for (auto& i : this->eigv_)
		rval.push_back(i);
	trc.dprint(rval.size(), " all parameters");
	return rval;
}

vector<Par*>
pset::
get_nostate () {
	Trace trc(3,"pset::get_nostate ", this->desc());
	vector<Par*> rval;

	for (auto& i : this->pmap_) {
		if (i.second->is_nostate())
			rval.push_back(i.second);
	}
	trc.dprint(rval.size(), " no-state parameters");
	return rval;
}

vector<Par*>
pset::
get_indep () {
// Make up a list of all indep variables in this pset
// Create the list from scratch each time in case the state
// of variables have changed
	Trace trc(3,"pset::get_indep ", this->desc());
	vector<Par*> indep;

	for (auto& i : this->pmap_) {
		if (i.second->is_indep())
			indep.push_back(i.second);
	}
	trc.dprint(indep.size(), " independent parameters");
	return indep;
}

vector<Par*>
pset::
get_derived () {
	Trace trc(3,"pset::get_derived ", this->desc());
	vector<Par*> rval;
	for (auto& i : this->pmap_) {
		trc.dprint(i.first," is ",i.second->state_desc());
		if (i.second->is_derived())
			rval.push_back(i.second);
	}
	trc.dprint(rval.size(), " derived parameters");
	return rval;
}

vector<Par*>
pset::
get_fixed () {
	Trace trc(3,"pset::get_fixed ", this->desc());
	vector<Par*> rval;
	for (auto& i : this->pmap_) {
		if (i.second->is_fixed())
			rval.push_back(i.second);
	}
	trc.dprint(rval.size(), " fixed parameters");
	return rval;
}

vector<Par*>
pset::
get_aux () {
	Trace trc(3,"pset::get_aux(const) ", this->desc());
	vector<Par*> rval;

	for (auto& i : this->pmap_) {
		if (i.second->is_aux())
			rval.push_back(i.second);
	}
	trc.dprint(rval.size(), " Aux parameters");
	return rval;
}

//------------------------------------------------------------------
// Eigenvector access functions
//------------------------------------------------------------------

void
pset::
getadev(vector<complex<Ad>>& rval) {
// Fill input "rval" with the current advalues of the
// eigenvector parameters.
	Trace trc(2,"pset::getadev");

	vector<Eigpar*>& evp = this->eigv();
	trc.dprint(evp.size(), " eigenvector components");
	if (evp.empty())
		return;

	size_t nrev = evp.size(); // # real ev
	size_t ncev = nrev/2;     // # complex ev

	// re-allocate rval if necessary
	if (rval.size() != ncev) {
		trc.dprint("re-allocated input from ",rval.size()," to ",ncev);
		rval = vector<complex<Ad>>(ncev);
	}
	// copy each separately - not contiguous
	for (size_t i=0; i<ncev; i++) {
		rval[i].real(evp[2*i]->advalue());
		rval[i].imag(evp[2*i+1]->advalue());
	}
}

std::vector<double>
pset::
getev() {
// Returns the current eigenvector for this pset as a vector of doubles;
// if the autmatic-differentiation representation is needed, use getadev()
// This will be the real representation of complex gc.
	vector<Eigpar*>& gcp = eigv();
	size_t nrgc = gcp.size();
	vector<double> rval(nrgc);

	// then fill rval
	for (size_t i=0; i<nrgc; i++) {
		rval[i] = gcp[i]->value();
	}
	// should use move()
	return rval;
}

vector<complex<double>>
pset::
get_eigenvector(size_t idx) {
// returns a complex vector, the idx (0b) solns of all eigenvector
// components
	// check that idx is in range
	size_t naltv = this->nsolns();
	if (idx >= naltv)
		throw runtime_error(vastr("attempt to extract eigenvector ",idx,
				" but there are only ",naltv," solns in this pset"));

	// get pointers to all eigenvector parameters
	vector<Eigpar*>& evp = eigv();
	size_t ncev = evp.size()/2;  // number of complex ev components
	if (ncev == 0)
		throw runtime_error(vastr(this->desc()," has no eigenvector components"));
	vector<complex<double>> rval(ncev);
	for (size_t i=0; i<ncev; i++) {
		double xr = evp[2*i]->solns.at(idx);
		double xi = evp[2*i+1]->solns.at(idx);
		rval[i] = complex<double>(xr,xi);
	}
	// should use move()
	return rval;
}

std::vector<std::vector<double> >
pset::
get_evmatrix() {
// XXX is this used?
// Gather *all* values of all eigenvector parameters and return them
// in a vector of vectors of doubles, so the return value
// is the real representation of complex eigenvector matrix. If the
// return value is rval, then rval[i][j] is the real (i%2=0)
// or imag (i%2=1) part of the the (i/2)th gc in the jth mode.
// XXX why not just return a vector<complex<double>>?
	Trace trc(1,"get_evmatrix");
	std::vector<std::vector<double> > rval;

	// get pointers to all eigenvector parameters
	vector<Eigpar*>& evp = eigv();

	size_t nrev = evp.size();  // number of real eigenvector components

	if (nrev == 0) {
		trc.dprint("returning empty: no eigenvector parameters");
		return rval;
	}
	size_t nc = this->nsolns();

	// rval is a set of nc nrev-vectors of doubles
	vector<double> col(nrev, 0.0);
	for (size_t j=0; j<nc; j++) {
		for (size_t i=0; i<nrev; i++) {
			col[i] = evp[i]->solns.at(j);
		}
		rval.push_back(col);
	}
	return rval;
}

void
pset::
updateEv(std::vector<double> const& ev) {
// Eigenvector components are kept as parameters: the real
// and imag parts of each component are kept in a separate parameter
// with a name like "ev[1].real" or "ev[22].imag".
// The values are changed like any other parameter, but note that
// because the input ev are doubles, not Ads, the updating is
// done with Par::value(x) instead of Par::advalue(x).
// This convenience function allows the ev parameters to be updated
// from a vector, e.g. an eigenvector from a start-point calculation
	updateEv(ev.size(), &ev[0]);
}

void
pset::
updateEv(size_t n, double const* ev) {
	Trace trc(2,"pset::updateEv", " n ", n);

	// just cycle through pset::eigv_
	// Ad::value(x) just sets the value, not derivatives
	for (size_t i=0; i<n; i++) {
		eigv_[i]->value(ev[i]);
	}
}

//------------------------------------------------------------------
// end of eigenvector access functions
//------------------------------------------------------------------

void
pset::
back_solns(pset& from) {
// push the current values of parameters in list "from"
// onto the back-end of my "solns"

	assert(this->verify(from));
	auto j = this->pmap_.begin();
	for (auto& i : from.pmap()) {
		j->second->solns.push_back(i.second->value());
		j++;
	}
	for (size_t i=0; i<eigv_.size(); i++)
		eigv_[i]->solns.push_back(from.eigv()[i]->value());
}

void
pset::
front_solns(pset& from) {
// push the current values of parameters in list "from"
// onto the front-end of my "solns"

	assert(this->verify(from));
	auto j = from.pmap().cbegin();
	for (auto& i : this->pmap_) {
		i.second->solns.push_front(j->second->value());
		j++;
	}
	for (size_t i=0; i<eigv_.size(); i++)
		eigv_[i]->solns.push_front(from.eigv()[i]->value());
}

// store and fetch for pset
pset*
pset::
fetch (std::string const& mid) {
// Fetch a set of parameters previously stored with pset::store.
// Throws a runtime_error exception if the pset does not exist
	Trace trc(1,"pset::fetch ",mid);
	pset* rval{nullptr};

	// a pset is stored in a file called "mid" serialized
	// with Receiver::serialize(map<string,Par>)
	Receiver file(mid);  // throws exception if it does not exist
	if (!file.good() || file.eof()) {
		string exc{vastr(mid," does not exist")};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	Fio* op = Fio::get(file);
	rval = dynamic_cast<pset*>(op);
	if (rval == nullptr) {
		string exc{vastr("reading ", mid, ", expected pset, got ",op->vid())};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	trc.dprint("returning ",rval->size()," parameters in ",rval->desc());
	return rval;
}

void
pset::
store (std::string const& mid) {
// Store a parameter set in a file named "mid"
	Trace trc(1,"pset::store ", desc());

	trc.dprint("mid ",mid,", desc ",desc(),", ",this->nsolns()," solns");

	if (this->empty()) {
		string exc{vastr("attempt to store empty pset \"", desc(),'\"')};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	Sender file(mid);
	Fio::put(file, *this);
}


// plot "toPlot" parameters to an "apf" plotfile
void
pset::
plot (std::string const& path, std::string const& title, std::string const& cid,
		vector<string> const& toPlot, bool append) {
	Trace trc(1,"pset::plot");
	size_t nv = nsolns();

	trc.dprint("path \"",path,"\", title \"",title,"\", cid \"",cid,"\"");
	trc.dprint(nv," data points");

	if (nv == 0) {
		return;
	}

	// only one curve:
	Curve* curve = new Curve(this->desc(), cid, "");
	curve->params.clear();

	size_t j;
	std::ostringstream os;

	// make up a list of parameters to plot; plot all parameters
	// (except eigenvector components) if toPlot is empty
	if (toPlot.empty()) {
		for (auto& par : this->pmap_) {
			if (!par.second->is_eigv()) {
				curve->params.add(par.second);
			}
		}
	} else {
		for (j=0; j < toPlot.size(); j++) {
			Par* pp = this->findp(toPlot[j]);
			if (!pp) {
				throw runtime_error(vastr("parameter \"",toPlot[j],"\" was requested ",
					"to plot but is not in the parameter list:\n",*this, '\n'));
			}
			curve->params.add(pp);
		}
	}

	// make sure the filename ends in .apf
	string file{path};
	if (rsubstr(path, 4) != ".apf")
		file += ".apf";
	// write the Apf file
	vector<Curve*> curves;
	curves.push_back(curve);
	vector<string> toplot;	// XXX just use toPlot? or re-write using Apf::exporter(vector)
	Apf::exporter(curves, toplot, file, append);
}

void
pset::
plot() {
// simple plotting: all parameters, use desc() as path, title, and cid
	string id = replace_char(this->desc(), ' ', '_');
	if (id.empty())
		id = "unknown";
	vector<string> toplot;
	this->plot(id, id, id, toplot, false);
}

Ad
pset::
absgcv(int gcno) {
// returns the absolute value of generalized coordinate gcno (1b)
// If the parameter absgcxxx does not exist, create it and add it to this
	Trace trc(2,"pset::absgcv ",gcno);
	Par* absgc_par = this->findp(vastr("absgc",gcno));
	Ad rval;
	if (absgc_par == nullptr) {
		// only create it if ev(gcno).real/imag exist - may not have been
		// created yet
		Par* er = this->findp(evname(gcno,false));
		Par* ei = this->findp(evname(gcno,true));
		if (er != nullptr && ei != nullptr) {
			string defn{vastr("absgc",gcno,"(abs(gc",gcno,"))")};
			string eqn{vastr("gcnorm*sqrt(ev",gcno,".real^2 + ev",
							gcno,".imag^2)")};
			Par newpar(defn, eqn);
			absgc_par = this->add(&newpar);
		} else
			return Ad(0.0);
	}
	rval = absgc_par->eval(*this);
	trc.dprint("returning ",rval);
	return rval;
}

Ad
pset::
absgcp(int gc1, int gc2) {
// returns the absolute value of the difference between 2 generalized
// coordinates, gc1 and gc2 (1b)
	Trace trc(2,"pset::absgcp ",gc1,", ",gc2);

	// findg will throw an exception if the parameter doesn't exist
	Par* ev1r = this->findg(evname(gc1,false));
	Par* ev1i = this->findg(evname(gc1,true));
	Par* ev2r = this->findg(evname(gc2,false));
	Par* ev2i = this->findg(evname(gc2,true));
	Ad gcnorm = this->parval("gcnorm");
	Ad diffr = ev1r->advalue() - ev2r->advalue();
	Ad diffi = ev1i->advalue() - ev2i->advalue();
	Ad rval = gcnorm*sqrt(diffr*diffr + diffi*diffi);
	trc.dprint("returning ",rval);
	return rval;
}

Ad
pset::
parval(std::string const& name, bool inrange) {
//------------------------------------------------------------------
// Returns the Ad value of parameter "name" from this pset, evaluating
// it first if necessary. If "inrange" is true the parameter limit with
// zero Ad derivatives is returned if the parameter value is out of limits.
// Eigenvector components are accessed with names like "ev32.real" and
// "ev32.imag"; generalized-coordinate components are accessed with names
// like "gc32.real" and "gc32.imag".
//------------------------------------------------------------------
	Trace trc(2,"parval ",desc(),'(',name,')');
	Ad rval{0.0};
	Par* pp{nullptr};
	// GC are not stored so if we are called with, e.g. "gc1.real"
	// we must get the eigenvector component and multiply it by gcnorm
	int rgcno = gcparse(name);
	if (rgcno >= 0) {
		pp = findp(revname(rgcno));
		if (pp) {
			Ad gcnorm = this->parval("gcnorm");
			rval = gcnorm*pp->advalue();
			trc.dprint("scaled eigenvector by ",gcnorm," to create gc");
		} else if (!this->monitoring()) {
		 	string exc{vastr(name," is not a parameter (parval)")};
		 	trc.dprint("throwing exception: ",exc);
		 	throw runtime_error(exc);
		}
	} else {
		pp = findp(name);
		if (pp != nullptr) {
			if (!pp->fresh) {
				pp->eval(*this);
				pp->fresh = true;
			}
			if (inrange)
				rval = pp->iradvalue();
			else
				rval = pp->advalue();
		} else if (this->monitoring()) {
			// absgcxxx parameters may not have been created yet: create them
			int cgcno;
			if (is_absgcname(name, cgcno))
				rval = this->absgcv(cgcno);
		} else {
			string exc{vastr(name," is not a parameter (parval)")};
			trc.dprint("throwing exception: ",exc);
			throw runtime_error(exc);
		}
	}
	trc.dprint("returning ",name," = ",rval);
	// assuming RVO is applied here
	return rval;
}

void
pset::
setpar(Ad const& x, std::string const& name) {
// set the value of a parameter from an Ad in presentation units;
// if the parameter does not exist it is created and a warning
// is printed so that you can just call this function from within
// a function without first defining the parameter.
// Intended to be called from user-written C++ functions
	Trace trc(2,"setpar");
	Par* pp = findp(name);

	// if the parameter does not exist create it, add it to this
	// pset and the global pset, and print a message
	if (pp == nullptr) {
		string lhs = vastr(name,"(user-defined)");
		string rhs = vastr(x.value());
		Par newp(lhs, rhs);
		newp.solns.clear();
		//!! this->pmap_.insert(make_pair(name, newp)); // add to this...
		pp = this->add(&newp);
		pp->set_aux();
		pp->advalue(x);
		gpset::get().add(pp);  // ...and the gpset
		flaps::info("creating new parameter \"", pp->summary(), "\"");
	}
	// convert the value to internal units
	Ad val = x/pp->convFactor();
	pp->advalue(val);
	trc.dprint("set ",pp->name," to ",val);
}

Ad
pset::
par_or_double(const string& str) {
// x is either a:
//   - parameter name: return its advalue from "this"
//   - a double: convert to Ad
	double t;
	if (str2double(str, t))
		return Ad(t);
	else
		return parval(str);
}

//------------------------------------------------------------------
// end of pset implementation
//------------------------------------------------------------------

//------------------------------------------------------------------
// gpset implementation
//------------------------------------------------------------------
void
make_gpset(pset& rval) {
// Make the standard builtin parameters with SI EU units and
// either SI or USCS presentation units (PU).
	Trace trc(1,"make_gpset");
	ostringstream os;
	// some conversion factors: check for "units" in the environment
	bool uscs{(getenv("units") != nullptr)};
	string pressure_conv{"<Pa>"};
	string vel_conv{"<m/s>"};
	string vel_limits{"[0:800]"};
	if (uscs) {
		pressure_conv = vastr("<",Pa2psf," psf/Pa>");
		vel_conv = vastr("<",mps2knot," knot/(m/s)>");
		vel_limits = "[0:1200]";
	}
	string pressure_limits{"[0:10e+7]"};
	// Note: make_pair is an STL function equivalent to pair<string,Par*>

	// alt: altitude is always either Fixed or Indep; default is zero
	os.str("");
	os << "alt(Altitude)[" << atmos::altmin() << ':' << atmos::altmax() << ']';
	if (uscs)
		os << "<" << m2ft << " ft/m>";
	else
		os << "<m>";
	Par* alt = new Par(os.str(), "0");
	rval.pmap().insert(make_pair("alt", alt));

	// cdpress: calibrated dynamic pressure (impact pressure)
	os.str("");
	os << "cdpress(Calib Dyn Press)" << pressure_limits << pressure_conv;
	Par* cdpress = new Par(os.str(), {"atmos::cdpress(vtas,alt)"});
	rval.pmap().insert(make_pair("cdpress", cdpress));

	// dpress: dynamic pressure
	os.str("");
	os << "dpress(Dynamic Pressure)" << pressure_limits << pressure_conv;
	Par* dpress = new Par(os.str(), {"rho*(vtas^2)/2","rhoref*(veas^2)/2"});
	rval.pmap().insert(make_pair("dpress", dpress));

	// freq: give it a small minimumn (not zero) to avoid, e.g.
	// problems computing sigma from growth
	os.str("");
	os << "freq(Frequency)[.01:100]<radps2Hz>";
	Par* freq = new Par(os.str(), {"rf*vtas","2.0*sigma/growth"});
	rval.pmap().insert(make_pair("freq", freq));

	// gcnorm: 2-norm of the gc - no equations, only Fixed or Indep
	Par* gcnorm = new Par("gcnorm(Gc 2-norm)", "0");
	rval.pmap().insert(make_pair("gcnorm", gcnorm));

	// growth factor: dimensionless
	Par* growth = new Par("growth(Growth Factor)",{"2.0*sigma/freq"});
	rval.pmap().insert(make_pair("growth", growth));

	// Mach number: I've left off the eqn in terms of vcas
	Par* mach = new Par("mach(Mach Number)[0:1.5]",{"vtas/vsound"});
	rval.pmap().insert(make_pair("mach", mach));

	// rf: the characteristic length is assumed to be 1 m 
	Par* rf = new Par("rf(Reduced Frequency)[0:]",{"freq/vtas"});
	rval.pmap().insert(make_pair("rf", rf));

	// rho: default density is atmos::rho(default_alt)
	os.str("");
	os << "rho(Std Atm Density)[0:10]";
	if (uscs)
		os << "<" << kgpm32slugpft3 << " (slug/ft^3)/(kg/m^3)>";
	else
		os << "<kg/m^3>";
	Par* rho = new Par(os.str(), {"atmos::rho(alt)"});
	rval.pmap().insert(make_pair("rho", rho));

	// rsf: the characteristic length is assumed to be 1 m 
	Par* rsf = new Par("rsf(Reduced Stability Factor)", {"sigma/vtas"});
	rval.pmap().insert(make_pair("rsf", rsf));

	// sdamp: structural damping - no eqn, always either Indep or Fixed
	Par* sdamp = new Par("sdamp(Structural Damping)", "0");
	rval.pmap().insert(make_pair("sdamp", sdamp));

	// sigma: stability factor
	Par* sigma = new Par("sigma(Stability Factor))",{"rsf*vtas","growth*freq/2.0"});
	rval.pmap().insert(make_pair("sigma", sigma));

	// spress: static pressure from tables
	os.str("");
	os << "spress(Static Pressure)" << pressure_limits << pressure_conv;
	Par* spress = new Par(os.str(), {"atmos::press(alt)"});
	rval.pmap().insert(make_pair("spress", spress));

	// temp: temperature from tables
	Par* temp = new Par("temp(Static Temperature)[0:690]<K>",
			{"atmos::temp(alt)"});
	rval.pmap().insert(make_pair("temp", temp));

	// tpress: total or stagnation pressure
	os.str("");
	os << "tpress(Stagnation Pressure)" << pressure_conv;
	Par* tpress = new Par(os.str(), {"dpress+spress"});
	rval.pmap().insert(make_pair("tpress", tpress));

	// vcas: calibrated airspeed
	os.str("");
	os << "vcas(Velocity (CAS))" << vel_limits << vel_conv;
	Par* vcas = new Par(os.str(), {"atmos::vcas(vtas,alt)"});
	rval.pmap().insert(make_pair("vcas", vcas));

	// veas: equivalent airspeed
	os.str("");
	os << "veas(Velocity (EAS))" << vel_limits << vel_conv;
	Par* veas = new Par(os.str(),{"sqrt(rho/rhoref)*vtas","sqrt(2.0*dpress/rhoref)"});
	rval.pmap().insert(make_pair("veas", veas));

	// vsound: sonic velocity in the std atmosphere
	os.str("");
	os << "vsound(Sonic Velocity)" << vel_limits << vel_conv;
	Par* vsound = new Par(os.str(), {"atmos::vsonic(alt)"});
	rval.pmap().insert(make_pair("vsound", vsound));

	// vtas: true airspeed
	os.str("");
	os << "vtas(Velocity (TAS))" << vel_limits << vel_conv;
	Par* vtas = new Par(os.str(), {"mach*vsound", "freq/rf", "sqrt(2*dpress/rho)",
			"sqrt(rhoref/rho)*veas"});
	rval.pmap().insert(make_pair("vtas", vtas));

	// all parameters are "fresh" until equations are made
	for (auto pp : rval.pmap())
		pp.second->fresh = true;
	// also each parameter must have its advalue member reallocated after
	// Ad::initialize() i.e. call pset::realloc()

	trc.dprint("returning ",rval.size()," global parameters");

} // make_gpset

gpset::
gpset() {
// gpset one-and-only constructor. First try fetching it; if this
// fails build it from scratch with make_gpset()
	Trace trc(1,"gpset constructor");
	
	try {
		pset* pl{nullptr};
		pl = pset::fetch(persistentName());
		if (pl == nullptr) {
			// fetch returned nullptr so make it from scratch:
			make_gpset(*this);
		} else {
			// set state to nostate, put the fetched parameters into this
			for (auto& pp : pl->pmap()) {
				if (!(pp.second->is_aux()))
					pp.second->set_nostate();
				this->pmap().insert(pp);
			}
			for (auto& ep : pl->eigv()) {
				ep->set_aux();
				this->eigv().push_back(ep);
			}
		}
	// in case fetch throws an exception:
	} catch (runtime_error& s) {
		trc.dprint(s.what(),": calling make_gpset");
		make_gpset(*this);
	} catch (const std::exception& e) {
		trc.dprint(e.what(),": calling make_gpset");
		make_gpset(*this);
	}

	// Give this one-and-only global pset a description
	this->desc("Global Parameters");

	trc.dprint("returning ",this->size()," global param");
} // gpset constructor

//------------------------------------------------------------------
// stuff for assigning equations for each parameter that has candidates in a pset
//------------------------------------------------------------------
void
pset::
equations () {
//------------------------------------------------------------------
// Set the "equation" member of each parameter in this list to
// one of it's candidates
// -----------------------------------------------------------------
	Trace trc(2,"pset::equations");

	// delete equations
	for (auto& pp : this->pmap())
		pp.second->equation.clear();

	for (auto& pi : pmap()) {
		Par* pp = pi.second;
		trc.dprint("\n", separator());
		trc.dprint("working on ",pp->name,", ",pp->candidates.size()," candidates");
		if (pp->is_indep() || pp->is_fixed() || pp->candidates.empty()
			|| !pp->equation.empty()) {
			trc.dprint("skip: state ",pp->state);
			continue;
		}
		try {
			this->choose_eqn(*pp);
		} catch (runtime_error& s) {
			throw runtime_error(vastr("cannot create equation for ",
				pp->name,": ",s.what()));
		}
		// did choose_eqn find an equation?
		if (pp->equation.empty())
			throw runtime_error(vastr(pp->name," is not derivable from the independent "
				"& fixed parameters. It must be set to a constant or "
				"derivable from the standard equation(s)"));
		trc.dprint("\n",separator());
	}
	this->has_eqns(true);
} // equations

int
pset::
choose_eqn (Par& par) {
// Find an acceptable equation for this parameter as a
// member of this parameter list
// Returns: the "work" to evaluate with the chosen evaluator or -1 if
//          no acceptable equation could be found
// Method:
//     for each candidate equation, check each of its dependent parameters;
//     if they are indep or fixed we are done; if they have candidates but no
//     equation we must call this function (recurse) to get an equation for that parameter.
// This function is recursive and not thread-safe!
	Trace trc(1,"choose_eqn ", par.name);
	string besteqn;
	bool best_constant{false};
	static Stack<string> par_stack;	// the names of all parameters this function has
												// been called with recursively

	if (!par.equation.empty()) {
		trc.dprint("quick return: already have equation: ",par.equation);
		return par.work;
	}

	if (par_stack.contains(par.name)) {
		trc.dprint("cycle in dep list for \"",par.name,"\", current stack: ",par_stack);
		return -2;
	}
	par_stack.push(par.name);

	// Check each candidate: check each parameter in the equation;
	// if it is derived and has no equation yet, call its choose_eqn()
	// to get one. If it is constant or independent we are done with that parameter.
	int work{10000};
	for (auto& eqn : par.candidates) {
		int worki{0}; // work for this eqn
		bool is_complex{false};
		bool constant{true};
		trc.dprint("testing ",par.name, " = ",eqn," ----------------");
		// parse the equation, get it's dependent parameters...
		vector<string> dep = adeqn::dependson(*this, eqn, is_complex);
		// ... then get each dependent's dependencies
		for (auto& dj : dep) {
			if (dj == par.name)  // skip this par
				continue;
			Par* pj = this->findp(dj);
			if (pj == nullptr) {
				std::string exc{vastr("unknown parameter: ",dj)};
				trc.dprint("rejecting equation \"",eqn,"\": ",exc);
				worki = -1;
				break;
			}
			trc.dprint(eqn," depends on ",pj->name,": ",pj->state);

			if (!(pj->is_indep() || pj->is_fixed() || pj->is_aux())
					&& !pj->candidates.empty()) {
				if (!pj->equation.empty()) {
					worki += pj->work;
					trc.dprint(pj->name," already has equation with work ", pj->work);
				} else {
					try {
						int workj = this->choose_eqn(*pj);
						if (workj < 0) {
							trc.dprint("rejecting equation \"",eqn,"\"");
							worki = -1;
							break;
						}
						worki += workj;
					} catch (runtime_error& s) {
						trc.dprint("rejecting equation \"",eqn,"\": ",s.what());
						worki = -1;
						break;
					}
				}
			} else {
				// something other than derived: no change to worki
				trc.dprint(pj->name," is ",pj->state_desc());
			}
			// if pj is not constant => this equation is not constant
			if (!pj->is_constant()) {
				constant = false;
				trc.dprint(eqn," is not constant: ",pj->name," is not constant");
			}
		}

		// eqn rejected?
		if (worki < 0) {
			trc.dprint("rejecting eqn \"",eqn,"\": worki ",worki);
			continue;
		}
		// if this equation is constant it's work is zero
		if (constant)
			worki = 0;

		// made it through all dependent parameters?
		if (worki >= 0 && worki < work) {
			work = worki;
			besteqn = eqn;
			best_constant = constant;
			trc.dprint("best eqn so far: ",eqn);
		}
	}

	// give it the best equation and evaluate it
	if (!besteqn.empty()) {
		trc.dprint("best equation: ",besteqn);
		par.equation = besteqn;
		par.work = work;
		// mark this parameter as constant if the equation but evaluate it first
		if (best_constant) {
			//!! par.evalf(*this);
			bool noexc{true};
			par.advalue(adeqn::eval(*this, par.equation, noexc));
			trc.dprint("constant (",par.name,") = ",par.stringPreValue()," ",par.presUnits());
		}
		par.constant = best_constant;
	}
	
	// remove this name from the stack
	assert(par_stack.top() == par.name);
	par_stack.pop();

	if (par.equation.empty()) {
		trc.dprint("no acceptable equation found");
		return -1;
	}
	// set this parameter's state to Derived
	par.set_derived();

	trc.dprint("returning work ",par.work," for ",par.name,", eqn ",besteqn, ", constant? ",par.constant);
	return par.work;
} // choose_eqn

ostream&
operator<<(ostream& s, const pset& pl) {
	s << "pset " << pl.desc() << "{\n";
	for (const auto& par : pl.pmap())
		s << par.second->longsummary() << endl;
	s << pl.eigv().size() << " eigenvector components\n";
	s << "}\n";
	return s;
}

//------------------------------------------------------------------
// gpset implementation
//------------------------------------------------------------------
std::string
gpset::
persistentName() {
// Returns the name the standard parameter list is to be stored
// and fetched with
	return "globalpset";
}

void
gpset::
store() {
// Store the global gpset in the Flaps database
// with the name given by persistentName()
	Trace trc(1,"gpset::store");

	// this is a static member function so we must use gpset::get()
	gpset& gps = gpset::get();

	trc.dprint(gps.size()," parameters");

	// mark all derived parameters "fresh"
	vector<Par*> derived = gps.get_derived();
	for (auto i : derived)
		i->fresh = true;

	// make a copy of the global pset to store
	pset tostore(gps);

	// delete all eigenvector components: they are stored separate
	// from other parameters, in a vector instead of a map
	tostore.eigv().clear();
		
	// delete all solns
	tostore.clear_solns();

	try {
		tostore.store(persistentName());
	} catch (runtime_error& s) {
		// error("cannot store std parameter list: ", s);
		trc.dprint("caught exception: ",s.what());
		return;
	}
	trc.dprint("stored ",tostore.size()," parameters");
}

gpset::
~gpset() {
	Trace trc(1,"gpset destructor");
	// store for subsequent programs
	this->store();
}

vector<string>
pset::
make_eigvk (int k, double init) {
// create real & imag component parameters for the real and imag
// parts of eigenvector component k, return a vector of the names.
// In general the units of an eigenvector are not known; frequently
// they are dimensionless. The user can add units with the "parameters"
// command.
	vector<string> rval;
	string initstr = vastr(init);
	string namer = evname(k, false);
	string namei = evname(k, true);
	// use brute force instead of findp to avoid recursion
	bool found{false};
	for(size_t i=0; i<eigv_.size(); i++) {
		if (eigv_[i]->name == namer) {
			found = true;
			break;
		}
	}

	if (!found) {
		Eigpar pr(vastr(namer,"(real(y_",k,"))"), initstr);
		pr.set_aux();
		this->add(&pr);
		rval.push_back(namer);
	}
	// imaginary part
	found = false;
	for(size_t i=0; i<eigv_.size(); i++) {
		if (eigv_[i]->name == namei) {
			found = true;
			break;
		}
	}
	if (!found) {
		Eigpar pi(vastr(namei,"(imag(y_",k,"))"), "0.0");
		pi.set_aux();
		this->add(&pi);
		rval.push_back(namei);
	}
	return rval;
}

void
pset::
make_eigv (int nev, double init) {
// create real & imag component parameters for nev
// complex eigenvector components (only the missing components)
	Trace trc(2,"make_eigv");
	for (int i=1; i<=nev; i++)
		vector<string> ei = make_eigvk(i, init);
}

void
gpset::
realloc() {
// a static function in gpset only - must be called
// after Ad::initialize() but before any clones XXX check this?
	Trace trc(1,"gpset::realloc");
	for (auto pp : pmap())
		pp.second->adrealloc();
	for (auto ep : eigv())
		ep->adrealloc();
	set_adpar_derivs();
}


Par*
gpset::
find(const std::string& name) {
// return a pointer to parameter named "name"
// If name refers to an eigenvector component that has not been
// created yet, it and all lesser components will be created
	Par* rval = gpset::get().findp(name);
	return rval;
}

string
tableHelper(const Par* pp, vector<int> colsize, char namefn, char valuefn) {
// helper function for operator<<(gpset)
// write:
//   name desc value+units limits conv eqn
// to the output string with widths for each in colsize
	Trace trc(3,"tableHelper",pp->summary());
	ostringstream os;
	double fac{pp->convFactor()};
	ios::fmtflags left{ios::left};
	os.setf(left);
	os << setw(colsize[0]) << pp->name << namefn << " ";
	// os << setw(colsize[1]) << '(' << pp->desc() << ')';
	os << setw(colsize[1]) << pp->desc() << " ";

	// value with units and multiple-Fixed values
	ostringstream valos;
	if (!pp->altval.empty())
		valos << valuefn << pp->altval[0]*fac;
	else
		valos << valuefn << pp->valuef()*fac;

	// units
	if (!pp->conv.punits.empty())
		valos << " " << pp->conv.punits;
	os << setw(colsize[2]) << valos.str() << " ";
	// limits
	os << setw(colsize[3]) << pp->stringLimits() << " ";
	// conversion factor
	if (!pp->conv.punits.empty())
		os << setw(colsize[4]) << pp->conv << " ";
	else
		os << setw(colsize[4]) << " " << " ";
	// equation
	if (!pp->equation.empty())
		os << setw(colsize[5]) << pp->equation;
	else
		os << setw(colsize[5]) << " ";
	// ignore state

	// multiple-fixed altval's
	if (!pp->altval.empty()) {
		os << "\n  " << pp->altval.size() << " altval:";
		for (auto ai : pp->altval)
			os << "  " << ai*fac;
	}
	return os.str();
}

ostream&
gpset::
table(ostream& s, const string& title) {
// Special op for gpset: print a full table
	Trace trc(3,"gpset::operator<<");
	std::string ti;
	std::map<char, std::string> footnote;
	gpset& pl = gpset::get();
	char valuefn{' '};
	char namefn{' '};
	vector<Par*> indep = pl.get_indep();
	vector<Par*> derived = pl.get_derived();
	vector<Par*> fixed = pl.get_fixed();
	vector<Par*> nostate = pl.get_nostate();
	vector<Par*> aux = pl.get_aux();
	vector<Eigpar*>& ev = pl.eigv();
	std::ostringstream os;
	int pagewidth = page_width();
	size_t n = pl.size();

	std::string minstr, maxstr;

	std::vector<std::string> names(n, "");
	std::vector<std::string> desc(n, "");
	std::vector<std::string> limits(n, "");
	std::vector<std::string> conv(n, "");
	std::vector<std::string> valunits(n, "");
	std::vector<std::string> eqns(n, "");
	std::vector<std::string> states(n, "");

	// grab:
	//   XXX name desc limits conv value+units eqn state
	//   name desc value+units limits conv eqn
	// for each parameter
	size_t i{0};
	for (const auto& par : pl.pmap()) {
		const Par* pp = par.second;

		names[i] = pp->name;
		limits[i] = pp->stringLimits();
		desc[i] = pp->desc();
		
		if (!pp->conv.eunits.empty())
			conv[i] = vastr(pp->conv);

		valunits[i] = pp->stringPreValue();
		// if (pp->hasDefault())
		// 	value[i] += " (default: " + pp->stringPreDefault() + ')';

		if (!pp->conv.punits.empty())
			valunits[i] +=  "  " + pp->conv.punits;

		if (!pp->equation.empty())
			eqns[i] = pp->equation;

		states[i] = pp->state_desc();
		i++;
	}

	// Columns are:
	//   XXX name desc limits conv value+units eqns states
	//   name desc value+units limits conv eqn
	// Figure out how wide to make each column
	std::vector<int> colsize(7, 0);
	for (i=0; i<n; i++) {
		colsize[0] = std::max(colsize[0], (int)names[i].size());
		colsize[1] = std::max(colsize[1], (int)desc[i].size());
		// colsize[4] = std::max(colsize[4], (int)valunits[i].size());
		colsize[2] = std::max(colsize[2], (int)valunits[i].size());
		colsize[3] = std::max(colsize[3], (int)limits[i].size());
		colsize[4] = std::max(colsize[4], (int)conv[i].size());
		colsize[5] = std::max(colsize[5], (int)eqns[i].size());
		colsize[6] = std::max(colsize[6], (int)states[i].size());
	}
	trc.dprintv(colsize, "column sizes");

	// title
	s << title << endl;

	// first print the independents...
	if (!indep.empty()) {
		os.str("");
		os << indep.size() << " Independent Parameters";
		ti =  os.str();
		s << stringCenter(pagewidth, ti) << "\n";
		s << stringCenter(pagewidth, string(ti.size(), '-')) << "\n";
		for (auto pp : indep)
			s << tableHelper(pp, colsize, namefn, valuefn) << endl;
		s << separator() << endl;
	}

	//... then the fixed parameters...
	if (!fixed.empty()) {
		os.str("");
		os << fixed.size() << " Fixed Parameters";
		ti = os.str();
		s << stringCenter(pagewidth, ti) << endl;
		s << stringCenter(pagewidth, string(ti.size(), '-')) << endl;
		for (auto pp : fixed)
			s << tableHelper(pp, colsize, namefn, valuefn) << endl;
		s << separator() << endl;
	}

	// ... then derived & constant parameters...
	// only print values of constant parameters
	if (!derived.empty()) {
		os.str("");
		os << " Derived Parameters";
		ti = os.str();
		s << stringCenter(pagewidth, ti) << endl;
		s << stringCenter(pagewidth, string(ti.size(), '-')) << endl;
		for (auto pp : derived) {
			namefn = ' ';
			valuefn = ' ';
			if (pp->is_constant()) {
				valuefn = '*';
				footnote[valuefn] = "constant";
			}
			s << tableHelper(pp, colsize, namefn, valuefn) << endl;
		}

		if (!footnote.empty()) {
			s << "--------------------\n";
			for (auto& j : footnote)
				s << j.first << " " << j.second << endl;
		} 
		s << separator() << endl;
	}

	// ... then Aux parameters...
	namefn = ' ';
	valuefn = ' ';
	if (!aux.empty()) {
		os.str("");
		os << aux.size() << " Auxiliary Parameters";
		ti = os.str();
		s << stringCenter(pagewidth, ti) << endl;
		s << stringCenter(pagewidth, string(ti.size(), '-')) << endl;
		// int lineargc = ev.size() - nlev.size();
		for (auto pp : aux) {
			if (!pp->is_eigv()) {
				s << tableHelper(pp, colsize, namefn, valuefn) << endl;
			}
		}
		// XXX here print nlev?
		// if (lineargc > 0)
			// s << lineargc << " linear eigenvector components\n";
			s << ev.size() << " eigenvector components\n";
		s << separator() << endl;
	}

	vector<string> const& adparvec = Ad::adnames();
	string adparnames = Ad::toString();
	s << adparvec.size() << " Automatic differentiation parameters: ";
	if (adparvec.size() < 10)
		s << adparnames;
	else
		s << roff(adparnames, 10, 100);
	// string sep("");
	// for (auto& nm : adparnames) {
	// 	s << sep << nm;
	// 	sep = ", ";
	// }
	s << endl;
	s << separator() << endl;

	// Parameters with nostate: unused
	if (!nostate.empty()) {
		ti = "Unused Parameters";
		s << stringCenter(pagewidth, ti) << endl;
		s << stringCenter(pagewidth, string(ti.size(), '-')) << endl;
		for (auto pp : nostate)
			s << tableHelper(pp, colsize, namefn, valuefn) << endl;
		s << separator() << endl;
	}
	return s;
}  // operator<<(gpset)
//------------------------------------------------------------------
// end gpset implementation
//------------------------------------------------------------------


// for finding an element of get_indep(), etc functions above
const Par*
find(const std::vector<const Par*>& list, const std::string& nm) {
	for (size_t i=0; i<list.size(); i++)
		if (list[i]->name == nm)
			return list[i];
	return nullptr;
}

Par*
find(std::vector<Par*>& list, const std::string& nm) {
	for (size_t i=0; i<list.size(); i++)
		if (list[i]->name == nm)
			return list[i];
	return nullptr;
}


#ifdef MAIN
#undef MAIN

using namespace std;

const bool Cmvd = false;
const double Mach = 0.9;  // for plotStdParam()

int
main (int argc, char *argv[]) {
	Ftmpdir ftmp;
	try {
		Par* altp = gpset::get().at("alt");
		altp->set_indep();
		Par* machp = gpset::get().at("mach");
		machp->value(0.8);
		machp->set_fixed();
		Par* freqp = gpset::get().at("freq");
		freqp->set_indep();
		Par* sigmap = gpset::get().at("sigma");
		sigmap->set_indep();
		Par* vsoundp = gpset::get().at("vsound");
		vsoundp->set_derived();
		gpset::get().make_equations();

		double altmin{atmos::altmin()};
		double altmax{atmos::altmax()};
		gpset::get().clear_solns();
		int nstep{101};
		for (int i=0; i<nstep; i++) {
			double alt = altmin + i*(altmax-altmin)/(double)nstep;
			altp->value(alt);
			vsoundp->fresh = false;
			// vsoundp->eval(gpset::get());
			// altp->solns.push_back(alt);
			// vsoundp->solns.push_back(vsoundp->value());
			// evaluate all par...
			gpset::get().eval();
			// and put current values on solns
			gpset::get().back_solns(gpset::get());
		}

		vector<string> toPlot;
		toPlot.push_back ("alt");
		toPlot.push_back ("vsound");
		gpset::get().plot ("pset.apf", "test", "runid", toPlot);
	} catch (runtime_error& s) {
		cerr << "Error: " << s << endl;
		exit(1);
	}
}
#endif // MAIN
