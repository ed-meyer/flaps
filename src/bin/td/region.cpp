//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "region.h"

using namespace std;

// classes Bndry and Region implementation


Bndry::
Bndry (Tok const& opt) {
// Parse a start "Bndry", e.g. "freq[1:10]". A Bndry consists of a
// Par (name) or a Curve id.
// Method: create a new parameter by parsing the description, then
// call the Bndry(Par*) constructor which will add the parameter to
// the gpset if it is new, and convert min and max to internal if
// a conversion factor is not part of lhs
	Trace trc(2,"Bndry(opt) constructor");

	trc.dprint("input opt: ",opt);

	closest = false;
	min = 0.0;
	max = 0.0;

	// special case: lhs is a curve id, e.g. mode_1.a and no rhs
	if (opt.lhs.size() > 4 && opt.lhs.substr(0,4) == "mode"
			&& opt.svec.empty()) {
		curveid = opt.lhs;
	} else if (compare(opt.lhs,"ord") && !opt.ivec.empty()) {
		ordinal = opt.ivec;
	} else if (compare(opt.lhs,"mode") && !opt.ivec.empty()) {
		modes = opt.ivec;
	} else if (opt.srhs == "closest") {
		closest = true;
		if (opt.ropt.empty()) {
			string exc("incomplete \"closest\" option: "
					"should have the form \"sigma=closest{0}\"");
			throw runtime_error(exc);
		}
		str2double(opt.ropt, min);
		max = min;
		parname = opt.lhs;
	} else {
		try {
			// Par parser constructor: only pass srhs, not svec,
			// let upstream deal with multiple rhs
			Par ap(opt.lhs, opt.srhs);
			trc.dprint("creating a Bndry with: ",ap.desc());

			// Grab the limits and make a note as to whether they are
			// in internal or external units so we can convert them to
			// internal later if necessary
			if (ap.has_min()) {
				min = ap.min();
			} else {
				min = ap.value();
			}
			if (ap.has_max()) {
				max = ap.max();
			} else {
				max = ap.value();
			}
			bool internal = false;
			if (ap.has_conv())
				internal = true;

			// ok if it is a new parameter - add it to the gpset
			// The Bndry limits are the parameter limits in this case
			Par* existing = gpset::find(ap.name);
			Par* par;
			if (existing != nullptr) {
				par = existing;
			} else {
				par = gpset::get().add(&ap);
			}
			// save the name only:
			parname = par->name;

			// Convert the limits to internal if necessary
			double fac = par->convFactor();
			if (!internal) {
				min /= fac;
				max /= fac;
				trc.dprint("converted min, max to internal: ",min,", ",max);
			}
		} catch (runtime_error& s) {
			throw runtime_error(vastr("bad start-region definition (",opt,"): ", s.what()));
		}
	}
	trc.dprint("constructed: ",*this);
}

string
Bndry::
summary(pset& plt) const {
	ostringstream os;
	double fac{1.0};

	// just a curveid...
	if (!curveid.empty()) {
		os << "curveid=" << curveid;
		return os.str();
	} else if (!modes.empty()) {
		return vastr("mode(s): ",modes);
	} else if (!ordinal.empty()) {
		return vastr("ordinal(s): ", ordinal);
	} else if (closest) {
		return vastr(parname, " closest to ", min);
	}

	const Par* par = plt.findp(parname);

	if (par)
	 	fac = par->convFactor();

	os << parname;

	if (min == max) {
		os << " = " << min*fac;
	} else {
		os << " [" << min*fac
		<< ':' << max*fac << ']';
	}
	return os.str();
}

Region::
Region (const string& options) {
// parse a set of options and create a list of Bndrys
// This constructor allows for multi-values options, e.g.
//   start{vtas=(100,200)}
	Trace trc(1,"Region(string) constructor");

	vector<Tok*> toks = flaps::lexer(options);
	for (auto tok : toks) {
		if (tok->svec.size() > 1) {
			for (size_t i=0; i<tok->svec.size(); i++) {
				Tok p(*tok);
				p.srhs = tok->svec[i];
				p.svec.clear();
				p.svec[0] = p.srhs;
				push_back(Bndry(p));
			}
		} else {
			push_back(Bndry(*tok));
		}
	}
}

string
Region::
summary(pset& plt) const {
	ostringstream os;

	if (this->empty())
		os << "empty";

	for (size_t j=0; j<size(); j++) {
		if (j > 0)
			os << ", ";
		os <<  (*this)[j].summary(plt);
	}
	return os.str();
}

bool
Bndry::
isPoint() const {
	if (!parname.empty())
		return min == max;
	return false;
}

Bndry*
Region::
pointBndry() {
// if this Region has a point Bndry (one where min==max)
// return a pointer to that Bndry, otherwise return nullptr
	Trace trc(1,"Region::pointBndry");
	for (auto& reg : *this) {
		if (reg.isPoint()) {
			trc.dprint("returning ",reg);
			return &reg;
		}
	}
	return nullptr;
}

vector<string>
Region::
curveidBndry() {
// if this Region has a Bndry which is simply a curveid
// return that curveid
	Trace trc(1,"Region::curveidBndry");
	vector<string> rval;
	for (auto& reg : *this) {
		if (!reg.curveid.empty()) {
			trc.dprint("returning ",reg.curveid);
			rval.push_back(reg.curveid);
		}
	}
	return rval;
}

int
getmodeno(const string& cid) {
// given a curveid (cid), what is the mode number; e.g.
//   curve,mode_11.c.3 -> 11
	string::size_type idx = cid.find("mode_");
	string::size_type start = idx + 5;
	string::size_type end = start;
	while (end != cid.size() && !(cid[end] == '.' || cid[end] == ','))
		end++;
	string mstr = cid.substr(start, end-start);
	return stoi(mstr);
}

bool
Region::
curveidOk(string& cid) {
// if this Region has curveid's does cid match any of them?
	Trace trc(1,"curveidOk");

	// watch out for uncertainty curves: cid like mode_1_hi or mode_1_lo
	string tail = rsubstr(cid, 2);
	if (tail == "hi" || tail == "lo")
		return false;

	int ntest{0};
	smatch mch;
	// check the curveid in each Bndry in this
	for (auto& bndry : *this) {
		if (!bndry.curveid.empty()) {
			ntest++;
			if (regex_match(cid, mch, regex(bndry.curveid))) {
				trc.dprint(cid," matches ",bndry.curveid);
				return true;
			}
		}
		if (!bndry.modes.empty()) {
			ntest++;
			int modeno = getmodeno(cid);
			if (this->contains_mode(modeno)) return true;
		}
	}
	// nothing matched: return false
	if (ntest > 0)
		return false;
	// no curvid or modes in this Region: return true
	return true;
}

bool
Region::
hasCurveid(string& id) {
// does this Region contain one or more curveid Bndrys,
// and if so do any match "id"?
	Trace trc(1,"Region::hasCurveid");
	int nid{0};
	for (auto& bndry : *this) {
		if (!bndry.curveid.empty()) {
			nid++;
			if (bndry.curveid == id) {
				trc.dprint("returning ",bndry.curveid);
				return true;
			}
		}
	}
	if (nid > 0)
		return false;
	else
		return true;
}

vector<int>
Region::
ordinal_numbers() const {
	Trace trc(1,"Region::ordinal_numbers");
	vector<int> rval;
	for (auto& bndry : *this) {
		if (!bndry.ordinal.empty()) {
			for (auto i : bndry.ordinal) {
				rval.push_back(i);
			}
		}
	}
	trc.dprint("returning ",rval.size()," ordinals: ",rval);
	return rval;
}

vector<int>
Region::
mode_numbers() const {
	Trace trc(1,"Region::mode_numbers");
	vector<int> rval;
	for (auto& bndry : *this) {
		if (!bndry.modes.empty()) {
			for (auto i : bndry.modes) {
				rval.push_back(i);
			}
		}
	}
	trc.dprint("returning ",rval.size()," modes: ",rval);
	return rval;
}

bool
Region::
contains_ordinal(int ordinal) const {
// does this Region contain "ordinal", i.e. did the
// option have something like "ordinal=1"
	Trace trc(1,"Region::contains_ordinal ",ordinal);
	vector<int> ordinals = this->ordinal_numbers();
	trc.dprint(ordinals.size()," ordinals");
	// if no ordinals were specified all are ok
	if (ordinals.empty()) {
		trc.dprint("quick (true) return: no ordinals specified");
		return true;
	}

	for (auto i : ordinals) {
		if (ordinal == i) {
			trc.dprint("returing true");
			return true;
		}
	}
	trc.dprint("returning false");
	return false;
}

bool
Region::
contains_mode(int mode) const {
// does this Region contain "mode", i.e. did the
// option have something like "mode=1"
	Trace trc(1,"Region::contains_mode ",mode);
	vector<int> modes = this->mode_numbers();
	trc.dprint(modes.size()," modes");
	// if no ordinals were specified all are ok
	if (modes.empty()) {
		trc.dprint("quick (true) return: no modes specified");
		return true;
	}

	for (auto i : modes) {
		if (mode == i) {
			trc.dprint("returing true");
			return true;
		}
	}
	trc.dprint("returning false");
	return false;
}


ostream&
operator<<(ostream& s, const Bndry& t) {
	s << t.summary(gpset::get());
	return s;
}

ostream&
operator<<(ostream& s, const Region& t) {
	s << t.summary(gpset::get());
	return s;
}

const Par*
Region::
find (pset& plt, string const& name) const {
// Look for a Bndry in "this" with "parname" == "name",
// return a pointer to that parameter in the global pl
	for (size_t i=0; i<size(); i++) {
		if ((*this)[i].parname == name) {
			const Par* existing = plt.findp(name);
			return existing;
		}
	}
	return 0;
}

Bndry*
Region::
get (string const& name) {
// return a pointer to the Bndry in this list with
// parameter name "name"
	Trace trc(1,"Region::get ",name);
	for (auto& bndry : *this) {
		if (bndry.parname == name) {
			trc.dprint("returning ",bndry);
			return &bndry;
		}
	}
	return nullptr;
}

bool
Region::
contains (string const& name, double val) const {
// Does this Region contain a Bndry with parameter
// "name" and if so does its limits contain "val"?
// Returns: false if this Region contains "name" and "val"
//          does not lie in the interval, true otherwise
	for (auto& bndry : *this) {
		if (bndry.parname == name) {
			if (bndry.min <= val && val <= bndry.max)
				return true;
			else
				return false;
		}
	}
	// Bndry not found - just return false
	return false;
}

const Par*
Region::
outOfRange (pset& pl, string& oormsg) const {
// Find a parameter in pl which is out of these Bndrys
	Trace<> trc("outOfRange");
	for (const auto& bndry : *this) {
		if (!bndry.isPoint()) {
			const Par* pp = pl.findp(bndry.parname);
			if (pp) {
				trc.dprint("check ",pp->summary()," against ",bndry);
				double val = pp->value();
				double fac = pp->convFactor();
				if (bndry.min > val || val >= bndry.max) {
					ostringstream os;
					os << bndry.parname << " is out of the range ["
						<< bndry.min*fac << ':' << bndry.max*fac << ']';
					oormsg = os.str();
					trc.dprint("out of range: returning ",pp->summary());
					return pp;
				}
			}
		}
	}
	return nullptr;
}
