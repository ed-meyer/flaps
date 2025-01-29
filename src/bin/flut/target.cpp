//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Target implementation
// Functions for treating targets in flut
// Public interface:
//		Target(const string& lhs, const string& rhs, bool closest);
//       parses "options", create targets from them

#include <string>
#include <vector>

#include "Par.h"
#include "pset.h"
#include "target.h"
#include "text.h"
#include "trace.h"

using namespace std;

Target::
Target (string const& lhs, string const& rhs, bool cls) {
// create a target by parsing a token lhs=rhs
	Trace trc(1,"Target constructor");

	trc.dprint(lhs," = ",rhs);

	parname = lhs;
	double prevalue{0.0}; // default: 0
	if (!rhs.empty())
		str2double(rhs, prevalue); // presentation value
	// convert to internal units
	Par* pp = gpset::find(lhs);
	if (pp == nullptr) {
		string exc = vastr("target parameter \"",lhs,"\" has not been defined");
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	convfactor = pp->conv.factor;
	value = prevalue/convfactor;
	closest = cls;
}

Target::
Target (string const& par, double val, bool cls, string const& limit) :
		parname(par), value(val), closest(cls) {
	Trace trc(1,"Target double constructor");
	// lower or upper limit?
	if (limit == "lower")
		is_lowerlimit = true;
	else if (limit == "upper")
		is_upperlimit = true;
	// get conversion factor
	Par* pp = gpset::find(par);
	if (pp == nullptr)
		throw runtime_error(vastr("target parameter \"",par,
			"\" has not been defined"));

	convfactor = pp->conv.factor;
}

ostream&
operator<<(ostream& s, const Window& t) {
	s << t.name << "[" << t.min*t.convfactor << ":" << t.max*t.convfactor << "]";
	return s;
}

ostream&
operator<<(ostream& s, const Target& t) {
	if (t.closest) {
		s << t.parname << " closest to " << t.value*t.convfactor;
	} else if (t.is_lowerlimit) {
		s << "lower limit for " << t.parname << " = " << t.value*t.convfactor;
	} else if (t.is_upperlimit) {
		s << "upper limit for " << t.parname << " = " << t.value*t.convfactor;
	} else {
		s << t.parname << " = " << t.value*t.convfactor;
	}
	// any windows?
	if (!t.windows.empty())
		s << " within " << t.windows;
	return s;
}

ostream&
operator<<(ostream& s, const vector<Target>& t) {
	// separate limit-targets from targets
	vector<string> lim;
	for (auto& targ : t) {
		if (targ.is_lowerlimit || targ.is_upperlimit) {
			if (find(lim.begin(),lim.end(),targ.parname) == lim.end())
				lim.push_back(targ.parname);
		} else
			s << targ << endl;
	}
	if (!lim.empty())
		s << "limits for " << lim << endl;
	return s;
}

