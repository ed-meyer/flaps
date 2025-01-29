//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <sys/types.h>
#include <sys/stat.h>  // stat()

#include "config.h"
#include "lexer.h"
#include "pset.h"
#include "specs.h"
#include "text.h"
#undef Complex

using namespace std;

// return a reference to the Specs: an automatic Singleton
Specs& specs() {
	static Specs thespecs;
	return thespecs;
}

bool
parse_specs (string const& options) {
	Trace trc(1,"parse_specs");
	bool rval = true;
	Specs& sp = specs();
	ostringstream os;
	vector<string> ignore;

	// lambdas for parsing Toks
	vector<Tok*> unrec = flaps::lexer(options, {
		{"i", [&](const Tok& p) {
			sp.input.insert(sp.input.end(), p.svec.begin(), p.svec.end());
			// all rhs options are assumed to be Par defn
			for (auto i : p.roptvec) {
				if (i.empty())
					continue;
				trc.dprint("treating ropt \"",i,"\"");
				vector<Tok*> unrec = flaps::lexer(i, {
					{".*", [&](const Tok& s) {
						sp.params.push_back(new Par(s.lhs, s.srhs));
						return true;
					}}});
			}
			return true;}},
		{"o", [&](const Tok& p) { sp.results = p.svec; return true;}},
		{"octave", [&](const Tok& p) { return sp.use_octave = true; }},
		{"matlab", [&](const Tok& p) { return sp.use_matlab = true; }},
	});

	if (!unrec.empty()) {
		if (unrec.size() == 1)
			flaps::warning("the following option was not recognized: ",*unrec[0]);
		else
			flaps::warning("the following options were not recognized: ",unrec);
	}

	// update gpset with the sp.params
	pset& gp = gpset::get();
	for (auto pi : sp.params) {
		Par* pp = gp.findp(pi->name);
		if (pp == nullptr)
			throw runtime_error(vastr("parameter \"",pi->name,"\" is not defined"));
		pp->value(pi->value());
		trc.dprint("set ",*pp);
	}

/*----------------  end of parse_specs ----------------------------------------*/
	return rval;
}

std::ostream&
operator<<(std::ostream& s, const Specs& t) {
	s << "octave? " << t.use_octave << endl;
	s << "matlab? " << t.use_matlab << endl;
	s << "parameters: " << t.params << endl;
	s << "inputs: " << t.input << endl;
	s << "result: " << t.results << endl;
	return s;
}
