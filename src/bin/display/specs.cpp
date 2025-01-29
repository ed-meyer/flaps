//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <string>
#include <tuple>
#include <vector>

#include "config.h"
#include "Ad.h"
#include "exim.h"
#include "fio.h"
#include "lexer.h"
#include "matrix.h"
#include "specs.h"
#include "pset.h"

using namespace std;

Specs&
specs() {
	static Specs thespecs;
	return thespecs;
}

ostream&
operator<<(ostream& s, const Specs& t) {
#ifdef NEVER // new syntax
	if (t.matrices.empty()) {
		s << "no matrices specified\n";
		return s;
	}
#endif // NEVER // new syntax
	s << "filter " << t.filter << endl;
	s << "format " << t.format << endl;
	s << "output file " << t.outfile << endl;
	if (!t.params.empty()) {
		s << "parameters:\n";
		for (auto pi : t.params)
			s << "  " << get<0>(pi) << " = " << get<1>(pi) << endl;
	}
	if (!t.derivs.empty()) {
		s << "derivatives:\n";
		for (auto di : t.derivs)
			s << "  " << di << endl;
	}
#ifdef NEVER // new syntax
	s << "matrices:\n";
	for (auto di : t.matrices)
		s << "  " << di << endl;
#endif // NEVER // new syntax
	return s;
}

vector<pair<string,Specs>>
parser () {
/*------------------------------------------------------------------
 * parse the options, return a vector of (mid,specs) pair
 * Syntax:
 *    mid{specs}, ...
 * where specs are
 *   filter
 *   matview
 *   o(ut)?
 *   summary
 *   deriv(s)?=(0,par_name,...)
 *------------------------------------------------------------------*/
	Trace trc(2,"parser");
	vector<pair<string,Specs>> rval;

	// each Tok is a matrix name with optional specs
	vector<Tok*> unrec = flaps::lexer("", {
		{".*", [&](const Tok& p) {	// matrix name?
			if (!p.svec.empty())
				return false;
			Specs sp;
			if (p.lopt.empty()) {
				rval.push_back({p.lhs,sp});
				return true;
			}
			// treat the matrix options
			vector<Tok*> unrec = flaps::lexer(p.lopt, {
				{"filter", [&](const Tok& p) { sp.filter = p.rvec[0]; return true; }},
				{"matview", [&](const Tok& p) { return sp.matview = true; }},
				{"o(ut)?", [&](const Tok& p) { sp.outfile = p.srhs; return true; }},
				{"format", [&](const Tok& p) { sp.format = p.srhs; return true; }},
				{"summary", [&](const Tok& p) { sp.format = "summary"; return true; }},
				{"deriv(s)?", [&](const Tok& p) {
					sp.derivs.insert(sp.derivs.end(),p.svec.begin(),p.svec.end());
					return true;
				}},
				{".*", [&](const Tok& p) {	// parameter value?
					if (p.svec.empty()) return false;	// matrix name?
					sp.params.push_back({p.lhs,p.rvec[0]});
					return true;
				}},
			});
			if (!unrec.empty())
				flaps::warning(unrec.size()," options were not recognized: ",unrec);
			rval.push_back({p.lhs,sp});
			trc.dprint("specs for ",p.lhs,":\n",sp);
			return true;
		}},
	});
	if (!unrec.empty())
		flaps::warning(unrec.size()," options were not recognized: ",unrec);

	return rval;
}
