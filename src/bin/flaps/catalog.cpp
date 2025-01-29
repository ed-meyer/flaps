//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "config.h"
#include "fio.h"
#include "functions.h"
#include "lexer.h"
#include "settings.h"

using namespace std;

static void doCat (string const& pat, bool sum);

//--------------------------------------------------------------------
// list catalog entries on cout
//--------------------------------------------------------------------
bool
catalog (string const& options) {
	Trace trc(2,"catalog");
	bool summary{false};	// give a summary, not just names
	string pat{".*"};	// default: all


	// lambdas for handling options
	// if options is empty do not call lexer: it will try
	// to read options from cin
	vector<Tok*> unrec;
	if (!options.empty()) {
		unrec = flaps::lexer (options, {
			{"all", [&](const Tok& p) { pat = ".*"; return true; } },
			{"sum(mary)?", [&](const Tok& p) { return summary = true; } },
			{".*", [&](const Tok& p) {
					pat = p.lhs;
					if (!p.svec.empty())
						pat += string("=") + p.svec[0];
					return true;
				}
			}
		});
	}

	if (!unrec.empty())
		throw runtime_error(vastr("unrecognized specs: ",unrec));
	
	doCat(pat, summary);
	return true;
}

static void
doCat (string const& pat, bool sum) {
	Trace trc(2,"doCat");
	vector<string> cat;

	// get the catalog from the data manager...
	try {
		cat = fio::catalog(pat, sum);
	} catch (runtime_error& s) {
		cerr << "catalog: " << s << endl;
	}

	// ... and print it
	if (cat.empty()) {
		if (pat.empty())
			cout << "   there are no matrices in the Flaps database\n";
		else
			cout << pat << " is not in the Flaps database\n";
	} else {
		if (!pat.empty() && pat != ".*") {
			cout << "     Flaps Data Inventory: " << pat << endl;
			cout << "     ---------------------\n";
		} else {
			cout << "     Flaps Data Inventory\n";
			cout << "     --------------------\n";
		}
		for (auto& s : cat)
			cout << stringIndent(3,s) << endl;
	}
}

