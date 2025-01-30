//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// parameters: allow the user to define new parameters or modify
// existing parameters

#include "config.h"

#include <iostream>
#include <stdexcept>
#include <string>

#include "conv.h"
#include "lexer.h"
#include "main.h"  // setup signal, debugger handling
#include "matrix.h"
#include "trace.h"
#include "version.cpp"

using namespace std;

bool parser ();

int
main (int argc, char** argv) {
	string progname{argv[0]};
	T_(Trace trc(1,progname,":main");)

	// print program header
	(void)version(progname);

	try {
		// parse options & create parameters
		parser ();
	} catch (runtime_error& s) {
		flaps::error(s.what());
		return 1;
	} catch (const std::exception& t) {
		flaps::error(t.what());
		return 1;
	} catch (...) {
		flaps::error("caught unknown exception");
		return 1;
	}

	return (flaps::nerrors());
}

//--------------------------------------------------------------------
// parameters { "lnsbf(Left nacelle side-bending freq)[0:10]<radps2Hz>=1",...}
//-------------------------------------------------------------------
bool
parser () {
	T_(Trace trc(2,"parser");)

	// no Parsers - just get all Tok*
	vector<Tok*> toks = flaps::lexer("");
	for (auto pi : toks) {
		// all Toks are assumed to be Par definitions
		Par np(pi->lhs, pi->svec);
		np.pref = 1;  // defined as a user preference
		if (np.candidates.empty())	// no equations? must be auxilliary
			np.set_aux();
		// check to see if it exists already
		Par* existing = gpset::find(np.name);
		if (existing != nullptr) {
			string oldpar{vastr(*existing)};
			// upgrade the existing parameter with stuff given here
			existing->upgrade(&np);
			flaps::info("changed \"",oldpar,"\" to \"",*existing,"\"");
		} else {
			gpset::get().add(&np);
			flaps::info("created new parameter: ",np);
		}
	}
	return true;
}
