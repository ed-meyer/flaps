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

#include <iostream>
#include <sys/types.h>
#include <unistd.h>

#include "custom.h"
#include "lexer.h"
#include "main.h"  // setup signal, debugger handling
#include "trace.h"
#include "version.cpp"

using namespace std;

int
main (int argc, char** argv) {
	T_(Trace trc(1,argv[0],":main");)
	string progname{argv[0]};

	try {
		string entrypt;	// the command
		string options;		// options for the command
		vector<Tok*> unrec = flaps::lexer ("", {
			{".*", [&](const Tok& p) {
				entrypt = p.lhs;
				options=p.lopt;
				return true;
			}}
		});
		if (!unrec.empty())
			throw runtime_error(vastr("unrecognized option(s): ",unrec));

		// get a pointer to the function assuming it was added to custom.so
		// in flaps
		int (*fcn)(const string& options);
		if (!Custom::is_defined(entrypt)) {
			string exc{vastr("\"call\" function ",entrypt," has not been defined")};
			string fcns;
			string sep;
			for (auto fi : Custom::map()) {
				fcns += sep + fi.first;
				sep = ", ";
			}
			exc += vastr(", defined functions: ",fcns);
			throw runtime_error(exc);
		}
		Custom::get(entrypt, fcn);
		fcn(options);

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

	// print a summary of errors/warnings...
	flaps::error_summary();

	return (flaps::nerrors());
}

