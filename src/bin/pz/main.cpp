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
#include "main.h"
#include "settings.h"
#include "specs.h"
#include "version.cpp"

using namespace std;

int
main (int argc, char** argv) {
	Trace trc(1,argv[0]);
	string progname{argv[0]};

	// if there is an argument it is the output matrix name
	if (argc > 1)
		specs().outputname = string(argv[1]);

	version(progname);
	try {
		// most of the work is done in parser
		Matrix* output = parser();

		output->store();

		flaps::info("Output:\n", *output);
	} catch (runtime_error& s) {
		flaps::error(s.what());
	} catch (exception& t) {
		flaps::error("caught exception: ",t.what());
	}

	return flaps::nerrors();
}

