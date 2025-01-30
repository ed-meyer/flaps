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

#include "concurrent.h"
#include "Exception.h"
#include "flutcurve.h"
#include "startpts.h"
#include "main.h"  // setup signal, debugger handling
#include "matrix.h"
#include "specs.h"
#include "process.h"
#include "trace.h"
#include "version.cpp"

using namespace std;

int
main (int argc, char** argv) {
	T_(Trace trc(1,argv[0],":main");)
	string progname{argv[0]};

	try {
		// parse preferences
		parser (progname);

		// print program header
		(void)version(progname);

		// do some set up
		pre_process(progname);

		// get start points
		vector<Flutcurve*> startpoints = startpts();
		if (startpoints.empty()) {
			flaps::warning("no start points found");
			return flaps::nerrors();
		}
		print_start(startpoints);

		// All threads track modes from the stack of start points
		// until there are no more to track
		Fstack fstack(startpoints);
		concurrent(&fstack, trace_curves);

		// the root thread prints and plots results
		post_process(startpoints);

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

