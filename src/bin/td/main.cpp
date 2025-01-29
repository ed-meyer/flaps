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

#include "main.h"  // setup signal, debugger handling
#include "hopf.h"
#include "matrix.h"
#include "specs.h"
#include "trace.h"
#include "version.cpp"

using namespace std;

Specs& specs() {
	static Specs thespecs;
	return thespecs;
}

int
main (int argc, char** argv) {
	Trace trc(1,argv[0],":main");
	string progname{argv[0]};
	Specs& sp = specs();


	try {
		// parse specs
		// n: order of the mass, stiffness, etc
		// nbeta: number of RFA betas
		parser (progname, sp);

		// print program header
		(void)version(progname);

		// 2) make equations
		gpset::get().equations();
		// 3) print a table of parameters
		gpset::get().table(cout, "Parameters");

		hopf(sp);

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

