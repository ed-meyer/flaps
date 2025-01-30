//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Main program for a standalone amviz

#include <cstdio>
#include <iostream>
#include <string>
#include <unistd.h>	// isatty

#include "amviz.h"
#include "fio.h"
#include "specs.h"
#include "trace.h"
#include "util.h"

// this turns on floating-point exception handling which catches
// one in gtk somewhere
// #include "main.h"

using namespace std;

IMPLEMENT_APP_NO_MAIN(AmvizApp)

bool
AmvizApp::
OnInit() {
// this OnInit will only be executed if AmvizApp is IMPLEMENT_APPed
	T_(Trace trc(1,"AmvizApp::OnInit");)
	// amviz_run creates the AmvizFrame then runs animate in a thread
	amviz_run();
	return true;
}

static void
usage() {
	cout << "Usage: amviz file.uf\n";
	exit(1);
}

int
main(int argc, char** argv) {
// Main function for the standalone program amviz

	// we can't include main.h like with other Flaps programs
	// because wxWidgets and OpenGL get floating-point exceptions
	char* dbg = getenv("DEBUG");
	if (dbg != nullptr) Trace::debug(std::stoi(dbg));
//!! #endif // NEVER // not needed if main.h included
//!! feclearexcept(FE_DIVBYZERO|FE_OVERFLOW);	// turn off fp exceptions
	// Trace::logfile("amviz.err");
	T_(Trace trc(1,"amviz main");)
	// this will not create a new ftmp if FTMP is in the environment and
	// the file exists and is accessible; it should not go inside the try{}
	// space or it will get deleted prematurely
	Ftmpdir ftmp;

	try {
		// parse the specs: some may come on the cmdline - convert them
		// to a string of options
		Specs& sp = specs();
		string options = cmdline(argc, argv);
		// XXX allow specifying the uf file in the File menu?
		if (options.empty() && isatty(fileno(stdin)))
			usage();
		if (!parse_specs(options, sp))
			return 1;
		// This is the standalone program amviz so the command line
		// should have a uf file name: either sp.ufname or sp.plotfiles[0]
		if (sp.ufname.empty()) {
			if (sp.plotfiles.empty() && isatty(fileno(stdin)))
				usage();
			if (!sp.plotfiles.empty())
				sp.ufname = sp.plotfiles[0];
		}
				
		// call Amviz::instance() to create the singleton
		Amviz& amviz = Amviz::instance();
		// make a note that we are running standalone before amviz_kit
		// because it will set the gct to the modes for visualization
		amviz.standalone = true;
		// load the data
		amviz_kit();
		// and create the plot: wxEntry will call OnInit
		int argc{1};
		char* argv[2];
		char* prog = strdup("amviz");
		argv[0] = prog;
		argv[1] = nullptr;
		argc = 1;
		wxEntry(argc, argv);

	} catch (exception& s) {
		string msg = vastr("caught exception: ",s.what());
		flaps::error(msg);
	} catch (...) {
		string msg{vastr("caught unknown exception")};
		flaps::error(msg);
	}
}	// main
