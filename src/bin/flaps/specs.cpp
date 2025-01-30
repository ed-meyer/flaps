//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// get Specs for flaps from settings in the flaps script
#include "config.h"
#include "conv.h"
#include "functions.h"
#include "lexer.h"
#include "settings.h"
#include "specs.h"
#include "trace.h"

using namespace std;

Specs&
flaps_specs() {
	static Specs thespecs;
	return thespecs;
}

bool
parser (string const& settings) {
	T_(Trace trc(1,"parser");)

	Specs& sp = flaps_specs();

	// parse these two options taking the input from cin (empty string)
	vector<Tok*> unrec = flaps::lexer(settings, {
		{"^d(ebug)?", [&](const Tok& p) {
			//!! Settings::defaults.global(settings);
			sp.debug = p.ivec[0];
			flaps::info("setting debug = ",sp.debug," for subsequent processes");
			putEnv(vastr("DEBUG=",sp.debug));		// XXX for now
			putEnv("KEEPFTMP=1");
			Trace::debug(sp.debug);
			return true; }},
		{"^debugger", [&](const Tok& p) {
			sp.debugger = p.srhs;
			return true; }},
		{"^callgrind", [&](const Tok& p) {
			sp.debugger = "valgrind --tool=callgrind";
			return true; }},
		{"^time", [&](const Tok& p) {
			sp.debugger = "valgrind --tool=callgrind";
			flaps::info("timings will be taken for subsequent processes "
				"and stored in a file named callgrind.out.pid");
			return true; }},
		{"^valgrind", [&](const Tok& p) {
			sp.debugger = "valgrind";
			return true; }},
		{"^units", [&](const Tok& p) {
			Settings::defaults.global(settings);	// XXX for now
			return true;
		}},
		{"^wait", [&](const Tok& p) {
			Settings::defaults.global("wait");	// XXX for now
			return true;
		}},
	});
		
	if (!unrec.empty()) {
		string exc{vastr(unrec.size()," options were not recognized: ",unrec)};
		flaps::warning(exc);
	}

	T_(trc.dprint("specs: {\n",sp,"}");)
	return true;
}

ostream&
operator<<(ostream& s, const Specs& t) {
	s << "debug " << t.debug << endl;
	s << "debugger " << t.debugger << endl;
	return s;
}
