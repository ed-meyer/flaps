//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <sstream>
#include <string>
#include <vector>

#include "config.h"
#include "cmds.h"
#include "fio.h"
#include "lexer.h"
#include "settings.h"
#include "trace.h"

using namespace std;

bool
save (string const& options) {
	T_(Trace trc(1,"save");)
	string outfile{"savefile"};	// default: savefile
	vector<string> tosave{".*"};	// default: save all

	vector<Tok*> unrec = flaps::lexer(options, {
      {"^o(utput)?$", [&](const Tok& p) { outfile = p.srhs; return true; }},
		{".*", [&](const Tok& p) {
			if (p.svec.empty()) {	// no rhs? matrix name to save
				tosave.push_back(p.lhs);
				return true;
			}
			return false;		// unrecognized
		}}
	});

	if (!unrec.empty())
		throw runtime_error(vastr(unrec.size()," unrecognized options: ",unrec));

	// default: save all
	if (tosave.empty())
		tosave.push_back(".*");

	vector<string> saved = fio::save (tosave, outfile);

	if (saved.empty()) {
		return false;
	}
	flaps::info("saved ",saved.size()," matrices:");
	for (size_t i=0; i<saved.size(); i++) {
		cout << setw(20) << saved[i];
		if ((i+1)%4 == 0)
			cout << endl;
	}
	cout << endl;
	return true;
}

bool
restore (string const& options) {
	T_(Trace trc(1,"restore");)
	string infile;

	vector<Tok*> unrec = flaps::lexer(options, {
		{"^i(nput)?", [&](const Tok& p) {infile = p.svec[0]; return true; }},
		{".*", [&](const Tok& p) {
			if (p.svec.empty()) {
				infile = p.lhs;
				return true;
			}
			return false;
		}}
	});

	if (!unrec.empty())
		throw runtime_error(vastr(unrec.size()," unrecognized options: ",unrec));

	vector<string> res = fio::restore (infile);

	if (res.empty()) {
		return false;
	}
	flaps::info("restored ",res.size()," matrices:");
	for (size_t i=0; i<res.size(); i++) {
		cout << setw(20) << res[i];
		if ((i+1)%4 == 0)
			cout << endl;
	}
	return true;
}

#ifdef MAIN

int
main(int argc, char* argv[]) {
	Ftmpdir ftmp;

	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " pf_file_to_restore\n";
		exit(1);
	}

	string file(argv[1]);

	try {
		fio::restore(file);
		string empty;
		vector<string> cat = fio::catalog(empty);
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
}
#endif
