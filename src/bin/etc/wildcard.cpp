//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <iostream>
#include <fnmatch.h>
#include <string>
#include <vector>

#include "trace.h"

using namespace std;

void usage() {
	cerr << "Usage: wildcard pattern string ...\n";
	exit(1);
}

int
main(int argc, char** argv) {
	Trace trc(1, "wildcard");
	vector<string> args = cstr2tok(argc, argv, 1);
	if (args.empty())
		usage();
	string pattern = args[0];
	trc.dprint("checking wildcard for ",args);
	int flags{FNM_EXTMATCH};
	for (size_t i=1; i<args.size(); i++) {
		int res = fnmatch(pattern.c_str(), args[i].c_str(), flags);
		switch (res) {
			case 0:
				cout << args[i] << " matches " << pattern << endl;
				break;
			case FNM_NOMATCH:
				cout << args[i] << " does not match " << pattern << endl;
				break;
			default:
				cout << "\"" << pattern << "\" is invalid\n";
		}
	}
	return 0;
}
