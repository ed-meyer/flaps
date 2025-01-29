//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "trace.h"
#include "Regex.h"

using namespace std;

regex
make_regex(const string& pat) {
// convenience function for std::regex, catches exceptions
// from bad patterns and throws an exception with a better
// error message
	Trace trc(1,"make_regex ",pat);
	regex rval;
	try {
		rval = regex(pat);
	} catch (const regex_error& e) {
		string exc{vastr("bad regular expression (",pat,"): ",regex_msg(e.code()))};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	return rval;
}

regex
make_regex_ic(const string& pat) {
// create a regex ignoring case
// convenience function for std::regex, catches exceptions
// from bad patterns and throws an exception with a better
// error message
	Trace trc(1,"make_regex_ic ",pat);
	regex rval;
	try {
		rval = regex(pat, regex_constants::icase);
	} catch (const regex_error& e) {
		string exc{vastr("bad regular expression (",pat,"): ",regex_msg(e.code()))};
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}
	return rval;
}


#ifdef MAIN


int
main (int argc, char **argv) {
	Trace trc(1,argv[0]);
	string rx;
	string test;

	// If no arguments given just test a compilcated re...
	if (argc < 2) {
		flaps::info ("Usage: ",argv[0]," regular_expression item_to_test "
				"[match_number] [match_number]...  (matches ignoring case)");
//		exit(1);
      rx = "^[ \t]*(if[ \t]*\\([ \t]*!?)?(.*:)?(atlas|elfini)[ \t]*\\{.*$";
		test = "atlas { junk\nmorejunk";
		//test = "atlas { junkmorejunk";
	} else {
		rx = argv[1];
		if (argc > 2) {
			test = argv[2];
		} else {
			ostringstream os;
			while(std::cin >> test)
				cout << test;
			trc.dprint("got <",test,"> from cin");
		}
	}

	// try matching using POSIX basic syntax, case sensitive,
	// no matching of newlines...
	smatch mch;
	try {
		if (regex_match(test, mch, regexMake(rx))) {
			cout << test << " matches " << rx << endl;
			for (auto& m : mch) {
				cout << "sub-match: " << m << endl;
			}
		} else {
			cout << test << " does not match " << rx << endl;
		}
	} catch (runtime_error& s) {
		std::cerr << s << std::endl;
	}
}
#endif
