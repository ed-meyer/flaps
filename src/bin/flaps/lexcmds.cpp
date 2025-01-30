//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-


#include <cmath>
#include <iostream>
#include <regex>
#include <string>
#include <utility>
#include <vector>

#include "cmds.h"
#include "text.h"
#include "trace.h"

using namespace std;
using stridx = string::size_type;


vector<pair<string,string> >
lexcmds(const string& in) {
// split a Flaps script "in" into commands and their arguments
	T_(Trace trc(1,"lexcmds");)
	char obrace{'{'};
	// all lower-case char for commands, :0123... for aids
	string lcchar("abcdefghijklmnopqrstuvwxyz:0123456789");
	vector<pair<string,string> > rval;

	stridx pos{0};
	string cmd;
	while (pos != string::npos) {
		// look for the next lower-case char starting at pos
		stridx start = in.find_first_of(lcchar, pos);
		if (start == string::npos || in.substr(start,3) == "end")
			break;
		// the next char (besides whitespace) will be an open brace
		// if this is a command with preferences
		pos = in.find_first_not_of(lcchar, start);
		if (pos == string::npos)
			cmd = in.substr(start);
		else
			cmd = in.substr(start, pos-start);
		if (cmd == "end")
			break;

		// preferences? search for next non-white char, if it is
		// an open brace => start of preferences
		string white{" \t\n"};
		pos = in.find_first_not_of(white, pos);

		// then search for it's closing brace
		string preferences;
		if (in[pos] == obrace) {
			start = pos+1;
			pos = find_delim(in, pos);
			if (pos == string::npos)
				throw runtime_error(vastr("no closing brace for ",cmd));
			preferences = in.substr(start, pos-start);
			pos++;  // move 1 past closing brace, ready for next cmd
			// make preferences one long string: replace cr-lf with commas
			preferences = replace_char(preferences, '\n', ',');
		}
			
		rval.push_back(make_pair(cmd,preferences));
	}
	return rval;
}
