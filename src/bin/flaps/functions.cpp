//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <map>
#include <string>

#include "config.h"
#include "functions.h"
#include "pset.h"
#include "settings.h"
#include "specs.h"
#include "trace.h"

using namespace std;

using mp = std::pair<std::string,Fcnptr>;

// all the functions which may be called in a flaps script file
static map<string,Fcnptr> Functions {
	mp("catalog", catalog),
	mp("export", exporter),
	mp("gyro", gyro),
	mp("import", importer),
	mp("lti", lti),
	mp("settings", parser),	// XXX change parser to settings?
	mp("restore", restore),
	mp("save", save),
	mp("symmetrize", symmetrize),
	mp("unix_cmd", unix_cmd)
};

vector<string>&
fcn_names() {
// returns a vector of all defined function names
	static vector<string> rval;
	if (rval.empty()) {
		for (auto& fp : Functions) {
			rval.push_back(fp.first);
		}
	}
	return rval;
}

Fcnptr
fcnptr(const string& name) {
// Returns a pointer (Fcnptr) to function "name"; throws an
// "out_of_range" exception if unrecognized
	// at() will throw an exception if "name" is not in fcns
	for(auto& fun : Functions)
		if (regex_match(name, regex(fun.first)))
			return fun.second;
	throw runtime_error(vastr(name, " is not a valid function name"));
	return nullptr;
}

bool
unix_cmd (const string& cmd) {
	Trace trc(1,"unix_cmd");
	bool rval = true;
	ostringstream os;

	// Expand environment variables in the command line
	string cmdline(cmd);
	expenv(cmdline);

	// unix_preamble (cmdline);
	int status = system(cmdline.c_str());

	if (status != 0) {
		string exc = vastr(cmd, " returned status ",status);
		trc.dprint("unix cmd error: : ",exc);
		flaps::error(exc);
		return false;
	}

	// colophon(cmdline, errorCount(), true, 0);

	return rval;
}
