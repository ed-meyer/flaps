//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Public interface:
//   string load_df (vector<string>& blocks)		user describing functions

#include <iostream>
#include <regex>
#include <string>
#include <vector>
#include <unistd.h>

#include "fio.h"
#include "load.h"
#include "message.h"
#include "text.h"
#include "trace.h"

using namespace std;

// Functions for loading user-written DF code
// Public interface:
//   bool load_df (vector<string>& blocks);

bool make_userdf(vector<string>& blocks);
static string userdf_string();

string
load_df (vector<string>& blocks) {
// Build a file, userdf.cpp, on FTMP, containing the user's 
// DF functions; compile it and build a dynamic-load library, userdf.so,
// made available to all subsequent processes via LD_PRELOAD
// Returns: the name of the library on FTMP
	Trace trc(1,"load_df");
	string rval;

	if (blocks.empty()) {
		trc.dprint("quick return: no blocks");
		return rval;
	}
	// change to temp directory FTMP using RAII class Chdir
	string ftmp = getftmp();
	Chdir tmp(ftmp);

	trc.dprint(blocks.size()," blocks\n",blocks);
	
	// check read access for each block
	for (auto& bi : blocks) {
		if (access(bi.c_str(), R_OK) == -1)
			throw std::system_error(errno, std::generic_category(),
				vastr("user-function (",bi,") is not available"));
	}

	// create or add to ftmp/userdf.cpp
	make_userdf(blocks);

	// compile the new userdf.cpp into userdf.so
	string froot = getEnv("FROOT");
	string out = vastr(ftmp,"/userdf.out");
	rval = "userdf.so";
	// XXX the commands to compile should come from lib/Makefile
	ostringstream os;
	os << "g++ -std=c++17 -ggdb3 -O0 -Wall -I" << froot
		<< "/include -shared -fPIC -o " << rval
		<< " userdf.cpp >" << out << " 2>&1";
	trc.dprint("running \'", os.str(), "\' to include user code");
	int stat = system(os.str().c_str());

	// if compilation failed throw an exception with it's output
	if (stat != 0) {
		os.str("");
		os << "building userdf.so failed: status " << stat;
		os << " (" << exitstatus(stat) << ")";
		ifstream ifs(out);
		if (ifs) {
			os << "\ncompilation output:\n";
			char c;
			while(ifs.get(c))
				os << c;
		}
		trc.dprint("throwing exception: ",os.str());
		throw runtime_error(os.str());
	}

	flaps::info("added user-written function(s) \"",blocks,
			"\" into subsequent processes");
	return rval;
}

bool
make_userdf(vector<string>& blocks) {
// add user's source code to userdf.cpp by adding a line to it like
// "#include block". userdf_string() will return an existing userdf.cpp
// if there is one on the cwd; otherwise it will return a new one.
	Trace trc(1,"make_userdf");
	bool rval{true};
	ostringstream os;

	// get the prototype userdf.cpp
	string input{userdf_string()};
	trc.dprint("input: {\n",input,"\n}");

	// now write the new userdf.cpp including the source code
	// with #include statements
	ofstream output("userdf.cpp", std::ios::trunc);

	// convert input string to a stringstream so we can use getline;
	// then we do these:
	// 1) look for a line containing the string "include user source"
	//    and put the #include's right after it
	// 3) look for a line containing 'insert Df_fcns here' - put a Df_fcn
	//    for each block after it; for example
	//       Df_fcns("bilinear",bilinear::df_qs,bilinear::df_fq,bilinear::df_interp),
	//       Df_fcns("gap", gap::df_qs, gap::df_fq, gap::df_interp),
	//           ...
	stringstream in(input);
	string line;
	while (!in.eof()) {
		getline(in, line);
		if (line.find("include user source") != string::npos) {
			output << line << endl;
			for (auto& bi : blocks) {
				output << "namespace " << stringBasename(bi) << " {\n";
				output << "#include \"" << bi << "\"\n";
				output << "}\n";
			}
		} else if (line.find("insert Df") != string::npos) {
			for (auto& bi : blocks) {
				string base{stringBasename(bi)};
				output << "Df_fcns(\"" << base << "\", " << base << "::df_qs, "
					<< base << "::df_fq, " << base << "::df_interp)," << endl;
			}
		} else
			output << line << endl;
	}
	return rval;
}

static string
userdf_string() {
// If there is a userdf.cpp in the current directory, read it,
// otherwise return the prototype userdf.cpp
	Trace trc(1,"userdf_string");

	ifstream in("userdf.cpp");
	if (!in) {
		// userdf.cpp has not been created yet - return a raw string
		return R"(
#include <cmath>
#include <string>
#include <vector>
#include "df.h"
#include "trace.h"
#include "interp.h"
using namespace std;
// #include user source here

bool
userdf (const string& dfid, Qfcn& qs_fcn, Ffcn& fq_fcn, Interpfcn& interp_fcn) {
// Initial userdf.cpp - this will get added to in src/bin/flaps:load_df,
// compiled into userdf.so and loaded into all flaps programs via the
// LD_PRELOAD environment variable
	static vector<Df_fcns> entrypts{
	// insert Df_fcns here
	};
	for (auto ei : entrypts) {
		if (ei.dfid == dfid) {
			qs_fcn = ei.qs_fcn;
			fq_fcn = ei.fq_fcn;
			interp_fcn = ei.interp_fcn;
			return true;
		}
	}
	return false;
})";
	}

	// read a userdf.cpp from a previous call to load_df()
	string line;
	ostringstream os;
	while (!in.eof()) {
		getline(in, line);
		os << line << endl;
	}
	trc.dprint("returning existing userdf.cpp");
	return os.str();
}

#ifdef MAIN

int
main(int argc, char** argv) {
	debug(2);
	putEnv("FTMP=/tmp");
	vector<string> entrypts{"junk", "trash"};
	load_fcns(entrypts);
	return 0;
}

#endif // MAIN
