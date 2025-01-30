//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// functions for printing warning and error messages

#include <iostream>
#include <cerrno>
#include <cstring>
#include <unistd.h>  // for getpid()

#include "fio.h"
#include "message.h"
#include "text.h"

using namespace std;

namespace flaps {

static int ErrorCount = 0;
static int WarningCount = 0;

void
incr_nwarnings() {
	WarningCount++;
}
void
incr_nerrors() { ErrorCount++; }

// This function is here to allow setting a break point in ddd
// it would not be necessary (?) if ddd had _cxa_throw()
void
throwOutOfRange(std::string const& s) {
	throw std::out_of_range(s);
}

/*--------------------------------------------------------------------
 * error/warningCount:      returns the number of errors/warnings
 *------------------------------------------------------------------*/
int
nerrors (void) { return ErrorCount; }
int
nwarnings (void) { return WarningCount; }

static
string
helper(const string& msg, int count) {
	ostringstream os;
	os << "+++++++++++++++  " << count << " " << msg << " message";
	if (count > 1)
		os << 's';
	os << " printed  ++++++++++++++++++\n";
	return os.str();
}

void
error_summary() {
// write to stdout a summary of the number of warning/error messages
// that have been printed
	if (ErrorCount > 0 || WarningCount > 0) {
		cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		if (ErrorCount > 0)
			cout << helper("error", ErrorCount);
		if (WarningCount > 0)
			cout << helper("warning", WarningCount);
	}
}

} // namespace flaps

vector<string>
backtrace() {
	vector<string> rval;

	// use RAII Fpipe
	Fpipe p(vastr("gdb -batch -ex bt -p ",getpid()));
	size_t bufsize{1024};
	char* buf = new char[bufsize];
	while (getline(&buf, &bufsize, p()) != -1) {
		string str(buf);
		if (str.find("LWP") == string::npos)
			rval.push_back(str);
	}
	return rval;
}

#ifdef MAIN

void
testbacktrace() {
	cout << "Test of backtrace:\n";
	cout << backtrace() << endl;
}

void
throw_excep() {
	// throw excep("test of excep");
	string exc("test of excep");
	throw runtime_error(exc);
}

int
main(int argc, char** argv) {

	testbacktrace();

	flaps::error ("Now is the time for all good"
				" men to come to the aid of their country; "
				"Four-score and thirty years ago our fathers did something"
				" which I cannot remember");
	// errorSummary();

	try {
		throw_excep();
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
   exit(0);
}

#endif
