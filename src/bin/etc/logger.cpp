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
#include <fstream>
#include <mutex>
#include <string>

#include "trace.h"

using namespace std;

void
usage(string const& prog) {
	cerr << "Usage: " << prog << " path tmplog\n";
	cerr << "   path    the log file to modify\n";
	cerr << "   tmplog  file with the data to append\n";
	exit(1);
}

int
main(int argc, char** argv) {
	Trace trc(1,"logger");
	string prog(argv[0]);

	if (argc < 3)
		usage(prog);
	// lock this bit of code to prevent another instance of logger 
	// XXX no - must lock the file with fcntl
	// std::mutex mtx;
	// std::lock_guard<std::mutex> guard(mtx);

	// the first arg is the file to append to
	string path(argv[1]);

	// the second arg is the temporary log - one line
	string tmplog(argv[2]);
	ifstream ifs(tmplog);
	string repline;
	if (ifs) {
		getline(ifs, repline);
	} else {
		cerr << "could not open " << tmplog << endl;
		exit(1);
	}

	// open the log file in append mode
	ofstream ofs(path, std::ios::app);
	if (ofs) {
		ofs << repline << endl;
		trc.dprint("appended ",repline);
	} else {
		cerr << "could not open " << path << endl;
		exit(1);
	}
}
