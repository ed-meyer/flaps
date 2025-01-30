//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// catalog file: print an inventory of objects on "file"

#include <regex>
#include <unistd.h> // access

#include "Curve.h"
#include "trace.h"
#include "fio.h"

using namespace std;

int
dofile (string& filename);

void
usage(string const& prog) {
	cerr << "Usage: " << prog << " file1 [file2]...\n";
	cerr << "   file*    Flaps savefiles or plotfiles (.pf) paths\n";
	cerr << "Inventories each file\n";
	exit(1);
}


int
main(int argc, char** argv) {
	Ftmpdir ftmp;

	vector<string> args = cstr2tok(argc, argv);

	if (argc < 2)
		usage(args[0]);

	for (size_t i=1; i<args.size(); i++) {
		try {
			dofile(args[i]);
		} catch (runtime_error& s) {
			cout << s << endl;
		} catch (std::exception& e) {
			cout << e.what() << endl;
		}
	}
}

int
dofile (string& filename) {
// inventory "filename" by restoring & cataloging
	T_(Trace trc(1,"dofile ",filename);)
	int rval{0};

	// check existence of the file, may need to append .pf
	if (access(filename.c_str(), R_OK) == -1) {
		string plotfile = filename + ".pf";
		if (access(plotfile.c_str(), R_OK) == -1) {
			string exc = vastr("file \"",filename,"\" is not available, ",
					"also tried ",plotfile);
			throw runtime_error(exc);
		}
		filename = plotfile;
	}

	// restore the file
	vector<string> cat = fio::restore(filename);
	// cout << cat << endl;

	size_t namewidth{0};
	for (auto name : cat)
		namewidth = std::max(namewidth, name.size());
	namewidth += 3;

	// read each object & print it's name & datatype
	string title{vastr(filename, " has ", cat.size(), " objects:")};
	cout << endl << stringCenter(40, title) << endl;
	cout << stringCenter(40, string(title.size(),'-')) << endl;
	cout << setw(namewidth) << std::left << "   name" << "         data type\n";
	cout << setw(namewidth) << std::left << "   ----" << "         ---------\n";
	for (auto name : cat) {
		Receiver s(name);
		string type;
		s.serialize(type);
		cout << setw(namewidth) << std::left << name << "       " << type << endl;
	}


	// if there are curves print the parameter names
	for (auto name : cat) {
		if (name.substr(0,5) == "curve") {
			Curve* cp = Curve::fetch(name);
			cout << "\n   " << name << " parameters:\n";
			int np{0};
			for (auto pp : cp->params.pmap()) {
				cout << std::setw(14) << pp.first;
				if (++np%6 == 0)
					cout << endl;
			}
			cout << endl;
			break;
		}
	}
	return rval;
}
