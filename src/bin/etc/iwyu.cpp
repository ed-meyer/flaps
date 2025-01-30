//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// iwyu: includes 

#include <iostream>
#include <fstream>
#include <regex>

#include "trace.h"

using namespace std;

int
dofile (const string& filename, bool sysh);

void
usage(string const& prog) {
	cerr << "Usage: " << prog << " file1 [file2]...\n";
	cerr << "   file*    C/C++ file path\n";
	cerr << "Checks if #include files are needed by brute-force:\n";
	cerr << "comment each #include in turn and try compiling; if successful\n";
	cerr << "keep it commented\n";
	exit(1);
}


int
main(int argc, char** argv) {
	string prog(argv[0]);

	if (argc < 2)
		usage(prog);

	for (int i=1; i<argc; i++) {
		string file(argv[i]);
		try {
			// flaps includes
			dofile(file, false);
			// system includes
			dofile(file, true);
		} catch (runtime_error& s) {
			cout << s.what() << endl;
		} catch (std::exception& e) {
			cout << e.what() << endl;
		}
	}
}

int
dofile (const string& filename, bool sysh) {
	T_(Trace trc(1,"dofile ",filename);)
	int rval{0};
	vector<string> lines;
	string line;
	ostringstream os;
	ifstream ifs(filename);
	if (!ifs) {
		os << filename << " is not available";
		T_(trc.dprint("throwing exception: ",os.str());)
		throw runtime_error(os.str());
	}
	// read the file, put each line into lines
	while (!ifs.eof()) {
		istream& strm = getline(ifs, line);
		if (strm)
			lines.push_back(line);
	}
	ifs.close();
	T_(trc.dprint("read ",lines.size(),"-line file");)
	T_(trc.dprint("last line: \"",lines[lines.size()-1],"\"");)

	// the make command
	os.str("");
	os << "make " << filename;
	string make(os.str());
	make[make.size()-1] = 'o';
	T_(trc.dprint("compiling with \"",make,"\"");)
	// rm .o command
	os.str("");
	os << "rm -f " << filename;
	string rmo(os.str());
	rmo[rmo.size()-1] = 'o';
	T_(trc.dprint("rm .o file with \"",rmo,"\"");)

	// regex re("#include [<\"]\([^\">]*\)[\">]");
	//	re = regex(R"(#include[ \t]*[<\"]([^\">]*)[\">])");
	regex re;
	if (sysh) {
		re = regex(R"(#include[ \t]*<([^>]*)>)");
	} else {
		re = regex(R"(#include[ \t]*\"([^\"]*)\")");
	}

	size_t i;
	size_t k{0};
	while (true) {
		ofstream ofs(filename, std::ios::out|std::ios::trunc);
		bool goth = false;
		size_t lineno{0};
		string hfile;
		for (i=0; i<lines.size(); i++) {
			smatch mch;
			if (hfile.empty() && i >= k && regex_match(lines[i], mch, re)) {
				hfile = mch[1];
				T_(trc.dprint("got include file ",hfile);)
				k = i+1;
				lineno = i;
				goth = true;
			} else {
				ofs << lines[i] << endl;
			}
		}
		if (!goth)
			break;
		system(rmo.c_str());
		int stat = system (make.c_str());
		T_(trc.dprint("make returned ",stat);)
		if (stat == 0) {
			os.str("");
			// compilation successful: comment out this include
			lines[lineno] = os.str();
			rval++;
			T_(trc.dprint(hfile," unused: replace line: ",lines[lineno]);)
		}
	}
	T_(trc.dprint(rval," unused include files");)
	return rval;
}

