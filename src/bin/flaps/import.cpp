//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <string>
#include <utility>

#include "config.h"
#include "exim.h"
#include "fma.h"
#include "functions.h"
#include "lexer.h"
#include "matrix.h"
#include "settings.h"
#include "trace.h"

using namespace std;

static void rename(vector<Matrix*> ml, vector<pair<string,string>> newnames);
static string get_format(string const& path);

static bool
input_f (const Tok& opt) {
// handle a file name: import the file according to its extension
// Optional output file syntax:
//   name = path
	string output;	// name of the matrix output
	string path;
	if (opt.srhs.empty()) {
		path = opt.lhs;
	} else {
		output = opt.lhs;
		path = opt.srhs;
	}
	// check for path{oldname=newname, id=...}
	vector<std::pair<string,string>> newnames;
	string id;
	string vzid;
	vector<int> dof;
	if (!opt.lopt.empty()) {
		vector<Tok*> unrec = flaps::lexer(opt.lopt, {
			{ "id", [&](const Tok& p) { id = p.srhs; return true; }},
			{ "vzid", [&](const Tok& p) { vzid = p.srhs; return true; }},
			{ "dof", [&](const Tok& p) { dof = p.ivec; return true; }},
			{ "^.*", [&](const Tok& p) {
				newnames.push_back(make_pair(p.lhs,p.srhs)); return true;}}
		});
	}
	string fileformat = get_format(path);
	if (fileformat.empty()) {
		flaps::error("could not determine file format for ",path);
		return flaps::nerrors();
	}

	flaps::info ("Import of ",fileformat," file ",path);

	vector<Matrix*>  ml;
	if (fileformat == "output4") {
		ml = Output4::importer(path, output);
	} else if (fileformat == "matlab") {
		ml = Matlab::importer (path, output);
	} else if (fileformat == "matrixmarket") {
		Matrix* mp = MM::importer(path, output);
		ml.push_back(mp);
	} else if (fileformat == "universal") {
		// the Fma constructor will import a UF and store the Fma
		// If the user did not give a vzid use the file name without ext
		if (vzid.empty()) {
			vzid = stringBasename(path);
			flaps::info("defaulting the visualization id to \"",
				vzid,"\"");
		}
		Fma fma(path, vzid);
	} else {
		flaps::error("unrecognized file format: ", fileformat);
		return true;
	}

	// rename the matrices if requested
	rename(ml, newnames);
	// tag the matrix names with "id"
	if (!id.empty())
		for (auto m : ml)
			m->mid(m->mid() + '.' + id);

	// Store the matrices
	for (auto mp : ml) {
		// extract dof XXX ???
		if (!dof.empty() && mp->rsize() == mp->csize()) {
			size_t n = dof.size();
			Matrix* newmp = new Matrix(mp->mid(),mp->desc(),n,n,mp->is_complex());
			vector<int> akey;
			int nr = mp->rsize();
			for (int i=1; i<=nr; i++)
				akey.push_back(i);
			if (mp->is_complex()) {
				complex<double> one(1.0);
				complex<double> zero(0.0);
				blas::copy(akey,akey,one,mp->celem(),dof,dof,zero,newmp->celem());
			} else {
				blas::copy(akey,akey,1.0,mp->elem(),dof,dof,0.0,newmp->elem());
			}
			delete mp;
			mp = newmp;
		}

		flaps::info ("    ", *mp);
		mp->store();
	}
	return true;
}

bool
importer (string const& options) {
// Syntax
//   path
//   output=path
//   vzid=path
	T_(Trace trc(2,"importer");)

try {
	vector<Tok*> unrec = flaps::lexer(options, {
		{"^.*", [&](const Tok& p) { return input_f(p); }}
	});
} catch (runtime_error& s) {
	flaps::error(s.what());
}

	return flaps::nerrors();
}

static string
get_format(string const& path) {
// Try to determine the format of the ASCII file in "path":
// 1) if it has a file extension
//    op4:  Nastran output4
//    mat:  Matlab
//    uf:  Universal file
//    mm:  Matrix Market
// 2) read the first line and try to determine the format from it
	T_(Trace trc(1,"get_format");)

	string::size_type idx = path.rfind('.');
	if (idx != string::npos) {
		string ext = path.substr(idx+1);
		if (ext == "op4")
			return "output4";
		else if (ext == "mat")
			return "matlab";
		else if (ext == "uf")
			return "universal";
		else if (ext == "mm")
			return "matrixmarket";
	}
	T_(trc.dprint("no extension - try reading the first line");)

	// Next try  reading the first line and testing it
	regex matrixMarketRe{"^%%MatrixMarket"};
	// The first line of a nastran output4 file is
	//    NC NR NF NTYPE NAME
	// where NR is negative if BIGMAT=TRUE
	char const* output4Pat = "^ +[0-9]+ +-?[0-9]+ +[0-9]+";
	regex output4Re(output4Pat);

#ifdef NEVER // needs work
	// try Matlab (binary file) ...
	if (isMatlabFile(path)) {
		return "matlab";
	}
#endif // NEVER : needs work

	// ASCII file formats - open the file and read the
	// first line, compare against various regular expressions
	ifstream file = ifstream(path);
	if (!file) {
		string fpath = vastr(getftmp(), '/', path);
		file = ifstream(fpath);
		if (!file)
			throw runtime_error(vastr("import file \"",path,"\" is not available"));
	}
	string line;
	getline(file, line);

	T_(trc.dprint("testing \"",line," to determine file format");)

	string filetype;
	if (regex_match(line, output4Re))
		filetype = "output4";
	else if (regex_match(line, matrixMarketRe))
		filetype = "matrixmarket";
	else if (line.substr(0,6) == "    -1")
		filetype = "universal";
	return filetype;  // may be empty
}

void
rename(vector<Matrix*> ml, vector<pair<string,string>> newnames) {
// rename a set of matrices from newnames[i].first to newname[i].second
	T_(Trace trc(1,"rename");)

	for (auto m : ml) {
		for (size_t i=0; i<newnames.size(); i++) {
			T_(trc.dprint("testing ",m->mid()," =? ",newnames[i].first," or ",newnames[i].second);)
			// old=new
			if (m->mid() == newnames[i].first) {
				m->mid(newnames[i].second);
				flaps::info("renamed ",newnames[i].first," to ",m->mid());
				break;
			}
			// new = old
			if (m->mid() == newnames[i].second) {
				m->mid(newnames[i].first);
				flaps::info("renamed ",newnames[i].second," to ",m->mid());
				break;
			}
		}
	}
}
