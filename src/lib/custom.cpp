//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// class Custom: maintain a map of entry-point names and pointers to the functions

// NOTE: this file must be compiled with custom code included,
//       so lots of headers are required
#include "config.h"  // included here for HAVE_BOOST...

#include <any>
// #include <boost/any.hpp>
#ifdef HAVE_BOOST_CORE_DEMANGLE_HPP
#include <boost/core/demangle.hpp>
#endif
#include <exception>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <unistd.h>
#include <vector>

#include "Ad.h"
#include "conv.h"
#include "custom.h"
#include "df.h"
#include "exim.h"
#include "fio.h"
#include "lapack.h"
#include "lexer.h"
#include "message.h"
#include "pset.h"
#include "Pz.h"
#include "settings.h"
#include "trace.h"

using namespace std;

// Custom map: access with Custom::map()
//!! static std::map<std::string, std::any> the_map;

// include custom source code here
string
test1(const string& str) { return vastr("this is test1(",str,")"); }
double
test2(double t) { cerr << "test2: " << t << endl; return t; }

static std::map<std::string, std::any> the_map {
// add mapitems here
		{"test1", test1},
		{"test2", test2} };

std::map<std::string, std::any>&
Custom::
map() { return the_map; }


// Custom function pointers are kept in a map of string keys and std::any data;
// note that std::any requires C++17 (or Boost)
string
Custom::
add (string const& path) {
// just make a vector and call the vector version
	vector<string> paths;
	paths.push_back(path);
	vector<string> entrypts = Custom::add(paths);
	return entrypts[0];
}

vector<string>
Custom::
add (vector<string> const& names) {
// Build a file, custom.cpp, on FTMP, containing custom C++ functions;
// compile it and build a dynamic-load library, custom.so, made available to
// all subsequent processes via LD_PRELOAD. "names" may be empty in which case
// custom.so is built anyway so that LD_PRELOAD can be set in the flaps program.
// Input: a list of strings which are one of:
//   1) entry point: either the entry point name or the name with a .cpp
//                   extension, already compiled into custom.so (find() returns ptr)
//   2) file names:  with a .cpp extension, and has not been compiled into
//                   custom.so (find() returns nullptr)
//   3) paths:       relative or absolute paths of a files to be compiled
//                   (find() returns nullptr)
// Returns: a list of entry point names
	Trace trc(1,"Custom::add");

	trc.dprint("input names: ",names);

	vector<string> rval;		// entry point names to return
	vector<string> files;	// those of "names" which need compiling

	// 1) check each item with find() to see if it is already in custom.so
	for (auto& fi : names) {
		string entrypt = stringBasename(fi);
		if (Custom::is_defined(entrypt))
			rval.push_back(entrypt);
		else
			files.push_back(fi);
	}

	// lambda to test for read access
	auto readok = [] (const string& file) {
		if (access(file.c_str(), R_OK) == -1)
			return false;
		return true;
	};
	
	// save the cwd in case some files are there instead of ftmp
	string cwd = get_cwd() + '/';

	// change to temp directory using RAII class Chdir
	string ftmp = getftmp();
	Chdir tmp(ftmp);

	// check read access for each entry point file; save the
	// paths in "paths" - some may be on cwd
	vector<string> paths;
	for (auto& ei : files) {
		if (ei[0] == '/' && readok(ei))	// full path
			paths.push_back(ei);
		else if (readok(cwd+ei))
			paths.push_back(cwd+ei);
		else if (readok(ftmp+'/'+ei))
			paths.push_back(ftmp+'/'+ei);
		else
			throw std::system_error(errno, std::generic_category(),
				vastr("custom function (",ei,") is not available"));
	}
			
	trc.dprint(paths.size()," paths:",paths);

	// create or add to ftmp/custom.cpp
	vector<string> entrypts = make_custom(paths);

	// compile the new custom.cpp into custom.so
	string froot = getEnv("FROOT");
	string out = vastr(ftmp,"/custom.out");
	// XXX the commands to compile should come from lib/Makefile
	ostringstream os;
	os << "g++ -std=c++17 -ggdb3 -O0 -Wall -I" << froot
		<< "/include -I" << froot << "/share/libflaps -shared -fPIC -o custom.so "
		<< "custom.cpp >" << out << " 2>&1";
	trc.dprint("running \'", os.str(), "\' to include custom code");
	int stat = system(os.str().c_str());

	if (stat != 0) {
		os.str("");
		os << "building custom.so failed: status " << stat;
		os << " (" << exitstatus(stat) << ")";
		ifstream adds(out);
		if (adds) {
			os << "\ncompilation output:\n";
			char c;
			while(adds.get(c))
				os << c;
		}
		trc.dprint("throwing exception: ",os.str());
		throw runtime_error(os.str());
	}

	if (!entrypts.empty())
		flaps::info("added custom function(s)\n   ",entrypts,
			"\ninto subsequent processes");

	for (auto& ei : entrypts)
		rval.push_back(ei);

	return rval;
}

vector<string>
Custom::
make_custom(const vector<string>& paths) {
// add user's source code to custom.cpp by adding a line to it
// like "#include entrypt", where "entrypt" is the base name (path
// without the directory) of the path and *must* be the name of the
// function which is the primary entry point.
// Function "Custom::open()" will open an existing custom.cpp if
// there is one on the cwd - this will contain custom code from previous
// calls to Custom::add(); otherwise it will return the prototype code
// as a raw string.
	Trace trc(1,"make_custom");

	// create a vector of entry point names from the paths to return
	vector<string> entrypts;
	for (auto& pi : paths)
		entrypts.push_back(stringBasename(pi));

	// open a stream to either an existing custom.cpp or the prototype
	ifstream input = Custom::open();

	// now write the new custom.cpp including the source code
	// with #include statements
	ofstream output("custom.cpp", std::ios::trunc);

	// read each line with getline, then we do these:
	// 1) look for a line containing the string "// include custom ..."
	//    and put the #include's right after it
	// 2) look for a line containing "// insert, e.g.	map_item", put a line like
	//        make_pair("entrypt",entrypt),
	//    where "entrypt" is the name of the entry point and entrypt is a pointer
	//    to the function; in other words two instances of the entry point name
	//    with the first in quotes.
	string line;
	bool gotinc{false};
	bool gotmap{false};
	while (!input.eof()) {
		getline(input, line);
		//!! input.getline(line);
		output << line << endl;
		//!! if (line.find("include custom source here") != string::npos)
		if (!gotinc && line == "// include custom source code here") {
			gotinc = true;
			for (auto& pi : paths)
				output << "#include \"" << pi << "\"\n";
		} else if (!gotmap && line.find("add mapitems here") != string::npos) {
			gotmap = true;
			for (auto& ei : entrypts)
				output << "{\"" << ei << "\", " << ei << "},\n";
		}
	}
	return entrypts;
}

std::ifstream
Custom::
open() {
// Returns an input stream from which custom.cpp can be read; it may
// be either an existing custom.cpp or a prototype
// Note: this function assumes C++11 or later, which has rvalue and
// move semantics for file streams - this means the ifstream created
// here is moved to the calling function so when that function goes
// out of scope the file is closed
	Trace trc(1,"Custom::open");

	// try cwd first, then the prototype string
	string path{"custom.cpp"};
	if (access(path.c_str(), R_OK|W_OK) == -1) {
		trc.dprint("returning prototype custom.cpp");
		return Custom::prototype();
	}

	trc.dprint("returning existing custom.cpp");

	// rename the existing custom.cpp to "tmp", open it and return its ifstream ...
	if (rename("custom.cpp", "tmp") == -1)
		throw std::system_error(errno, std::generic_category(),"rename failed");
	// ... open the renamed file
	ifstream rval("tmp");
	if (!rval) {
		string exc = vastr("cannot open ",get_cwd(),"/tmp");
		throw runtime_error(exc);
	}
	// return the ifstream
	return rval;
}

ifstream
Custom::
prototype() {
// open the initial (prototype) froot/share/libflaps/custom.cpp
	string path = vastr(getfroot(),"/share/libflaps/custom.cpp");
	ifstream rval(path);
	if (!rval)
		throw runtime_error(vastr("cannot open ",path));
	return rval;
}


bool
Custom::
is_defined(const string& name) {
	Trace trc(2,"Custom::is_defined");
	trc.dprint("searching for \"",name,"\"");
	for (auto& ei : Custom::map()) {
		trc.dprint("entrypt: ",ei.first);
		if (ei.first == name)
			return true;
	}
	return false;
// is "name" in the list of entry points?
	auto pos = Custom::map().find(name);
	if (pos == Custom::map().end()) {
		return false;
	}
	return true;
}

std::string
Custom::
datatype(const std::string& name) {
	auto pos = the_map.find(name);
	if (pos == the_map.end())
		return "undefined";
#ifdef HAVE_BOOST_CORE_DEMANGLE_HPP
	return boost::core::demangle(pos->second.type().name());
#else
	return pos->second.type().name();
#endif
}

#ifdef MAIN
//!! string
//!! test1(const string& str) { return vastr("this is test1(",str,")"); }
//!! double
//!! test2(double t) { cerr << "test2: " << t << endl; return t; }
int
main(int argc, char **argv) {
	Trace trc(1,argv[0]);

	//!! Custom::add("test1", test1);
	//!! Custom::add("test2", test2);
	try {
		string (*t)(const string&);
		if (Custom::is_defined("test1")) {
			cerr << "test1 is defined: datatype \"" << Custom::datatype("test1") << "\"\n";
			if (Custom::get("test1", t))
				cerr << "got test1: " << t("test1") << endl;
		}
		double (*q)(double);
		if (Custom::get("test2", q))
			cerr << "got test2: " << q(3.14) << endl;
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}
	return 0;
}
#endif // MAIN
