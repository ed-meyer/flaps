//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "config.h"
#include <unistd.h>

#include "cmds.h"
#include "fio.h"
#include "functions.h"
#include "load.h"
#include "main.h" // setup signal and fp exception handling
#include "Pz.h"
#include "Regex.h"
#include "trace.h"
#include "unistd.h"

using namespace std;


static string cmdline(int argc, char *argv[]);
string makePattern (vector<string> const& fcns_or_progs);

bool NoExe{false};
// run the commands, return total cpu time
// valid programs are returned by prog_names, valid functions
// are returned by fcn_names() in functions.cpp.
double do_cmds(vector<pair<string,string> >& cmds);
vector<string>& prog_names();

bool ignoreError(string const& fcn);

bool isBlank(string const& line);
string split (const string& file);

int
main (int argc, char *argv[]) {
// Flaps executive program: read the input file, split blocks,
// and run each function or program
	T_(Trace trc(2,"flaps");)
	int rval{0};
	Ftmpdir ftmp;		// create a temporary directory for data, etc

	// parse the command line
	// may not return if, e.g. "flaps help", or error
	string file;
	try {
		file = cmdline(argc, argv);
		double cputime{0.0};
		string vz("vz");
		string viz("viz");
		if (getenv("NOVZ") != nullptr) {
			vz = "";
			viz = "";
		}

		// message of the day
		motd();

		// set process group to my pid for "spy"
		if (setpgid(getpid(), getpid()) < 0)
			throw std::system_error(errno, std::generic_category(),
				"cannot set process group id");

		// break the file into blocks, return the program in "input"
		string input = split(file);
		// expand environment variables in input
		input = expenv(input);
		// split the input into commands and their options
		vector<pair<string,string> > cmds = lexcmds(input);
		// run each command
		cputime = do_cmds(cmds);
		cerr << "cpu time: " << cputime << endl;
	} catch (std::exception const& se) {
		cout << "terminating " << file << ": " << se.what() << endl;
		rval = 1;
	}
	// wait for background processes to finish
	bgprocess();
	return rval;
}


static void
usage(const string& prog) {
    cerr << "Usage: " << prog << " [options] [path]\n";
    cerr << "   catalog  - list all objects on flaps savefile or plot file \"path\"\n";
    cerr << "   clean    - delete flaps temporary directories left over from\n"
		 <<   "              previous runs\n";
    cerr << "  -d [n]    - set debug level to n, a number from 1-3 (default: 1)\n";
    cerr << "   help     - open the Flaps manual\n";
    cerr << "  -k        - keep the temporary directory\n";
    cerr << "   kill [n] - kill flaps job \"n\" where \"n\" is the job number from the\n"
		   << "              spy command. Default: all flaps jobs\n";
	 cerr << "   path     - input file\n";
}

static string
cmdline(int argc, char *argv[]) {
// parse the command line, return the input file name
	T_(Trace trc(1,"cmdline");)
	string rval;

	// First convert arguments to a vector of strings skipping argv[0]
	vector<string> args = cstr2tok(argc, argv, 1);

	// check each arg
	for (size_t i=0; i<args.size(); i++) {
		string arg = args[i];
		string next;
		if (i+1 < args.size())
			next = args[i+1];
		if (arg == "catalog") {
			if (next.empty()) {
				cerr << "Syntax: flaps catalog file\n"
					<< "  where \"file\" is a flaps savefile or plot file\n";
				exit(1);
			}
			string cmd = vastr("inv ",next);
			system(cmd.c_str());
			exit(0);
		} else if (arg == "clean") {
			string cmd{"fclean"};
			if (next == "-r")
				cmd += " -r";
			system(cmd.c_str());
			exit(0);
		} else if (arg == "-d") {
			int level{1};  // default
			if (next[0] != '-')
				level = stoi(next);
			//!! Trace::debug(level);	// XXX where?
			putEnv(vastr("DEBUG=",level));
		} else if (arg.substr(0,2) == "-h" || arg == "help") {
			string cmd = vastr("evince ",getfroot(),"/doc/manual.pdf");
			system(cmd.c_str());
			exit(0);
		} else if (arg == "-k") {
			putEnv("KEEPFTMP=1");
		} else if (arg == "kill") {
			string cmd{"fkill"};
			for (size_t j=i+1; j<args.size(); j++)
				cmd += " " + args[j];
			system(cmd.c_str());
			exit(0);
		} else if (arg[0] == '-') {
			flaps::error ("unrecognized cmd line arg (",arg,")");
			usage("flaps");
			exit(1);
		} else {
			rval = arg;
			if (rsubstr(rval, 3) != ".fp")
				rval += ".fp";
			// check for existance and read access
			if (access(rval.c_str(), R_OK) == -1)
				throw std::system_error(errno, std::generic_category(),
					vastr("input file \"",rval,"\" is not available"));
		}
	}

	return rval;
}

string
makePattern (vector<string> const& cmds) {
// Create an ECMAscript regex pattern for finding command
// names using the list of names in "cmds".
// Allow for expressions like
//   flut:vso {...}
// and explicit paths:
// if (!/home/flut {mass=gmass001, ... }
// note that the portion of the regex that treats the directories
// matches the longest string starting with a slash and ending
// with a slash - that is, it does not have the usual [^/]* character
// class.
// Note the path must be a full path: it must start with /
// the pattern returned has a few sub-matches, e.g. for the above example:
// sub-match 0:  entire string
// sub-match 1:  if(!
// sub-match 2:  /home/
// sub-match 3:  flut
	T_(Trace trc(1,"makePattern");)
	// string start(R"(^\s*(if\s*\(\s*!?)?(/.*/)?()");
	string start(R"(()");
	// string end(R"()\s*\{.*$)");
	// string end(R"())");
	string end(R"()(:.*)?)");
	string rval(start);
	for (size_t i=0; i<cmds.size(); i++) {
		if (i > 0)
			rval.push_back('|');
		rval += cmds[i];
	}
	rval += end;
	T_(trc.dprint("returning ",rval);)
	return rval;
}

double
do_cmds(vector<pair<string,string> >& cmds) {
// call a set of functions or programs depending on what "cmd" is
// Returns the total cpu time for all the cmds
	T_(Trace trc(1,"do_cmds");)
	double rval{0.0};

	regex programRe = make_regex (makePattern(prog_names()));
	regex functionRe = make_regex (makePattern(fcn_names()));

	// check that each cmd is recognized
	for (auto& cp : cmds) {
		if (!regex_match(cp.first, programRe) &&
				!regex_match(cp.first, functionRe)) {
			flaps::error(cp.first," is not a valid Flaps command");
			return rval;
		}
	}
	
	// if the user included the -n option we are done
	if (NoExe)
		return rval;

	for (auto& cp : cmds) {
		string cmd = cp.first;
		// program? ...
		if (regex_match(cmd, programRe)) {
			bool ignoreErr = ignoreError(cmd);
			// run in background? only vz, viz, and display for now
			bool wait{true};
			if (cmd == "vz" || cmd == "viz" || cmd == "display") {
				wait = false;
				if (getenv("NOVZ") != nullptr)
					continue;
			}

			double cputime;
			int status = exec_prog (cmd, cp.second, wait, ignoreErr, cputime);
			rval += cputime;
			if (status != 0)
				throw runtime_error(vastr(cmd," had ",status," errors"));

		// ... or function?
		} else if (regex_match(cmd, functionRe)) {
			bool ignore_err{false};
			call_function (cmd, cp.second, ignore_err);
		} else {
			string exc = vastr("\"",cmd,"\" is not a valid command");
			throw runtime_error(exc);
		}
	}

	return rval;
}

vector<string>&
prog_names() {
// Returns a vector of the names of all flaps programs which
// may be called from a flaps script
	static vector<string> rval{"call", "display", "fem", "flut",
		"octlab", "parameter(s)?", "pz", "td", "viz", "vz"};
	return rval;
}

bool
ignoreError(string const& fcn) {
// Errors in certain programs and functions  should always be ignored
	vector<string> ign{"amviz", "display", "viz", "vz"};
	if (std::find(ign.begin(), ign.end(), fcn) != ign.end())
		return true;
	return false;
}

bool
isBlank(string const& line) {
	string::size_type idx = line.find_first_not_of(" \t\n");
	return (idx == string::npos);
}

string
split (const string& file) {
// split the input flaps file into blocks enclosed by, e.g.
//    stiffcn {{
//       ...
//    }}
// Returns the first block which is the control program
// All other blocks are written to FTMP with filenames the same
// as the block name
	T_(Trace trc(1,"split ",file);)
	ostringstream os;
	
	// if "file" is empty take the input from cin
	// The actual file MUST have a .fp extension, but the
	// user only needs to supply the base name
	istream* inp = &cin;
	if (!file.empty()) {
		inp = new ifstream(file);
		if (!inp)
			throw std::system_error(errno, std::generic_category(),
				vastr("cannot open ",file));
	}

	regex endrx = make_regex("^[ \t]*end[ \t]*$");
	regex cpprx = make_regex("^#(include|define|ifdef|endif)");
	// line consisting of only a comment: #
	regex commentrx = make_regex("\\s*#.*$");
	// #ifdef env_var
	regex ifdefrx = make_regex("^\\s*#ifdef\\s*(\\w*)");
	regex endifrx = make_regex("^#endif.*$");
	string line;
	while (inp->good()) {
		getline(*inp, line);
		// strip leading whitespace...
		string::size_type first{0};
		while (first < line.size() && std::isblank(line[first])) first++;
		if (first > 0)
			line = line.substr(first);
		// strip comments - watch out for #ifdef...
		if (regex_match(line, commentrx)) {
			T_(trc.dprint("got comment: ",line);)
			// if this is a #ifdef look in the environment
			// for the macro - if not there skip lines until #endif
			smatch m;
			if (regex_match(line, m, ifdefrx)) {
				string var = m[1];
				T_(trc.dprint("got macro ",var);)
				// macro not defined - skip lines until #endif
				if (getenv(var.c_str()) == nullptr) {
					T_(trc.dprint("env var ",var," not found: skip line until endif");)
					while(inp->good()) {
						getline(*inp, line);
						if (regex_match(line, endifrx)) {
							T_(trc.dprint("got endif");)
							break;
						}
					}
				}
			}
			continue;
		}
		// ... or trailing comment
		string::size_type idx = line.find('#');
		if (idx != string::npos) {
			// drop the trailing comment
			line = line.substr(0, idx);
		}
		// end line?
		if (regex_match(line, endrx)) {
			break;
		}
		// getline returns the line without the newline -  add it
		if (!isBlank(line)) {
			if (line[0] == '!')
				os << "unix_cmd{" << line.substr(1) << "}\n";
			else
				os << " " << line << endl;
		}
	}
	// the input file is now in os.str()
	string rval = os.str();
	T_(trc.dprint("input file {\n",rval,"\n}\n");)
	if (rval.empty())
		throw runtime_error("input file is empty");


	// Regular expressions for finding blocks
	// ECMAscript regex, case sensitive, no newline
	regex beginblock = make_regex(R"(^\s*([_[:alpha:]][.\w]*)\s*\{\{[\s\n]*$)");
	regex endblock = make_regex(R"(^\s*}}[\s\n]*$)");
	smatch mch;

	// Read the rest of the file searching for blocks of
	// code delimited like this:
	// mysubb {{
	// ...
	// }}
	// Each block will be written into the FTMP/data directory with
	// the block name. Recognized block types are
	// 1) C++ Functions for creating or modifying matrices: .cpp extension,
	//    and the name must match the main entry point (see pz 'code' option)
	// 2) Matlab m files: .m extension, used in the octlab command
	// 3) Function for create a describing-function: no extension required
	vector<string> entrypts;
	vector<string> df_blocks;
	while (!inp->eof()) {
		std::getline(*inp, line);
		if (regex_match(line, mch, beginblock)) {
			string filename = vastr(getftmp(),"/",mch[1]);
			T_(trc.dprint("got block name <",filename,">");)
			// does the block name have an extension?
			string ext;
			string::size_type idx = filename.rfind('.');
			if (idx != string::npos)
				ext = filename.substr(idx+1);
			ofstream ofs(filename);
			string block_type;
			while (!inp->eof()) {
				std::getline(*inp, line);
				if (line.find("df_qs") != string::npos)
					block_type = "df";
				line = expenv(line);
				if (regex_match(line, endblock)) {
					break;
				}
				ofs << line << endl;
			}
			// if the extension is .cpp the entry pt is the name
			// without the extension; otherwise if it is not a df block
			// leave it on FTMP with the block name
			if (block_type == "df")
				df_blocks.push_back(filename);
			else if (ext == "cpp")
				entrypts.push_back(filename);
				//!! entrypts.push_back(stringBasename(filename));
		}
	}

	// accumulate all preload library paths
	vector<string> libs;

	// compile all entrypts into ftmp/custom.so
	// do this even if entrypts is empty to build an empty custom.so,
	// and have LD_PRELOAD set for child processes
	Custom::add(entrypts);
	libs.push_back(vastr(getftmp(), "/custom.so"));

	// ... and df blocks into dflib
	string dflib = load_df(df_blocks);
	if (!dflib.empty())
		libs.push_back(vastr(getftmp(),"/",dflib));

	// Set LD_PRELOAD to the libraries so they will be available to all
	// child processes; the shared objects created may be missing
	// symbols so add libflaps.so also
	if (!libs.empty()) {
		string sep{""};
		os.str("");
		os << "LD_PRELOAD=";
		for (auto& lib : libs) {
			os << sep << lib;
			sep = ":";
		}
		os << sep << getfroot() << "/lib/libflaps.so";
		putEnv(os.str());
	}

	return rval;
}
