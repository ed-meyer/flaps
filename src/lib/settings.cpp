//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// namespace Settings: functions to maintain a collection
// of settings for this process and child processes

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
#include <vector>

#include "fptype.h"
#include "lexer.h"
#include "message.h"
#include "trace.h"
#include "settings.h"

using namespace std;

#ifndef MAIN

// settings are kept in a map of string keys and std::any data;
// note that std::any requires C++17 (or Boost)
// static std::map<std::string, std::any> list;
// static std::map<std::string, boost::any> Settings::list;
// map<string, std::any>&
// Settings::
// the_map() { return list; }


Settings Settings::defaults;

Settings::
Settings() {
}

Settings::
~Settings() {
}

bool
Settings::
is_defined(const string& name) const {
	auto pos = the_map.find(name);
	if (pos == the_map.end()) {
		return false;
	}
	return true;
}

bool
Settings::
remove(const string& name) {
// delete a setting, return true if successful,
// false if the setting does not exist
	auto pos = the_map.find(name);
	if (pos != the_map.end()) {
		the_map.erase(pos);
		return true;
	}
	return false;
}

bool
is_int(const std::any& t) {
	int s;
	return t.type() == typeid(s);
}

bool
is_string(const std::any& t) {
	std::string s;
	return t.type() == typeid(s);
}

bool
is_vstring(const std::any& t) {
	std::vector<std::string> s;
	return t.type() == typeid( std::vector<std::string>);
}

bool
is_double(const std::any& t) {
	double s;
	return t.type() == typeid(s);
}

std::string
Settings::
datatype(const std::string name) {
	auto pos = the_map.find(name);
	if (pos == the_map.end())
		return "undefined";
#ifdef HAVE_BOOST_CORE_DEMANGLE_HPP
	return boost::core::demangle(pos->second.type().name());
#else
	return pos->second.type().name();
#endif
}

std::ostream&
operator<<(std::ostream& s, const Settings& t) {
	for (auto& ti : t.the_map)
		s << ti.first << " = " << boolalpha << ti.second;
	return s;
}

std::ostream&
operator<<(std::ostream& s, const std::any& t) {
		if (t.type() == typeid(bool))
			s << boolalpha << std::any_cast<bool>(t);
		else if (t.type() == typeid(std::string))
			s << std::any_cast<std::string>(t);
		else if (t.type() == typeid(double))
			s << std::any_cast<double>(t);
		else if (t.type() == typeid(int))
			s << std::any_cast<int>(t);
		else if (t.type() == typeid(std::vector<int>))
			s << std::any_cast<vector<int>>(t);
		else if (t.type() == typeid(std::vector<string>))
			s << std::any_cast<vector<string>>(t);
		else if (t.type() == typeid(std::vector<double>))
			s << std::any_cast<vector<double>>(t);
#ifdef HAVE_BOOST_CORE_DEMANGLE_HPP
		s << " (" << boost::core::demangle(t.type().name()) << ")\n";
#else
		s << " (" << t.type().name() << ")\n";
#endif
	return s;
}

/*--------------------------------------------------------------------
 * Settings::global(string& pref):  set various global settings
 * These same options can be set for an entire job by setting
 * the corresponding environment variable; for example
 *   DEBUGGER=valgrind flaps pv >out 2>err
 * will run each process under the memory debugger "valgrind"
 *
 * Valid options:
 *    units=(si|uscs) use International System of Units (SI) or United
 *                    States Customary System (USCS) for output. Default: SI
 *    env="environment variable"
 *                   set an environment variable, e.g. env="DISPLAY=:0"
 *    d(ebug)?       debug level (0, 1, 2, or 3)
 *    debugger       the name of a debugger, memory checker, or profiler:
 *                   ddd
 *                   valgrind
 *    pagewidth      max width for subsequent printout
 *    o(utput)?      output file name
 *    wait           the next process will print it's process id then
 *                   sleep until it receives a SIGINT signal, usually from
 *                   a debugger after it has been attached to that process.
 *                   For example, to run gdb on flut the input file has:
 *                      settings{wait}
 *                      flut { ... }
 *                   and in another window type
 *                      gdb flut -> attach pid
 *    
 * Debuggers, memory checkers, profilers
 * -------------------------------------
 *
 * valgrind
 * --------
 *    include options in the command, e.g.
 *        settings{memorycheck="valgrind --logfile=valgrind --leak-check=yes
 *              --freelist-vol=10000000 --show-reachable=yes"}
 *
 *    Or create a .valgrindrc file either on HOME or local directory;
 *    my $HOME/.valgrindrc contains:
 *
 *     --tool=memcheck
 *     --memcheck:leak-check=yes
 *     --logfile=valgrind.%p
 *     --memcheck:num-callers=20
 *
 *    Then all you need is
 *       settings{debugger=valgrind}
 *    or
 *       settings{valgrind}
 *------------------------------------------------------------------*/

bool
Settings::
global (string const& preferences) {
// Set some preferences so they will be accessible to all processes
// The allowable settings:
//   output=path     redirect stdout to path
//   env=var         set 'var' in the environment, e.g. "TMP=/tmp"
//   d(ebug)=n       set DEBUG=n in the environment
//   units=SI|USCS   set system of units
//   valgrind        run subsequent processes under valgrind
//   callgrind       run subsequent processes under callgrind
//   time            visualize the output with kcachegrind
//   massif          run subsequent processes under massif
	Trace trc(1,"global");
	size_t i;
	ostringstream os;
	string si_units{"SI"};
	string uscs_units{"USCS"};

	// split the preferences string into Toks
	vector<Tok*> opts = flaps::lexer(preferences);
	trc.dprint("preferences: ",preferences);

	for (i=0; i<opts.size(); i++) {
		trc.dprint("got <", *opts[i]);
		// split the option into tokens
		vector<string> toks = string2tok(opts[i]->lhs, " ");
      if (opts[i]->compare("^o(utput)?$")) {
         fflush (stdout);
			if (opts[i]->svec.empty()) {
				os << "no rhs with the output preference: " << opts[i];
				throw runtime_error(os.str());
			}
         freopen (opts[i]->svec[0].c_str(), "w", stdout);
         trc.dprint("set Output = <",opts[i]->svec[0],">");
      } else if (opts[i]->compare("^env(ironment)?$")) {
			string env{opts[i]->svec[0]};
			// strip quotes: might have PATH='$PATH:/tmp'
			stripquotes(env);
			// ... then expand environment variables
			env = expenv(env);
			putEnv(env);
			flaps::info("setting ",env," in the environment");
      } else if (opts[i]->compare("^d(ebug)?$")) {
			if (opts[i]->rhscompare("^eq.*")) {
				putEnv("DEBUGEQN=1");
				flaps::info("Debugging equations");
			} else {
				int lvl = opts[i]->ivec[0];
				if (lvl < 0) lvl = 0;
				if (lvl > 0) flaps::info ("Debug has been set to ",lvl);
				putEnv (vastr("DEBUG=",lvl));
				// also set it for the current process (flaps)
				Trace::debug(lvl);
			}
      } else if (opts[i]->compare("^unit(s)?")) {
			if (compare(opts[i]->srhs, si_units)) {
				Settings::defaults.set("units", si_units);
			} else if (compare(opts[i]->srhs, uscs_units)) {
				Settings::defaults.set("units", uscs_units);
				putEnv(vastr("units=",uscs_units));
				flaps::info("presentation units will be in USCS");
			} else {
				flaps::warning("unrecognized units: ",opts[i]->srhs);
			}
      } else if (opts[i]->compare("^debugger$")) {
			if (opts[i]->svec.empty()) {
				throw runtime_error("no rhs for debugger preference");
			}
			os.str("");
			os << "DEBUGGER=" << opts[i]->svec[0];
			if (opts[i]->svec[0] == "0" || opts[i]->svec[0] == "none") {
				flaps::info("Subsequent processes will not be run in a debugger");
			}
			flaps::info ("Put \"",os.str(),"\" into the environment");
			putEnv(os.str());
      } else if (toks[0] == "valgrind") {
			string env = vastr("DEBUGGER=",*opts[i]);
			flaps::info("Subsequent processes will be run under ", *opts[i]);
			putEnv(env);
      } else if (opts[i]->compare("^callgrind$")) {
			string env("DEBUGGER=valgrind --tool=callgrind");
			flaps::info("running under callgrind: use kcachegrind to analyze callgrind.*");
			putEnv(env);
      } else if (opts[i]->compare("^massif$")) {
			string env("DEBUGGER=valgrind --tool=massif");
			flaps::info("Subsequent processes will be run under massif");
			putEnv(env);
      } else if (opts[i]->compare("^page[_]?width$")) {
			if (opts[i]->is_int()) {
				os.str("");
				os << "PAGEWIDTH=" << opts[i]->svec[0];
				putEnv (os.str());
			} else {
				flaps::warning("page width must be a positive integer: ", opts[i]);
			}
      } else if (opts[i]->compare("^timer$")) {
			if (opts[i]->is_int()) {
				os.str("");
				os << "TIMER=" << opts[i]->svec[0];
				putEnv (os.str());
			} else {
				flaps::warning("timer must be a positive integer: ", opts[i]);
			}
      } else if (opts[i]->compare("^profile$") || opts[i]->compare("^perf")) {
			// Linux profiler "pref" needs to use a version matching the kernel
			string release = get_uname("release");
			string::size_type idx = release.find('.', release.find('.')+1);
			string cmd = vastr("DEBUGGER=perf_",release.substr(0,idx)," record -g -F 100");
			putEnv(cmd);
			os.str("");
			os << "profiling subsequent programs with \"" << cmd << "\"\n"
				<< "it may be necessary to do the following:\n"
				<< "   sudo sh -c 'echo -1 >/proc/sys/kernel/perf_event_paranoid'\n"
				<< "Run the following when the job completes:\n"
				<< "   perf report -g none|graph -i perf.data\n";

			flaps::info (os.str());
      } else if (opts[i]->compare("^malloc$")) {
			string file("malloc_trace.out");
			if (!opts[i]->svec.empty())
				file = opts[i]->svec[0];
			os.str("");
			os << "MALLOC_TRACE=" << file;
			putEnv (os.str());
			flaps::info("malloc tracing turned on");
			flaps::info("when finished run:\n    mtrace ",file);
      } else if (opts[i]->compare("^wait$")) {
			putEnv ("WAITFORDEBUG=1");
			flaps::info("flaps will pause to wait for gdb to attach");
      } else {
			flaps::error("unrecognized preference (", *opts[i], ')');
		}
   }
	return true;
}

std::string
Settings::
toString() const {
// write all defined settings to a string
	ostringstream os;
	os << the_map;
	// for (auto si : the_map)
	// 	os << si << endl;
	return os.str();
}

#else // MAIN

int
main(int argc, char **argv) {
// Usage:
// settings "preferences_string"
//   the first arg is taken as a preferences string, e.g.
//      settings "a=b, c=d"
	Trace trc(1,argv[0]);

	// test flaps::lexer
	try {
		string preferences;
		if (argc == 2) {
			preferences = argv[1];
			vector<Tok*> ol = flaps::lexer(preferences);
			trc.dprint("got preferences<",ol,">");
			for (auto oi : ol) {
				if (oi->is_char()) {
					if (oi->svec.size() > 1) {
						Settings::defaults.set(oi->lhs, oi->svec);
					} else {
						if (oi->srhs == "true") {
							bool set{true};
							Settings::defaults.set(oi->lhs,set);
						} else if (oi->srhs == "false") {
							bool set{false};
							Settings::defaults.set(oi->lhs,set);
						} else {
							Settings::defaults.set(oi->lhs,oi->srhs);
						}
					}
				} else if (oi->is_int()) {
					if (oi->ivec.size() > 1)
						Settings::defaults.set(oi->lhs, oi->ivec);
					else
						Settings::defaults.set(oi->lhs, oi->ivec[0]);
				} else if (oi->is_real()) {
					Settings::defaults.set(oi->lhs, oi->rvec[0]);
				} else {
					flaps::error("unrecognized pref: ",*oi);
				}
			}
			cout << "currently defined settings:\n" << Settings::defaults << endl;
		} else {
		// Test Settings
			double cf;
			if (Settings::defaults.get("minstepsize", cf))
				cout << "minstepsize: " << cf << endl;
			else
				cout << "minstepsize is not defined\n";
			bool b;
			if (Settings::defaults.get("constrained",b))
				cout << "constrained: " << b << endl;
			else
				cout << "constrained is not defined\n";
			Settings::defaults.set("constrained",true);
			Settings::defaults.get("constrained",b);
			cout << "new constrained: " << b << endl;
		}
	} catch (runtime_error& s) {
		cerr << "caught exception: " << s << endl;
	}

	// test Settings::global
	Settings::global("pagewidth=200, env='path=/tmp'");
}
#endif // MAIN
