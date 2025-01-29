//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Functions for executing commands (programs and library functions)
// Public interface:
// double exec_prog: execute a program and wait for it or register it
//                 as a background process
// bool call_function: call a library or local function
// bool unix_command:  run a unix command
//
// References
// 1) this is known as APUI:
//   @book{stevens2008advanced,
//     title={Advanced programming in the UNIX environment},
//     author={Stevens, W Richard and Rago, Stephen A},
//     year={2008},
//     publisher={Addison-Wesley}
//   }

#include <chrono>
#include <sys/types.h>
#include <sys/wait.h>

#include "config.h"
#include "cmds.h"
#include "fio.h"
#include "functions.h"
#include "lexer.h"
#include "run_program.h"
#include "settings.h"
#include "specs.h"
#include "trace.h"

using namespace std;

static void
preamble (const std::string& task, const std::string& options, bool waiton);
static void
colophon (const std::string& task, int nerror, bool waiton,
		pid_t pid, const double& cputime);
static void
unix_preamble (std::string const& commandline);


void
bgprocess(pid_t pid) {
// Register a background process so they can be waited on when
// flaps is finished. If pid == 0 each registered process is
// waited on
	Trace trc(1,"bgprocess ", pid);
	static vector<pid_t> bgprocs;
	if (pid > 0)
		bgprocs.push_back(pid);
	else {
		for (auto p : bgprocs) {
			int stat;
			trc.dprint("waiting on bg proc ",p);
			waitpid(p, &stat, 0);
			trc.dprint("pid ",p," terminated with status ",stat);
			if (stat != 0) {
				cout << p << " terminated abnormally: "
					<< exitstatus(stat) << " (see stderr)\n";
			}
		}
	}
}

int
exec_prog (const string& prog, const string& options,
		bool waiton, bool ignoreErr, double& cputime) {
/*--------------------------------------------------------------------
 * Print Flaps cmd header, exec cmd, pipe option list, and wait.
 * Input
 *   prog  program name, possibly the full path
 *
 * Returns:		the exit status of the program
 *------------------------------------------------------------------*/
	Trace trc(1,"exec_prog");
	string name;
	string exc;
	ostringstream os;
	Specs& sp = flaps_specs();

	// get program name without the directory
	string::size_type idx = prog.rfind('/');
	if (idx == string::npos)
		name = prog;
	else
		name = prog.substr(idx+1);

	// Expand environment variables in the options
	string optstring = expenv(options);

	// Programs and scripts are on $froot/bin
	// XXX no - use the users PATH
	string progpath = vastr(getEnv("FROOT"),"/bin/",prog);

	// Print a header with the execution options, etc...
	preamble (prog, optstring, waiton);

	// ... and finally run the program
	pid_t childpid;
	int status = run_program (prog, optstring, waiton, childpid, cputime,
		sp.debug, sp.debugger);
	// ... and add the pid to the cerr message if it is >0 (background)
	if (childpid > 0)
		cerr << " (" << childpid << ")\n";

	// Pause if it is running in the background so its output
	// doesnt get mingled with the next process's
	// Also, register the child's pid for later waiting
	if (!waiton) {
		bgprocess(childpid);
		sleep(2);
	}

	// ... and the program trailer
	colophon (name, status, waiton, childpid, cputime);

	// Throw an exception if there was an error and
	// ignoreErr is false
	if (!ignoreErr && status != 0) {
		string err = prog + string(" failed");
		trc.dprint("throwing exception: ",err);
		throw (runtime_error(err));
	}

	// if this prog left a gmon.out run it through gprof
	if (access("gmon.out", R_OK) != -1) {
		os.str("");
		os << "gprof " << progpath << " >" << prog << '.' << childpid << ".gprof";
		system(os.str().c_str());
	}
	return status;
}

bool
call_function (const string& name, const string& options, bool ignoreErr) {
	Trace trc(1,"call_function");
	Fcnptr fun = fcnptr(name);

	// Expand environment variables in the options
	expenv(options);

	preamble (name, options, true);
	// call the function
	bool rval = fun(options);

	int nerr = flaps::nerrors();
	colophon(name, nerr, true, 0, 0.0);
	if (nerr > 0 && !ignoreErr)
		throw runtime_error(vastr(name, " failed"));

	return rval;
}

double
callRealFunction (const char *funName, double (*fun)(char const*), char const* opt, bool ignoreErr) {
	Trace trc(1,"callRealFunction");
	preamble (funName, opt, true);
	double rval = fun(opt);
	int nerr = flaps::nerrors();
	colophon(funName, nerr, true, 0, 0.0);
	if (nerr > 0 && !ignoreErr)
		throw runtime_error(vastr(funName, " failed"));
	return rval;
}

bool
unix_command (const string& cmd, bool ignoreErr) {
	Trace trc(1,"unix_command");
	bool rval = true;
	ostringstream os;

	// Expand environment variables in the command line
	string cmdline(cmd);
	expenv(cmdline);

	unix_preamble (cmdline);
	int status = system(cmdline.c_str());

	if (status != 0) {
		string exc = vastr(cmd, " returned status ",status);
		trc.dprint("unix cmd error: : ",exc);
		flaps::error(exc);
		return false;
	}

	colophon(cmdline, flaps::nerrors(), true, 0, 0.0);

	return rval;
}

static
void
preamble (const string& progpath, const string& options, bool waiton) {
/*------------------------------------------------------------------
 * print a message of the form:
 *
 *    ------------------  Flaps Command path ----------
 *    cmd {
 *      option1=rhs1
 *         ...
 *    ------------------------------------------------------
 *
 *------------------------------------------------------------------*/
	Trace trc(1,"preamble");
	string task;
	string line;

	trc.dprint("prog path ",progpath,", options: ",options);

	// task is the name of the command without path info
	string::size_type idx = progpath.rfind("/");
	if (idx != string::npos)
		task = progpath.substr(idx+1);
	else
		task = progpath;

	cout << "\n\n ------------  Flaps command " <<
				progpath << " -----------------\n";
	cout << "   " << task;

	cout << " {\n";
	if (!options.empty()) {
		vector<Tok*> ol = flaps::lexer(options);
		for (auto opt : ol)
			cout << "       " << *opt << endl;
	}

   cout << " --------------------------------------------------------\n\n";
	cout.flush();

	// Print a message in stderr
	cerr << "starting " << task;
	if (waiton)
		cerr << "... ";

	return;
}

static
void
unix_preamble (string const& commandline) {
/*------------------------------------------------------------------
 * print a message of the form:
 *
 *    ------------------  Unix Command commandline  ---------- {
 *
 *------------------------------------------------------------------*/
	Trace trc(1,"unix_preamble");

	cout << "\n\n ------------  Unix Command " <<
				commandline << "   -----------------{\n";

	cout.flush();

	// Print a message in stderr
	cerr << "starting " << commandline << "... ";

	return;
}

static void
colophon (const string& task, int status, bool waiton,
		pid_t pid, const double& cputime) {
// print messages in cout, cerr regarding the completion of a task
	Trace trc(1,"colophon");

	if (task.empty())
		return;

	if (waiton) {
		cout << "} --------------- Flaps command \"" << task
			<< "\" (" << pid << ") finished";
		if (status == 0) {
			cout << " (" << cputime << " sec)\n";
			cerr << " (" << cputime << ")\n";
		} else {
			cout << " with " << exitstatus(status) << "\n";
			cerr << " " << exitstatus(status) << endl;
		}
	} else
		cout << "} --------------- Flaps process " << pid << " \"" << task
			<< "\" is running in the background\n";
	fflush (stdout);

	return;
}

void
wait_for_group() {
// Wait for all members of my process group to finish...
// cannot use waitpid() if we have grandchildren
// Note: this function only works when called from a process
// named "control"
	Trace trc(1,"wait_for group");
	char buf[256];
	pid_t pid = getpid();
	ostringstream os;

	// Open a pipe to the Unix (POSIX) command "ps -eo pgid,comm"
	// then read each line, check it's process group (pgid)
	const char* ps = "ps -eo pgid,comm";
	// XXX I don't understand why ps -g pid does not work
	// os << "ps -g " << pid;
	// const char* ps = os.str().c_str();
	//
	// os << "ps -e -o \"%r %c\"";

	trc.dprint("searching for pgid ",pid);
	// keep running ps until the only processes left in
	// my group are me (control), ps, or sh
	while(1) {
		// FILE* pipe = popen(os.str().c_str(), "r");
		trc.dprint("popen(",ps,")");
		FILE* pipe = popen(ps, "r");

		while (1) {
			char* cp = fgets(buf, sizeof(buf), pipe);
			if (!cp) {
				trc.dprint("no more processes: done");
				pclose(pipe);
				return;
			}
			string line(buf);
			vector<string> toks = string2tok(line, " \t\n");
			trc.dprint("read ",buf);
			int pgid;
			// more than 2 items on a line? probably <defunct>,
			// or if the 1st item is not an int, ignore this line
			if (toks.size() != 2 || !str2int(toks[0], pgid))
				continue;
			// if one of my group is still running break, pause, and re-try
			if (pgid == (int)pid) {
				trc.dprint("found a proc with my gpid: ",buf);
				if (toks[1] != "control" && toks[1] != "ps" && toks[1] != "sh")
					break;
			}
		}
		pclose(pipe);
		sleep(1);
	}
}
