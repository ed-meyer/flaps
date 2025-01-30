//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// a C++-style interface to fork/execvp

#include <string>
#include <unistd.h>
#include <sys/wait.h>
#include <sys/types.h> // open
#include <sys/stat.h> // open
#include <fcntl.h> // open
// #include <arpa/inet.h>	/* for inet_ntoa() */
#include <iostream>
#include <fstream>

#include "csignals.h"
#include "exe.h"
#include "fio.h"  // get_cwd()
#include "settings.h"
#include "trace.h"

using namespace std;

void
redirect(int fd, const string& file);
void
addToEnvironment(char* const *envp);

/*--------------------------------------------------------------------
 * exec - spawn a process, pipe data to its stdin & (optionally) wait
 *        Since this function uses execvp to run the process, it
 *        works either with executables or shell scripts (ref Stevens,
 *        APUI p 208)
 * Input
 *    arg arg[0] is the path to the program, or the name of a debugger
 * Returns
 *    -1 if the spawn was not successful
 *    exit status of the child if wait=true
 *    child pid if wait=false
 *	Discusion
 *		To get a string with an explanation of the exit status use
 *		exitstatus(int status) (text.h)
 *------------------------------------------------------------------*/
int
exec (const vector<string>& arg, const vector<string>& env,
	const string& data, const string& outfile, const string& errfile, bool wait) {

	T_(Trace trc(2,"exec");)
   int pfd[2];
	int rval{0};
	string exc;

	T_(trc.dprint("my pid ",getpid()," path ",arg[0],", ",env.size()," extra env");)
	T_(trc.dprint(", data <",data, "> out file<",outfile,"> errfile<",errfile);)
	T_(trc.dprint("> wait? ",wait, " args ",arg);)

	pfd[0] = pfd[1] = -1;

	string path = arg[0];

	// set up a pipe to send data to child's stdin
   if (!data.empty()) {
		rval = pipe(pfd);
		if (rval == -1)
			throw std::system_error(errno, std::generic_category(),
				vastr("cannot open pipe to ",path));

		T_(trc.dprint("opened fd ",pfd[0]," and ",pfd[1]," for piping data");)
	}

	// flush stdout/err prior to fork to avoid duplicate outputs
	fflush(stdout);
	fflush(stderr);

	pid_t childpid = fork();

	rval = childpid;

	// XXX throw exception or return -1?
   if (childpid == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("cannot fork process ",path));

	// Child process: redirect stdin to the pipe if there is data to read,
	// close stdin otherwise. Redirect stdout, stderr if requested
	// If DEBUGGER is in the environment (or the additional
	// environment envp), start the process under the debugger
	// specified by DEBUGGER.
   if (childpid == 0) {
		//  Ignore interupts if running in background
		if (!wait)
			Signal::ignore(SIGINT);
		// put extra env in the environment
		for (auto ei : env)
			putEnv(ei);
		cout.flush();
		// redirect file descriptor 1 (stdout) to "outfile"...
		redirect(STDOUT_FILENO, outfile);
		// ... and fd 2 to errfile
		redirect(STDERR_FILENO, errfile);

		// attach stdin to the pipe
      if (!data.empty()) {
			if (dup2 (pfd[0], STDIN_FILENO) == -1)
				throw std::system_error(errno, std::generic_category(),
					vastr("cannot open pipe to ",path));

			if (pfd[0] != STDIN_FILENO)
				if (close(pfd[0]) == -1)
					throw std::system_error(errno, std::generic_category(),
						vastr("cannot close pipe to ",path));

			if (pfd[1] != STDIN_FILENO)
				if (close(pfd[1]) == -1)
					throw std::system_error(errno, std::generic_category(),
						vastr("cannot close pipe to ",path));
      } else {
			T_(trc.dprint("no data to send: redirect stdin to /dev/null");)
			if (freopen("/dev/null", "r", stdin) == 0)
				flaps::warning("cannot redirect cin for ",path,": ",strerror(errno));
		}
		T_(trc.dprint("calling execvp argv<",arg,">");)

		// execvp needs char** argv
		vector<const char*> argv;
		for (auto& ai : arg)
			argv.push_back(ai.c_str());
		argv.push_back(nullptr);

      execvp (argv[0], (char* const*)&argv[0]);		// should not return

		throw std::system_error(errno, std::generic_category(),
			vastr("cannot execute \"", path, "\""));
   }

   T_(trc.dprint("Process id for <",path,"> is ",childpid);)

	// Parent process (childpid != 0): pipe data to child
   if (!data.empty()) {
		if (close(pfd[0]) == -1)
			throw std::system_error(errno, std::generic_category(),
				vastr("cannot close pfd[0] to \"", path, "\""));

		T_(trc.dprint("writing ",data.size(),"-char data to childs input pipe");)
      if (write (pfd[1], data.c_str(), data.size()) == -1)
			throw std::system_error(errno, std::generic_category(),
				vastr("cannot send data to \"", path, "\""));

      if (close(pfd[1]) == -1)
			throw std::system_error(errno, std::generic_category(),
				vastr("cannot close pfd[1] to \"", path, "\""));
   }

	if (wait) {
		int status;
		rval = 0;
		while (waitpid(childpid, &status, 0) < 0) {
			if (errno != EINTR)
				throw std::system_error(errno, std::generic_category(),
					vastr("waiting for \"",path,"\""));
		}
		if (rval != -1) {
			T_(trc.dprint(getpid(),": waitpid(",childpid,") exit status 0x",status);)
			rval = status;
		}
	} else {
		rval = childpid;
	}

   T_(trc.dprint("returning ",rval);)
   return rval;
}

void
redirect(int fd, const string& file) {
	T_(Trace trc(1,"redirect");)
	string exc;

	if (file.empty())
		return;

	// open the file with low-level C fcn: returns a file descriptor
	int newfd = open(file.c_str(),O_WRONLY|O_CREAT|O_TRUNC, S_IRWXU);
	if (newfd == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("cannot open ",file));

	// copy the new file descriptor to "fd"
	if (dup2(newfd, fd) == -1)
		throw std::system_error(errno, std::generic_category(),
			vastr("cannot redirect fd ",fd," to ",file));
	// close unused file descriptor
	close(newfd);
}

void
addToEnvironment(char* const *envp) {
	T_(Trace trc(1,"addToEnvironment");)
	int i;
	char* ep;
	if (envp == 0)
		return;
	for (i=0; envp[i]; i++) {
		ep = strdup(envp[i]);
		T_(trc.dprint("added \"",ep,"\" to environment");)
		putenv(ep);
	}
}

string
which (const string& cmd, int* place) {
// Returns the full path of a command, or an empty string if not found
	T_(Trace trc(3,"which ",cmd);)
   char *path = getenv("PATH");
	string rval;
	string empty;

	errno = 0;  // in errno.h

	if (place != nullptr)  // default place: not found
		*place = -1;

	// if the command is already a full path,
	// return a copy of it
	if (cmd[0] ==  '/') {
      if (access(cmd.c_str(), X_OK) == 0) {
			T_(trc.dprint("returning ",cmd);)
			if (place != nullptr)
				*place = 0;
			return cmd;
		} else {
		// the file may exist but not have execute permission: how to
		// tell the user?
			if (access(cmd.c_str(), R_OK) == 0)
				flaps::warning(cmd," does not have execute permission");
			return empty;
		}
	}

	// search current directory if no PATH, return full path in
	// case we can't exec files on local directory
   if (path == nullptr) {
		if (access(cmd.c_str(), X_OK) == 0) {
			rval = vastr(get_cwd(),'/',cmd);
			T_(trc.dprint("returning ",rval);)
			if (place != nullptr)
				*place = 0;
			return rval;
		} else {
			if (access(cmd.c_str(), R_OK) == 0)
				flaps::warning(cmd," does not have execute permission");
			return rval;
		}
	}

	string pathstr(path);
	int k{0};
   T_(trc.dprint("searching for ",cmd," on path ",pathstr);)
	vector<string> paths = string2tok(pathstr, ":");
	for (auto& pi : paths) {
		string testpath = vastr(pi,'/',cmd);
		struct stat statbuf;
		if (access(testpath.c_str(), X_OK) == 0) {
			if (stat(testpath.c_str(), &statbuf) == -1)
				throw std::system_error(errno, std::generic_category(),
					vastr("cannot obtain status of ", testpath));

			if (S_ISREG(statbuf.st_mode)) {
				T_(trc.dprint("returning <%s>",testpath);)
				if (place != nullptr)
					*place = k;
				return testpath;
			} else {
				string exc = vastr(testpath, " is not a regular file");
				T_(trc.dprint("throwing exception: ",exc);)
				throw runtime_error(exc);
			}
		}
		k++;
	}

	// not found - return empty
	T_(trc.dprint("returning empty");)
	return rval;
}


#ifdef MAIN

int
main(int argc, char **argv) {
	if (argc < 2) {
		cerr << "Usage: exe prog [arg ...]\n";
		exit(1);
	}
	vector<string> arg;
	for (int i=1; i<argc; i++) {
		arg.push_back(string(argv[i]));
	}
	vector<string> env;

	string data("");

	// test which
	int place;
	cout << "executing " << which(arg[0], &place) << "(" << place << ")\n";
	// test new c++ exec...
	exec (arg, env, data, "exe.out", "exe.err", true);

}

#endif // MAIN
