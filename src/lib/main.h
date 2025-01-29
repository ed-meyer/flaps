//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// A header meant to be included in each flaps main program,
// i.e. each main.cpp
// Set up some signal-trapping, and pause if WAITFORDEBUG is in the
// environment
#ifndef MAIN_H
#define MAIN_H

#include "config.h"
#include <cfenv> // for feenableexcept
#include <cstdlib> // for getenv
#include <iostream> // for cout, cerr
#include <stdexcept> // for runtime_error
#include <string>
#include <sys/types.h> // for getpid
#include <unistd.h>  // for getpid, sleep
#include "csignals.h"	// my cover for signals.h
#include "trace.h"

class Main {
public:
	Main() {
		// if DEBUG is in the environment set Trace::debug()
		char* dbg = getenv("DEBUG");
		if (dbg != nullptr) Trace::debug(std::stoi(dbg));
		// if WAITFORDEBUG is in the environment pause until the user
		// starts gdb with this program, attaches my pid, then types
		// 'signal SIGINT' in gdb. Use a lambda for the signal handler
		// which just returns
		if (getenv("WAITFORDEBUG")) {
			Signal::set(SIGINT, [](int sig){return;});
			std::cerr << "waiting for gdb to\n";
			std::cerr << " attach " << getpid() << std::endl;
			std::cerr << "then type \"signal SIGINT\"";
			sleep(14400);  // 4 hours
		}
		// use a lambda to handle signals, throw an exception with a
		// description of the signal
#ifdef NEVER // cannot return from the handler: infinite loop
		Signal::catchall([](int sig) {
				// std::cerr << "got signal " << sig << ": " << Signal::desc(sig) << std::endl;
				std::cerr << "got signal " << sig  << std::endl;
			// throw std::runtime_error(vastr("got signal ",sig,": ",Signal::desc(sig)));
			// XXX don't exit - maybe we are already in exit(1);
		});
#endif // NEVER : cannot return from the handler: infinite loop

		// turn on floating-point exception handling for a few floating-point
		// errors: divide-by-zero, overflow; assuming underflow just goes to
		// zero, and inexact is too common. Note: lapack fcn drscl.f:134
		// underflows easily. Even though these are called *exceptions* they
		// are NOT C++ exceptions, just signals, so the signal-trapping lambda
		// above will catch them and throw a C++ exception.
		feenableexcept(FE_DIVBYZERO | FE_OVERFLOW);
	}
};

// one-and-only instance, assuming this header is included
// only once per program
Main setup;

#endif // MAIN_H
