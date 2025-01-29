//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Some functions to make handling Posix signals easier in C++
#ifndef CSIGNALS_H
#define CSIGNALS_H

#include <csignal>
#include <string>

#include "message.h"

class Signal {
public:
	// text to signal number and vice-versa:
	static std::string desc(int signo);
	static std::string macro(int signo);
	static int number(const std::string& macro);

	// set signal actions:
	static void set(int sig, void (*handler)(int));
	static void ignore(int sig);
	static void def(int sig);
	static void block(int sig);
	static void unblock(int sig);

	// set most signals to "handler"
	static void catchall(void (*handler)(int));

	// test the action for "sig"
	static bool istrapped(int sig);
	static bool ispending(int sig);
};

#endif	/* CSIGNALS_H */
