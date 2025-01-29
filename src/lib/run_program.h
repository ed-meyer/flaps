//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef RUN_PROGRAM_H
#define RUN_PROGRAM_H

#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

// Run a program, pipe data to it's stdin, and either wait for it to
// finish or run it in the background.
int
run_program (const std::string& program, const std::string& options, bool waiton,
		pid_t& childpid, double& cputime, int debug=0, std::string const& debugger="");

// create a gdb.in containing options for gdb startup
void
write_gdbin(const std::string& options, const std::string& debugger);

// C++ style interface to fork/execvp:
int
exec (const std::vector<std::string>& arg,
		const std::vector<std::string>& env,
		const std::string& data, const std::string& outfile,
		const std::string& errfile, bool wait);

std::string which(const std::string& cmd, int* place=nullptr);

#endif	// RUN_PROGRAM_H
