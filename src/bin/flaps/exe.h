//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef EXE_H
#define EXE_H

#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// register a background process to be waited on when
// bgprocess is called with pid=0 or no arg
void
bgprocess(pid_t pid);

// C++ style interface to fork/execvp:
int
exec (const std::vector<std::string>& arg,
		const std::vector<std::string>& env,
		const std::string& data, const std::string& outfile,
		const std::string& errfile, bool wait);

std::string which(const std::string& cmd, int* place=nullptr);

#endif	/* EXE_H */
