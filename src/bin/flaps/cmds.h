//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef Commands_H
#define Commands_H

#include <string>
#include <vector>

// break the input file into commands+options
std::vector<std::pair<std::string,std::string> > lexcmds(const std::string& in);

int
exec_prog (const std::string& prog, const std::string& options,
		bool waiton, bool ignoreErr, double& cputime);
bool
call_function (const std::string& name,
		const std::string& options, bool ignoreErr);
bool
unix_command (const std::string& cmd, bool ignoreErr);

// register the pid of a background process (pid > 0) so it can
// be waited on, or wait for all registered processes to finish (pid=0
// or no arg)
void
bgprocess(pid_t pid=0);

#endif // COMMANDS_H
