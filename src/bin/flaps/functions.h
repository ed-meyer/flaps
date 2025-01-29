//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <string>
#include <vector>

// To add a function:
// 1) add a decl to this list:
bool catalog (const std::string& specs);
bool call (const std::string& specs);
bool exporter (const std::string& specs);
bool lti (const std::string& specs);
bool gyro (const std::string& specs);
bool importer (const std::string& specs);
bool parser (const std::string& specs);
bool restore (const std::string& specs);
bool save (const std::string& specs);
bool symmetrize (const std::string& specs);
bool unix_cmd (const std::string& specs);

// 2) add to the map in functions.cpp
using Fcnptr = bool (*)(const std::string& specs);

// return a vector of all defined function names
std::vector<std::string>& fcn_names();
// return a pointer to function "name"
Fcnptr fcnptr(const std::string& name);

// message of the day
void motd();

#endif // FUNCTIONS_H
