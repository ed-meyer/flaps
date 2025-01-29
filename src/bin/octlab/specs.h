//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Specs for octlab
#ifndef Specs_h
#define Specs_h

#include <string>
#include <map>
#include <vector>

#include "Par.h"

class Specs {
public:
	bool use_octave{false};
	bool use_matlab{false};
	std::vector<Par*> params;
	std::vector<std::string> input;
	std::vector<std::string> results;
};

std::ostream& operator<<(std::ostream& s, const Specs& t);

// get specs, set various policies
bool parse_specs(std::string const& optionList);

// return a reference to the Specs: an automatic Singleton
Specs& specs();

#endif // Specs_h
