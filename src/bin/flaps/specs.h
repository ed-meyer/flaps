//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef specs_h
#define specs_h

#include <string>
#include <vector>

class Specs {
public:
	int debug{0};
	std::string debugger;
};

Specs& flaps_specs();	// ref to the global specs
								//
// parser is declared in functions.h:
//!! void parser(const std::string&);

std::ostream& operator<<(std::ostream& s, const Specs& t);

//!! void parser (std::string const& prog);

//!! void print_prefs(std::string const& prog);
#endif // specs_h
