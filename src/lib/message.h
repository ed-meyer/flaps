//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef MESSAGE_H
#define MESSAGE_H

#include <string>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "text.h"
#include "vastr.h"


namespace flaps {

int nerrors (void);
void incr_nerrors();
int nwarnings (void);
void incr_nwarnings();
void error_summary();

// variadic-template functions for printing warning and error messages
template<typename T, typename... Types>
void
warning (const T& first, const Types&... args) {
	std::ostringstream os;
	int indent{3};
	os << "Warning: " << first;
	os << vastr(args...);
	// format the string before printing
	std::string rval = roff(os.str(), indent, page_width()-indent);
	std::cout << rval << "\n\n";
	incr_nwarnings();
}

template<typename T, typename... Types>
void
error (const T& first, const Types&... args) {
	std::ostringstream os;
	int indent{3};
	os << "Error: " << first;
	os << vastr(args...);
	// format the string before printing
	std::string rval = roff(os.str(), indent, page_width()-indent);
	std::cout << rval << "\n\n";
	incr_nerrors();
}

template<typename T, typename... Types>
void
info (const T& first, const Types&... args) {
// variadic-template functions for short messages to the user
	std::ostringstream os;
	os << first;
	os << vastr(args...);
	std::string rval = stringIndent(3, os.str());
	std::cout << rval << "\n\n";
	std::cout.flush();
}
} // namespace flaps


inline
std::ostream&
operator<<(std::ostream& s, const std::runtime_error& t) {
	s << t.what();
	return s;
}

// backtrace attempt to print a stack trace on cerr
std::vector<std::string> backtrace();

#endif // MESSAGE_H
