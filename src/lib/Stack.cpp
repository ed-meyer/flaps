//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <iostream>
#include <string>

#include "Stack.h"

using namespace std;

// specialization of Stack<string>
std::ostream&
operator<<(std::ostream& s, const Stack<std::string>& t) {
	if (t.empty())
		return s;
	std::string sep;
	for (const auto& ti : t.c) {
		s << sep << ti;
		sep = ", ";
	}
	return s;
}
