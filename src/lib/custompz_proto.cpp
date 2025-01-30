//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

// Beginning prototype for building a custom.cpp containing custom
// code for creating and modifying matrices
// This file should be installed in $FROOT/share/libflaps

#include <map>
#include <string>
#include <vector>
#include "Ad.h"
#include "conv.h"
#include "df.h"
#include "exim.h"
#include "lapack.h"
#include "pset.h"
#include "Pz.h"
#include "settings.h"
#include "trace.h"

using namespace std;

using Real = Ad;
using Complex = complex<Ad>;
using CMatrix = vector<complex<Ad>>;

// #include custom source here

using map_item = std::pair<std::string, CustomEval>;
CustomEval
CustomPz::
find (const string& entrypt) {
	T_(Trace trc(1,"CustomPz::find", ": ",entrypt);)
      static std::map<string, CustomEval> functions {
	};
	if (functions.empty()) {
		T_(trc.dprint("returning nullptr: no map");)
		return nullptr;
	}
	auto k = functions.find(entrypt);
	if (k == functions.end()) {
		// throw runtime_error(vastr("function \"",entrypt,"\" has not been loaded"));
		T_(trc.dprint("returning nullptr: not in map");)
		return nullptr;
	}
	T_(trc.dprint("returning fcn pointer ",k->second," for ",entrypt);)
	return k->second;
}
