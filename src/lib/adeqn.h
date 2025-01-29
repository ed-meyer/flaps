//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Public interface to adeqn.yy, the automatic differentiation
// equation parser. The interface consists of
//   adeqn::eval(pset& plt, string eqn, bool& constant, bool parse_only)

#ifndef Adeqn_h
#define Adeqn_h 1

#include <string>
#include <vector>

#include "Ad.h"
#include "pset.h"

void
setLinearGC(bool to);
bool
linearGC();

namespace adeqn {

// evaluate an equation, return the Ad value; optionally do
// not throw exceptions on error if noexc=true
Ad
eval (pset& plt, std::string const& eqn, bool noexc=true);

std::vector<std::string>
dependson(pset& plt, std::string const& eqn, bool& isComplex);
} // namespace adeqn


#endif // Adeqn_h
