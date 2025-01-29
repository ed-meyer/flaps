//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef STDATM_H
#define STDATM_H


#include "Ad.h"

namespace atmos {
Ad temp (Ad const& z);
Ad press (Ad const& z);
Ad rho (Ad const& z);
Ad vsonic (Ad const& z);
double altmin();
double altmax();
double rhoref();
double spressref();
Ad vcas(const Ad& vtas, const Ad& alt);
Ad cdpress(const Ad& vtas, const Ad& alt);
Ad delta(const Ad& alt);	// pressure ratio
}

#endif /* STDATM_H */
