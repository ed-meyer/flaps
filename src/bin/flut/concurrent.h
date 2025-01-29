//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef CONCURRENT_H
#define CONCURRENT_H

#include <vector>
#include "flutcurve.h"

// trace modes concurrently?
int concurrently(int act=-1);

// see "C++ Concurrency in Action" p. 24 for how to pass curves
// as a reference instead of a pointer
void
concurrent(Fstack* curves, void (*fcn)(Fstack*,int));

#endif // CONCURRENT_H
