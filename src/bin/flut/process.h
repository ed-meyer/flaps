//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef PROCESS_H
#define PROCESS_H

#include "flutcurve.h"

void pre_process(std::string const& prog);

void print_start(std::vector<Flutcurve*>& startpts);

void post_process(std::vector<Flutcurve*>& curves);

// return the order of the largest matrix
size_t get_order (size_t* basic=nullptr);

#endif // PROCESS_H
