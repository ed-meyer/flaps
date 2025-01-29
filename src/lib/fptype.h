//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef FPTYPE_H
#define FPTYPE_H

#include <cmath>

// divide num/den accounting for tiny den
double safe_divide (double num, double den);

namespace flaps {

// random numbers between -1 and 1:
double xrand();

bool
is_valid(double const& a);

} // namespace flaps


double ulp_frac (double a, double b);
int ulp_double (double a, double b);
bool is_equal (double const& a, double const& b, int sigfig);
bool is_lessthan (double a, double b, int sigfig);
bool is_greaterthan (double a, double b, int sigfig);

int
fcmp(double x1, double x2, double epsilon);

double closest_10(double x);

#endif // FPTYPE_H
