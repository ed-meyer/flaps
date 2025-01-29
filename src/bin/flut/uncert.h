//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef UNCERT_H
#define UNCERT_H 1

#include <string>
#include <vector>

class Uncert;

#include "flutcurve.h"
#include "interval.h"

class Uncert {
	Flutcurve* unc;  // unc par are indep
	std::vector<std::string> uncpar; // names of the uncertain parameters
	std::vector<double> uncfac; // uncertainty factors: +/-fac*value
public:
	Flutcurve* lo;
	Flutcurve* hi;
	// constructor
	Uncert(const Flutcurve& curve, const std::vector<std::string>& uncnames,
		const std::vector<double>& factors, const std::vector<std::string>& constpar);
	// set hi & lo
	void eval(Point& pt);
	void update(std::vector<Interval>& xu, double coord);
	std::vector<Interval> compute_f(Point& pt);
};

#endif // UNCERT_H
