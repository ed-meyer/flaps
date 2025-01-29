//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef TARGET_H
#define TARGET_H
// Functions for treating targets in flut
// Public interface:
//   vector<Target*> targetlist(string& options)
//       parses "options", create targets from them

#include <string>
#include <iostream>

#include "pset.h"

class Window {
	std::string name;
	double min;		// in EU
	double max;		// in EU
	double convfactor{1.0};	// EU -> PU
public:
	Window() noexcept = default;
	Window(std::string nm, double mi, double ma, double cf) :
		name(nm), min(mi/cf), max(ma/cf), convfactor(cf) {}
	Window(const Window&) noexcept = default;		// copy constructor
	Window& operator=(const Window&) = default;	// assignment op
	Window& operator=(Window&&) = default;	// assignment op
	bool inrange(pset& p) {
		double v = p.parval(name).value();
		if (v > min && v < max)
			return true;
		return false;
	}
	friend std::ostream& operator<<(std::ostream& s, const Window& t);
};


class Target {
public:
	std::string parname;  // parameter name
	double value{0.0};        // in "equation" units (EU)
	double convfactor{1.0};   // mult value to get "presentation" units (PU)
	std::vector<Window> windows;
	std::string sortpar;		// sort parameter name
	bool closest{false};
	bool is_lowerlimit{false};	// is this just a parameter limit?
	bool is_upperlimit{false};	// is this just a parameter limit?

	// constructors:
	Target() noexcept = default;
	Target(std::string const& par, std::string const& val, bool cls=false);
	Target(std::string const& par, double val, bool cls=false,
			std::string const& limit="");
	Target(const Target&) noexcept = default;		// copy constructor
	Target& operator=(const Target&) = default;	// assignment op
	Target& operator=(Target&&) = default;	// assignment op
};

std::ostream&
operator<<(std::ostream& s, const Target& t);

std::ostream&
operator<<(std::ostream& s, const std::vector<Target>& t);

#endif // TARGET_H
