//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef REGION_H
#define REGION_H

#include <string>
#include <vector>

#include "lexer.h"
#include "pset.h"

// Classes for describing regions, including
// - intervals of parameter values, e.g. freq[0:10]
// - fixed values of parameters, e.g. vtas=(100,200,300)
// - curve ids, e.g. mode_1
// Class Region is just a vector of Bndry with supporting functions.
//
// A Bndry can be an interval of a parameter (parname[l:h]),
// a point value for a parameter (parname=value), or just
// a curveid, e.g. mode_1, mode_2, etc.
class Bndry {
public:
	std::string parname;
	double min;   // internal units
	double max;   // internal units
	std::string curveid;
	// only take ordinal places where the fixed parameter has the
	// desired value
	std::vector<int> ordinal;
	std::vector<int> modes;
	bool closest;

	// constructors
	Bndry (Tok const& opt);

	std::string summary(pset&) const;
	// isPoint() return true if this is a point Bndry
	bool isPoint() const;
};

class Region : public std::vector<Bndry> {
public:
	// constructors:
	Region() {}
	Region(std::vector<Bndry>& reg);
	Region(std::string const& options);

	std::string summary(pset&) const;
	const Par* find(pset& plt, std::string const& name) const;
	Bndry* get(std::string const& name);
	// is there a parameter in "pl" that is out of the region Bndrys?
	const Par* outOfRange (pset& pl, std::string& oormsg) const;
	Bndry* pointBndry();  // is there a point Bndry?
	// curveid Bndrys?
	std::vector<std::string> curveidBndry();
	bool hasCurveid(std::string& id);
	bool curveidOk(std::string& id);
	// ordinal Bndrys?
	std::vector<int> ordinal_numbers() const;
	std::vector<int> mode_numbers() const;
	bool contains_ordinal(int ordinal) const;
	bool contains_mode(int mode) const;
	bool contains (std::string const& name, double val) const;
};

std::ostream&
operator<<(std::ostream& s, const Bndry& t);

std::ostream&
operator<<(std::ostream& s, const Region& t);

#endif // REGION_H
