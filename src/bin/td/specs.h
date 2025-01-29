//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef specs_h
#define specs_h

// time-domain analysis of LCO
// References
//    seydel1988equilibrium
//    doedel2007lecture
#include <string>
#include <vector>

#include "matrix.h"
#include "Par.h"
#include "region.h"

class Specs {
public:
	std::vector<double> mass;		// (n,n) real matrix
	Matrix* stif{nullptr};
	std::vector<double> beta;
	std::vector<std::vector<double>> R;	// 2+nbeta (n,n) real matrices
	std::vector<double> P;
	int ireal{1};	// 0-based index of the ev comp to hold zero for normalization
	int maxstep;
};

// return a reference to *the* specs
Specs& specs();

// parse the user's options from cin
void
parser (std::string const& prog, Specs& specs);

void print_specs(std::string const& prog);

// is "desc" a matrix in the flutter equation (mass, stif, gaf, etc)?
bool is_flutmat (std::string const& desc);


#endif // specs_h
