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

#include <string>
#include <tuple>
#include <vector>
#include "matrix.h"
#include "Par.h"

class Specs {
public:
	bool units{false};	// default: already in SI
	std::vector<std::tuple<std::string,std::string>> input;
	std::string outputname;
	std::vector<double> betas;
	Matrix* output_matrix{nullptr};
	bool forcer0{false};
	bool kscale{false};
	bool optbeta{false};
	bool extrap{false};
	double smooth{0.0};
	std::vector<Matrix*> matrices;
};

Specs& specs();
Matrix* parser();

#endif // specs_h
