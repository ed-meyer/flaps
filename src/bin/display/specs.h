//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef SPECS_H
#define SPECS_H

class Specs {
public:
	double filter{0.0};
	std::string format;
	bool matview{false};
	std::string outfile;
	std::vector<std::tuple<std::string,double>> params;
	std::vector<std::string> derivs;
	//!! std::vector<std::string> matrices;
};

std::vector<std::pair<std::string,Specs>> parser();

Specs& specs();

std::ostream&
operator<<(std::ostream& s, const Specs& t);

#endif // SPECS_H
