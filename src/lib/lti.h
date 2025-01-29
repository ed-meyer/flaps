//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef LTI_H
#define LTI_H

#include "Pz.h"

// Internal time delays
class Itd {
public:
	std::vector<Elem> elem;
	double deltat;

	Itd() : deltat(0.0) {}
	Itd(std::vector<Elem> e, double dt) : elem(e), deltat(dt) {}
};

namespace LTI {

int
eval (pset& plt, int nr, int nc, std::vector<std::complex<Ad>>& result,
	Matrix& A,
   Matrix& B,
   Matrix& C,
   Matrix& D,
   std::vector<Itd>& Atd,
   std::vector<double>& itd,
   std::vector<double>& otd,
   std::vector<int>& Sidx,
	// user-specified:
	const std::string& psiname,
	const std::string& stifname,
	std::vector<std::pair<int,double>>& ecolk,
   std::vector<std::string>& igains,
   std::vector<std::string>& iphases,
   std::vector<std::string>& ogains,
   std::vector<std::string>& ophases);
}	// namespace LTI

#endif // LTI_H
