//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Some numerical algorithms taken from "Numerical Recipes in C":
// @article{teukolsky1992numerical,
//   title={Numerical recipes in C},
//   author={Teukolsky, Saul A and Flannery, Brian P and Press, WH and Vetterling, WT},
//   journal={SMR},
//   volume={693},
//   number={1},
//   pages={59--70},
//   year={1992}
// }
#ifndef NRC_H
#define NRC_H

#include <algorithm>
#include <functional>
#include <limits>
#include <string>
#include <vector>

using Fdfcn = std::function<void(int n, double x, double* f)>;
using Rootfcn = std::function<double(double)>;

namespace NRC {
double root(double x0, double x1, double tol, Rootfcn fcn);
bool fdapprox (int n, double x, double xmin, double xmax,
		double* y, Fdfcn fcn, double& err);
}	// namespace NRC

#endif  // NRC_H
