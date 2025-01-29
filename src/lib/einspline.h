//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef EINSPLINE_H
#define EINSPLINE_H

#include "Ad.h"

// ------------  #include "nubasis.h"

class NUBasis {
public:
	// constructor
	NUBasis() {}
	NUBasis (const std::vector<double>& grid_pts);

	std::vector<double> grid;
  // xVals is just the grid points, augmented by two extra points on
  // either side.  These are necessary to generate enough basis functions. 
  std::vector<double> xVals;

  // dxInv[3*i+j] = 1.0/(grid(i+j-1)-grid(i-2))
  std::vector<double> dxInv;
};

// XXX make this a member. Note: x may get change if out-of-range so
// it's a non-const reference
int
get_NUBasis_funcs_d (NUBasis& basis, double& x, double bfuncs[4]);
int
get_NUBasis_funcs_d (NUBasis& basis, Ad& x, Ad bfuncs[4]);


// -----------------  end of  NUBASIS_H

// -----------   #include "bspline_base.h"

using NUBspline_1d_d = struct {
  std::vector<double>  coefs;
  NUBasis x_basis;
};

class NUBspline_2d {
public:
  std::vector<double>  coefs;
  int x_stride{0};
  // NUBasis x_basis;
  // NUBasis y_basis;

	// constructor
	NUBspline_2d() {}
	// NUBspline_2d(const std::vector<double>& x_grid,
	// 		const std::vector<double>& y_grid, const std::vector<double>& data);
	NUBspline_2d(const std::vector<double>& data, NUBasis& x_basis, NUBasis& y_basis);
	// evaluators
	Ad eval(const Ad a[4], const Ad b[4], Ad cb[4], int ix, int iy);
	double eval(double a[4], double b[4], int ix, int iy);

};

class NUBspline_3d {
public:
  std::vector<double>  coefs;
  int x_stride{0};
  int y_stride{0};
  NUBasis x_basis;
  NUBasis y_basis;
  NUBasis z_basis;
	// constructors
	NUBspline_3d() {}
	NUBspline_3d(const std::vector<double>& data,
			NUBasis& x_basis, NUBasis& y_basis, NUBasis& z_basis);

	// evaluators
	double eval(double a[4], double b[4], double c[4], int ix, int iy, int iz);
	Ad eval(Ad a[4], Ad b[4], int ix, int iy, int iz);
	Ad eval(Ad a[4], Ad b[4], Ad c[4],
			int ix, int iy, int iz);
};


// evaluation
void
eval_NUBspline_2d (NUBspline_2d* spline, 
		    double x, double y, double& val);

void
eval_NUBspline_3d (NUBspline_3d* spline, 
		    double x, double y, double z, double& val);

#endif // EINSPLINE_H
