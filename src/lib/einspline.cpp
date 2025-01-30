//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <vector>
#include "blas.h"
#include "einspline.h"
#include "trace.h"

using namespace std;

static int
general_grid_reverse_map (const vector<double>& grid, double& x);
static int
general_grid_reverse_map (const vector<double>& grid, Ad& x);
static void
get_NUBasis_funcs_di (NUBasis&  basis, int i, double bfuncs[4]);
static void
get_NUBasis_d2funcs_di (NUBasis& basis, int i,
			double bfuncs[4], double dbfuncs[4], double d2bfuncs[4]);

NUBasis::
NUBasis (const vector<double>& grid_pts) {
// NUBasis constructor
	T_(Trace trc(2,"NUBasis constructor");)
	T_(trc.dprint(grid_pts.size()," grid points");)
  this->grid = grid_pts;
  int N = grid_pts.size();
 
  this->xVals = vector<double>(N+5, 0.0);

  this->dxInv = vector<double>(3*(N+2), 0.0);
  for (int i=0; i<N; i++)
    this->xVals[i+2] = grid_pts[i];
  const double* g = &grid_pts[0];
  // Extend grid points on either end to provide enough points to
  // construct a full basis set
 this->xVals[0]   = g[ 0 ] - 2.0*(g[1]-g[0]);
 this->xVals[1]   = g[ 0 ] - 1.0*(g[1]-g[0]);
 this->xVals[N+2] = g[N-1] + 1.0*(g[N-1]-g[N-2]);
 this->xVals[N+3] = g[N-1] + 2.0*(g[N-1]-g[N-2]);
 this->xVals[N+4] = g[N-1] + 3.0*(g[N-1]-g[N-2]);
  for (int i=0; i<N+2; i++) 
    for (int j=0; j<3; j++) 
      this->dxInv[3*i+j] = 
	1.0/(this->xVals[i+j+1]-this->xVals[i]);
}

void
solve_NUB_deriv_interp_1d_d (NUBasis& basis, const double* data, int datastride,
	  double* p, int pstride, double abcdInitial[4], double abcdFinal[4]) {
// compute spline coefficients
  int M = basis.grid.size();
  int N = M+2;
  // Banded matrix storage.  The first three elements in the
  // tinyvector store the tridiagonal coefficients.  The last element
  // stores the RHS data.
  vector<double> bands(4*N, 0.0);

  // Fill up bands
  for (int i=0; i<4; i++) {
    bands[i]         = abcdInitial[i];
    bands[4*(N-1)+i] = abcdFinal[i];
  }
  for (int i=0; i<M; i++) {
    get_NUBasis_funcs_di (basis, i, &(bands[4*(i+1)]));
    bands[4*(i+1)+3] = data[datastride*i];
  }
    
  // Now solve:
  // First and last rows are different
  bands[4*0+1] /= bands[4*0+0];
  bands[4*0+2] /= bands[4*0+0];
  bands[4*0+3] /= bands[4*0+0];
  bands[4*0+0] = 1.0;
  bands[4*1+1] -= bands[4*1+0]*bands[4*0+1];
  bands[4*1+2] -= bands[4*1+0]*bands[4*0+2];
  bands[4*1+3] -= bands[4*1+0]*bands[4*0+3];
  bands[4*0+0] = 0.0;
  bands[4*1+2] /= bands[4*1+1];
  bands[4*1+3] /= bands[4*1+1];
  bands[4*1+1] = 1.0;
  
  // Now do rows 2 through M+1
  for (int row=2; row < N-1; row++) {
    bands[4*(row)+1] -= bands[4*(row)+0]*bands[4*(row-1)+2];
    bands[4*(row)+3] -= bands[4*(row)+0]*bands[4*(row-1)+3];
    bands[4*(row)+2] /= bands[4*(row)+1];
    bands[4*(row)+3] /= bands[4*(row)+1];
    bands[4*(row)+0] = 0.0;
    bands[4*(row)+1] = 1.0;
  }

  // Do last row
  bands[4*(M+1)+1] -= bands[4*(M+1)+0]*bands[4*(M-1)+2];
  bands[4*(M+1)+3] -= bands[4*(M+1)+0]*bands[4*(M-1)+3];
  bands[4*(M+1)+2] -= bands[4*(M+1)+1]*bands[4*(M)+2];
  bands[4*(M+1)+3] -= bands[4*(M+1)+1]*bands[4*(M)+3];
  bands[4*(M+1)+3] /= bands[4*(M+1)+2];
  bands[4*(M+1)+2] = 1.0;

  p[pstride*(M+1)] = bands[4*(M+1)+3];
  // Now back substitute up
  for (int row=M; row>0; row--)
    p[pstride*(row)] = bands[4*(row)+3] - bands[4*(row)+2]*p[pstride*(row+1)];
  
  // Finish with first row
  p[0] = bands[4*(0)+3] - bands[4*(0)+1]*p[pstride*1] - bands[4*(0)+2]*p[pstride*2];
}

static void
find_NUBcoefs_1d_d (NUBasis& basis,
		    const double* data,  int dstride,
		    double *coefs, int cstride) {
	int M = basis.grid.size();
	// Setup boundary conditions
	double bfuncs[4], dbfuncs[4], abcd_left[4], abcd_right[4];
	// Left boundary
	get_NUBasis_d2funcs_di (basis, 0, bfuncs, dbfuncs, abcd_left);
	// natural b-spline:
	abcd_left[3] = 0.0;

	// Right boundary
	get_NUBasis_d2funcs_di (basis, M-1, bfuncs, dbfuncs, abcd_right);
	// natural b-spline:
	abcd_right[3] = 0.0;

	// Now, solve for coefficients
	solve_NUB_deriv_interp_1d_d (basis, data, dstride, coefs, cstride,
			 abcd_left, abcd_right);
}

NUBspline_2d::
NUBspline_2d (const vector<double>& data,NUBasis& x_basis, NUBasis& y_basis) {
// NUBspline_2d constructor

  // create the bases
  // this->x_basis = new NUBasis (x_grid);
  // this->y_basis = new NUBasis (y_grid);

  // int My = y_grid.size();
  int ny = y_basis.grid.size();
  int My = ny;
  // int Nx = x_grid.size() + 2;
  int nx = x_basis.grid.size();
  int Nx = nx + 2;
  // int Ny = y_grid.size() + 2;
  int Ny = ny + 2;
  
  this->x_stride = Ny;

  this->coefs = vector<double>(Nx*Ny, 0.0);

  // First, solve in the X-direction 
  for (int iy=0; iy<My; iy++) {
    // int doffset = iy;
    // int coffset = iy;
    find_NUBcoefs_1d_d (x_basis, &data[iy], My,
			&this->coefs[iy], Ny);
  }
  
  // Now, solve in the Y-direction
  for (int ix=0; ix<Nx; ix++) {
    // int doffset = ix*Ny;
    // int coffset = ix*Ny;
    find_NUBcoefs_1d_d (y_basis, &this->coefs[ix*Ny], 1, 
			&this->coefs[ix*Ny], 1);
  }
}


NUBspline_3d::
NUBspline_3d (const vector<double>& data,
		NUBasis& x_basis, NUBasis& y_basis, NUBasis& z_basis) {
// NUBspline_3d constructor: modification from the original - passing
// the basis' to avoid computing them more than once

  int nx = x_basis.grid.size();
  int ny = y_basis.grid.size();
  int My = ny;
  int nz = z_basis.grid.size();
  int Mz = nz;
  int Nx = nx + 2;
  int Ny = ny + 2;
  int Nz = nz + 2;
  
  this->x_stride = Ny*Nz;
  this->y_stride = Nz;

  this->coefs = vector<double>(Nx*Ny*Nz, 0.0);

	// First, solve in the X-direction 
	for (int iy=0; iy<My; iy++) {
		for (int iz=0; iz<Mz; iz++) {
			int doffset = iy*Mz+iz;
			int coffset = iy*Nz+iz;
			find_NUBcoefs_1d_d (x_basis, &data[doffset], My*Mz,
				&this->coefs[coffset], Ny*Nz);
		}
	}
  
  // Now, solve in the Y-direction
	for (int ix=0; ix<Nx; ix++) {
		for (int iz=0; iz<Nz; iz++) {
			int doffset = ix*Ny*Nz + iz;
			int coffset = ix*Ny*Nz + iz;
			find_NUBcoefs_1d_d (y_basis, &this->coefs[doffset], Nz, 
				&this->coefs[coffset], Nz);
		}
	}
  
	// Now, solve in the z-direction
	for (int ix=0; ix<Nx; ix++) {
		for (int iy=0; iy<Ny; iy++) {
			int doffset = (ix*Ny + iy)*Nz;
			int coffset = (ix*Ny + iy)*Nz;
			find_NUBcoefs_1d_d (z_basis, &this->coefs[doffset], 1, 
				&this->coefs[coffset], 1);
		}
  }
}

int
get_NUBasis_funcs_d (NUBasis& basis, double& x, double bfuncs[4]) {
  double b1[2], b2[3];
  // int i = (*basis->grid->reverse_map)(basis->grid, x);
  int i = general_grid_reverse_map(basis.grid, x);
  int i2 = i+2;
  double*  dxInv = &basis.dxInv[0];
  double*  xVals = &basis.xVals[0];

  b1[0]     = (xVals[i2+1]-x)  * dxInv[3*(i+2)+0];
  b1[1]     = (x-xVals[i2])    * dxInv[3*(i+2)+0];
  b2[0]     = (xVals[i2+1]-x)  * dxInv[3*(i+1)+1] * b1[0];
  b2[1]     = ((x-xVals[i2-1]) * dxInv[3*(i+1)+1] * b1[0]+
	       (xVals[i2+2]-x) * dxInv[3*(i+2)+1] * b1[1]);
  b2[2]     = (x-xVals[i2])    * dxInv[3*(i+2)+1] * b1[1];
  bfuncs[0] = (xVals[i2+1]-x)  * dxInv[3*(i  )+2] * b2[0];
  bfuncs[1] = ((x-xVals[i2-2]) * dxInv[3*(i  )+2] * b2[0] +
	       (xVals[i2+2]-x) * dxInv[3*(i+1)+2] * b2[1]);
  bfuncs[2] = ((x-xVals[i2-1]) * dxInv[3*(i+1)+2] * b2[1] +
	       (xVals[i2+3]-x) * dxInv[3*(i+2)+2] * b2[2]);
  bfuncs[3] = (x-xVals[i2])    * dxInv[3*(i+2)+2] * b2[2];
  return i;
}

int
get_NUBasis_funcs_d (NUBasis& basis, Ad& x, Ad bfuncs[4]) {
  Ad b1[2], b2[3];
  // int i = (*basis->grid->reverse_map)(basis->grid, x);
  int i = general_grid_reverse_map(basis.grid, x);
  int i2 = i+2;
  double*  dxInv = &basis.dxInv[0];
  double*  xVals = &basis.xVals[0];

	int nxval = basis.xVals.size();
	int ndxinv = basis.dxInv.size();
 	if (i2+3 >= nxval)
		throw runtime_error(vastr("i2+3 ",i2+3," > xVals size ",nxval));
 	if (3*(i+2)+2 >= ndxinv)
		throw runtime_error(vastr("3*(i+2)+2 ",3*(i+2)+2,"  > dxInv size ",ndxinv));
  b1[0]     = (xVals[i2+1]-x)  * dxInv[3*(i+2)+0];
  b1[1]     = (x-xVals[i2])    * dxInv[3*(i+2)+0];
  b2[0]     = (xVals[i2+1]-x)  * dxInv[3*(i+1)+1] * b1[0];
  b2[1]     = ((x-xVals[i2-1]) * dxInv[3*(i+1)+1] * b1[0]+
	       (xVals[i2+2]-x) * dxInv[3*(i+2)+1] * b1[1]);
  b2[2]     = (x-xVals[i2])    * dxInv[3*(i+2)+1] * b1[1];
  bfuncs[0] = (xVals[i2+1]-x)  * dxInv[3*(i  )+2] * b2[0];
  bfuncs[1] = ((x-xVals[i2-2]) * dxInv[3*(i  )+2] * b2[0] +
	       (xVals[i2+2]-x) * dxInv[3*(i+1)+2] * b2[1]);
  bfuncs[2] = ((x-xVals[i2-1]) * dxInv[3*(i+1)+2] * b2[1] +
	       (xVals[i2+3]-x) * dxInv[3*(i+2)+2] * b2[2]);
  bfuncs[3] = (x-xVals[i2])    * dxInv[3*(i+2)+2] * b2[2];
  return i;
}

void
get_NUBasis_funcs_di (NUBasis&  basis, int i, double bfuncs[4]) {
  int i2 = i+2;
  double b1[2], b2[3];
  double x = basis.grid[i];
  double*  dxInv = &basis.dxInv[0];
  double*  xVals = &basis.xVals[0]; 

  b1[0]     = (xVals[i2+1]-x)  * dxInv[3*(i+2)+0];
  b1[1]     = (x-xVals[i2])    * dxInv[3*(i+2)+0];
  b2[0]     = (xVals[i2+1]-x)  * dxInv[3*(i+1)+1] * b1[0];
  b2[1]     = ((x-xVals[i2-1]) * dxInv[3*(i+1)+1] * b1[0]+
	       (xVals[i2+2]-x) * dxInv[3*(i+2)+1] * b1[1]);
  b2[2]     = (x-xVals[i2])    * dxInv[3*(i+2)+1] * b1[1];
  bfuncs[0] = (xVals[i2+1]-x)  * dxInv[3*(i  )+2] * b2[0];
  bfuncs[1] = ((x-xVals[i2-2]) * dxInv[3*(i  )+2] * b2[0] +
	       (xVals[i2+2]-x) * dxInv[3*(i+1)+2] * b2[1]);
  bfuncs[2] = ((x-xVals[i2-1]) * dxInv[3*(i+1)+2] * b2[1] +
	       (xVals[i2+3]-x) * dxInv[3*(i+2)+2] * b2[2]);
  bfuncs[3] = (x-xVals[i2])    * dxInv[3*(i+2)+2] * b2[2];
}


void
get_NUBasis_d2funcs_di (NUBasis& basis, int i,
			double bfuncs[4], double dbfuncs[4], double d2bfuncs[4]) {
  double b1[2], b2[3];
  double x = basis.grid[i];
  int i2 = i+2;
  double*  dxInv = &basis.dxInv[0];
  double*  xVals = &basis.xVals[0];

  b1[0]     = (xVals[i2+1]-x)  * dxInv[3*(i+2)+0];
  b1[1]     = (x-xVals[i2])    * dxInv[3*(i+2)+0];
  b2[0]     = (xVals[i2+1]-x)  * dxInv[3*(i+1)+1] * b1[0];
  b2[1]     = ((x-xVals[i2-1]) * dxInv[3*(i+1)+1] * b1[0]+
	       (xVals[i2+2]-x) * dxInv[3*(i+2)+1] * b1[1]);
  b2[2]     = (x-xVals[i2])    * dxInv[3*(i+2)+1] * b1[1];
  bfuncs[0] = (xVals[i2+1]-x)  * dxInv[3*(i  )+2] * b2[0];
  bfuncs[1] = ((x-xVals[i2-2]) * dxInv[3*(i  )+2] * b2[0] +
	       (xVals[i2+2]-x) * dxInv[3*(i+1)+2] * b2[1]);
  bfuncs[2] = ((x-xVals[i2-1]) * dxInv[3*(i+1)+2] * b2[1] +
	       (xVals[i2+3]-x) * dxInv[3*(i+2)+2] * b2[2]);
  bfuncs[3] = (x-xVals[i2])    * dxInv[3*(i+2)+2] * b2[2]; 

  dbfuncs[0] = -3.0 * (dxInv[3*(i  )+2] * b2[0]);
  dbfuncs[1] =  3.0 * (dxInv[3*(i  )+2] * b2[0] - dxInv[3*(i+1)+2] * b2[1]);
  dbfuncs[2] =  3.0 * (dxInv[3*(i+1)+2] * b2[1] - dxInv[3*(i+2)+2] * b2[2]);
  dbfuncs[3] =  3.0 * (dxInv[3*(i+2)+2] * b2[2]);

  d2bfuncs[0] = 6.0 * (+dxInv[3*(i+0)+2]* dxInv[3*(i+1)+1]*b1[0]);
  d2bfuncs[1] = 6.0 * (-dxInv[3*(i+1)+1]*(dxInv[3*(i+0)+2]+dxInv[3*(i+1)+2])*b1[0] +
		        dxInv[3*(i+1)+2]* dxInv[3*(i+2)+1]*b1[1]);
  d2bfuncs[2] = 6.0 * (+dxInv[3*(i+1)+2]* dxInv[3*(i+1)+1]*b1[0] -
		        dxInv[3*(i+2)+1]*(dxInv[3*(i+1)+2] + dxInv[3*(i+2)+2])*b1[1]);
  d2bfuncs[3] = 6.0 * (+dxInv[3*(i+2)+2]* dxInv[3*(i+2)+1]*b1[1]);
}


Ad
NUBspline_2d::
eval(const Ad a[4], const Ad b[4], Ad cb[4], int ix, int iy) {
  // Ad a[4], b[4];
  // int ix = get_NUBasis_funcs_d (this->x_basis, x, a);
  // int iy = get_NUBasis_funcs_d (this->y_basis, y, b);
  
  int xs = x_stride;
	// lambda for checking for range errors
	auto C = [&] (int i, int j) {
		int idx = (ix+i)*xs + iy+j;
	  // rf-rsf interp in test case gaf.fp has ix == 4 causing this error
	  // but there seems to be no harm because it only triggers with C(3,*)
	  // but a[3] is zero in this case so returning 0.0 is ok and prevents
	  // an out-of-range access of coeffs
		if (idx >= (int)coefs.size()) {
			//!! string exc{vastr("(",i,',',j,") ix ",ix,", iy ",iy,
				//!! 	", xs ",xs," > coef.size ",coefs.size())};
			//!! throw runtime_error(exc);
			//!! cerr << exc << endl;
			return 0.0;
		}
	  return coefs[idx];
	 };

//	Ad rval =
//			a[0]*(C(0,0)*b[0] + C(0,1)*b[1] + C(0,2)*b[2] + C(0,3)*b[3]) +
//			a[1]*(C(1,0)*b[0] + C(1,1)*b[1] + C(1,2)*b[2] + C(1,3)*b[3])+
//			a[2]*(C(2,0)*b[0] + C(2,1)*b[1] + C(2,2)*b[2] + C(2,3)*b[3])+
//			a[3]*(C(3,0)*b[0] + C(3,1)*b[1] + C(3,2)*b[2] + C(3,3)*b[3]);
	 // the return value is a triple product
	 //   a^t C b
	vector<double> c(16);
	for (int i=0; i<4; i++)
		for (int j=0; j<4; j++)
			c[i+4*j] = C(i,j);
	//!! vector<Ad> cb(4);
	blas::gemv("n", 4, 4, 1.0, &c[0], 4, &b[0], 1, 0.0, &cb[0], 1);
	Ad rval;
	blas::dot(4, &a[0], 1, &cb[0], 1, rval);
  return rval;
}

double
NUBspline_2d::
eval(double a[4], double b[4], int ix, int iy) {
  // double a[4], b[4];
  // int ix = get_NUBasis_funcs_d (x_basis, x, a);
  // int iy = get_NUBasis_funcs_d (y_basis, y, b);
  
  int xs = x_stride;
  // lambda for checking for range errors
  auto C = [&] (int i, int j) {
		int idx = (ix+i)*xs + iy+j;
		if (idx >= (int)coefs.size()) {
			string exc{vastr("(",i,',',j,") ix ",ix,", iy ",iy,
					", xs ",xs," > coef.size ",coefs.size())};
			//!! throw runtime_error(exc);
			cerr << exc << endl;
			return 0.0;
		}
	  return coefs[idx];
	 };
  double rval = (a[0]*(C(0,0)*b[0]+C(0,1)*b[1]+C(0,2)*b[2]+C(0,3)*b[3])+
	  a[1]*(C(1,0)*b[0]+C(1,1)*b[1]+C(1,2)*b[2]+C(1,3)*b[3])+
	  a[2]*(C(2,0)*b[0]+C(2,1)*b[1]+C(2,2)*b[2]+C(2,3)*b[3])+
	  a[3]*(C(3,0)*b[0]+C(3,1)*b[1]+C(3,2)*b[2]+C(3,3)*b[3]));
	return rval;
}

// inline void
// eval_NUBspline_3d (NUBspline_3d*  spline,
// 		double x, double y, double z, double&  val) {
double
NUBspline_3d::
eval(double a[4], double b[4], double c[4], int ix, int iy, int iz) {

  // double a[4], b[4], c[4];
  // int ix = get_NUBasis_funcs_d (spline->x_basis, x, a);
  // int iy = get_NUBasis_funcs_d (spline->y_basis, y, b);
  // int iz = get_NUBasis_funcs_d (spline->z_basis, z, c);
  double*  coefs = &this->coefs[0];
  
  int xs = x_stride;
  int ys = y_stride;
#define P(i,j,k) coefs[(ix+(i))*xs+(iy+(j))*ys+(iz+(k))]
  double rval =
	  (a[0]*(b[0]*(P(0,0,0)*c[0]+P(0,0,1)*c[1]+P(0,0,2)*c[2]+P(0,0,3)*c[3])+
		b[1]*(P(0,1,0)*c[0]+P(0,1,1)*c[1]+P(0,1,2)*c[2]+P(0,1,3)*c[3])+
		b[2]*(P(0,2,0)*c[0]+P(0,2,1)*c[1]+P(0,2,2)*c[2]+P(0,2,3)*c[3])+
		b[3]*(P(0,3,0)*c[0]+P(0,3,1)*c[1]+P(0,3,2)*c[2]+P(0,3,3)*c[3]))+
	  a[1]*(b[0]*(P(1,0,0)*c[0]+P(1,0,1)*c[1]+P(1,0,2)*c[2]+P(1,0,3)*c[3])+
		b[1]*(P(1,1,0)*c[0]+P(1,1,1)*c[1]+P(1,1,2)*c[2]+P(1,1,3)*c[3])+
		b[2]*(P(1,2,0)*c[0]+P(1,2,1)*c[1]+P(1,2,2)*c[2]+P(1,2,3)*c[3])+
		b[3]*(P(1,3,0)*c[0]+P(1,3,1)*c[1]+P(1,3,2)*c[2]+P(1,3,3)*c[3]))+
	  a[2]*(b[0]*(P(2,0,0)*c[0]+P(2,0,1)*c[1]+P(2,0,2)*c[2]+P(2,0,3)*c[3])+
		b[1]*(P(2,1,0)*c[0]+P(2,1,1)*c[1]+P(2,1,2)*c[2]+P(2,1,3)*c[3])+
		b[2]*(P(2,2,0)*c[0]+P(2,2,1)*c[1]+P(2,2,2)*c[2]+P(2,2,3)*c[3])+
		b[3]*(P(2,3,0)*c[0]+P(2,3,1)*c[1]+P(2,3,2)*c[2]+P(2,3,3)*c[3]))+
	  a[3]*(b[0]*(P(3,0,0)*c[0]+P(3,0,1)*c[1]+P(3,0,2)*c[2]+P(3,0,3)*c[3])+
		b[1]*(P(3,1,0)*c[0]+P(3,1,1)*c[1]+P(3,1,2)*c[2]+P(3,1,3)*c[3])+
		b[2]*(P(3,2,0)*c[0]+P(3,2,1)*c[1]+P(3,2,2)*c[2]+P(3,2,3)*c[3])+
		b[3]*(P(3,3,0)*c[0]+P(3,3,1)*c[1]+P(3,3,2)*c[2]+P(3,3,3)*c[3])));
#undef P
	return rval;
}

Ad
NUBspline_3d::
eval(Ad a[4], Ad b[4], Ad c[4], int ix, int iy, int iz) {

  // double a[4], b[4], c[4];
  // int ix = get_NUBasis_funcs_d (spline->x_basis, x, a);
  // int iy = get_NUBasis_funcs_d (spline->y_basis, y, b);
  // int iz = get_NUBasis_funcs_d (spline->z_basis, z, c);
  // double*  coefs = &this->coefs[0];
  
  int xs = x_stride;
  int ys = y_stride;
// #define P(i,j,k) coefs[(ix+(i))*xs+(iy+(j))*ys+(iz+(k))]
  // lambda for checking for range errors
  auto P = [&] (int i, int j, int k) {
		int idx = (ix+i)*xs + (iy+j)*ys + iz+k;
		if (idx >= (int)this->coefs.size()) {
			string exc{vastr("(",i,',',j,',',k,") ix ",ix,", iy ",iy,", iz ",iz,
					", xs ",xs,", ys ",ys," > coef.size ",coefs.size())};
			//!! throw runtime_error(vastr("(",i,',',j,") ix ",ix,", iy ",iy,", stride ",stride," > coef.size ",coefs.size()));
			cerr << exc << endl;
			return 0.0;
		}
	  return coefs[idx];
	 };

//  Ad rval =
//	  (a[0]*(
//			b[0]*(P(0,0,0)*c[0]+P(0,0,1)*c[1]+P(0,0,2)*c[2]+P(0,0,3)*c[3])+
//			b[1]*(P(0,1,0)*c[0]+P(0,1,1)*c[1]+P(0,1,2)*c[2]+P(0,1,3)*c[3])+
//			b[2]*(P(0,2,0)*c[0]+P(0,2,1)*c[1]+P(0,2,2)*c[2]+P(0,2,3)*c[3])+
//			b[3]*(P(0,3,0)*c[0]+P(0,3,1)*c[1]+P(0,3,2)*c[2]+P(0,3,3)*c[3]))+
//	  a[1]*(
//			b[0]*(P(1,0,0)*c[0]+P(1,0,1)*c[1]+P(1,0,2)*c[2]+P(1,0,3)*c[3])+
//			b[1]*(P(1,1,0)*c[0]+P(1,1,1)*c[1]+P(1,1,2)*c[2]+P(1,1,3)*c[3])+
//			b[2]*(P(1,2,0)*c[0]+P(1,2,1)*c[1]+P(1,2,2)*c[2]+P(1,2,3)*c[3])+
//			b[3]*(P(1,3,0)*c[0]+P(1,3,1)*c[1]+P(1,3,2)*c[2]+P(1,3,3)*c[3]))+
//	  a[2]*(
//			b[0]*(P(2,0,0)*c[0]+P(2,0,1)*c[1]+P(2,0,2)*c[2]+P(2,0,3)*c[3])+
//			b[1]*(P(2,1,0)*c[0]+P(2,1,1)*c[1]+P(2,1,2)*c[2]+P(2,1,3)*c[3])+
//			b[2]*(P(2,2,0)*c[0]+P(2,2,1)*c[1]+P(2,2,2)*c[2]+P(2,2,3)*c[3])+
//			b[3]*(P(2,3,0)*c[0]+P(2,3,1)*c[1]+P(2,3,2)*c[2]+P(2,3,3)*c[3]))+
//	  a[3]*(
//			b[0]*(P(3,0,0)*c[0]+P(3,0,1)*c[1]+P(3,0,2)*c[2]+P(3,0,3)*c[3])+
//			b[1]*(P(3,1,0)*c[0]+P(3,1,1)*c[1]+P(3,1,2)*c[2]+P(3,1,3)*c[3])+
//			b[2]*(P(3,2,0)*c[0]+P(3,2,1)*c[1]+P(3,2,2)*c[2]+P(3,2,3)*c[3])+
//			b[3]*(P(3,3,0)*c[0]+P(3,3,1)*c[1]+P(3,3,2)*c[2]+P(3,3,3)*c[3])));

	// pc[i] = b[j]*P[i,j,k]*c[k]
	// rval = a^t * pc
	vector<Ad> pc(4);
	Ad pcb;
	Ad rval;
	for (int i=0; i<4; i++) {
		vector<double> Pi(16);
		for (int j=0; j<4; j++) {
			for (int k=0; k<4; k++) {
				Pi[j+4*k] = P(i,j,k);
			}
		}
		blas::gemv("n", 4, 4, 1.0, &Pi[0], 4, &c[0], 1, 0.0, &pc[0], 1);
		blas::dot(4, &pc[0], 1, &b[0], 1, pcb);
		rval += a[i]*pcb;
	}

	//!! XXX test the 2 methods
//!!	if (!is_equal(rval, rval2, 8))
//!!		throw runtime_error(vastr("2 method do not agree: ",rval," and ",rval2));

// #undef P
  return rval;
}

int
general_grid_reverse_map (const vector<double>& grid, double& x) {
  // int N = grid->points.size();
  int N = grid.size();
	if (x <= grid[0]) {
		x = grid[0];
		return (0);
	} else if (x >= grid[N-1]) {
		x = grid[N-1];
		return (N-1);
  } else {
    int hi = N-1;
    int lo = 0;
    bool done = false;
    while (!done) {
      int i = (hi+lo)>>1;
      if (grid[i] > x)
	hi = i;
      else
	lo = i;
      done = (hi-lo)<2;
    }
    return (lo);
  }
}

int
general_grid_reverse_map (const vector<double>& grid, Ad& x) {
  // int N = grid->points.size();
  int N = grid.size();
	if (x.value() <= grid[0]) {
		x.value(grid[0]);  // retain derivatives
		return (0);
	} else if (x.value() >= grid[N-1]) {
		x.value(grid[N-1]);	// retain derivatives
		return (N-1);
  } else {
    int hi = N-1;
    int lo = 0;
    bool done = false;
    while (!done) {
      int i = (hi+lo)>>1;
      if (grid[i] > x.value())
         hi = i;
      else
         lo = i;
      done = (hi-lo)<2;
    }
    return (lo);
  }
}

#ifdef MAIN
#undef MAIN

#include <iostream>

int
main(int argc, char **argv) {
	vector<double> xs{0.0, 1.0};
	vector<double> ys{0.0, 1.0, 2.0};
	vector<double> data;
	//  0  1  2
	//  3  4  5
	double t{0.0};
	for (size_t i=0; i<xs.size(); i++) {
		for (size_t j=0; j<ys.size(); j++) {
			data.push_back(t++);
		}
	}
	NUBspline_2d spline(xs, ys, data);

	double x{0.5};
	double y{1.5};
	double z = spline.eval(x,y);
	cout << "at (" << x << ',' << y << ") = " << z << endl;
}
#endif // MAIN
