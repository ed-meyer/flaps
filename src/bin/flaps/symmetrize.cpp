//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "config.h"
#include "functions.h"
#include "fio.h"
#include "lexer.h"
#include "matrix.h"
#include "settings.h"
#include "trace.h"

using namespace std;

static void dosym (Matrix& mat, vector<int>& dofs);


//--------------------------------------------------------------------
// 
//-------------------------------------------------------------------
bool
symmetrize (string const& options) {
	T_(Trace trc(2,"symmetrize");)
	vector<int> dofs;
	string mid;
	string newmid;

	// parse the options with lambdas
	vector<Tok*> unrec = flaps::lexer (options, {
		{"dof", [&](const Tok& p) { dofs = p.ivec; return true; }},
		{"i", [&](const Tok& p) { mid = p.srhs; return true; }},
		{"o", [&](const Tok& p) { newmid = p.srhs; return true; }},
	});

	// fetch constructor
	Matrix mat(mid);
	// symmetrize it...
	dosym(mat, dofs);
	// and store with the new name
	mat.mid(newmid);
	mat.store();

	return true;
}

static void
dosym (Matrix& mat, vector<int>& dof) {
	double* rp = mat.elem();
	complex<double>* cp = mat.celem();
	// dof must be in pairs
	assert(dof.size()%2 == 0);
	int nr = mat.rsize();
	int nc = mat.csize();
	assert(nr == nc);
	// first do the pairs of dof rows
	for (size_t ip = 0; ip<dof.size(); ip += 2) {
		int i1b = dof[ip];
		int j1b = dof[ip+1];
		int i0b = i1b - 1;
		int j0b = j1b - 1;
		for (int k=0; k<nr; k++) {
			if (k == i0b || k == j0b)
				continue;
			if (mat.is_complex()) {
				double cir = cp[IJ(i0b,k,nr)].real();
				double cii = cp[IJ(i0b,k,nr)].imag();
				double cjr = cp[IJ(j0b,k,nr)].real();
				double cji = cp[IJ(j0b,k,nr)].imag();
				double tr = (abs(cir)+abs(cjr))/2.0;
				double ti = (abs(cii)+abs(cji))/2.0;
				double dir = (cir < 0.0 ? -tr : tr);
				double dii = (cii < 0.0 ? -ti : ti);
				double djr = (cjr < 0.0 ? -tr : tr);
				double dji = (cji < 0.0 ? -ti : ti);
				cp[IJ(i0b,k,nr)] = complex<double>(dir,dii);
				cp[IJ(j0b,k,nr)] = complex<double>(djr,dji);
				// repeat for the columns
				cir = cp[IJ(k,i0b,nr)].real();
				cii = cp[IJ(k,i0b,nr)].imag();
				cjr = cp[IJ(k,j0b,nr)].real();
				cji = cp[IJ(k,j0b,nr)].imag();
				tr = (abs(cir)+abs(cjr))/2.0;
				ti = (abs(cii)+abs(cji))/2.0;
				dir = (cir < 0.0 ? -tr : tr);
				dii = (cii < 0.0 ? -ti : ti);
				djr = (cjr < 0.0 ? -tr : tr);
				dji = (cji < 0.0 ? -ti : ti);
				cp[IJ(k,i0b,nr)] = complex<double>(dir,dii);
				cp[IJ(k,j0b,nr)] = complex<double>(djr,dji);
			} else {
				double ti = abs(rp[IJ(i0b,k,nr)]);
				double tj = abs(rp[IJ(j0b,k,nr)]);
				double t = (ti+tj)/2.0;
				double ri = (rp[IJ(i0b,k,nr)] < 0.0 ? -t : t);
				double rj = (rp[IJ(j0b,k,nr)] < 0.0 ? -t : t);
				rp[IJ(i0b,k,nr)] = ri;
				rp[IJ(j0b,k,nr)] = rj;
				rp[IJ(k,i0b,nr)] = ri;
				rp[IJ(k,j0b,nr)] = rj;
			}
		}
	}
	// then symmetrize the whole matrix if real
	if (mat.is_complex())
		return;
	for (int i=0; i<nr; i++) {
		for (int j=i+1; j<nr; j++) {
			double t = (rp[IJ(i,j,nr)]+rp[IJ(j,i,nr)])/2.0;
			rp[IJ(i,j,nr)] = t;
			rp[IJ(j,i,nr)] = t;
		}
	}
}

