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

#include <assert.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "config.h"
#include "Ad.h"
#include "blas.h"
#include "conv.h"
#include "exim.h"
#include "fem.h"
#include "lapack.h"
#include "matrix.h"
#include "Pz.h"
#include "settings.h"
#include "trace.h"

using namespace std;

vector<double>
Fem::
freevibration() {
// Do a free-vibration solution with the output nodal mass,
// stiffness, and possibly gyro; return the (nr,nmodes) modes matrix
	Trace trc(1,"freevibration");

	int nr = retained().size();

	// create a free-vibration gct based on the
	// statically-merged mass, stiffness, and possibly gyro
	vector<double> m = mass_->data(); // copy
	vector<double> k = stif_->data(); // copy
	if ((int)mass_->rsize() != nr)
		throw runtime_error(vastr(nr," retained dof, but mass is order ",mass_->rsize()));
	if ((int)stif_->rsize() != nr)
		throw runtime_error(vastr(nr," retained dof, but stif is order ",stif_->rsize()));
	vector<double> w(nr, 0.0);
	vector<double> rval(nr*nmodes()); // real eigenvectors
	int info = lapack::dsygv(1,"v","u",nr,&k[0],nr,&m[0],nr,&w[0]);
	if (info != 0) {
		ostringstream os;
		os << "could not compute the free-vibration solution";
		if (info > nr)
			os << ": singular mass matrix";
		else
			os << ": dsygv returned " << info;
		flaps::warning(os.str());
	} else {
		// save nmodes eigenvectors
		blas_copy(nr*nmodes(), &k[0], 1, &rval[0], 1);
		// ... and print the branch frequencies and eigenvector summary
		vector<double> freqs;
		for (int i=0; i<nmodes(); i++)
			freqs.push_back(sqrt(w[i])/(2.0*flaps::pi));
		flaps::info("free-vibration frequencies (Hz):",freqs);
	}
	return rval;
}

static
void
filter(vector<double>& T) {
	double meps = sqrt(std::numeric_limits<double>::epsilon());
	double eps = blas_snrm2(T.size(), &T[0], 1)*meps;
	for (auto& t : T) {
		if (std::abs(t) < eps)
			t = 0.0;
	}
}
void
Fem::
modal_reduction(const vector<double>& T, const string& id) {
// Reduce all output matrices with a triple-product: T'AT
// where T is a modes matrix from, e.g. freevibration(), BMA, CMS
// Output: all Fem matrices with ".id" appended to the mid
	Trace trc(1,"modal_reduction");

	size_t n = this->retained().size();
	size_t nmodes = T.size()/n;
	// lambda to create a matrix id with input id
	auto makemid = [&](const string& base) {
		string rval{base};
		if (!id.empty())
			rval += "." + id;
		return rval;
	};

#ifdef NEVER // new args
	const vector<double>& T = this->modal();
	nmodes(T.size()/n);
	if (nmodes()*n != T.size())
#else // NEVER // new args
	if (nmodes*n != T.size())
#endif // NEVER // new args
		throw runtime_error(vastr("dimension error in modal transform: should be (",
			n, ",",nmodes,"), but T size = ",T.size()));
#ifdef NEVER // new args
	if (nvz() == 0)
		throw runtime_error("vzdata has not been called before modal_reduction");
#endif // NEVER // new args

	// save the modes matrix for LTI psi
	Matrix modesmatrix(makemid("modes"), "modes", n, nmodes, false);
	modesmatrix.data() = T;
	modesmatrix.store();

	// save data necessary for modal animation
	// post-multiply gct_nodal by T: gct = gct_nodal*T (nvz,n)*(n,nmodes)
	vector<double> gct(nvz()*nmodes);
	blas_sgemm("n","n",nvz(),nmodes,n,1.0,&gct_nodal_[0],nvz(),&T[0],n,0.0,&gct[0],nvz());
	// create & save matrices without vzid
	vzmatrices(id, gct);

	// triple-product on all output matrices
	vector<double> tat(nmodes*nmodes, 0.0);
	if (stif_ != nullptr) {
		string mid{makemid("stif")};
		Matrix modalstif(mid, "stif", nmodes, nmodes, false);
		lapack::triprod(n, nmodes, &T[0], stif_->elem(), modalstif.elem());
		trc.dprintm(nmodes,nmodes,nmodes,modalstif.elem(), "modal stif");
		MM::exporter(mid, "modal stif", modalstif.elem(), nmodes,nmodes);
		modalstif.store();
	}
	if (mass_ != nullptr) {
		string mid{makemid("mass")};
		Matrix modalmass(mid, "mass", nmodes, nmodes, false);
		lapack::triprod(n, nmodes, &T[0], mass_->elem(), modalmass.elem());
		trc.dprintm(nmodes,nmodes,nmodes,modalmass.elem(),"modal mass");
		MM::exporter(mid, "modal mass", modalmass.elem(), nmodes,nmodes);
		modalmass.store();
	}

	Matrix* modalgaf{nullptr};
	if (dlm != nullptr) {
		vector<Matrix*> modaldlms;
		for (size_t i=0; i<dlm->nodal_gafs.size(); i++) {
			string mid = makemid(vastr("gaf",i+1));
			Matrix* modaldlm = new Matrix(mid, "gaf", nmodes, nmodes, true);
			lapack::triprod(n, nmodes, &T[0], dlm->nodal_gafs[i]->celem(), modaldlm->celem());
			// give the matrix a Pz_const pz with dlms_[i]->pz.params
			Pz* pz = dlm->nodal_gafs[i]->pz[0]->clone();
			pz->rsize() = nmodes;
			pz->csize() = nmodes;
			modaldlm->pz.push_back(pz);
			modaldlms.push_back(modaldlm);
			trc.dprintm(nmodes,nmodes,nmodes,modaldlm->celem(), mid);
			modaldlm->store();
		}
		// interpolate dlm aero
		modalgaf = new Matrix(makemid("gaf"), "gaf (dlm)", nmodes, nmodes, true);
		// interpolate modaldlms_ and add an IntPz
		modalgaf->pz.push_back(new IntPz(modaldlms));
		modalgaf->store();

#ifdef NEVER // experimental
		// Plotting requested?
		vector<Elem> elements;
		for (int i=0; i<nmodes; i++)
			elements.push_back(Elem(i,i));
		string plotpar{"rf"};
		string fixedpar{"ac"};
		Par* pp = gpset::find(fixedpar);
		pp->value(-0.152);
		modalgaf->plot_apf("modalgaf.apf", elements, 101, plotpar);
#endif // NEVER // experimental
	}

	// transforms for gyro and pua: rows of T at their nodes
	vector<int> ckey(nmodes,0);
	for (size_t i=0; i<nmodes; i++)
		ckey[i] = i+1;			// all columns

	if (!this->gyro_elements.empty()) {
		// does gyro need a custom function?
		bool needpz{false};
		for (auto& ge : gyro_elements) {
			if (!ge.amom_name().empty())
				needpz = true;
			vector<Freedom> rotations = get_rotations(ge.node());
			ge.modal() = vector<double>(3*nmodes, 0.0);
			// mergex(this->retained(), ckey, 1.0, T, rotations, ckey, 0.0, ge.modal());
			blas::copy(this->retained(), ckey, 1.0, T, rotations, ckey, 0.0, ge.modal());
			trc.dprintm(3,nmodes,3,ge.modal(),vastr("node",ge.node(),"modal"));
		}
		// do a triple-product on the full gyro
		string mid{makemid("gyro")};
		Matrix modalgyro(mid, "gyro", nmodes, nmodes, false);
		lapack::triprod(n, nmodes, &T[0], gyro_->elem(), modalgyro.elem());
		trc.dprintm(nmodes,nmodes,nmodes,modalgyro.elem(),"modal gyro");
		MM::exporter(mid, "modal gyro", modalgyro.elem(), nmodes,nmodes);

		// Now that each gyro element has it's modal transform we must
		// re-create gyrofcn:
		if (needpz) {
			string entrypt;
			string path = make_gyrofcn(entrypt, id);
			Custom::add(path);
			CustomPz* pz = new CustomPz(entrypt, nmodes, nmodes, 0);
			modalgyro.pz.push_back(pz);
		}
		modalgyro.store();
	}

	// for pgaf replace pgaf.T with pgaf.T*T XXX use a copy of pgaf?
	for (auto& ge : this->pgaf_elements) {
		vector<Freedom> rotations = get_rotations(ge.node());
		vector<double> TT(3*nmodes, 0.0);
		vector<double> tmp(3*nmodes, 0.0);
		blas::copy(this->retained(), ckey, 1.0, T, rotations, ckey, 0.0, tmp);
		blas_sgemm("n","n",3,nmodes,3,1.0,&ge.T[0],3,&tmp[0],3,0.0,&TT[0],3);
		// set relatively small elements in T to zero to reduce number of
		// Ad necessary for nonlinear pgaf
		filter(TT);
		MM::exporter("pgafT.mm", "modal pgaf T", &TT[0], 3, nmodes);
		ge.T = TT;
	}
	// for pgaf
	if (!pgaf_elements.empty()) {
		string entrypt;
		if (!pgaf_custom.empty()) {
			entrypt = Custom::add(pgaf_custom);
		} else {
			string path = make_pgaffcn(entrypt, id);
			// ... add it to custom.so ...
			Custom::add(path);
		}
		// ... and include it in the modal gaf parameterizations
		CustomPz* pz = new CustomPz(entrypt, nmodes, nmodes, 0);
		Matrix modalpgaf("pgaf","pgaf", nmodes,nmodes,true);
		modalpgaf.pz.push_back(pz);
		modalpgaf.store();
		// and the dlm matrix parameterizations
		if (modalgaf != nullptr) {
			modalgaf->pz.push_back(pz);
			modalgaf->store();
		}
	}

	// if this is a BMA model create CMS uf file for comparison
	//!! cms();
}
