//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// substructuring and branch-modes model reduction

#include "config.h"

#include <assert.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "blas.h"
#include "conv.h"
#include "exim.h"
#include "fem.h"
#include "lapack.h"
#include "lexer.h"
#include "matrix.h"
#include "trace.h"

using namespace std;

std::ostream&
operator<<(std::ostream& s, const Comp& t) {
	s << t.ss->id() << ": " << t.nmodes << " modes";
	return s;
}
std::ostream& operator<<(std::ostream& s, const BMA& t) {
	s << "root: " << t.root << std::endl;
	for (auto& bi : t.branches)
		s << "branch: " << bi << std::endl;
	return s;
}

std::ostream& operator<<(std::ostream& s, const CMS& t) {
	s << "CMS substructures: " << t.ss << std::endl;
	return s;
}

// parsing functions (_f)
bool
Fem::
ss_f (const Tok& tok) {
// ss{ id=s, beams{...}, mass{...},   }
	T_(Trace trc(2,"ss_f");)
	string id;
	vector<int> node_numbers;
	Material matl;
	vector<Beam*> beams;
	vector<Mass*> masses;
	Freev freev;
	vector<int> proj;

	// lambda handlers
	vector<Tok*> unrec = flaps::lexer(tok.lopt, {
		{"id", [&](const Tok& p) { id = p.srhs; return true;}},
		//!! {"node(s)?", [&](const Tok& p) { node_numbers = p.ivec; return true; }},
      {"material", [&](const Tok& p) {matl = parse_material(p); return true;}},
      {"beam(s)?", [&](const Tok& p) {beams = parse_beam(p, matl, freev); return true; }},
      {"conmass", [&](const Tok& p) {return parse_conmass(p,freev,beams,masses); }},
		{"project", [&](const Tok& p) {proj = p.ivec; return true;}},
	});

	if (!unrec.empty())
		throw runtime_error(vastr("unrecognized preferences: ",unrec));

	SS* ss = new SS(id,beams,masses,freev,proj);
	T_(trc.dprint("got ss: ",*ss);)

	substructures_.push_back(ss);
	return true;
}

SS*
Fem::
find_ss(const string& id) {
	T_(Trace trc(2, "find_ss");)
	for (auto& si : substructures_) {
		if (id == si->id()) {
			T_(trc.dprint("returning \"",si->id(),"\"");)
			return si;
		}
	}
	throw runtime_error(vastr("substructure \"",id,"\" is undefined"));
}

bool
Fem::
bma_f(const Tok& p) {
// definition of a BMA: root=ss{n}, branches=(ss1{n}, ss2,..),
	T_(Trace trc(2,"Fem::bma_f");)
	this->bma = new BMA;
	vector<Tok*> unrec = flaps::lexer(p.lopt, {
		{"root",[&](const Tok& p) { 
			bma->root.ss = find_ss(p.srhs);
			if (!p.ropt.empty())
				if (!str2int(p.ropt, bma->root.nmodes))
					throw runtime_error(vastr("unrecognized nmodes option in ",p));
			return true;
		}},
		{"branch(es)?", [&](const Tok& p) {
			for (size_t i=0; i<p.svec.size(); i++) {
				Comp comp;
				comp.ss = find_ss(p.svec[i]);
				if (!p.roptvec[i].empty())
					str2int(p.roptvec[i], comp.nmodes);
				bma->branches.push_back(comp);
			}
			return true;
		}}
	});
	T_(trc.dprint("got BMA: ",*bma);)
	return true;
}

bool
Fem::
cms_f(const Tok& tok) {
// definition of a CMS: ss=(ss1{n}, ss2{n}, ..),
	T_(Trace trc(2,"Fem::cms_f");)
	this->cms = new CMS;
	vector<Tok*> unrec = flaps::lexer(tok.lopt, {
		{"ss", [&](const Tok& p) {
			for (size_t i=0; i<p.svec.size(); i++) {
				Comp comp;
				comp.ss = find_ss(p.svec[i]);
				if (!p.roptvec[i].empty())
					str2int(p.roptvec[i], comp.nmodes);
				cms->ss.push_back(comp);
			}
			return true;
		}}
	});
	T_(trc.dprint("got CMS: ",*cms);)
	return true;
}


SS::
SS (const string& id, const vector<Beam*>& beams,
	const vector<Mass*>& masses, Freev& freev, const vector<int>& project) {
// SS constructor
	id_ = id;
	beams_ = beams;
	retained_ = freev;
	masses_ = masses;
	proj_ = project;
}

void
Freev::
add(const Freedom& f) {
// Add freedom "f" to this set of freedoms if it is not already there
   if (std::find(begin(),end(),f) == end())
         push_back(f);
}


void
SS::
assemble(const vector<SS*> substructures) {
// For this substructure:
// - Split freedoms into interior and attachment
// - assemble the stiffness and mass matrices for a substructure
// - compute Gia = K_{ii}^{-1}K_{ia}
	T_(Trace trc(1,"SS::assemble");)

	size_t n = retained_.size();
	// divide retained_ into interior and attachment freedoms
	for (auto& fi : retained_) {
		for (auto& si : substructures) {
			if (si->id() == id()) continue;
			auto begin = si->retained().begin();
			auto end = si->retained().end();
			if (std::find(begin,end,fi) != end)
				attach_.add(fi);
		}
	}
	// ... then if a freedom is not in attach_ it must be interior
	auto b = attach_.begin();
	auto e = attach_.end();
	for (auto& fi : retained_) {
		if (std::find(b, e, fi) == e)
			interior_.add(fi);
	}

	stif_ = vector<double>(n*n, 0.0);
	mass_ = vector<double>(n*n, 0.0);
	// Merge the stiffness matrix
	for (auto bp : beams_)
		blas::copy(bp->freedoms(), bp->freedoms(), 1.0,
			bp->data(), this->retained_, this->retained_, 1.0, stif_);
	T_(trc.dprintm(n,n,n,stif_,vastr(id()," stif matrix"));)

	// Merge the mass matrix
	// if any of the mass elements has a cg parameter (cg_p),
	// create a custom function XXX downstream
	// create each (6,6) elemental mass matrix
	for (auto& mi : this->masses_) {
		vector<double> elem(36, 0.0);
		vector<double> cgv = mi->orient();
		vector<Freedom> freedoms;
		for (int i=0; i<6; i++)
			freedoms.push_back(Freedom(mi->node(),i+1));
		if (!cgv.empty()) {
			assert(cgv.size() == 3);
			blas::scal(3, mi->cg(), &cgv[0], 1);
		}
		conmass(mi->mass(), mi->moi(), cgv, elem);
		T_(trc.dprintm(6,6,6,&elem[0],vastr("node ",mi->node()," mass ",mi->mass()));)
		// insert into the gross matrix
		blas::copy(freedoms, freedoms, 1.0, elem, retained_, retained_, 1.0, mass_);
	}
	T_(trc.dprintm(n,n,n,mass_,vastr(id()," mass matrix"));)

	// compute the transformation from attach to interior dof
	Gia_ = gia();
	// cantilevered vibration solution
}


vector<double>
SS::
gia() {
// Compute the Guyan transformation to condense the interior dof
// to the interace:
//   G_{ia} = -K^{-1}_{ii} K_{ia}
// where "a" subscripts are attachment dof, and "i" are interior
	T_(Trace trc(1,"gia");)

	vector<Freedom> ret = this->retained();	// all freedoms for this branch
	vector<Freedom> att = this->attach();		// attachment freedoms
	vector<Freedom> inter = this->interior();	// interior freedoms
	T_(trc.dprint("retained: ",ret);)
	T_(trc.dprint("attach: ",att);)
	T_(trc.dprint("interior: ",inter);)

	// extract K_{ii} and K_{ia}
	int na = att.size();
	int ni = inter.size();
	vector<double> kia(ni*na, 0.0);
	kii = vector<double>(ni*ni, 0.0);	// kii is a member
	blas::copy(ret, ret, 1.0, stif(), inter, inter, 0.0, kii);
	blas::copy(ret, ret, 1.0, stif(), inter, att, 0.0, kia);
	T_(trc.dprintm(ni,ni,ni,kii,vastr(id(),"_kii"));)
	T_(trc.dprintm(ni,na,ni,kia,vastr(id(),"_kia"));)

	int info = lapack::dposv("u",ni,na,&kii[0],ni,&kia[0],ni);
	if (info != 0) {
		blas::copy(ret, ret, 1.0, stif(), inter, inter, 0.0, kii);
		string file{vastr(id(),"_kii.mm")};
		MM::exporter(file,"non-positive-definite Kii",kii.data(),ni,ni);
		throw runtime_error(vastr("Kii is not positive-definite: see ",file)); 
	}

	// Gia is returned from dposv in kia
	// change the sign
	blas::scal(ni*na, -1.0, kia.data(), 1);
	T_(trc.dprintm(ni,na,ni,kia,vastr(id(),"_gia"));)
	return kia;
}

vector<double>
SS::
lump(vector<double>& kaa) {
// Condense this branches mass onto the attachment dof using the
// rigid transformation T_{ia} = rigid_T:
//   M_{aa} = - [M_{ai} T_{ia} + (M_{ai} T_{ia})'] + T_{ia}' M_{ii} T_{ia}
// XXX or should we return T_{ra} where r = a+i, all freedoms in this branch
// this would require inserting Gia and I_a into the correct locations
//   M_{aa} = T' M T
//   T = |      I_a          |
//       | -K^{-1}_ii K_{ia} |
// where "a" subscripts are attachment dof, and "i" are interior
// The return mass matrix must be *added* to the root substructure attachment
	T_(Trace trc(1,"lump");)

	vector<Freedom> ret = this->retained();
	vector<Freedom> att = this->attach();
	vector<Freedom> inter = this->interior();
	T_(trc.dprint("retained: ",ret);)
	T_(trc.dprint("attach: ",att);)
	T_(trc.dprint("interior: ",inter);)

	int na = att.size();
	int ni = inter.size();

	// Gia = -K_{ii}^{-1} K_{ia}
	vector<double> gia = Gia();
	T_(trc.dprintm(ni,na,ni,gia,"gia");)

	vector<double> rval(na*na, 0.0);

	// insert I_a and G_{ia} into an (nr,na) matrix T
	int nr = ret.size();
	vector<double> Ia(na*na, 0.0);
	vector<double> T(nr*na, 0.0);
	for (int i=0; i<na; i++)
		Ia[i*(1+na)] = 1.0;
	blas::copy(att, att, 1.0, Ia, ret, att, 0.0, T);
	blas::copy(inter, att, 1.0, gia, ret, att, 0.0, T);
	// triple products for kaa, maa
	lapack::triprod(nr, na, T.data(), mass().data(), rval.data());
	kaa = vector<double>(na*na, 0.0);
	lapack::triprod(nr, na, T.data(), stif().data(), kaa.data());

	T_(trc.dprintm(na,na,na,rval,vastr(id(),"_maa"));)
	T_(trc.dprintm(na,na,na,kaa,vastr(id(),"_kaa"));)

	// kaa should have 1-6 rb modes
	vector<double> kaas{kaa};
	vector<double> w(na, 0.0);
	lapack::dsyev("v", "u", na, kaas.data(), na, w.data());
	T_(trc.dprint("eigenvalues of kaa: ",w);)
	return rval;
}

ostream&
operator<<(ostream& s, const vector<SS*>& t) {
	if (!t.empty()) {
		s << "substructures:\n";
		for (auto ss : t)
			s << *ss << endl;
	}
	return s;
}

ostream&
operator<<(ostream& s, const SS& ss) {
	s << "substructure " << ss.id() << endl;
	// for (auto i : ss.node_numbers())
		// s << i << endl;
	// s << "   beams:\n" << ss.beams() << endl;
	s << "   beam nodes:\n";
	for (auto i : ss.beams())
		s << "     " << i->node1() << ", " << i->node2() << endl;
	s << "   freedoms: " << ss.retained() << endl;
	
	return s;
}

void
SS::
projection(vector<double>& v) {
// project each column of v onto the corresponding 
// axis (1-3) in member "proj_"
	T_(Trace trc(2,"SS::projection");)
	T_(trc.dprint("ss ",id()," proj: ",proj_);)
	// v is (m,n) where n = proj_.size()
	int n = proj_.size();
	int m = v.size()/n;
	vector<double> p(n);
	vector<double> rval(m*n);
	for (int j=0; j<n; j++) {
		int i = proj_[j]-1; // 0b
		// extract row i of v
		for (int k=0; k<n; k++)
			p[k] = v[i+k*m];
		// normalize p to keep unit mass normalization
		blas::normalize2(n, &p[0], 1, true);
		// mult v*p -> column j of rval
		blas::gemv("n",m,n,1.0,&v[0],m,&p[0],1,0.0,&rval[j*m],1);
		// normalize column j of rval
		// blas::normalize_inf(m, &rval[j*m], 1, true);
	}
	// replace v with rval
	v = rval;
}

vector<double>
BMA::
transform (const Freev& gross) {
// Branch-mode reduction:
// 1) merge lumped_mass from each branch with root mass
// 2) do an eigensolution with root mass and stif for nmodes
// 3) extract dof rows of eigenvectors at attachment node,
//    premultiply by rigid transform, insert in branch dof rows
//    to give rigid motion of branch
// 4) do an eigensolution on the branch cantilevered at
//    the attachment nodes; insert these eigenvectors in the
//    branch dof rows, trailing columns
	T_(Trace trc(1,"BMA::transform");)

	assert(root.ss != nullptr);
	// 1) add lumped_masses, lumped_stif to root attachments
	vector<double> M = root.ss->mass(); // copy
	vector<double> K = root.ss->stif();  // copy
	vector<Freedom> ret{root.ss->retained()};
	size_t nr{root.ss->retained().size()};
	for (auto bp : branches) {
		vector<double> kaa;
		vector<double> maa = bp.ss->lump(kaa);
		vector<Freedom> att{bp.ss->attach()};
		blas::copy(att, att, 1.0, maa, ret, ret, 1.0, M);
		blas::copy(att, att, 1.0, kaa, ret, ret, 1.0, K);
	}

	// 2) eigensolution on root: vectors are returned in K
	vector<double> w(nr, 0.0);
	int info = lapack::dsygv(1,"v","u",nr,&K[0],nr,&M[0],nr,&w[0]);
	if (info != 0) {
		ostringstream os;
		os << "could not compute the branch-modes root eigensolution";
		if (info > (int)nr) {
			os << ": singular mass matrix (" << info-nr << ") - see root_mass.mm";
			MM::exporter("root_mass.mm","Branch-mode root mass",&M[0],nr,nr);
		} else
			os << ": dsygv returned " << info;
		throw runtime_error(os.str());
	}

	// we return an (n,m) matrix of branch modes, where n = Fem retained(),
	// and m = root_p->nmodes() + sum(branches->nmodes())
	int n = gross.size();
	int mr = root.nmodes;
	int m{mr};
	for (auto bp : branches)
		m += bp.nmodes;
	vector<double> rval(n*m, 0.0);
	// the eigenvectors were returned in K - save only the first mr
	vector<double> v(nr*mr, 0.0);
	blas::copy(nr*mr, &K[0], 1, &v[0], 1);
	T_(trc.dprintm(nr,mr,nr,&K[0],"root eigenvectors");)
	// the eigenvectors go in columns 1-mr; ckey is a vector of
	// ints from 1:m, next is the 1b column # to start putting vectors
	vector<int> ckey;		// rval column keys
	for (int i=0; i<m; i++)
		ckey.push_back(i+1);
	vector<int> ckeyr;   // root eigenvectors columns
	for (int i=0; i<mr; i++)
		ckeyr.push_back(i+1);
	// copy the root eigenvectors to rval
	blas::copy(root.ss->retained(),ckeyr,1.0,v,gross,ckey,0.0,rval);
	int next = mr+1;

	// 3) for each branch: extract attachment rows, multiply by Gia & merge into rval
	for (auto bp : branches) {
		int na = bp.ss->attach().size();
		vector<double> att(na*mr, 0.0);
		blas::copy(root.ss->retained(),ckeyr,1.0,v,bp.ss->attach(),ckeyr,0.0,att);
		T_(trc.dprintm(na,mr,na,att,vastr(bp.ss->id()," attachment motion"));)
		// rigid: (ni,na) # interior freedoms in the branch by # attach freedoms
		// att:   (na,mr) # attach freedoms by root nmodes (mr)
		// vint = rigid*att: (ni,mr)
		int ni = bp.ss->interior().size();
		vector<double> vint(ni*mr, 0.0);
		blas::gemm("n","n",ni,mr,na,1.0,bp.ss->Gia().data(),ni,&att[0],
			na, 0.0, vint.data(), ni);
		blas::copy(bp.ss->interior(),ckeyr,1.0,vint,gross,ckey,0.0,rval);
		// do a triple-product vint' Kii vint to check determinate-ness
		int nbp = bp.ss->retained().size();
		vector<double> vbp(nbp*mr, 0.0);
		blas::copy(bp.ss->interior(),ckeyr,1.0,vint,bp.ss->retained(),ckeyr,0.0,vbp);
		blas::copy(bp.ss->attach(),ckeyr,1.0,att,bp.ss->retained(),ckeyr,0.0,vbp);
		vector<double> vkv(mr, 0.0);
		for (int j=0; j<mr; j++) {
			double* vp = &vbp[j*nbp];
			//!! lapack::triprod(nbp, mr, vbp.data(), bp->stif().data(), vkv.data());
			lapack::triprod(nbp, 1, vp, bp.ss->stif().data(), &vkv[j]);
			vkv[j] /= std::max(w[j],1.0);
		}
		T_(trc.dprint("vkv scaled = ",vkv);)
		// eigenvectors of the branch cantilevered at the attachment
		// start in column "next"
		vector<int> ckeyi;
		int nmodes = bp.nmodes;
		assert(nmodes > 0);
		for (int i=0; i<nmodes; i++)
			ckeyi.push_back(next+i);
		next += nmodes;

		// 4) eigensolution on cantilevered branch XXX do this in SS ctor?
		vector<double> K(ni*ni, 0.0);
		vector<double> M(ni*ni, 0.0);
		blas::copy(bp.ss->retained(), bp.ss->retained(), 1.0, bp.ss->stif(),
			bp.ss->interior(), bp.ss->interior(), 0.0, K);
		blas::copy(bp.ss->retained(), bp.ss->retained(), 1.0, bp.ss->mass(),
			bp.ss->interior(), bp.ss->interior(), 0.0, M);
		T_(trc.dprintm(ni,ni,ni,K,vastr(bp.ss->id()," cantilevered stif"));)
		T_(trc.dprintm(ni,ni,ni,M,vastr(bp.ss->id()," cantilevered mass"));)
		w = vector<double>(ni, 0.0);
		int info = lapack::dsygv(1,"v","u",ni,&K[0],ni,&M[0],ni,&w[0]);
		if (info != 0)
			throw runtime_error(vastr("eigensolution on cantilevered ",bp.ss->id(),
				" failed: info = ",info));
		T_(trc.dprint("cantilevered eigenvalues: ",w);)
		// save the first nmodes eigenvectors
		vector<double> v(ni*nmodes, 0.0);
		blas::copy(ni*nmodes, &K[0], 1, &v[0], 1);
		if ((int)bp.ss->proj().size() == nmodes)
			bp.ss->projection(v);
		string title = vastr(bp.ss->id()," cantilevered modes");
		T_(trc.dprintm(ni,nmodes,ni,v,title);)
		MM::exporter(title,title,&v[0],ni,nmodes);
		// the eigenvectors go in columns ckeyi
		blas::copy(bp.ss->interior(),ckeyi,1.0,v,gross,ckey,0.0,rval);
		// print the cantilevered frequencies
		vector<string> sum;
		for (int i=0; i<nmodes; i++)
			sum.push_back(lapack::eigstring(sqrt(w[i])/(2.0*flaps::pi),ni,&v[i*ni]));
		flaps::info(bp.ss->id()," frequencies (Hz): ",sum);
	}
	T_(trc.dprintm(n,m,n,rval,"branch-mode transform");)
	// export to a matrix-market file for viewing
	MM::exporter("branch_modes.mm", "Branch modes", &rval[0], n, m);

	return rval;
}
	
vector<double>
CMS::
transform (const Freev& gross) {
// Create a CMS transformation from nodal to CMS generalized coordinates
// - for each substructure:
//   - extract mii, kii
//   - eigensolution on mii, kii, save phi_im
//   - compute G_{ia} = -K^{-1}_{ii} K_{ia}
	T_(Trace trc(2,"CMS::transform");)

	// the output matrix (rval) has 2 sets of columns
	// - 0-(nm-1): interior modes for each ss (nm = sum(nmodes))
	// - nm-(nm+na-1) where na is the number of unique attachment dof
	// Column keys are Freedoms: the first set are the first
	// nmodes interior freedoms for each ss, the second are the unique
	// attachment dof
	vector<Freedom> rcol;		// columns of rval: 0-nm+na-1
	int n = gross.size();		// number of rows in rval
	// determine the number of columns in rval
	vector<Freedom> att;		// all unique attachment freedoms
	for (auto& comp : this->ss) {
		SS* ss = comp.ss;
		int ni = comp.nmodes;
		const Freev& ssi = ss->interior();
		// the first nmodes interior freedoms:
		for (int i=0; i<ni; i++)
			rcol.push_back(ssi[i]);
		for (auto& fi : ss->attach())
			add2vector(fi, att);
	}
	// 2nd set: attachment dof
	for (auto& fi : att)
		rcol.push_back(fi);

	int nc = rcol.size();		// number of columns in rval
	T_(trc.dprint(nc," modes total, all attachment freedoms: ",att);)
	T_(trc.dprint("rval columns: ",rcol);)

	vector<double> rval(n*nc);

	for (auto& comp : this->ss) {
		// extract mii, kii
		SS* ss = comp.ss;
		const vector<Freedom>& ret = ss->retained();
		const vector<Freedom>& it = ss->interior();
		int ni = it.size();
		vector<double> mii(ni*ni);
		vector<double> kii(ni*ni);
		blas::copy(ret,ret,1.0,ss->mass(),it,it,0.0,mii);
		blas::copy(ret,ret,1.0,ss->stif(),it,it,0.0,kii);
		// eigensolution: eigenvectors are returned in K
		vector<double> w(ni, 0.0);
		int info = lapack::dsygv(1,"v","u",ni,kii.data(),ni,mii.data(),ni,&w[0]);
		if (info != 0) {
			ostringstream os;
			os << "could not compute the cantilevered eigensolution for " << ss->id();
			if (info > (int)ni) {
				string mmfile = vastr(ss->id(),"_mii.mm");
				os << ": singular mass matrix (" << info-ni << ") - see " << mmfile;
				MM::exporter(mmfile,"singular mass",mii.data(),ni,ni);
			} else
				os << ": dsygv returned " << info;
			throw runtime_error(os.str());
		}
		T_(trc.dprint(ni," eigenvalues of ",ss->id()," interior: ",w);)
		int m = comp.nmodes;
		if (m <= 0)
			throw runtime_error(vastr("nmodes was not given for ss ",ss->id()));
		vector<Freedom> imcol(m);
		// the column keys for kii are the first m interior freedoms
		for (int i=0; i<m; i++)
			imcol[i] = it[i];
		blas::copy(it,imcol,1.0,kii.data(),gross,rcol,0.0,rval.data());
	}
	// ... then insert all the Gia matrices (aka constraint modes), and
	// the identity in the attachment rows/cols
	for (auto& comp : this->ss) {
		SS* ss = comp.ss;
		blas::copy(ss->interior(),ss->attach(),1.0,ss->Gia(), gross, rcol, 0.0, rval);
		int na = ss->attach().size();
		vector<double> identity(na*na, 0.0);
		for (int i=0; i<na; i++)
			identity[(na+1)*i] = 1.0;
		blas::copy(ss->attach(),ss->attach(),1.0,identity, gross, rcol, 0.0, rval);
	}
	return rval;
}

