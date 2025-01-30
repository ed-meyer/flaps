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
#include "lti.h"
#include "trace.h"

using namespace std;

namespace LTI {

void
evalgains(pset& plt, const vector<string>& gains, const vector<string>& phases,
		vector<double> td, vector<complex<Ad>>& g) {
	T_(Trace trc(2,"evalgains");)
	Ad sigma = plt.parval("sigma");
	Ad freq = plt.parval("freq");
	for (size_t i=0; i<gains.size(); i++) {
		Ad gain = plt.par_or_double(gains[i]);
		Ad phase = plt.par_or_double(phases[i]);
		double tdi = td[i];
		Ad ar(-sigma*tdi);
		Ad ai(phase - freq*tdi);
		complex<Ad> ex{ar,ai};
		g[i] = gain*exp(ex);
		T_(trc.dprint("g[",i,"] = ",g[i]);)
	}
}

int
eval (pset& plt, int nr, int nc, vector<complex<Ad>>& result,
	Matrix& A,
   Matrix& B,
   Matrix& C,
   Matrix& D,
   vector<Itd>& Atd,
   vector<double>& itd,
   vector<double>& otd,
   vector<int>& Sidx,
	// user-specified:
	const string& psiname,
	const string& stifname,
	vector<pair<int,double>>& ecolk,
   vector<string>& igains,
   vector<string>& iphases,
   vector<string>& ogains,
   vector<string>& ophases) {
// Evaluate an LTI, putting the matrix into "result"
	complex<Ad> zero(Ad(0.0),Ad(0.0));
	complex<Ad> t;
	T_(Trace trc(1,"LTI::eval");)

	// read K and psi
	Matrix K(stifname);
	Matrix psi(psiname);
	if (psi.is_complex())
		throw runtime_error("psi must be real");
	double* psip = psi.elem();
	
	int ne = K.rsize();
	int ns = A.rsize();
	if (nr != (int)(ns + ne))
		throw runtime_error(vastr("LTI::eval got nr ",nr,", expected ",ns+ne));
	int ni = psi.rsize();
	if (B.rsize() != A.rsize())
		throw runtime_error(vastr(B," does not match ",A));
	if (B.csize() != psi.rsize())
		throw runtime_error(vastr(B," does not match ",psi));
	if (K.rsize() != psi.csize())
		throw runtime_error(vastr(K," does not match ",psi));
	int no = C.rsize();

	Ad sigma = plt.parval("sigma");
	Ad freq = plt.parval("freq");
	complex<Ad> s(sigma,freq);

	T_(trc.dprint("s = ",s);)

	// (2,2) block: \prod_{k=1}^{n_t} A - sI
	// Do this first since if psi and KE were not specified we are done
	// XXX what does this mean?
	// CAdvector work(ns*ns);
	vector<complex<Ad>> work(ns*ns);
	A.eval(plt, work);
	T_(trc.dprintm(ns,ns,ns,work,"A");)
	// first mult by internal time delays...
	for (auto& td : Atd) {
		complex<Ad> ct = exp(-s*td.deltat);
		T_(trc.dprint("internal td: ",ct);)
		for (size_t k=0; k<td.elem.size(); k++) {
			int i = td.elem[k].row;
			int j = td.elem[k].col;
			//!! result[IJ(ne+i,ne+j,nr)] *= ct;
			work[IJ(i,j,ns)] *= ct;
		}
	}
	// ... then subtract sI...
	for (int j=0; j<ns; j++)
		work[IJ(j,j,ns)] -= s;
		//!! result[IJ(ne+j,ne+j,nr)] -= s;
	T_(trc.dprintm(ns,ns,ns,work,"td*A - sI");)
	// ... and insert into result
	for (int j=0; j<ns; j++) {
		for (int i=0; i<ns; i++) {
			result[IJ(ne+i,ne+j,nr)] = work[IJ(i,j,ns)];
		}
	}

	// if neither psi nor KE were specified we are done XXX???
	//!! if (psi.empty() && KE == nullptr)
		//!! return true;

	// S:
	//CAdvector S(ni);
	vector<complex<Ad>> S(ni);
	for (int i=0; i<ni; i++) {
		if (Sidx[i] == 0)
			S[i] = complex<Ad>(1.0);
		else if (Sidx[i] == 1)
			S[i] = s;
		else
			S[i] = s*s;
	}
	T_(trc.dprintm(ni,1,ni,S,"S");)
	// G, H: input and output gain and time delays
	// CAdvector G(ni);
	vector<complex<Ad>> G(ni);
	evalgains(plt, igains, iphases, itd, G);
	// CAdvector H(no);
	vector<complex<Ad>> H(no);
	evalgains(plt, ogains, ophases, otd, H);
	// change sign on H XXX or do this when (1,1) and (1,2) blocks created?
	for (int i=0; i<no; i++)
		H[i] *= -1.0;
	T_(trc.dprintm(ni,1,ni,G,"G");)
	T_(trc.dprintm(no,1,no,H,"H");)

	// evaluate KE -> ke
	//!! KE->eval(plt, ke);
	// CAdvector ke(ne*no);
	vector<complex<Ad>> ke(ne*no);
	// CAdvector wk(ne*ne);
	vector<complex<Ad>> wk(ne*ne);
	K.eval(plt, wk);
	for (int j=0; j<no; j++) {
		int col = ecolk[j].first;
		double sf = ecolk[j].second;
		for (int i=0; i<ne; i++)
			ke[IJ(i,j,ne)] = sf*wk[IJ(i,col,ne)];
	}
	T_(trc.dprintm(ne,no,ne,ke,"KE");)

	double dmax{0.0};
	for (auto& d : D.data())
		dmax = std::max(abs(d), dmax);
	if (dmax > 0.0) {
		// (1,1) block: -KE*H*D*G*S*Psi
		// CAdvector HDGS(no*ni);
		vector<complex<Ad>> HDGS(no*ni);
		double* dp = D.elem();
		for (int i=0; i<no; i++) {
			for (int j=0; j<ni; j++) {
				HDGS[IJ(i,j,no)] = H[i]*dp[IJ(i,j,no)]*G[j]*S[j];
			}
		}
		// CAdvector KEHDGS(ne*ni);
		vector<complex<Ad>> KEHDGS(ne*ni);
		// ftn_cgemm("n", "n", ne, ni, no, Complex(1.0), ke->gm, ne,
		// 		&HDGS[0], no, Complex(0.0), &KEHDGS[0], ne);
		for (int i=0; i<ne; i++) {
			for (int j=0; j<ni; j++) {
				t = zero;
				for (int k=0; k<no; k++) {
					t += ke[IJ(i,k,ne)]*HDGS[IJ(k,j,no)];
				}
				KEHDGS[IJ(i,j,ne)] = t;
			}
		}
		// vector<ADs> KEHDGSPsi(ne*ne, zero);
		// multCR (ne, ne, ni, &KEHDGS[0], psi->gm, &KEHDGSPsi[0]);
		// for (j=0; j<ne; j++) {
		// 	ftn_ccopy (ne, &KEHDGSPsi[IJ(0,j,ne)], 1, &result[IJ(0,j,nr)], 1);
		// }
		for (int i=0; i<ne; i++) {
			for (int j=0; j<ne; j++) {
				t = zero;
				for (int k=0; k<ni; k++) {
					t += KEHDGS[IJ(i,k,ne)]*psip[IJ(k,j,ni)];
				}
				result[IJ(i,j,nr)] = t;
			}
		}
	}

	// (1,2) block: -KE*H*C
	// CAdvector HC(no*ns);
	vector<complex<Ad>> HC(no*ns);
	double* cp = C.elem();
	for (int i=0; i<no; i++) {
		for (int j=0; j<ns; j++) {
			HC[IJ(i,j,no)] = H[i]*cp[IJ(i,j,no)];
		}
	}
	T_(trc.dprintm(no,ns,no,HC,"HC");)
	// vector<ADs> KEHC(ne*ns, zero);
	//ftn_cgemm("n", "n", ne, ns, no, Complex(1.0), ke->gm, ne,
	//		&HC[0], no, Complex(0.0), &KEHC[0], ne);
	//for (j=0; j<ns; j++) {
	//	ftn_ccopy (ne, &KEHC[IJ(0,j,ne)], 1, &result[IJ(0,ne+j,nr)], 1);
	//}
	for (int i=0; i<ne; i++) {
		for (int j=0; j<ns; j++) {
			t = zero;
			for (int k=0; k<no; k++) {
				t += ke[IJ(i,k,ne)]*HC[IJ(k,j,no)];
			}
			result[IJ(i,j+ne,nr)] = t;
		}
	}

	T_(trc.dprintm(ni,ne,ni,psip,"psi");)

	// (2,1) block: B*G*S*Psi: (ns,ni)(ni,ni)(ni,ni)(ni,ne) = (ns,ne)
	// CAdvector GSPsi(ni*ne);	// (ni,ni)(ni,ni)(ni,ne)
	vector<complex<Ad>> GSPsi(ni*ne);
	for (int i=0; i<ni; i++) {
		T_(trc.dprint("S[",i,"] = ",S[i]);)
		for (int j=0; j<ne; j++) {
			GSPsi[IJ(i,j,ni)] = G[i]*S[i]*psip[IJ(i,j,ni)];
		}
	}
	T_(trc.dprintm(ni,ne,ni,GSPsi,"GSPsi");)
	//vector<ADs> BGSPsi(ns*ne, zero);
	// multRC(ns, ne, ni, &B[0], &GSPsi[0], &BGSPsi[0]);
	// for (j=0; j<ne; j++) {
	// 	ftn_ccopy (ns, &BGSPsi[IJ(0,j,ns)], 1, &result[IJ(ne,j,nr)], 1);
	// }
	double* bp = B.elem();
	// CAdvector BGSPsi(ns*ne);
	vector<complex<Ad>> BGSPsi(ns*ne);
	for (int i=0; i<ns; i++) {
		for (int j=0; j<ne; j++) {
			t = zero;
			for (int k=0; k<ni; k++) {
				t += bp[IJ(i,k,ns)]*GSPsi[IJ(k,j,ni)];
			}
			//!! result[IJ(i+ne,j,nr)] = t;
			BGSPsi[IJ(i,j,ns)] = t;
		}
	}
	T_(trc.dprintm(ns,ne,ns,BGSPsi,"BGSPsi");)
	// insert into result
	for (int i=0; i<ns; i++) {
		for (int j=0; j<ne; j++)
			result[IJ(i+ne,j,nr)] = BGSPsi[IJ(i,j,ns)];
	}

	// Scale the last ns columns
	// scaleX(ne, ns, result);
	T_(trc.dprintm(nr,nc,nr,result,"LTI::result");)

	return true;
}
}	// namespace LTI
