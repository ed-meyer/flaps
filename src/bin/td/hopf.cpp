//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "exim.h"
#include "hopf.h"
#include "lapack.h"
#include "matrix.h"
#include "pac.h"
#include "pset.h"
#include "Pz.h"
#include "trace.h"

using namespace std;

void
Fmatrix(vector<Ad>& F);
int
fjac(vector<double>& x, vector<double>& f, vector<double>& jac);
int
bif(complex<double> w, vector<complex<double>> v, double vtas);
string
run_octlab(Specs& sp);

#ifdef NEVER // new fjac
static void
invert(vector<double> P) {
// invert an (n,n) matrix, return it in the input P
	T_(Trace trc(1,"invert");)

	int n = sqrt(P.size());
	T_(trc.dprintm(n,n,n,P,"inverting");)
	vector<double> I(n*n, 0.0);
	for (int i=0; i<n; i++)
		I[i*(n+1)] = 1.0;
	vector<int> ipiv(n,0);

	//!! int info = lapack::dgesv(n,n,P.data(),n,ipiv.data(),I.data(),n);
	char equed[2];
	vector<double> r(n, 0.0);
	vector<double> c(n, 0.0);
	vector<double> ferr(n, 0.0);
	vector<double> berr(n, 0.0);
	vector<double> pf(n*n, 0.0);
	vector<double> x(n*n, 0.0);
	vector<double> psave = P;
	double rcond;
	int info = lapack::dgesvx("E", "N", n, n, P.data(), n, pf.data(), n,
		ipiv.data(), equed, r.data(), c.data(), I.data(), n, x.data(), n,
		&rcond, ferr.data(), berr.data());
	if (info != 0)
		throw runtime_error(vastr("failed to invert M+rho/2*R2: dgesv info ",info));
	T_(trc.dprint("rcond ",rcond);)
	T_(trc.dprintm(n,n,n,x,"inverted");)
	// check by multiplying x*P
	blas_sgemm("n", "n", n, n, n, 1.0, psave.data(), n,
			x.data(), n, 0.0, pf.data(), n);
	T_(trc.dprintm(n,n,n,pf,"P*x");)
	//!! P = I;
	P = x;
}
#endif // NEVER // new fjac

template<typename Dta, typename Dtb, typename Dtc>
void
block_insert(const vector<Dta>& A, int ic, int jc,
		Dtb b, vector<Dtc>& C) {
// Insert (na,na) matrix A times b into 0-based block (ic,jc) of
// the (nc,nc) C; nc must be a multiple of na!
	int na = sqrt(A.size());
	int nc = sqrt(C.size());
	if (nc%na > 0)
		throw runtime_error(vastr("attempt to insert an order ", na,
			" matrix into an order ", nc, " matrix"));

	if ((ic+1)*na > nc || (jc+1)*na > nc)
		throw runtime_error(vastr("attempt to insert an order ", na,
			" matrix into block (", ic, ',', jc, ") of an order ", nc, " matrix"));

	for (int i=0; i<na; i++) {
		for (int j=0; j<na; j++) {
			C[IJ(i+na*ic,j+na*jc,nc)] = b*A[IJ(i,j,na)];
		}
	}
}

void
categorize(vector<double> const& vr, vector<double> const& wr, vector<double> const& wi,
		vector<complex<double>>& w, vector<vector<complex<double>>>& v) {
// create complex eigenvalues/eigenvectors from the vr/(wr,wi) from
// lapack::dgeevx. 2 cases to consider:
// 1) wi[j] = 0:  the eigenvector is real, vr[:j]
// 2) wi[j] = -wi[j+1]:  the eigenvector is complex, v[j] = vr[j] + i*vr[j+1]
//                       only return wi[j], ignore wi[j+1]
	T_(Trace trc(1,"categorize");)

	v.clear();
	w.clear();
	int n = wr.size();
	for (int j=0; j<n; j++) {
		vector<complex<double>> vj;
		if (wi[j] == 0.0) {
			w.push_back(wr[j]);
			for (int i=0; i<n; i++)
				vj.push_back(complex<double>(vr[i+j*n]));
			normalize_inf(n, vj.data(), 1, true);
			T_(trc.dprint("real: ",wr[j],", v = ",flaps::summarize(n,vj.data()));)
		} else if (is_equal(wr[j],wr[j+1],8) && is_equal(wi[j],-wi[j+1],8)) {
			complex<double> wj(wr[j],wi[j]);
			w.push_back(complex<double>(wr[j],wi[j]));
			for (int i=0; i<n; i++)
				vj.push_back(complex<double>(vr[i+j*n], vr[i+(j+1)*n]));
			normalize_inf(n, vj.data(), 1, true);
			v.push_back(vj);
			j++;		// skip the conjugate
			T_(trc.dprint(w.size(),") ",wj,", v = ",flaps::summarize(n,vj.data()));)
		} else {
			complex<double> wj(wr[j],wi[j]);
			w.push_back(complex<double>(wr[j],wi[j]));
			for (int i=0; i<n; i++)
				vj.push_back(vr[i+j*n]);
			normalize_inf(n, vj.data(), 1, true);
			v.push_back(vj);
			T_(trc.dprint("no conjugate: ",w.size(),") ",wj,", v = ",flaps::summarize(n,vj.data()));)
		}
			
	}
	T_(trc.dprint("returning ",w.size()," eigenvalues/vectors");)
}
int
hopf(Specs& sp) {
	T_(Trace trc(1,"hopf");)
	int n = sp.stif->rsize();
	int nbeta = sp.beta.size();
	int nt = n*(2 + nbeta);
	vector<double> x(nt, 0.0);
	vector<double> f(nt, 0.0);
	vector<double> J(nt*nt, 0.0);
	// use the global pset
	pset& plt = gpset::get();

	Par* vtasp = plt.findp("vtas");
	double vtas = vtasp->value();
	Par* dpressp = plt.findp("dpress");
	double dpress = dpressp->value();
	Par* rhop = plt.findp("rho");
	double rho = rhop->value();

	// check the Jacobian
	vector<string> pnames;
	string reason{"check jacobian"};
	//!! check_jac(x, f, jac, fjac, nullptr, nullptr, reason, pnames);

	// compute eigenvalues of J = df/dy at y=0 over a range of velocities
	// where f = Fy so J = F at y=0
	double v0 = vtasp->min();
	double v1 = vtasp->max();
	double vdel{0.0};
	if (sp.maxstep > 1)
		vdel = (v1 - v0)/(sp.maxstep-1);
	vector<double> results;
	//!! int nres = n + 2*nbeta;
	int nres = nt+3;
	vector<double> res(nres, 0.0);
	vector<Ad> work(nt*nt);
	vector<complex<double>> w;
	vector<vector<complex<double>>> v;
	vector<double> wr(nt, 0.0);
	vector<double> wi(nt, 0.0);
	vector<double> vl(nt*nt, 0.0);
	vector<double> vr(nt*nt, 0.0);
	bool biffound{false};
	for (int j = 0; j<sp.maxstep; j++) {
		vtas = v0 + j*vdel;
		vtasp->value(vtas);
		plt.eval();
		vtas = vtasp->value();
		dpress = dpressp->value();
		rho = rhop->value();

		// J = F when y=0
		Fmatrix(work);
		extract(work, "", J.data());

		//!! lapack::dgeev("N","V",nt,J.data(),nt,wr.data(),wi.data(),
			//!! vl, nt, vr.data(), nt);
		int ihi, ilo;
		vector<double> scale(nt, 0.0);
		double abnrm;
		vector<double> rconde(nt, 0.0);
		vector<double> rcondv(nt, 0.0);
		lapack::dgeevx("B","V","V","E",nt,J.data(),nt,wr.data(),wi.data(),
			vl.data(), nt, vr.data(), nt, ilo, ihi, scale.data(), abnrm,
			rconde.data(),rcondv.data());

		// create w and v: only pos wi
		categorize(vr,wr,wi,w,v);

		T_(trc.dprint("real(w)",wr);)
		// search for wr >= 0
		cout << "v " << vtas << " max real: " <<
			*std::max_element(wr.cbegin(),wr.cend()) << endl;
		for (size_t i=0; i<w.size(); i++) {
			if (w[i].real() >= 0) {
				cout << "Hopf bifurcation found at v = " << vtas << ", w: " << w[i]
					<< ", vector: " << flaps::summarize(nt,v[i].data()) << endl;
				biffound = true;
				bif(w[i], v[i], vtas);
			}
		}
		// save the results
		int k{0};
		std::fill(res.begin(), res.end(), 0.0);
		res[k++] = vtas;
		res[k++] = dpress;
		res[k++] = rho;
		for (size_t i=0; i<w.size(); i++) {
			//!! if (k == nres)
				//!! throw runtime_error("more than n+nbeta pos freq");
			//!! if (wi[i] >= 0.0)
				res[k++] = w[i].real();
				//!! res[k++] = wi[i];
		}
		vappend(results, nres, res.data());
		if (biffound)
			break;
	}
	// plot the results - results is (2*nt+1, nc) real
	vector<tuple<string,string>> params;
	params.push_back(make_tuple("vtas", "Velocity (tas)"));
	params.push_back(make_tuple("dpress", "Dynamic Pressure (Pa)"));
	params.push_back(make_tuple("rho", "Density (kg/m^3)"));
	for (int i=0; i<nt; i++) {
		string wri = vastr("wr",i+1);
		params.push_back(make_tuple(wri,wri));
	}
	Apf::exporter("hopf", params, results, "hopf.apf", false);
	return 0;
}

void
Fmatrix(vector<Ad>& F) {
// compute F as an (nt,nt) Ad matrix
	T_(Trace trc(1,"Fmatrix");)
	Specs& sp = specs();

	Ad dpress = gpset::find("dpress")->advalue();
	Ad vtas = gpset::find("vtas")->advalue();
	Ad rho = gpset::find("rho")->advalue();

	int nbeta = sp.beta.size();
	int nt = sqrt(F.size());
	int n = sqrt(sp.mass.size());
	// nt = (2*nbeta)*n
	if (nt != (2+nbeta)*n)
		throw runtime_error(vastr("dimensions are wrong: n ",n,", nbeta ",nbeta,", nt ",nt));

	// insert identity matrices...
	vector<Ad> I(n*n, 0.0);
	for (int i=0; i<n; i++)
		I[i*(n+1)] = 1.0;
	block_insert(I, 0, 1, 1.0, F);
	// ... and the R matrices premultiplied by Pinverse
	vector<Ad> PR(n*n, 0.0);
	Ad one(1.0);
	for (int k=0; k<nbeta; k++) {
		blas::gemm("n", "n", n, n, n, 1.0, sp.P.data(), n,
				sp.R[k+3].data(), n, 0.0, PR.data(), n);
		block_insert(PR, 1, k+2, dpress, F);
		block_insert(I, k+2, 1, one, F);
		Ad scl = -sp.beta[k]*vtas;
		block_insert(I, k+2, k+2, scl, F);
	}
	// -P(K - qR0)
	vector<Ad> KR(n*n);
	vector<complex<Ad>> work(n*n);
	sp.stif->eval(gpset::get(), work);
	for (int i=0; i<n*n; i++)
		KR[i] = work[i].real();

	blas::axpy(n*n, -dpress, sp.R[0].data(), 1, KR.data(), 1);

	vector<Ad> PKR(n*n, 0.0);
	blas::gemm("n", "n", n, n, n, 1.0, sp.P.data(), n,
				KR.data(), n, 0.0, PKR.data(), n);
	block_insert(PKR, 1, 0, -1.0, F);

	// -P(G + V - q/v_t R_1) = -P(G + V - rho v_t/2 R_1)
	vector<Ad> GR(n*n);
	Ad t = -rho*vtas/2.0;
	for (int i=0; i<n*n; i++)
		GR[i] = t*sp.R[1][i];
#ifdef NEVER // no gyro or viscous damping yet
	Matrix* gyro = Matrix::find_desc("gyro");
	if (gyro != nullptr)
		GR = gyro->data();
	Matrix* vdamp = Matrix::find_desc("vdamp");
	if (vdamp != nullptr)
		blas_axpy(n*n, 1.0, vdamp->elem(), 1, GR.data(), 1);
#endif // NEVER // no gyro or viscous damping yet
	vector<Ad> PGR(n*n, 0.0);
	blas::gemm("n", "n", n, n, n, 1.0, sp.P.data(), n,
				GR.data(), n, 0.0, PGR.data(), n);
	block_insert(PGR, 1, 1, -1.0, F);

	T_(trc.dprintm(nt,nt,nt,F,"F matrix");)
}

void
Jmatrix(vector<Ad> const& yad, vector<double>& f, vector<double>& F, vector<double>& J) {
// compute J(y) = df/dy = F + dF/dy*y
// where y' = F(y)y
// Compute Fy in Ad, then extract f, F, and J
	T_(Trace trc(1,"Jmatrix");)
	Specs& sp = specs();
	int n = sp.stif->rsize();
	int nbeta = sp.beta.size();
	int nt = (2+nbeta)*n;
	int ntsq = nt*nt;

	// compute F as (nt,nt) Ad
	vector<Ad> Fad(ntsq);
	Fmatrix(Fad);
	// extract the values -> F
	extract(Fad, "", F.data());
	// multiply Fad*yad = fad
	vector<Ad> fad(nt);
	blas::gemv("n",nt,nt,1.0,F.data(),nt,yad.data(), 1, 0.0, fad.data(), 1);

	// J = d(Fy)/dy
	// copy F to J...
	J = F;
	// ...then overwrite columns with the Ad derivatives of Fy that
	// are Ad components of y
	vector<double> df(nt);
	vector<string> const& adpar = Ad::adnames();
	for (size_t k=0; k<adpar.size(); k++) {
		int j = evparse(adpar[k]);
		if (j >= 0) {
			extract(fad, adpar[k], df.data());
			blas_copy(nt, df.data(), 1, &J[j*nt], 1);
		}
	}
}

int
bif(complex<double> w, vector<complex<double>> v, double vtas) {
// trace a bifurcation from gcnorm=0 to gmax
	T_(Trace trc(1,"bif");)
	Specs& sp = specs();

	int n = sp.stif->rsize();
	int nbeta = sp.beta.size();
	int nt = (2+nbeta)*n;
	int nx = 2*nt + 3;
	int nf = 2*nt + 2;
	// make freq, gcnorm & vtas indep
	Par* freq = gpset::find("freq");
	freq->set_indep();
	Par* gcnorm = gpset::find("gcnorm");
	gcnorm->set_indep();
	Par* vtasp = gpset::find("vtas");
	vtasp->set_indep();
	// create eigenvector components (nt complex)
	gpset::get().make_eigv(nt);
	// initialize Ad
	vector<string> adparnames{"freq", "gcnorm", "vtas"};
	Ad::initialize(adparnames);
	gpset::get().realloc();
	// x: {freq, gcnorm, vtas, y}
	vector<double> x(nx);
	x[0] = w.imag();
	x[1] = 0.0;		// gcnorm
	x[2] = vtas;
	blas_copy(2*nt, (double*)&v[0], 1, &x[3], 1);
	vector<double> f(nf);
	vector<double> jac(nf*nx);
	fjac(x, f, jac);
	T_(trc.dprintm(nx,1,nx,x.data(), "bif x");)
	T_(trc.dprintm(nf,1,nf,f.data(), "bif f");)
	T_(trc.dprintm(nf,nx,nf,jac.data(), "bif jac");)
	return 0;
}

int
fjac(vector<double>& x, vector<double>& f, vector<double>& jac) {
	Specs& sp = specs();

	int n = sp.stif->rsize();
	int nbeta = sp.beta.size();
	int nt = (2+nbeta)*n;	// size of F
	int nx = x.size();
	// sanity checks:
	// x should be 2*nt + 3: y + (freq, gcnorm, vtas)
	if (nx != 2*nt + 3)
		throw runtime_error(vastr("nx = ",nx,", should be ",2*nt+3));
	// f should be 2*nt + 2: Fx + normalization
	int nf = f.size();
	if (nf != 2*nt + 2)
		throw runtime_error(vastr("nf = ",nf,", should be ",2*nt+2));
	// and jac should be (nf,nx)
	int nfx = jac.size();
	if (nfx != nf*nx)
		throw runtime_error(vastr("jac size = ",nfx,", should be ",nf*nx));

	// x layout: freq, gcnorm, vtas, y
	int ny = 2*nt;
	vector<double> y(ny, 0.0);
	for (int i=0; i<ny; i++)
		y[i] = x[i+3];

	// update gpset
	Par* freq = gpset::find("freq");
	freq->value(x[0]);
	Par* gcnorm = gpset::find("gcnorm");
	gcnorm->value(x[1]);
	Par* vtas = gpset::find("vtas");
	vtas->value(x[2]);
	gpset::get().updateEv(y);
	gpset::get().eval();
	// set Ad parameter derivatives to 1
	gpset::get().set_adpar_derivs();
	// get the eigenvector as Ad vector
	vector<complex<Ad>> adev(nt);
	gpset::get().getadev(adev);
	
	// compute F as an (nt,nt) Ad matrix
	vector<Ad> F(nt*nt);
	Fmatrix(F);

	// multiply F*x in Ad
	vector<complex<Ad>> Fx(nt);
	blas::gemv("n",nt,nt,1.0,F.data(),nt,adev.data(), 1, 0.0, Fx.data(), 1);

	// f:
	for (int i=0; i<nt; i++) {
		f[2*i] = Ad::real(Fx[i]).data()[0];
		f[2*i+1] = Ad::imag(Fx[i]).data()[0];
	}

	// extract the values of F into a real (nt,nt)
	vector<double> Fr(nt*nt);
	extract(F, "", Fr.data());
	// cast it to complex
	vector<complex<double>> Fc(nt*nt);
	for (int i=0; i<nt*nt; i++)
		Fc[i] = Fr[i];
	// insert it into the (nf,nx) jac
	blas::real_rep(nt,nt,Fc.data(),nt,&jac[3*nf],nf);

	// insert partials wrt freq, gcnorm, vtas
	vector<string> const& adpar = Ad::adnames();
	for (size_t k=0; k<adpar.size(); k++) {
		int j;
		if (adpar[k] == "freq")
			j = 0;
		else if (adpar[k] == "gcnorm")
			j = 1;
		else if (adpar[k] == "vtas")
			j = 2;
		else {
			j = evparse(adpar[k]) + 3;
		}
		for (int i=0; i<nt; i++) {
			jac[2*i+j*nf] = Ad::real(Fx[i]).data()[k+1];
			jac[2*i+1+j*nf] = Ad::real(Fx[i]).data()[k+1];
		}
	}

	// last 2 eqns: y'y-1 and y[ireal]
	for (int j=0; j<ny; j++)
		jac[nt+(j+3)*nf] = 2.0*x[j];
	jac[nt+1+(sp.ireal+3)*nf] = 1.0;

	return 0;
}	// fjac

void
biffcn (vector<double> const& z, vector<double>& g) {
// function to integrate eqn 7.24 (Seydel) to locate a Hopf bifurcation
// given an initial estimate of y, T, V
// compute g(z) = z' where
// z = [y T v h]
// g = [ Tf(y,v) 0 0 Tf_y(y,v)h ] = [TFy 0 0 TJh]
// f = Fy, J = df/dy = (dF/dy)y + F
// bcs: [ y(0)-y(1) h(0)-h(1) f_y[1:]h(0)  h_1(0)-1 ]
// y' = f(y,v) = F(y,v)y
	Specs& sp = specs();
	int n = sp.stif->rsize();
	int nbeta = sp.beta.size();
	int nt = (2+nbeta)*n;	// size of F
	int ntsq = nt*nt;
	vector<double> y(nt, 0.0);
	vector<double> h(nt, 0.0);
	int ih = nt + 2;
	for (int i=0; i<nt; i++) {
		y[i] = z[i];
		h[i] = z[ih++];
	}
	double T = z[nt];
	//!! double V = z[nt+1];	// z[nt+1] is velocity
	vector<double> f(nt, 0.0);
	vector<double> F(ntsq, 0.0);
	vector<double> J(ntsq, 0.0);
	// create yad from y: set Ad parameters
	vector<Ad> yad(nt);
	vector<string> const& adpar = Ad::adnames();
	for (int i=0; i<nt; i++)
		yad[i] = y[i];
	for (size_t k=0; k<adpar.size(); k++) {
		int j = evparse(adpar[k]);
		if (j >= 0)
			yad[j].deriv(k,1.0);
	}
	// get f, F and J
	Jmatrix(yad, f, F, J);

	// now create g = [ Tf(y,v) 0 0 Tf_y(y,v)h ] = [TFy 0 0 TJh]
	for (int i=0; i<nt; i++)
		g[i] = T*f[i];
	g[nt] = 0.0;
	g[nt+1] = 0.0;
	vector<double> Jh(nt);
	blas_sgemv("n", nt, nt, T, J.data(), nt, h.data(), 1, 0.0, Jh.data(), 1);
	for (int i=0; i<nt; i++)
		g[nt+2+i] = Jh[i];
}

void
locate_bif(double V, double T, vector<double>& y) {
// z = [y T v h]
	Specs& sp = specs();

	int nt = y.size();
	int nz = 2*nt + 2;
	vector<double> z(nz, 0.0);
	for (int i=0; i<nt; i++) {
		z[i] = y[i];
	}
	z[nt] = T;
	z[nt+1] = V;

   // open memory-mapped files to send/receive data to matlab/octave
   double* enviar = open_octlab("g", nz);
   double* recibir = open_octlab("z", nz);

   // start matlab/octave with the initial betas
   string file = run_octlab(sp);

}

	
string
run_octlab(Specs& sp) {
// create a script to run in matlab/octave
// XXX run it in a different thread?
	int maxfunevals{60};

	string cmd("matlab");
	string script("locate_bif");
	string file = vastr(script,".m");
	ofstream ofs(file, std::ios::trunc);
	ofs << "cd " << getftmp() << endl;
	ofs << "x = [";
	string sep;
	for (auto bi : sp.beta) {
		ofs << sep << bi;
		sep = "; ";
	}
	ofs << "];\n";
	ofs << "g = ode\n";
	ofs << "g.ODEFcn = hopfpt\n";
	ofs << "sol = solve(g,0.0,1.0);\n";
	ofs << "quit\n";
	// create the command line to run: only matlab works - octave doesn't have memmap
	string cmdline{vastr(cmd," -batch \"",script,"\"&")};
	//!! if (matlab_cmd == "octave")
		//!! cmdline = vastr(matlab_cmd, " optbeta.m");

 // unset LD_PRELOAD while running matlab (libstdc++ conflict)
   string ld_preload = getEnv("LD_PRELOAD");
   if (!ld_preload.empty())
      unsetenv("LD_PRELOAD");

	// run matlab/octave in the background
	int status = system(cmdline.c_str());

	return cmdline;
}
