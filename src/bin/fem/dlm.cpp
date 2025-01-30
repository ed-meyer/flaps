//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <string>
#include <vector>

#include "config.h"
#include "conv.h"
#include "fem.h"
#include "lapack.h"

using namespace std;

vector<complex<double>>
Dsolve(int np, int nc, vector<complex<double>> D, vector<complex<double>> w);

void computeInputMeshSegment(int n, int m, const vector<double>& x0, double span,
	double dihedral, double sweep, double cr, double tr, vector<double>& Xi,
	vector<double>& Xo, vector<double>& Xr, vector<double>& dXav);
void computeSurfaceSegment(int n, int m, const vector<double>& x0, double span,
	double dihedral, double sweep, double cr, double tr,
	vector<double>& X, vector<int>& conn);
void computeInfluenceMatrix(vector<complex<double>>& D, double omega,
	double U, double M, int np, const vector<double>& Xi, const vector<double>& Xo,
	const vector<double>& Xr, const vector<double>& dXav, int symmetric,
	bool steadykernel, double epstol);
void getModeBCs(int np, int n, double* mode, double* X, int* conn,
	double* vwash, double* dwash);
void washes(int np, int n, int nc, const double* mode, const double* X,
	const int* conn, double* vwash, double* dwash);
void nodal_forces(int np, int n, int nc, const complex<double>* Cp,
	double* X, const int* conn, complex<double>* forces);
void add_reflen(double reflen);

ostream&
operator<<(ostream& s, const Dlm& t) {
	s << t.panels[1] << " spanwise panels, " << t.panels[0] << " chordwise panels\n";
	s << "aero center(s) " << t.acs << endl;
	s << "reference length " << t.reflen << endl;
	s << "span " << t.semi_span << endl;
	s << "chord " << t.chord << endl;
	s << "taper " << t.taper << endl;
	s << "sweep " << t.sweep << endl;
	s << "dihedral " << t.dihedral << endl;
	s << "mach " << t.machs << endl;
	s << "reduced frequenies " << t.rfs << endl;
	s << "reduced stabilty factors " << t.rsfs << endl;
	return s;
}

const Node*
findnode(int n, const vector<Node>& grid) {
// find a Node number "n" in a vector<Node> (grid)
	for (auto& ni : grid)
		if (ni.number() == n)
			return &ni;
	throw runtime_error(vastr("DLM node ",n," is undefined"));
}

void
Dlm::
assemble (const Grid& grid, const Freev& freev) {
// Create DLM gaf matrices from this Dlm
	T_(Trace trc(1,"Dlm::assemble");)

	// default dlm nodes: all fem grid nodes
	if (this->nodes.empty()) {
		for (auto ni : grid.nodes())
			this->nodes.push_back(ni.number());
	}
	
	int m = this->panels[0];  // # chordwise panels
	// the spanwise panels are by default determined by the specified Fem nodes
	int n = this->nodes.size()-1;  // number of spanwise panels
	if (this->panels.size() == 1) {
		this->panels.push_back(n);  // spanwise panels
		const Node* node0 = findnode(this->nodes[0], grid.nodes());
		const Node* noden = findnode(this->nodes[n], grid.nodes());
		vector<double> x0 = node0->coord();
		vector<double> ds = noden->coord();
		blas::axpy(3, -1.0, &x0[0], 1, &ds[0], 1);
		this->semi_span = blas::snrm2(3, &ds[0], 1);
	} else {
		n = this->panels[1];  // # spanwise panels
		throw runtime_error("spanwise panels not implemented yet");
	}

	// The output matrix will be parameterized with at least rf, and possibly
	// rsf, mach, and ac. Check for multiple values specified and set defaults.
	// not all parameters effect the paneling (e.g. rf), but it's cheap so
	// for the sake of parameters that do (e.g. ac) just do it for all parameters
	vector<Par*> params;
	// reduced frequency: must have at least 1 (no defaults)
	if (this->rfs.empty())
		throw runtime_error("no reduced frequencies specified");
	Par* rf_p{gpset::find("rf")};
	if (rf_p != nullptr)
		params.push_back(rf_p);
	// reduced stability factor sigma/vtas
	if (this->rsfs.empty())
		this->rsfs.push_back(0.0);		// default: rsf=0
	Par* rsf_p{gpset::find("rsf")};
	if (this->rsfs.size() > 1 && rsf_p != nullptr)
		params.push_back(rsf_p);
	// Mach number
	if (this->machs.empty()) {
		this->machs.push_back(0.5);		// default: Mach=0.5
		flaps::info("default Mach number is ",this->machs[0]);
	}
	Par* mach_p{gpset::find("mach")};
	if (this->machs.size() > 1 && mach_p != nullptr)
		params.push_back(mach_p);
	// aerodynamic center
	if (this->acs.empty()) {
		this->acs.push_back(0.0);		// default: ac=0
		flaps::info("default aerodynamic center is ",this->acs[0]);
	}
	Par* ac_p{gpset::find("ac")};
	if (this->acs.size() > 1 && ac_p != nullptr)
			params.push_back(ac_p);

	flaps::info("creating doublet-lattice matrices with:\n",*this);

	int nc = freev.size();

	// 4 nested loops ac's, rf's, rsf's, mach's
	int ordinal{0};
	for (double ac : this->acs) {
		// set parameter values
		if (ac_p != nullptr)
			ac_p->value(ac);
		vector<double> root{ac, 0.0, 0.0};
		// create the panels
		vector<double> Xi(3*n*m);
		vector<double> Xo(3*n*m);
		vector<double> Xr(3*n*m);
		vector<double> dXav(n*m);
		computeInputMeshSegment(n, m, root, this->semi_span,
			this->dihedral, this->sweep, this->chord, this->taper, Xi, Xo, Xr, dXav);
		T_(trc.dprint("Xi: ",Xi,", Xo: ",Xo,", Xr: ",Xr);)

		// compute conn (panel node numbers) and X (nodal coord)
		this->conn = vector<int>(4*n*m, 0);
		this->X = vector<double>(3*(n+1)*(m+1), 0.0);
		computeSurfaceSegment(n, m, root, this->semi_span, this->dihedral,
			this->sweep, this->chord, this->taper, this->X, this->conn);
		T_(trc.dprint("X: ",this->X,", conn: ",this->conn);)

		// compute the structure-to-aero transform (nr,nc) where
		// nr = 3*(m+1)*(n+1) and nc = retained().size(). It depends on
		// ac: root
		int na = (m+1)*(n+1);
		int nr = 3*na;
		this->T = this->transform(m, n, this->X, grid.nodes(), freev);
		int nt = this->T.size();
		assert(nt == nr*nc);

		// check that reflen has been set...
		// XXX why not just set reflen=1?
		if (this->reflen == 0.0)
			this->reflen = 1.0; // default b: 1.0 so that rf = \omega/vtas
		// ... then add reflen to the equations for rf, rsf
		if (this->reflen != 1.0)
			add_reflen(this->reflen);

		//compute the D matrix for each k-value
		int np{n*m};  // number of panels
		vector<complex<double>> D(np*np);
		for (auto mach : this->machs) {
			if (mach_p != nullptr)
				mach_p->value(mach);
			for (auto rf : this->rfs) {
				rf_p->value(rf);
				double U{1.0};
				double omega{rf*U/this->reflen};
				T_(trc.dprint("rf = ",rf,", omega = ",omega);)
				int symm{1};
				bool steadykernel{false};
				double epstol{std::numeric_limits<double>::epsilon()};
				computeInfluenceMatrix(D, omega, U, mach, np, Xi, Xo, Xr,
					dXav, symm, steadykernel, epstol);
				T_(trc.dprintm(np,np,np,&D[0],"D matrix at rf=",rf);)
				vector<complex<double>> w(np, complex<double>(0.0));

				// get the downwash and normal wash
				vector<double> vwash(np*nc, 0.0);
				vector<double> dwash(np*nc, 0.0);
				washes(np, na, nc, &this->T[0], &this->X[0],
					&this->conn[0], &vwash[0], &dwash[0]);

				for (auto rsf : this->rsfs) {
					if (rsf_p != nullptr)
						rsf_p->value(rsf);
					
					// add w = i\omega/V*vwash + dwash
					// XXX w = (p/U)*vwash + dwash
					complex<double> ik(rsf, omega/U);
					w = vector<complex<double>>(np*nc, complex<double>(0.0,0.0));
					blas::copy(np*nc, &dwash[0], 1, (double*)&w[0], 2);
					for (size_t i=0; i<w.size(); i++)
						w[i] += ik*vwash[i];

					// compute the Cp = D^{-1}w
					vector<complex<double>> Cp = Dsolve(np, nc, D, w);
					T_(trc.dprint("Cp = ",Cp);)

					// compute the forces on each aero grid node: F (3,na,nc)
					// where na = # aero nodes, nc = # gc
					vector<complex<double>> forces(nr*nc, complex<double>(0.0));
					nodal_forces(np, na, nc, &Cp[0], &this->X[0], &this->conn[0], &forces[0]);
					T_(trc.dprintm(nr,nc,nr,forces,"forces");)

					// transform the forces to the FEM grid - first cast T to complex
					//   Q = T' F:   (nc,nr)*(nr,nc)
					vector<complex<double>> CT(nt, complex<double>(0.0));
					blas::copy(nt, &this->T[0], 1, reinterpret_cast<double*>(&CT[0]), 2);
					Matrix* qi = new Matrix(vastr("gaf",++ordinal,".nodal"), "gaf", nc, nc, true);
					complex<double>* Q = qi->celem();
					complex<double> alpha(1.0), beta(0.0);
					blas::gemm("t","n",nc,nc,nr,alpha,&CT[0],nr,&forces[0],nr,beta,Q,nc);
					T_(trc.dprintm(nc,nc,nc,Q,"Q at rf =",rf);)

					// give the Matrix (qi) a Pz_const parameterization with "params"
					qi->pz.push_back(new Pz_const(params,nc,nc));
					// save it for modal reduction
					//!! this->dlms_.push_back(qi);
					nodal_gafs.push_back(qi);
					// and store it
					qi->store();
				} // rsfs
			} // rfs
		} // machs
	} // acs

	// create the output dlm matrix...
	gaf_nodal = new Matrix("gaf.nodal","dlm gaf", nc, nc, true);

	// ... and add an IntPz to the matrix
	gaf_nodal->pz.push_back(new IntPz(nodal_gafs));

#ifdef NEVER // moved to class Dlm
	// copy this Dlm to the Fem
	this->dlm_ = new Dlm(de);
#endif // NEVER // moved to class Dlm

#ifdef NEVER // experimental
	// Plotting requested?
	vector<Elem> elements{ {33,33} };
	string plotpar{"rf"};
	string fixedpar{"ac"};
	Par* pp = gpset::find(fixedpar);
	pp->value(-0.152);
	this->gaf_->plot_apf("gaf.apf", elements, 101, plotpar);
#endif // NEVER // experimental

}

vector<complex<double>>
Dsolve(int np, int nc, vector<complex<double>> D, vector<complex<double>> w) {
// Input:
//   D     (np,np) complex
//   w     (np,nc) complex
//
	int nrhs{nc};
	vector<complex<double>> af(np*np, 0.0);
	vector<int> ipiv(np);
	char equed[2];
	vector<double> r(np, 0.0);
	vector<double> c(np, 0.0);
	vector<complex<double>> rval(np*nc, 0.0);
	double rcond;
	vector<double> ferr(nrhs,0.0);
	vector<double> berr(nrhs,0.0);

	lapack::zgesvx("e", "n", np, nrhs, &D[0], np, &af[0], np, &ipiv[0], equed,
		&r[0], &c[0], &w[0], np, &rval[0], np, &rcond, &ferr[0], &berr[0]);
	return rval;
}

vector<double>
Dlm::
transform(int m, int n, const vector<double>& X,
	const vector<Node>& grid, const vector<Freedom>& freedoms) {
// create the transformation matrix from the structure (grid) to
// the dlm grid (X).
// Input:
//   m         number of chordwise panels
//   n         number of spanwise panels
//   X         (3,m+1,n+1) matrix of coordinates
//   grid      FEM grid
//   freedoms  FEM freedoms
// Output: T (nr,nc) nr: number of displacements at aero nodes: 3*(m+1)*(n+1)
//                   nc: number of structural freedoms
	T_(Trace trc(1,"transform");)

	assert(!freedoms.empty());

	// for now we only treat spanwise grid == structural grid
	// each chordwise node must have the same y coord as X
	// X: (3,m+1,n+1)
	int nr = 3*(m+1)*(n+1);
	int nc = freedoms.size();
	this->T = vector<double>(nr*nc, 0.0);
	// for each Fem freedom (columns of T)
	for (int k=0; k<nc; k++) {
		int dof = freedoms[k].dof();
		int node = freedoms[k].node();
		// treat T as a (3,m+1,n+1,nc) matrix, X as (3,m+1,n+1)
		// for each chordwise node i at span j
		auto idx = find(this->nodes.begin(), this->nodes.end(), node);
		if (idx != this->nodes.end()) {
			int j = idx - this->nodes.begin();
			// tx, ty, tz
			if (dof < 4) {
				for (int i=0; i<m+1; i++)
					T[IJKL(dof-1,i,j,k,3,m+1,n+1)] = 1.0;		
			// ry: pitch
			} else if (dof == 5) {
				//!! const Node* ni = grid.find(node);
				const Node* ni = findnode(node, grid);
				double x = ni->x();
				for (int i=0; i<m+1; i++)
					T[IJKL(2,i,j,k,3,m+1,n+1)] = x - X[IJK(0,i,j,3,m+1)];
			// rz: yaw - not necessary for dlm but improves visualization
			} else if (dof == 6) {
				const Node* ni = findnode(node, grid);
				double x = ni->x();
				for (int i=0; i<m+1; i++)
					T[IJKL(1,i,j,k,3,m+1,n+1)] = X[IJK(0,i,j,3,m+1)] - x;
			}
		}
	}

	T_(trc.dprintm(nr,nc,nr,T,"aero transform");)
	return T;
}

void
add_reflen(double reflen) {
// add a reference length (aka "b") to the global equations for
// reduced frequency (rf) and reduced stability factor (rsf)

	string neweqn{vastr(reflen,"*freq/vtas")};
	Par* existing = gpset::find("rf");
	if (existing != nullptr) {
		existing->equation = neweqn;
		existing->candidates[0] = neweqn;
		flaps::info("changed \"rf = freq/vtas to \"rf = ",neweqn,"\"");
	}
	// rsf
	neweqn = vastr(reflen,"*sigma/vtas");
	Par newrsf{"rsf(Reduced stability factor)", {neweqn} };
	newrsf.pref = 1;
	existing = gpset::find("rsf");
	if (existing != nullptr) {
		string oldpar{vastr(*existing)};
		existing->upgrade(&newrsf);
		flaps::info("changed \"rsf = sigma/vtas to \"rsf = ",neweqn,"\"");
	}
}


// computeInfluenceMatrix
//     |
// computeQuadDoubletCoeff
//     |
// evalKernelNumerator
//     |
// evalK1K2Coeff
//     |
// approxKernelIntegrals

// computeInputMeshSegment
// computeSurfaceSegment
// addCpForces
// getModeBCs


// dlm.f90: A simple fortran code for oscillating aerodynamic analysis
//
// Copyright (c) 2015 Graeme Kennedy. All rights reserved. 
//
// The following code is a basic Doublet Lattice method implementation
// that is designed to be run from python. The computationally
// expensive computations are performed at the fortran level for
// efficiency. This code is intended to be used for flutter
// computations but could also be used for general subsonic oscilatory
// flows.

void
approxKernelIntegrals(complex<double>& I0, complex<double>& J0, double u1, double k1) {
  // Compute the approximate values of the integrals I0 and J0. These
  // integrals are required for the computation of the kernel function
  // at points along the bound vortex line. The integrals are
  // approximated based on the following expression from Desmarais:
  //
  // 1 - u/sqrt(1 + u^2) \approx \sum_{n} a_{n} exp(-p_{n}*u)
  //
  // Where p_{n} = b*2**n. The output I0 and J0 are defined as follows:
  //
  // I0 = int_{u1}^{infty} (1 - u/sqrt(1 + u^2)) du
  // J0 = int_{u1}^{infty} u (1 - u/sqrt(1 + u^2)) du
  //
  // The function input/output are defined as follows:
  // 
  // Input:
  // u1:  (M*R - x0)/(beta^2*x0)
  // k1:  omega*r1/U
  //
  // Output:
  // I0:  Approximate value of the integral I0
  // J0:  Approximate value of the integral J0 
  
  //!! integer :: n
  //!! real(kind=dtype) :: pn, a(12)
  vector<double> a{0.000319759140, -0.000055461471, 0.002726074362,
	  0.005749551566, 0.031455895072, 0.106031126212,
	  0.406838011567, 0.798112357155, -0.417749229098,
	  0.077480713894, -0.012677284771, 0.001787032960};
  //!! complex(kind=dtype) :: expn, invn, kval

  // Set the parameter b - the exponential in the kernel integral
  double b{0.009054814793};

  // Evaluate the integral for I0 and J0
  I0 = complex<double>(0.0);
  J0 = complex<double>(0.0);

  // Evaluate the integral for positive values of u1
  for (int n=0; n<12; n++) {
     double pn = b*pow(2.0, n+1);
     complex<double> kval(pn, k1);
     complex<double> expn = exp(-kval*u1);
     complex<double> invn = 1.0/kval;
     I0 = I0 + a[n]*expn*invn;
     J0 = J0 + a[n]*expn*(kval*u1 + 1.0)*invn*invn;
  }
}

void
evalK1K2Coeff(complex<double>& Kf1, complex<double>& Kf2, double r1,
	double u1, double k1, double beta, double R, double M) {
  // Compute the value of the K1 and K2 functions given the values of
  // the local panel variables. This code calls the function
  // approxKernelIntegrals to obtain the values of I0 and J0 which are
  // used in the evaluation of the kernel coefficients.
  //
  // Input:
  // u1:   (M*R - x0)/(beta^2*x0)
  // k1:   omega*r1/U
  // beta: sqrt(1.0 - M**2)
  // R:    sqrt(x0**2 + (beta*r1)**2)
  // M:    Mach number
  //
  // Output:
  // K1:   first kernel function
  // K2:   second kernel function
  
  // Local temporary variables
  complex<double> expk, I0, J0, I1, I2;
  complex<double> I10, I20, I11, I21;
  // real(kind=dtype) :: invsqrt, invR, u1pos
  // complex(kind=dtype) :: expk, I0, J0, I1, I2 
  // complex(kind=dtype) :: I10, I20, I11, I21

  // Constant definitions
  //complex(kind=dtype), parameter :: I = cmplx(0.0, 1.0, kind=dtype)
  complex<double> I{0.0, 1.0};
  double zero{0.0};
  double one{1.0};
  double two{2.0};

  double invR = one/R;
  double invsqrt = one/sqrt(one + u1*u1);

  if (u1 < zero) {
     // Use separate logic when the argument u1 is negative. This is
     // required since the approximate integrals for I0 and J0 are not
     // defined for negative values of u1. 
     approxKernelIntegrals(I0, J0, zero, k1);

     // Evaluate I1
     complex<double> I10 = one - I*k1*I0;
     
     // Evaluate 3*I2
     complex<double> I20 = two - I*k1*I0 + J0*k1*k1;

     // Evaluate the approximate integrals I0 and J0
     double u1pos = -u1;
     approxKernelIntegrals(I0, J0, u1pos, k1);

     // Compute the temporary variable values that will be used below
     expk = exp(-I*k1*u1pos);

     // Evaluate I1
     complex<double> I11 = (one - u1pos*invsqrt)*expk - I*k1*I0;
     
     // Evaluate 3*I2
     complex<double> I21 = ((two + I*k1*u1pos)*(one - u1pos*invsqrt)
          - u1pos*pow(invsqrt,3.0))*expk - I*k1*I0 + J0*k1*k1;

     I1 = complex<double>(2.0*I10.real() - I11.real(), I11.imag());
     I2 = complex<double>(2.0*I20.real() - I21.real(), I21.imag());

     // Recompute expk
     expk = exp(-I*k1*u1);
  } else {
     // Compute the temporary variable values that will be used below
     expk = exp(-I*k1*u1);

     // Evaluate the approximate integrals I0 and J0
     approxKernelIntegrals(I0, J0, u1, k1);

     // Evaluate I1
     I1 = (one - u1*invsqrt)*expk - I*k1*I0;
     
     // Evaluate 3*I2
     I2 = ((two + I*k1*u1)*(one - u1*invsqrt)
          - u1*pow(invsqrt,3.0))*expk - I*k1*I0 + J0*k1*k1;
	}

  // Compute the first component of the kernel function
  Kf1 = I1 + M*r1*invR*invsqrt*expk;
  
  // Compute the second component of the kernel function
  Kf2 = -I2 - I*k1*invsqrt*expk*pow(M*r1*invR, 2.0)
       - (M*r1*invR)*((one + u1*u1)*pow(beta*r1*invR, 2.0)
       + two + M*r1*u1*invR)*expk*pow(invsqrt,3.0);
}

void
evalKernelNumerator(complex<double>& Kf1, complex<double>& Kf2, double omega,
	double U, double beta, double M, double x0, double r1, double R, double T1,
	double T2, bool steadykernel, double epstol) {
  // Evaluate the two components of the kernel function which are
  // required for the evaluation of the influence coeffficients. These
  // coefficients are the difference between the oscillator and
  // zero-frequency components of the influence coefficients.
  // 
  // Input:
  // omega:      the frequency of oscillation
  // U:          the free-stream velocity
  // beta:       sqrt(1 - M**2)
  // M :         the free-stream Mach number
  // x0, y0, z0: the distances from the current panel location
  // coss, sins: the cos/sin of the dihedral angle of the sending panel
  // cosr, sinr: the cos/sin of the dihedral angle of the receiving panel
  //
  // Output
  // Kf1:        the influence

  double u1;

  // Constants used in this function
  //!! real(kind=dtype), parameter :: zero = 0.0_dtype
  //!! real(kind=dtype), parameter :: one = 1.0_dtype
  //!! real(kind=dtype), parameter :: two = 2.0_dtype
  double zero{0.0}, one{1.0}, two{2.0};
  complex<double> I{0.0, 1.0};

  // Compute the k1 and u1 coefficients used elsewhere
  double k1{omega*r1/U};
  if (r1 <= epstol)
     u1 = (M*R - x0)/(epstol*beta*beta);
  else
     u1 = (M*R - x0)/(r1*beta*beta);

  evalK1K2Coeff(Kf1, Kf2, r1, u1, k1, beta, R, M);

  // Compute the zero-frequency contributions from the coefficients
  double Kf10{zero};
  double Kf20{zero};
  if (steadykernel) {
     Kf10 = one + x0/R;
     Kf20 = -two - (x0/R)*(two + pow((beta*r1/R),2.0));
	}

  // Complete the values of the kernel function
  complex<double> expk = exp(-I*omega*x0/U);
  Kf1 = (Kf1*expk - Kf10)*T1;
  Kf2 = (Kf2*expk - Kf20)*T2;
}

void
computeQuadDoubletCoeff(complex<double>& dinf, double omega, double U,
	double beta, double M, double dxav, const double* xr,
	const double* xi, const double* xo, double e, double cosr, double sinr,
	double coss, double sins, bool steadykernel, double epstol) {
// Evaluate the influence coefficient between a sending panel and a
// recieving point using a quadratic approximation across a panel.
  double x0, y0, z0;
  double eta, zeta, F;

  // The influence coefficients for the different terms
  complex<double> dinf0{0.0}, dinf1{0.0}, dinf2{0.0};

  // The kernel functions evaluate at the different points
  complex<double> Ki1, Ki2, Km1, Km2, Ko1, Ko2;
  double r1, R;
  complex<double> A1, B1, C1, A2, B2, C2, alpha;

  // The coefficients for the horseshoe vortex computation
  double vy, vz, anrm, bnrm, ainv, binv;
  vector<double> a(3,0.0);
  vector<double> b(3,0.0);
  double fact;

  // Set a constant for later useage
  double zero{0.0}, half{0.5}, one{1.0};

	double PI{flaps::pi};
  fact = dxav/(8.0*PI);

  if (omega > 0.0) {
     // T1 = cos(gr - gs)
     double T1 = cosr*coss + sinr*sins;

     // Compute the kernel function at the inboard point
     x0 = xr[0] - xi[0];
     y0 = xr[1] - xi[1];
     z0 = xr[2] - xi[2];
     double T2 = (z0*coss - y0*sins)*(z0*cosr - y0*sinr);

     // Conmpute the distances
     r1 = sqrt(y0*y0 + z0*z0);
     R = sqrt(x0*x0 + beta*beta*(y0*y0 + z0*z0));
     evalKernelNumerator(Ki1, Ki2, omega, U, beta, M,
          x0, r1, R, T1, T2, steadykernel, epstol);
     
     // Evaluate the kernel function at the outboard point
     x0 = xr[0] - xo[0];
     y0 = xr[1] - xo[1];
     z0 = xr[2] - xo[2];
     T2 = (z0*coss - y0*sins)*(z0*cosr - y0*sinr);

     // Conmpute the distances
     r1 = sqrt(y0*y0 + z0*z0);
     R = sqrt(x0*x0 + beta*beta*(y0*y0 + z0*z0));
     evalKernelNumerator(Ko1, Ko2, omega, U, beta, M,
          x0, r1, R, T1, T2, steadykernel, epstol);
     
     // Evaluate the kennel function at the mid-point
     x0 = xr[0] - half*(xi[0] + xo[0]);
     y0 = xr[1] - half*(xi[1] + xo[1]);
     z0 = xr[2] - half*(xi[2] + xo[2]);
     T2 = (z0*coss - y0*sins)*(z0*cosr - y0*sinr);

     // Conmpute the distances
     r1 = sqrt(y0*y0 + z0*z0);
     R = sqrt(x0*x0 + beta*beta*(y0*y0 + z0*z0));
     evalKernelNumerator(Km1, Km2, omega, U, beta, M,
          x0, r1, R, T1, T2, steadykernel, epstol);

     // Compute the A, B and C coefficients for the first term
     A1 = (Ki1 - 2.0*Km1 + Ko1)/(2.0*e*e);
     B1 = (Ko1 - Ki1)/(2.0*e);
     C1 = Km1;

     // Compute the A, B and C coefficients for the second term
     A2 = (Ki2 - 2.0*Km2 + Ko2)/(2.0*e*e);
     B2 = (Ko2 - Ki2)/(2.0*e);
     C2 = Km2;

     // Compute horizontal and vertical distances from the origin in
     // the local ref. frame
     eta = y0*coss + z0*sins;
     zeta = -y0*sins + z0*coss;
     
     // First compute the F-integral
     if (abs(zeta) < epstol*e)
        F = 2*e/(eta*eta - e*e);
     else 
        F = atan(2*e*abs(zeta)/(eta*eta + zeta*zeta - e*e))/abs(zeta);

     // Compute the contribution from the integral of 
     // (A1*y**2 + B1*y + C1)/((eta - y)**2 + zeta**2)
     dinf1 = (((eta*eta - zeta*zeta)*A1 + eta*B1 + C1)*F
          + (0.5*B1 + eta*A1)*log((pow((eta - e),2.0) + zeta*zeta)/
          (pow((eta + e),2.0) + zeta*zeta)) + 2.0*e*A1);

     if (abs(zeta) < epstol*e)
        dinf2 = zero;
     else {
        // Compute the contribution from the integral of 
        // (A2*y**2 + B2*y + C2)/((eta - y)**2 + zeta**2)**2
        alpha = pow(e/zeta, 2.0)*(one - (eta*eta + zeta*zeta - e*e)/(2*e)*F);

        dinf2 = e/(eta*eta + zeta*zeta - e*e)*((
             (2.0*(eta*eta + zeta*zeta + e*e)*(e*e*A2 + C2) + 4.0*eta*e*e*B2))/
             ((pow(eta + e,2.0) + zeta*zeta)*(pow(eta - e,2.0) + zeta*zeta))
             - (alpha/(e*e))*((eta*eta + zeta*zeta)*A2 + eta*B2 + C2));
     }
	}

  if (steadykernel) {
     // Compute the term dinf0 from a horseshoe vortex method. First add
     // the contribution from the inboard and outboard vorticies
     a[0] = (xr[0] - xi[0])/beta;
     a[1] = (xr[1] - xi[1]);
     a[2] = (xr[2] - xi[2]);
     anrm = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
     ainv = one/(anrm*(anrm - a[0]));
     
     b[0] = (xr[0] - xo[0])/beta;
     b[1] = (xr[1] - xo[1]);
     b[2] = (xr[2] - xo[2]);
     //!! bnrm = sqrt(b(0)*b(0) + b(1)*b(1) + b(2)*b(2));
	  bnrm = blas::snrm2(3, &b[0], 1);
     binv = one/(bnrm*(bnrm - b[0]));
     
     vy =  a[2]*ainv - b[2]*binv;
     vz = -a[1]*ainv + b[1]*binv;
     
     // Now, add the contribution from the bound vortex
     ainv = one/(anrm*bnrm*(anrm*bnrm + a[0]*b[0] + a[1]*b[1] + a[2]*b[2]));
     vy = vy + (a[2]*b[0] - a[0]*b[2])*(anrm + bnrm)*ainv;
     vz = vz + (a[0]*b[1] - a[1]*b[0])*(anrm + bnrm)*ainv;
     
     // Compute the steady normalwash
     dinf0 = -(sinr*vy - cosr*vz);
	}

  // Add up all the contributions to the doublet
  dinf = fact*(dinf0 + dinf1 + dinf2);
  
}

void
computeInputMeshSegment(int n, int m, const vector<double>& x0, double span,
	double dihedral, double sweep, double cr, double tr, vector<double>& Xi,
	vector<double>& Xo, vector<double>& Xr, vector<double>& dXav) {
// This routine computes parts of the input mesh for a given lifting
// segment. The input consists of the root location, span, dihedral,
// sweep, root chord and taper ratio. This function can be called
// repeatedly to construct a model for a wing. Note that the sweep is
// the 1/4-chord sweep and the wing is always constructed such that
// it is parallel with the x-axis as required by the DLM theory.
// 
// Input:
// n:        the number of span-wise panels
// m:        the number of chord-wise panels
// span:     the semi-span of the segment
// dihedral: the wing dihedral
// sweep:    the wing sweep
// cr:       the root chord
// tr:       the taper ratio (tip chord = tr*cr)
//
// Output:
// Xi:   the inboard sending point (1/4 box chord in board)
// Xo:   the outboard sending point (1/4 box chord outboard)
// Xr:   the receiving point (3/4 box chord)
// dXav: the average panel length in the x-direction

  // integer, intent(in) :: n, m
  // integer :: i, j, counter
  // real(kind=dtype), intent(in) :: x0(3), span, dihedral, sweep, cr, tr
  // real(kind=dtype), intent(inout) :: Xi(3,n*m), Xo(3,n*m), Xr(3,n*m), dXav(n*m)
  // real(kind=dtype) :: yi, yo, yr, c

  int counter{0}; // 0b
  // do i = 1, n
  for (int i=1; i<=n; i++) {
     // do j = 1, m
	  for (int j=1; j<=m; j++) {
        // Compute the inboard doublet location at the 1/4 chord Note
        // that yp is the span-wise station and c is the chord position
        // relative to the 1/4 chord location of this segment. The
        // tan(sweep) takes care of the 1/4 chord sweep.
        double yi = (i-1.0)*span/n;
        double c = cr*(1.0 - (1.0 - tr)*yi/span)*((j - 0.75)/m - 0.25);
        Xi[IJ(0, counter,3)] = x0[0] + yi*tan(sweep) + c;
        Xi[IJ(1, counter,3)] = x0[1] + yi;
        Xi[IJ(2, counter,3)] = x0[2] + yi*tan(dihedral);

        // Compute the outboard doublet location at the 1/4 chord
        double yo = i*span/n;
        c = cr*(1.0 - (1.0 - tr)*yo/span)*((j - 0.75)/m - 0.25);
        Xo[IJ(0, counter,3)] = x0[0] + yo*tan(sweep) + c;
        Xo[IJ(1, counter,3)] = x0[1] + yo;
        Xo[IJ(2, counter,3)] = x0[2] + yo*tan(dihedral);

        // Compute the receiving point at the 3/4 chord
        double yr = (i-0.5)*span/n;
        c = cr*(1.0 - (1.0 - tr)*yr/span)*((j - 0.25)/m - 0.25);
        Xr[IJ(0, counter,3)] = x0[0] + yr*tan(sweep) + c;
        Xr[IJ(1, counter,3)] = x0[1] + yr;
        Xr[IJ(2, counter,3)] = x0[2] + yr*tan(dihedral);

        // Compute the average chord length of this panel
        dXav[counter] = 0.5*(cr/m)*(2.0 - (1.0 - tr)*(yi + yo)/span);

        // Update the counter location
        counter = counter + 1;
     }
  }
}

void
computeSurfaceSegment(int n, int m, const vector<double>& x0, double span,
	double dihedral, double sweep, double cr, double tr,
	vector<double>& X, vector<int>& conn) {
// This routine computes the surface points for a surface mesh
// corresponding to the given lifting segment. This routine computes
// all quadrilateral surface point locations and adds them to the
// vector X.
//
// Input:
// n:        the number of span-wise panels
// m:        the number of chord-wise panels
// span:     the semi-span of the segment
// dihedral: the wing dihedral
// sweep:    the wing sweep
// cr:       the root chord
// tr:       the taper ratio (tip chord = tr*cr)
//
// Output:
// X         the surface locations
  
  // Input/output specifications
  // integer, intent(in) :: n, m
  // integer, intent(inout) :: conn(4,n*m)
  // real(kind=dtype), intent(in) :: x0(3), span, dihedral, sweep, cr, tr
  // real(kind=dtype), intent(inout) :: X(3,(n+1)*(m+1))

  // Temporary variables
  // integer :: i, j, counter
  // real(kind=dtype) :: y, c

	// First, set the x, y, z locations in X:(3,m+1,n+1)
	for (int j=0; j<n+1; j++)
		for (int i=0; i<m+1; i++) {
			double y{j*span/n};
			double c{cr*(1.0 - (1.0 - tr)*y/span)*(double(i)/m - 0.25)};
			X[IJK(0,i,j,3,m+1)] = x0[0] + y*tan(sweep) + c;
			X[IJK(1,i,j,3,m+1)] = x0[1] + y;
			X[IJK(2,i,j,3,m+1)] = x0[2] + y*tan(dihedral);
		}

  // Now, set connectivity
  // Note: conn is an array of 4 0b node numbers for each panel where
  // the numbering is a 2D array increasing chordwise (x-direction),
  // then spanwise; the numbering for a panel goes clockwise. Each
  // conn is the index of a column of the (3,(n+1)*(m+1)) coordinate array X
  // do i = 1, n
     // do j = 1, m
  for (int i=0; i<n; i++) {
	  for (int j=0; j<m; j++) {
        // Compute the x/y/z locations of the connectivity mesh for
        // this lifting segment
        int panel = j + i*m;  // 0b
        conn[IJ(0, panel,4)] = j + i*(m+1);
        conn[IJ(1, panel,4)] = j + 1 + i*(m+1);
        conn[IJ(2, panel,4)] = j + 1 + (i+1)*(m+1);
        conn[IJ(3, panel,4)] = j + (i+1)*(m+1);
     }
  }
}

void
computeInfluenceMatrix(vector<complex<double>>& D, double omega, double U, double M,
	int np, const vector<double>& Xi, const vector<double>& Xo,
	const vector<double>& Xr, const vector<double>& dXav,
	int symmetric, bool steadykernel, double epstol) {
// This routine computes the complex influence coefficient
// matrix. The input consists of a number of post-processed
// connectivity and nodal locations are given and locations,
// determine
// 
// Input:
// omega: the frequency of oscillation
// U:     the velocity of the free-stream
// M:     the free-stream Mach number
// np:    number of panels
// Xi:    inboad sending point
// Xo:    outboard sending point
// Xr:    receiving point
// dXav:  average length in the x-direction of the panel
//
// Output:
// D:  complex coefficient matrix

  // Input/output types
  // logical, intent(in) :: steadykernel
  // integer, intent(in) :: np, symmetric
  // complex(kind=dtype), intent(inout) :: D(np, np)
  // real(kind=dtype), intent(in) :: omega, U, M, epstol
  // real(kind=dtype), intent(in) :: Xi(3,np), Xo(3,np), Xr(3,np), dXav(np)

  // Temporary data used internally
  //!! integer :: r, s
  //!! real(kind=dtype) :: beta, xrsymm(3), sinsymm
  //!! real(kind=dtype) :: pe(np), pcos(np), psin(np)
  //!! complex(kind=dtype) :: dtmp

  // Compute the compressibility factor
  double beta{sqrt(1.0 - M*M)};

  // Pre-processing step: Compute the sin/cos and length of all the
  // panels in the model
  vector<double> pe(np, 0.0);
  vector<double> pcos(np, 0.0);
  vector<double> psin(np, 0.0);
  //!! do r = 1, np
  for (int r=0; r<np; r++) {
     // Compute 1/2 the bound vortex length
     //!! pe(r) = 0.5*sqrt((Xo(2,r) - Xi(2,r))**2 + (Xo(3,r) - Xi(3,r))**2)
     double y = Xo[IJ(1,r,3)] - Xi[IJ(1,r,3)];
     double z = Xo[IJ(2,r,3)] - Xi[IJ(2,r,3)];
		pe[r] = 0.5*sqrt(y*y + z*z);
     // Compute the sin and cos of the dihedral
     pcos[r] = 0.5*y/pe[r];
     psin[r] = 0.5*z/pe[r];
  }

  if (symmetric == 0) {
     //!! do s = 1, np
	  for (int s = 0; s<np; s++) {
        //!! do r = 1, np
		  for (int r = 0; r<np; r++) {
           // Compute the panel influence coefficient
           //!! call computeQuadDoubletCoeff(D(r, s), omega, U, beta, M,
               //!!  dXav(s), Xr(:, r), Xi(:, s), Xo(:, s), pe(s), &
                //!! pcos(r), psin(r), pcos(s), psin(s), steadykernel, epstol)
           computeQuadDoubletCoeff(D[IJ(r,s,np)], omega, U, beta, M,
               dXav[s], &Xr[IJ(0,r,3)], &Xi[IJ(0,s,3)], &Xo[IJ(0,s,3)], pe[s],
               pcos[r], psin[r], pcos[s], psin[s], steadykernel, epstol);
        }
     }
  } else {
     //!! do s = 1, np
	  for (int s=0; s<np; s++) {
        //!! do r = 1, np
		  for (int r=0; r<np; r++) {
           // Compute the panel influence coefficient
           computeQuadDoubletCoeff(D[IJ(r,s,np)], omega, U, beta, M,
                dXav[s], &Xr[IJ(0,r,3)], &Xi[IJ(0,s,3)], &Xo[IJ(0,s,3)],
					 pe[s], pcos[r], psin[r], pcos[s], psin[s], steadykernel, epstol);

           // Compute the influence from the same panel, but the reflected point
           //!! xrsymm(1) =  Xr(1, r)
           //!! xrsymm(2) = -Xr(2, r)
           //!! xrsymm(3) =  Xr(3, r)
			  vector<double> xrsymm{Xr[IJ(0,r,3)], -Xr[IJ(1,r,3)], Xr[IJ(2,r,3)]};
           double sinsymm = -psin[r];

				complex<double> dtmp(0.0);
           computeQuadDoubletCoeff(dtmp, omega, U, beta, M,
                dXav[s], &xrsymm[0], &Xi[IJ(0,s,3)], &Xo[IJ(0,s,3)], pe[s],
                pcos[r], sinsymm, pcos[s], psin[s], steadykernel, epstol);
           D[IJ(r,s,np)] += dtmp;
        }
     }
  }
}

void
addCpForces(int np, int n, double qinf, complex<double>* Cp, double* X,
	const int* conn, complex<double>* forces) {
// Given the coefficient of pressure, compute the forces at each
// node. This distributes the force to each of the corresponding
// nodes. This function can be used to compute the forces that act on
// a finite-element mesh.
//
// Input:
// np:     the number of panels
// n:      the number of nodes
// qinf:   the dynamic pressure
// Cp:     the coefficient of pressure on each panel
// X:      the nodal locations
// conn:   the connectivity
//
// Output:
// forces: the force values at each node

  // Input/output types
  //!! integer, intent(in) :: np, n, conn(4,np)
  //!! complex(kind=dtype), intent(in) :: Cp(np)
  //!! real(kind=dtype), intent(in) :: qinf, X(3,n)
  //!! complex(kind=dtype), intent(inout) :: forces(3,n)

  // Temporary data used internally
  //!! integer :: i, j
  //!! real(kind=dtype) :: a(3), b(3), normal(3)
	vector<double> a(3,0.0);
	vector<double> b(3,0.0);
	vector<double> normal(3,0.0);
  //!! do i = 1, np
  for (int i=0; i<np; i++) {
     // Compute the a/b vectors
	  // (eem) note conn is 0b but the f90 code uses 1b indices,
	  //       so the +1 was deleted here
	  const int* co = &conn[i*4];
     //!! do j = 1, 3
	  for (int j=0; j<3; j++) {
	     // a(j) = X(j,conn(3,i)+1) - X(j,conn(1,i)+1)
        a[j] = X[IJ(j,co[2],3)]
		       - X[IJ(j,co[0],3)];
        b[j] = X[IJ(j,co[3],3)]
		       - X[IJ(j,co[1],3)];
     }

     // Compute n = a x b. Note that the normal is not normalized since
     // the area of the cell is equal to 1/2 the norm of a x b. This
     // accounts for the factor of .125 on each cell vertex (1/2*1/4).
     normal[0] = a[1]*b[2] - a[2]*b[1];
     normal[1] = a[2]*b[0] - a[0]*b[2];
     normal[2] = a[0]*b[1] - a[1]*b[0];

     // Add the forces
     //!! do j = 1, 4
	  for (int j=0; j<4; j++) {
        forces[IJ(0,co[j],3)] += 0.125*qinf*Cp[i]*normal[0];
        forces[IJ(1,co[j],3)] += 0.125*qinf*Cp[i]*normal[1];
        forces[IJ(2,co[j],3)] += 0.125*qinf*Cp[i]*normal[2];
     }
  }
}

void
nodal_forces(int np, int n, int nc, const complex<double>* Cp,
	double* X, const int* conn, complex<double>* forces) {
// Given the coefficient of pressure, compute the forces at each
// node. This distributes the force to each of the corresponding
// nodes. This function can be used to compute the forces that act on
// a finite-element mesh.
//
// Input:
// np:     the number of panels
// n:      the number of nodes
// nc      number of modes (columns of Cp)
// Cp:     (np,nc) complex coefficient of pressure on each panel
// X:      (3,n) the nodal locations
// conn:   (4,np) the connectivity: index of each of the 4 nodes on
//         each panel, clockwise from the lower (x,y) node
//
// Output:
// forces: (3,n,nc) complex force values at each node for each mode

	vector<double> a(3,0.0);
	vector<double> b(3,0.0);
	vector<double> normal(3*np,0.0);

	// compute the normal vector for each panel as the cross-product
	// of the 2 diagonals
  for (int j=0; j<np; j++) {
     // Compute the a/b vectors
	  // (eem) note conn is 0b but the f90 code uses 1b indices,
	  //       so the +1 was deleted here
	  const int* co = &conn[j*4];
	  for (int i=0; i<3; i++) {
	     // a(i) = X(i,conn(3,j)+1) - X(i,conn(1,j)+1)
        a[i] = X[IJ(i,co[2],3)]
		       - X[IJ(i,co[0],3)];
        b[i] = X[IJ(i,co[3],3)]
		       - X[IJ(i,co[1],3)];
     }

     // Compute n = a x b. Note that the normal is not normalized since
     // the area of the cell is equal to 1/2 the norm of a x b. This
     // accounts for the factor of .125 on each cell vertex (1/2*1/4).
		double* no = &normal[3*j];
		no[0] = a[1]*b[2] - a[2]*b[1];
		no[1] = a[2]*b[0] - a[0]*b[2];
		no[2] = a[0]*b[1] - a[1]*b[0];
	}

	// for each mode...
	for (int k=0; k<nc; k++) {
		// ...Add the forces XXX dlm4py *adds* A but we subtract: -qA so subtract
		for (int j=0; j<np; j++) {
			const int* co = &conn[j*4];
			double* no = &normal[3*j];
			complex<double> cp = Cp[IJ(j,k,np)];
			// ... for each node on this panel
			for (int i=0; i<4; i++) {
				complex<double>* f = &forces[IJK(0,co[i],k,3,n)];
				f[0] -= 0.125*cp*no[0];
				f[1] -= 0.125*cp*no[1];
				f[2] -= 0.125*cp*no[2];
			}
		}
	}
}

void
getModeBCs(int np, int n, const double* mode, const double* X,
	const int* conn, double* vwash, double* dwash) {
// Given the mode shape at the nodes of the DLM mesh, compute the
// contributions to the boundary conditions at the receiving
// points. Note that there are contributions from both the
// displacement and the normal component of the velocity.
// 
//
// Input:
// np:     the number of panels
// n:      the number of nodes
// mode:   the mode shape
// X:      the nodal locations
// conn:   the connectivity
//
// Output:
// vwash:  the normal wash due to velocities
// dwash:  the normal wash due to the change in orientation

  // Input/output types
  //!! integer, intent(in) :: np, n, conn(4,np)
  //!! real(kind=dtype), intent(in) :: X(3,n)
  //!! real(kind=dtype), intent(in) :: mode(3,n)
  //!! real(kind=dtype), intent(inout) :: vwash(np), dwash(np)

  // Temporary data used internally
  //!! integer :: i, j
  //!! real(kind=dtype) :: a(3), b(3), normal(3), nrm, dx
  //!! real(kind=dtype) :: up(3)


  //!! do i = 1, np
  for (int i=0; i<np; i++) {
     // Zero-out the components of the displacement (or velocity) at
     // the center point.
	  vector<double> a(3,0.0), b(3,0.0), normal(3,0.0);
	  vector<double> up(3, 0.0);

     // Compute the values of the displacement at the receiving point
	  // (eem) here again we delete the +1 to convert to 0b indexing
	  const int* co = &conn[i*4];
     //!! do j = 1, 3
	  for (int j=0; j<3; j++) {
        // up(j) = 0.5*( &
          //    0.25*mode(j, conn(1,i)+1) + 0.75*mode(j, conn(2,i)+1) + &
            //  0.25*mode(j, conn(4,i)+1) + 0.75*mode(j, conn(3,i)+1))
        up[j] = 0.5*(0.25*mode[IJ(j,co[0],3)] +
				         0.75*mode[IJ(j,co[1],3)] +
                     0.25*mode[IJ(j,co[3],3)] +
				         0.75*mode[IJ(j,co[2],3)]);
     }

     // Compute the a/b vectors
     //!! do j = 1, 3
	  for (int j=0; j<3; j++) {
        a[j] = X[IJ(j,co[2],3)] -
		         X[IJ(j,co[0],3)];
        b[j] = X[IJ(j,co[3],3)] -
		         X[IJ(j,co[1],3)];
     }

     // Compute n = a x b
     normal[0] = a[1]*b[2] - a[2]*b[1];
     normal[1] = a[2]*b[0] - a[0]*b[2];
     normal[2] = a[0]*b[1] - a[1]*b[0];

     // Compute the area in order to normalize the normal vector
	  // (eem) the area is 1/2 the norm - doesn't matter here
     //!!  nrm = sqrt(normal(1)**2 + normal(2)**2 + normal(3)**2)
	  double nrm = blas::snrm2(3, &normal[0], 1);

     // Set the normal displacement/velocity
     //!! vwash(i) = -(normal(1)*up(1) + normal(2)*up(2) + normal(3)*up(3))/nrm
	  blas::dot(3, &normal[0], 1, &up[0], 1, vwash[i]);
	  vwash[i] /= -nrm;

     // Compute d(mode)/dx
     // First, compute the average step in x
     double dx = 0.5*(
          (X[IJ(0, co[1],3)] -
			  X[IJ(0, co[0],3)]) +
          (X[IJ(0, co[2],3)] -
			  X[IJ(0, co[3],3)]));

     // Compute the derivative d(mode)/dx
     dwash[i] = -0.5*(
          (mode[IJ(2, co[1],3)] -
			  mode[IJ(2, co[0],3)]) +
          (mode[IJ(2, co[2],3)] -
			  mode[IJ(2, co[3],3)]))/dx;
  }
}

void
washes(int np, int n, int nc, const double* mode, const double* X,
	const int* conn, double* vwash, double* dwash) {
// Given the mode shape at the nodes of the DLM mesh, compute the
// contributions to the boundary conditions at the receiving
// points. Note that there are contributions from both the
// displacement and the normal component of the velocity.
// 
//
// Input:
// np:     the number of panels
// n:      the number of nodes = (m+1)*(n+1)
// nc:     the number of modes
// mode:   (3,n,nc) modes matrix
// X:      (3,n) the nodal locations
// conn:   (4,np) the connectivity
//
// Output:
// vwash:  (np,nc) the normal wash due to velocities
// dwash:  (np,nc) the normal wash due to the change in orientation


  // Temporary data used internally
  //!! integer :: i, j
  //!! real(kind=dtype) :: a(3), b(3), normal(3), nrm, dx
  //!! real(kind=dtype) :: up(3)

	T_(Trace trc(1,"washes");)

	// for each mode...
	for (int k=0; k<nc; k++) {
		// ... and each panel
		for (int j=0; j<np; j++) {
			// co is a vector of 4 ints: clockwise nodes for panel j
			const int* co = &conn[j*4];
			vector<double> a(3,0.0), b(3,0.0), normal(3,0.0);
			vector<double> up(3, 0.0);

			// Compute the values of the displacement at the receiving point
			// (eem) here again we delete the +1 to convert to 0b indexing
			for (int i=0; i<3; i++)
				up[i] = 0.5*(0.25*mode[IJK(i,co[0],k,3,n)] +
				             0.75*mode[IJK(i,co[1],k,3,n)] +
                         0.25*mode[IJK(i,co[3],k,3,n)] +
				             0.75*mode[IJK(i,co[2],k,3,n)]);
			T_(trc.dprint("up: ",up[0],", ",up[1],", ",up[2]);)

			// Compute the a/b vectors
			for (int i=0; i<3; i++) {
				a[i] = X[IJ(i,co[2],3)] - X[IJ(i,co[0],3)];
				b[i] = X[IJ(i,co[3],3)] - X[IJ(i,co[1],3)];
			}

			// Compute n = a x b
			normal[0] = a[1]*b[2] - a[2]*b[1];
			normal[1] = a[2]*b[0] - a[0]*b[2];
			normal[2] = a[0]*b[1] - a[1]*b[0];

			// Compute the area in order to normalize the normal vector
			// (eem) the area is 1/2 the norm - doesn't matter here
			double nrm = blas::snrm2(3, &normal[0], 1);
			T_(trc.dprint("nrm = ",nrm);)

			// Set the normal displacement/velocity
			vwash[IJ(j,k,np)] = -(normal[0]*up[0]+normal[1]*up[1]+normal[2]*up[2])/nrm;

			// Compute d(mode)/dx
			// First, compute the average step in x
			double dx = 0.5*((X[IJ(0,co[1],3)] - X[IJ(0,co[0],3)]) +
								  (X[IJ(0,co[2],3)] - X[IJ(0,co[3],3)]));
			T_(trc.dprint("dx = ",dx);)
			// Compute the derivative d(mode)/dx
			dwash[IJ(j,k,np)] =
				-((mode[IJK(2,co[1],k,3,n)] - mode[IJK(2,co[0],k,3,n)]) +
				  (mode[IJK(2,co[2],k,3,n)] - mode[IJK(2,co[3],k,3,n)]))/(2.0*dx);
		}
	}
	T_(trc.dprintm(np,nc,np,vwash,"vwash");)
	T_(trc.dprintm(np,nc,np,dwash,"dwash");)
}

void
Dlm::
vzaero(vector<int>& node_numbers, vector<double>& coords,
   vector<int>& segments, vector<double>& gct_nodal, int nc) {
// Add to vectors of node_numbers, coords, segments, and gct the
// aero grid
// conn   (4,n*m) node numbers for each (n*m) panel, clockwise
//        starting from the (upstream,inboard) corner
// X      (3,(n+1)*(m+1)) coordinates of each of the (m+1)*(n+1) nodes
// T      (3,(n+1),(m+1),nc) transformation from the Fem grid to the aero
//        grid, where nc is the number of Fem freedoms
	T_(Trace trc(1,"vzaero");)
	int m{panels[0]};		// # chordwise panels
	int n{panels[1]};		// # spanwise panels
	int np = n*m;
	int nnodes = (m+1)*(n+1);

	assert(3*nnodes*nc == (int)T.size());

	// find the largest node number, number the aero grid nodes
	// starting with maxnode+1
	int start{0};
	for (auto i : node_numbers)
		start = std::max(start, i);
	start++;
		
	// new nodes are numbered start to start+nnodes-1
	// X: (3,m+1,n+1)
	for (int j=0; j<nnodes; j++) {
		node_numbers.push_back(start+j);
		for (int i=0; i<3; i++)
			coords.push_back(X[IJ(i,j,3)]);
	}

	// add the 4 sides of each panel to "segments"
	// conn is 0b
	for (int j=0; j<np; j++) {
		const int* co = &conn[4*j];
		segments.push_back(co[0]+start);
		segments.push_back(co[1]+start);
		segments.push_back(co[2]+start);
		segments.push_back(co[3]+start);
		segments.push_back(co[0]+start);
		segments.push_back(0);
	}

	// append our T matrix to gct_nodal: treat it as a (3*nnodes,nc) array
	int oldnr = gct_nodal.size()/nc;
	int nrx = 3*nnodes;
	int newnr = oldnr + nrx;
	flaps::info("dlm grid nodes are rows ",oldnr+1," to ",newnr);
	if (gct_nodal.empty()) {
		gct_nodal = vector<double>(T.size(), 0.0);
		blas::copy(T.size(), &T[0], 1, (double*)&gct_nodal[0], 2);
	} else {
		vector<double> tmp(newnr*nc);
		for (int j=0; j<nc; j++) {
			blas::copy(oldnr, &gct_nodal[j*oldnr], 1, &tmp[j*newnr], 1);
			// T is double so skip factor for tmp is 2
			blas::copy(nrx, &T[j*nrx], 1, &tmp[IJ(oldnr,j,newnr)], 1);
		}
		gct_nodal = tmp;
		T_(trc.dprint("gct_nodal is now (",newnr,",",nc,")");)
	}
	T_(trc.dprint("coords are now ",coords.size());)
	T_(int nr = gct_nodal.size()/nc;)
	T_(trc.dprintm(nr,nc,nr,&gct_nodal[0],"vzaero new gct_nodal");)
}
