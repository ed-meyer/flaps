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
#include "Ad.h"
#include "matrix.h"
#include "Pz.h"

// elements of the E matrix
class Egc {
public:
	int gc;
	double factor;

	Egc(int g, double f) : gc(g), factor(f) {}
};

// Internal time delays
class Itd {
public:
	std::vector<Elem> elem;
	double deltat;

	Itd() : deltat(0.0) {}
	Itd(double dt) : deltat(dt) {}
};


class LTIPz : public Pz {
	void gainValues(std::vector<std::complex<double>>& igainval, std::vector<std::complex<double>>& ogainval);
	double interpA(double x, int eleno, int ider) const;
	void scaleX(int ne, int ns, std::vector<std::complex<double>>& a);

public:
	int ns;   // size of A
	int ni;   // number of inputs (col B, D)
	int no;   // number of outputs (rows C, D)
	// int nb;   // number of breakpoints in A interpolant
	int nix;  // number of elements of A that are variable
	std::vector<double> xb;    // nb paramA breakpoints for interpolating A
	std::vector<Elem> elem;		// nix interpolated A elements (zero-based)
	// int* row;
	// int* col;
	std::vector<double> deltati;  	// ni time delays
	std::vector<std::string> igain;	// input gains parameter names
	std::vector<std::string> iphase;	// input phase parameter names
	std::vector<double> deltato;  	// no time delays
	//!! std::vector<Par*> ogain;
	//!! std::vector<Par*> ophase;
	std::vector<std::string> ogain;
	std::vector<std::string> ophase;
	std::vector<int> sinput;	// (ni)
	std::vector<double> A;		// (ns,ns)
	std::vector<double> B;		// (ns,ni)
	std::vector<double> C;		// (no,ns)
	std::vector<double> D;		// (no,ni)
	std::vector<double> coef;	// (nb-1, 4, nix)
	std::vector<double> psi;	// (ni,ne)
	//!! Matrix* E{nullptr};			// (ne,no)
	Matrix* KE{nullptr}; 		// (ne,no) XXX custom Pz replaces makeKE
	std::vector<Egc> egc;      // no
	std::vector<Egc> sgc;      // ns
	std::vector<Itd> internalTimeDelays;
	std::string paramA;   // A is a function of this parameter
	std::string paramITD;   // internal time delays get multiplied by this parameter

	// constructors, desctructor
	LTIPz() {
		ns = ni = no = nix = 0;
		xb = deltati = deltato = A = B = C = D = 0;
	}
	// main constructor:
	LTIPz(std::string const& path, std::string const& psiName,
		std::string const& EName, std::string const& KEName,
		std::vector<Egc>& egc, std::vector<Egc>& sgc,
		std::vector<Par*> igain, std::vector<Par*> iphase,
		std::vector<Par*> ogain, std::vector<Par*> ophase, std::string const& nameITD);
	LTIPz(Receiver& s);
	LTIPz(const LTIPz&) noexcept = default;			// copy constructor
	LTIPz& operator=(const LTIPz& rhs) = default;	// assignment op

	~LTIPz() = default;		// dtor

	//!! void check() const;

	// clone or factory method ("virtual constructor")
	LTIPz* clone() const { return new LTIPz(*this); }

	// "Evaluate" the parameterization...
	bool eval (pset& plt, std::vector<std::complex<Ad>>& result, int nr, int nc);

	// evaluate using approximate exponentials for time-delays?
	static bool useApproxExp(bool yn);

	// printing stuff...
	virtual void osprint(std::ostream& s) const;
	//!! std::string desc() const;
	//!! void print(FILE* stream, const char* title) const;
	//!! std::string summary() const;
	//!! std::string table() const;
	// XXX eventually make virtual, use in all P14n:
	//!! std::string display(std::string const& path, std::string const& title) const;

	//!! std::vector<Real> origValues(std::string const& name) const;

	void smoothSplines();

	void plot(std::vector<Elem> const& elements, int nstep);
	void plotElem(int eleno, int nstep, std::string const& runid, bool append);
	void plotEigen();

	// Fio stuff
	static bool regd;		// register_class'd?
	static std::string id() { return std::string("LTIPz"); }
	virtual std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new LTIPz(s); }
	virtual void put (Sender& s) const;

	//!! void RFAABCD (std::vector<std::complex<double>> const& sval,
	//!!	std::vector<Real> beta);
	//!! void approx (std::vector<std::complex<double>> const& sval,
	//!!    std::vector<std::complex<double>>& beta, std::vector<Matrix*>& R);
};  // LTIPz
