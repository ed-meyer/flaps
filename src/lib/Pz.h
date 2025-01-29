//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef Pz_H
#define Pz_H

/*------------------------------------------------------------------
 * Matrix parameterization (Pz) classes
 * A Matrix stores the actual array in the "data()" member
 *
 *                          Pz (abstract base class)
 *                                      |
 *                 ----------------------------------------
 *                 |                                      |
 *            generative                               modifier
 *                 |                                      |
 *     -----------------------------------         ----------------
 *     |       |      |      |           |         |              |
 *  Pz_const  IntPz  RFAPz  CustomPz   LTIPz      LMPz          CustomPz
 *
 * Pz_const:   the matrix has been declared to be a function of a given
 *             set of parameters and the "data()" member is at these parameter
 *             values ("params" member)
 * IntPz       1, 2, or 3-dimensional interpolation using the "spline" member
 *             The "data()" member is not accessed
 * RFAPz       rational-function approximation;
 * LTIPz       Linear Time-Invariant state-space control equations
 * CustomPz    a custom function which either modifies a copy of the "data()"
 *             array or creates the result from scratch (generative)
 * LMPz        an equation that modifies a specific element of a copy of "data()"
 *
 * A parameterized matrix may have one (and only one) generative
 * Pz, however, note that a CustomPz can be either generative or a
 * modifier.
 *------------------------------------------------------------------*/

// Forward declarations
class Elem;
class Pz;
class Pz_const;
class IntPz;
class RFAPz;
class LTIPz;
class LMPz;
class CustomPz;

#include <cassert>
#include <functional>
#include <string>
#include <vector>
#include <algorithm>


#include "custom.h"
#include "interp.h"
#include "matrix.h"
#include "Ad.h"
#include "fio.h"
#include "pset.h"


// an Elem describes a matrix element, i.e. the row & column
class Elem {
public:
	size_t row;   // zero-based  XXX change to 1b?
	size_t col;   // zero-based
	Elem() : row(0), col(0) {}
	Elem (size_t i, size_t j) : row(i), col(j) {}
	Elem (std::string const&);  // parser constructor: "[i,j]"
	// default copy constructor, assignment operator ok
	bool operator==(Elem const& rhs) const { return (rhs.row == row && rhs.col == col); }
	// return the 0-based index of this elem in a 2D array with leading dimension ldim:
	size_t index(size_t ldim) const { return (row + col*ldim); }
	friend std::ostream& operator<<(std::ostream& s, const Elem& t);
};


/*------------------------------------------------------------------
 * Pz is the base class for all *Pz
 *------------------------------------------------------------------*/

class Pz : public Fio {
protected:
	size_t nrow, ncol;          // dimensions of the matrix
public:
	// constructors
	Pz() : nrow(0), ncol(0) {}
	Pz(size_t nr, size_t nc) : nrow(nr), ncol(nc) {}
	Pz(Receiver& s);
	// default copy constructor ok
	// default assignment operator ok
	virtual ~Pz();

	// clone or factory method ("virtual constructor")
	// XXX Pz should be an abstract class?
	virtual Pz* clone() const {return nullptr; }

	// Evaluate the parameterization and replace "result" (generative)
	// or modify "result" (modifier)
	virtual bool eval(pset& plt, std::vector<std::complex<Ad>>& result, size_t nr, size_t nc);
	
	// dependent parameter names: each type of pz is responsible for
	// determining it's dependent parameters
	virtual std::vector<std::string> dependson(pset& plt) {
		return std::vector<std::string>(); }
	std::vector<std::string> get_dependson(pset& plt);

	// access nrow, ncol:
	size_t& rsize() { return nrow; }
	size_t& csize() { return ncol; }

	// printing for all *Pz via this operator<< which calls the
	// virtual functions osprint so that the correct osprint gets called
	virtual void osprint(std::ostream& s) const {}
	friend std::ostream& operator<<(std::ostream& s, const Pz& p);

	// Fio required stuff
	static std::string id() { return std::string("Pz"); }
	virtual std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new Pz(s); }
	virtual void put (Sender& s) const;

}; // Pz


bool
copyAD (double const* from, size_t nrf, size_t ncf,
		std::vector<std::complex<Ad>>& to, size_t nrt, size_t nct);

bool
copyCAD (double const* from, size_t nrf, size_t ncf,
		std::vector<std::complex<Ad>>& to, size_t nrt, size_t nct);

// Pz_const is for a matrix that is a constant function of a
// set of parameters, so it just contains the names of the
// parameters and their values.
class Pz_const : public Pz {
public:
	std::vector<Par*> params;
	// std::vector<std::string> dep_par;  // names of dependent parameters
	// std::vector<double> dep_values;  // values of dependent parameters

	// main constructor:
	Pz_const(const std::vector<Par*>& par, size_t nr, size_t nc);
	Pz_const(Receiver& s);
	Pz_const& operator=(const Pz_const& rhs);
	virtual ~Pz_const() { /* necessary to avoid missing vtable? */ }

	virtual Pz* clone() const { return new Pz_const(*this); }

	// Evaluate the parameterization: do nothing: the Matrix contains the data
	virtual bool eval (pset& plt,std::vector<std::complex<Ad>>& result,size_t nr,size_t nc) {
		return true;
	}

	// printing: see comments in Pz
	virtual void osprint(std::ostream& s) const;
	// get dependent parameter names:
	virtual std::vector<std::string> dependson(pset& plt);

	// Fio/serialization stuff
	static bool regd;
	static std::string id() { return std::string("Pz_const"); }
	virtual std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new Pz_const(s); }
	virtual void put (Sender& s) const;
};  // Pz_const


//------------------------------------------------------------------
// IntPz: 1D Interpolation parameterization
//------------------------------------------------------------------

class IntPz : public Pz {
	std::vector<Ad> get_par_values(pset& plt);
public:
	Interp spline;
	int extrap{0};			// if > 0 => extrapolation ok
	double rho{0.0};       // smoothing factor
	std::string warning;		// ill-conditioned interpolation
	// constructors, desctructor
	IntPz() {}
	// main constructor: feed it a bunch of matrices
	IntPz(const std::vector<Matrix*>& m, int extrap=0, double rhop = 0.0);
	IntPz(Receiver& s);
	IntPz(const IntPz&) = default;						// copy constructor
	// IntPz(std::vector<IntPz*> list);  				// consolidation constructor
	IntPz& operator=(const IntPz& rhs) = default;	// assignment op

	~IntPz() {}

	// clone or factory method ("virtual constructor")
	IntPz* clone() const { return new IntPz(*this); }

	// Evaluate the parameterization
	bool eval (pset& plt, std::vector<std::complex<Ad>>& result, size_t nr, size_t nc);

	// printing: see comments in Pz
	virtual void osprint(std::ostream& s) const;
	// dependent parameter names:
	virtual std::vector<std::string> dependson(pset& plt) {
		return spline.dependson(); }
	// std::string desc() const;
	// void print(FILE* stream, const char* title) const;
	// std::string summary(const pset& plt) const;
	// std::string table(const pset& plt) const;
	// virtual std::string display(std::string const& path, std::string const& title) const;

	void plot(std::vector<int> elements, std::string const& var, int nstep);

	// Fio fio stuff
	static bool regd;        // register_class'd?
	static std::string id() { return std::string("IntPz"); }
	virtual std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new IntPz(s); }
	virtual void put (Sender& s) const;
};  // IntPz

// std::ostream&
// operator<<(std::ostream& s, const IntPz& t);

//------------------------------------------------------------------
// An LMPz (eLeMent parameterization) assigns an equation to
// (row,col) of a matrix
//------------------------------------------------------------------
class LMPz : public Pz {
public:
	std::string eqn;
	Elem elem;   // 0b
	char op;     // [row,col] op= adeqn::eval(eqn)
	std::vector<std::string> dep_par;  // dependent parameter names

	// constructors
	LMPz() {op = 0;}
	LMPz(size_t nr, size_t nc, Elem& el, std::string const& equation, char oper=0);
	LMPz(Receiver& s);
	// default copy constructor, assignment operator ok
	// LMPz(const LMPz&);	// copy constructor
	~LMPz();

	// clone or factory method ("virtual constructor")
	LMPz* clone() const { return new LMPz(*this); }


	// Evaluate the parameterization
	bool eval(pset& plt, std::vector<std::complex<Ad>>& result, size_t nr, size_t nc);

	// printing: see comments in Pz
	virtual void osprint(std::ostream& s) const;
	// dependent parameter names:
	virtual std::vector<std::string> dependson(pset& plt);
	// std::string desc() const;
	// std::string summary(const pset& plt) const;
	// std::string table(const pset& plt) const;
	// virtual std::string display(std::string const& path, std::string const& title) const;

		// Fio stuff
	static bool regd;
	static std::string id() { return std::string("LMPz"); }
	virtual std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new LMPz(s); }
	virtual void put (Sender& s) const;
}; // LMPz

// std::ostream&
// operator<<(std::ostream& s, const LMPz& t);

// RFA aero approximation
// Rodger's approximation:
// \mathbf{A} = \mathbf{R}_0 + p\mathbf{R}_1 + p^2 \mathbf{R}_2 +
// \sum_{i=1}^{m} \frac{p}{p + \beta_i} \mathbf{R}_{i+2}

class RFAPz : public Pz {
public:
	std::vector<double> beta;
	std::vector<std::vector<double>> R;
	double rsfmin, rsfmax;
	double rfmin, rfmax;

	// Constructor(s)
	RFAPz() {}
	RFAPz (std::vector<Matrix*> const& matrices, std::vector<double> betas,
			 bool kscale=false, bool forceR0=false, bool optbeta=false);
	RFAPz (size_t n, std::vector<std::complex<double>*> const& matrices,
			 std::vector<std::complex<double> > pval, std::vector<double> betas,
			 bool kscale, bool forceR0, bool optbeta);
	RFAPz(const RFAPz&) noexcept = default;			// copy constructor
	RFAPz& operator=(const RFAPz& rhs) = default;	// assignment op
	RFAPz (Receiver& s);   // serializing constructor
	~RFAPz() = default;

	// clone or factory method ("virtual constructor")
	RFAPz* clone() const { return new RFAPz(*this); }

	// compute an optimal set of beta's and R's
	static bool optimize (std::vector<std::complex<double>> const& pval,
		std::vector<double> const& beta,
		size_t n, std::vector<std::complex<double>*> const& matrices,
		std::vector<std::vector<double>>& R);

	// Evaluate the parameterization as a complex<Ad> matrix
	bool eval(pset& plt, std::vector<std::complex<Ad>>& result, size_t nr, size_t nc);
	// evaluate the parameterization as a complex<double> matrix
	static bool eval(std::complex<double> const& p, std::vector<double> const& beta,
		std::vector<std::vector<double>> const& R, std::vector<std::complex<double>>& result);

	// printing: see comments in Pz
	virtual void osprint(std::ostream& s) const;
	// dependent parameter names:
	virtual std::vector<std::string> dependson(pset& plt) {
		return std::vector<std::string>{"rsf","rf"};
	}
	// std::string desc() const;
	// std::string summary(const pset&) const;
	// std::string table(const pset&) const;
	// XXX obsolete void print(FILE* stream, const char* title) const;
	// virtual std::string display(std::string const& path, std::string const& title) const;

		// Fio stuff
	static bool regd;
	static std::string id() { return std::string("RFAPz"); }
	virtual std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new RFAPz(s); }
	virtual void put (Sender& s) const;
}; // RFAPz

// std::ostream&
// operator<<(std::ostream& s, const RFAPz& t);

/*------------------------------------------------------------------
 * CustomPz - custom function parameterization
 *
 * Support for custom functions describing or modifying matrices.
 * If a matrix has a user parameterization it contains
 * C++ code to evaluate certain terms of the matrix.
 *------------------------------------------------------------------*/

// this is the one-and-only definition of all user entry points:
using CustomEval = int (*)(pset& plt, int nr, int nc, std::vector<std::complex<Ad>>& a);
//!! using CustomEval = std::function<int(pset& plt, int nr, int nc, std::vector<std::complex<Ad>>& a)>;

class CustomPz : public Pz {
	// dependent parameter names:
	std::vector<std::string> dep_par;
public:
	std::string entrypt;	// name (entry pt) of the eval routine */
	// XXX n+nExtra is the size of the matrix
	size_t nExtra{0};			/* number of extra rows/columns */
	CustomEval function{nullptr};

	// Constructor(s)
	CustomPz() {}
		// XXX rm nExtra, let the caller figure it out?
	CustomPz (std::string const& filename, size_t nr, size_t nc, int nExtra=0);
	CustomPz (CustomEval fcn, size_t nr, size_t nc, int nExtra=0);
	CustomPz (Receiver& s);
	CustomPz(const CustomPz&);	// copy constructor
	CustomPz& operator=(const CustomPz& rhs);
	~CustomPz() { }

	// clone or factory method ("virtual constructor")
	CustomPz* clone() const { return new CustomPz(*this); }

	// Evaluate the parameterization as a complex AD array
	bool eval(pset& plt, std::vector<std::complex<Ad>>& result, size_t nr, size_t nc);

	// printing: see comments in Pz
	virtual void osprint(std::ostream& s) const;

	// dependent parameter names:
	virtual std::vector<std::string> dependson(pset& plt);

	// Fio required member & functions
	static bool regd;
	static std::string id() { return std::string("CustomPz"); }
	virtual std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new CustomPz(s); }
	virtual void put (Sender& s) const;
}; // class CustomPz

// unsteady aero derivatives for whirling motion of propellers
// as functions of blade angle at 3/4 span. Data is from figure 5 in
// Reed III, Wilmer H and Bland, Samuel R, "An analytical treatment of
//   aircraft propeller precession instability", NASA-TN D-659, 1961
void
whirl_deriv(const Ad& beta, Ad& cmq, Ad& czr,
		Ad& czpsi, Ad& cztheta, Ad& cmpsi);

// open a memory-mapped file and return a pointer to the memory
double*
open_octlab(std::string const& file, size_t n);

#ifdef NEVER // move to flaps
// Class LTIPz: Linear Time-Invariant Pz

// Internal time delays
class Itd {
public:
	std::vector<Elem> elem;
	double deltat;

	Itd() : deltat(0.0) {}
	Itd(double dt) : deltat(dt) {}
};


class LTIPz : public Pz {
	void gainValues(std::vector<std::complex<double>>& igainval,
			std::vector<std::complex<double>>& ogainval);
	void scaleX(int ne, int ns, std::vector<std::complex<double>>& a);

public:
	// constructors, desctructor
	LTIPz() = default;
	//!! LTIPz() : Pz() {}
	// main constructor:
	LTIPz(std::string const& path,
		std::vector<std::pair<int,double>> ecolk,
		std::vector<std::string> igains,
		std::vector<std::string> iphases,
		std::vector<std::string> ogains,
		std::vector<std::string> ophases,
		std::string psiname,
		std::string Kname);
	LTIPz(Receiver& s);
	LTIPz(const LTIPz&) noexcept = default;			// copy constructor
	LTIPz& operator=(const LTIPz& rhs) = default;	// assignment op
	~LTIPz() = default;		// dtor
	// clone or factory method ("virtual constructor")
	LTIPz* clone() const { return new LTIPz(*this); }

	// stuff that comes from lti file
#ifdef NEVER // not pointers
	Matrix* A{nullptr};		// (ns,ns)
	Matrix* B{nullptr};		// (ns,ni)
	Matrix* C{nullptr};		// (no,ns)
	Matrix* D{nullptr};		// (no,ni)
#else // NEVER // not pointers
	Matrix A;		// (ns,ns)
	Matrix B;		// (ns,ni)
	Matrix C;		// (no,ns)
	Matrix D;		// (no,ni)
#endif // NEVER // not pointers
	std::vector<Itd> Atd;
	std::vector<double> itd;	// ni input time delays
	std::vector<double> otd;	// no output time delays
	std::vector<int> Sidx;		// ni s exponents for the (ni,ni) diagonal S matrix

	// stuff supplied by the user to the constructor
	// no pair of columns of K and scale factors:
	std::vector<std::pair<int,double>> ecolk;
	std::vector<std::string> igains;	// input gain parameter names or values
	std::vector<std::string> iphases;	// input phase parameter names or values
	std::vector<std::string> ogains;	// output gain parameter names or values
	std::vector<std::string> ophases;	// output phase parameter names or values
	std::string psiname;					// psi Matrix name
	std::string Kname;					// K matrix name

	std::vector<std::string> dep_par;
	//!! Matrix* KE{nullptr}; 		// (ne,no) XXX custom Pz replaces makeKE

	// "Evaluate" the parameterization...
	bool eval (pset& plt, std::vector<std::complex<Ad>>& result, int nr, int nc);

	// evaluate using approximate exponentials for time-delays?
	//!! static bool useApproxExp(bool yn);

	// printing stuff...
	virtual void osprint(std::ostream& s) const;
	// dependent parameter names:
	virtual std::vector<std::string> dependson(pset& plt) { return dep_par; }

	// plotin stuff
	void plot(std::vector<Elem> const& elements, int nstep);
	void plotElem(int eleno, int nstep, std::string const& runid, bool append);
	void plotEigen();

	// Fio stuff
	static bool regd;		// register_class'd?
	static std::string id() { return std::string("LTIPz"); }
	virtual std::string vid() const { return id(); }
	static Fio* get (Receiver& s) { return new LTIPz(s); }
	virtual void put (Sender& s) const;
};  // LTIPz
#endif // NEVER // move to flaps

#endif // Pz_H
