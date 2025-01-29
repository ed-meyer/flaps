//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <cerrno>
#include <cstdio>	// mmap
#include <cstring>	// mmap
#include <fcntl.h>	// mmap
#include <regex>
#include <unistd.h>
#include <sys/mman.h>	// mmap

#include "config.h"
#include "adeqn.h"
#include "Ad.h"
#include "exim.h"
#include "extract.h"
#include "matrix.h"
#include "Pz.h"
#include "text.h"

using namespace std;
using namespace std::complex_literals;

vector<string>
deepdive(pset& plt, string& name) {
// get all of parameter "name"'s dependencies, and recursively
// call this function with each until parameters are reached
// that have no dependencies, then return only those.
	vector<string> rval;

	Par* pp = plt.findp(name);
	if (pp == nullptr)
		throw runtime_error(vastr("parameter \"",name,"\" is undefined"));
	// get its direct dependents
	vector<string> dep = pp->dependson(plt);
	if (dep.empty()) {
		rval.push_back(name);
		return rval;
	}
	// deep dive each dep
	for (auto& di : dep) {
		vector<string> depi = deepdive(plt, di);
		// put them on rval if not already there
		for (auto& dj : depi)
			if (find(rval.begin(),rval.end(), dj) == rval.end())
				rval.push_back(dj);
	}
	return rval;
}
	

vector<string>
Pz::
get_dependson(pset& plt) {
// returns a list of the immediate dependent parameter names this Pz
// is a function of (no deep-dives). Method: the evaluation function
// (eval) is called with monitoring turned on to capture the names of
// all dependent parameters.
	Trace trc(2,"Pz::get_dependson");

	size_t nr{rsize()};
	size_t nc{csize()};
	// get my immediate dependent parameters by evaluating with monitoring
	plt.monitor();
	vector<complex<Ad>> result(nr*nc);
	try {
		eval(plt, result, nr, nc);
	} catch (runtime_error& s) {
		trc.dprint("ignoring exception: ",s.what());
	}
	vector<string> rval = plt.rotinom();
	trc.dprint("returning ",rval.size()," dep param:\n",rval);
	return rval;
}

//------------------------------------------------------------------
// Elem implementation
//------------------------------------------------------------------
Elem::
Elem (string const& s) {
// Parsing constructor
// Given a string in the form "[i,j]", with i and j 1-based
// indices, parses the row (i) and column (j) number and constructs
// an Elem (with zero-base indices)
	int krow{0};
	int kcol{0};
	regex rx{R"(\[([^,]*),([^\]]*)\])"};
	smatch mch;
	if (!regex_match(s, mch, rx) || (mch.size() < 3) ||
			!str2int(mch[1], krow) || !str2int(mch[2], kcol)) {
		throw runtime_error(vastr("illegal matrix element definition (",s,")"));
	}
	row = krow - 1;
	col = kcol - 1;
}

std::ostream&
operator<<(std::ostream& s, const Elem& t) {
	s << '[' << t.row+1 << ',' << t.col+1 << ']';
	return s;
}

// Pz implementation

Pz::
~Pz() {
// virtual destructor: needs to be in Pz.c even though
// it is empty? causes undefined vtable?
}


bool
Pz::
eval(pset& plt, vector<complex<Ad>>& result, size_t nr, size_t nc) {
	// XXX nothing to do?
	// XXX should be pure virtual?
	return true;
}

Pz::
Pz(Receiver& s) {
	Trace trc(2,"Pz Receiver constructor");
	s.serialize(nrow);
	s.serialize(ncol);
}

void
Pz::
put (Sender& s) const {
	Trace trc(2,"Pz::put");
	s.serialize(nrow);
	s.serialize(ncol);
}

std::ostream&
operator<<(std::ostream& s, const Pz& pz) {
// printing is done by derived Pz member functions: osprint
	pz.osprint(s);
	return s;
}

// end Pz implementation

//------------------------------------------------------------------
// Pz_const is for a matrix that is a constant function of a
// set of parameters, so it just contains the names of the
// parameters and their values.
//------------------------------------------------------------------

// register class Pz_const for Fio::get/put serialization
bool Pz_const::regd = Fio::register_class(Pz_const::id(), Pz_const::get);

Pz_const::
Pz_const (const vector<Par*>& par, size_t nr, size_t nc) : Pz(nr,nc) {
	for (auto pp : par)
		params.push_back(pp->clone());
	// sort the parameters with a lambda as the criterion
	std::sort(params.begin(), params.end(),
			[] (const Par* ap, const Par* bp) { return ap->name < bp->name; } );
}

Pz_const::
Pz_const (Receiver& s) : Pz(s) {
	// s.serialize(dep_par);
	// s.serialize(dep_values);
	size_t n;
	s.serialize(n);
	for (size_t i=0; i<n; i++) {
		Fio* np = Fio::get(s);
		Par* pp = dynamic_cast<Par*>(np);
		if (pp == nullptr) {
			throw runtime_error(vastr("Pz_const::get: expecting Par, got ",np->vid()));
		}
		params.push_back(pp);
	}
}

void
Pz_const::
put (Sender& s) const {
	Trace trc(2,"Pz_const::put");
	Pz::put(s);
	s.serialize(params.size());
	for (auto pp : params) {
		Fio::put(s, *pp);
	}
	// s.serialize(dep_par);
	// s.serialize(dep_values);
}

Pz_const&
Pz_const::
operator=(const Pz_const& rhs) {
// assignment operator
	if (this != &rhs) {
		// invoke base class operator=
		Pz::operator=(rhs);
		params = rhs.params;
	}
	return *this;
}

void
Pz_const::
osprint(std::ostream& s) const {
	// s << "evaluated at ";
	string sep("");
	for (auto pp : params) {
		s << sep << pp->name << " = " << pp->value();
		sep = ", ";
	}
}

vector<string>
Pz_const::
dependson(pset& plt) {
// unlike the other Pz:: dependson() this one returns all
// the parameter names it was evaluated at, so it does
// not need the pset argument
	vector<string> rval;
	for (auto pp : params) {
		rval.push_back(pp->name);
	}
	return rval;
}

//------------------------------------------------------------------
// Interpolation
//------------------------------------------------------------------

// Function matrix_cmp is for sorting matrix pointers according to parameter values.
// Assumptions:
//   - each matrix has 1 Pz_const which contains all dependent parameters
//   - all matrices have the same dependent parameters and in the same order
bool
matrix_cmp(Matrix* ap, Matrix* bp) {
	int m{0};
	double eps = 512.0*std::numeric_limits<double>::epsilon();
	assert(ap->pz.size() == 1 && bp->pz.size() == 1);
	// Pz_const must have a set of Par* with values
	Pz_const* apz = dynamic_cast<Pz_const*>(ap->pz[0]);
	Pz_const* bpz = dynamic_cast<Pz_const*>(bp->pz[0]);
	size_t na{apz->params.size()};
	size_t nb{bpz->params.size()};
	assert(na == nb);
	// start with the first, assuming the last varies quickest (C-style row major)
	// check each parameter until a & b differ more than eps
	for (size_t i=0; i<na; i++) {
		if ((m = fcmp(apz->params[i]->value(), bpz->params[i]->value(), eps)) != 0)
			break;
	}
	if (m <= 0)
		return true;
	return false;
}

IntPz::
IntPz(const vector<Matrix*>& matrices, int extrapol, double rhop) {
/*------------------------------------------------------------------
 * IntPz constructor
 * Compute a set of tensor-product splines fitting a set of data
 * Boundary conditions are chosen to give a natural spline interpolant.
 * Input:
 *   matrices   a vector of pointers to Matrix; each
 *              matrix should have one Pz which tells the
 *              parameters it is a function of and the values of
 *              each parameter.
 *  extrapol    if > 0 it is ok to extrapolate; the default is 0, meaning
 *              if a variable is out if range, evaluate at the closest (min or max)
 *  rhop        smoothing parameter: 0->10: no smoothing -> max (default: 0)
 *------------------------------------------------------------------*/
	Trace trc(1,"IntPz constructor");
	vector<Matrix*> list;
	bool is_complex{false};

	trc.dprint(matrices.size()," matrices, extrap ",extrapol,", rho ",rhop);

	// extrapol means it is ok to eval out-of-range, not really extrapolating,
	// just hold the values constant
	extrap = extrapol;

	// save matrix dimensions: Pz::nrow/ncol
	nrow = matrices[0]->rsize();
	ncol = matrices[0]->csize();

	// rho determines smoothing: 0 means no smoothing, 1-10 more
	rho = rhop;

	// Copy the Matrix pointers to a new vector so we can sort them
	// (input list is const), and check that:
	//  - each matrix has 1 or more Pz_const
	//  - the parameters are the same
	//  Also check to see if any are complex - any that are real must
	//  be cast to complex
	vector<Par*> theparams;
	for(size_t i=0; i<matrices.size(); i++) {
		// dimensions must be the same
		if (matrices[i]->rsize() != nrow || matrices[i]->csize() != ncol)
			throw runtime_error(vastr(matrices[0]->mid()," and ",matrices[i]->mid(),
				" have different dimensions"));
		// 1 Pz_const pz...
		if (matrices[i]->pz.size() != 1)
			throw runtime_error(vastr(matrices[i]->mid()," has ",matrices[i]->pz.size(),
				" parameterizations: ",matrices[i]->pz));
		Pz_const* pzi = dynamic_cast<Pz_const*>(matrices[i]->pz[0]);
		assert(pzi != nullptr);
		// ... with theparams
		if (i == 0) {
			theparams = pzi->params;
			for (auto pp : theparams)
				pp->altval.clear();
		}
		assert(pzi->params.size() == theparams.size());
		trc.dprint("input ",i+1,": ",*matrices[i]);
		list.push_back(matrices[i]);
		is_complex = is_complex || matrices[i]->is_complex();
	}
	trc.dprint("function of ",theparams.size()," parameters");
		
	// sort matrices so that each parameter has increasing values
	trc.dprint("unsorted:\n",matrix_list_summary(list));
	sort(list.begin(), list.end(), matrix_cmp);
	trc.dprint("sorted:\n",matrix_list_summary(list));


	vector<vector<double> > data;
	// put each matrix into data, x, y, & z unique values
	for (size_t i=0; i<list.size(); i++) {
		// add this matrix's data to "data"
		data.push_back(list[i]->cvalues(""));  // cvalues casts to complex
		// may need to cast to complex
		if (is_complex && !list[i]->is_complex()) {
		}
		// put its parameter values (unique values only) into altval
		Pz_const* pzi = dynamic_cast<Pz_const*>(list[i]->pz[0]);
		for (size_t j=0; j<pzi->params.size(); j++) {
			double v = pzi->params[j]->value();
			if (find(theparams[j]->altval.begin(), theparams[j]->altval.end(), v)
					== theparams[j]->altval.end())
				theparams[j]->altval.push_back(v);
		}
	}

	// remove any with only 1 value
	for (auto& pp : theparams) {
		trc.dprint(pp->name," has ",pp->altval.size()," values: ",pp->altval);
		if (pp->altval.size() < 2)
			pp = nullptr;
	}
	theparams.erase(remove(theparams.begin(),theparams.end(),nullptr),
			theparams.end());

	// set the limits of theparams to the min/max of their altvals
	// and update gpset with these limits
	int ngpts{1};
	for (auto& pp : theparams) {
		size_t na = pp->altval.size();
		ngpts *= na;
		Par* existing = gpset::find(pp->name);
		if (existing == nullptr)
			throw runtime_error(vastr("parameter \"",pp->name,"\" is undefined"));
		existing->min(pp->altval[0]);
		existing->max(pp->altval[na-1]);
		//!! gpset::get().add(pp);
	}

	// 1 to 3 parameter interpolation: x, y, z
	vector<string> dep;
	for (auto pp : theparams)
		dep.push_back(pp->name);
	vector<double> x;
	vector<double> y;
	vector<double> z;
	for (auto xi : theparams[0]->altval)
		x.push_back(xi);
	if (theparams.size() > 1) {
		for (auto yi : theparams[1]->altval)
			y.push_back(yi);
	}
	if (theparams.size() > 2) {
		for (auto zi : theparams[2]->altval)
			z.push_back(zi);
	}
	// member "spline" holds any of the interpolation types
	spline = Interp(matrices[0]->mid(), dep, "data");
	// construct a vspline(1/2/3)
	interpolant* ip{nullptr};
	if (theparams.size() == 1) {
		ip = new vspline(x, data, extrap, rho);
	} else if (theparams.size() == 2) {
		ip = new vspline2(x,y,data);
	} else if (theparams.size() == 3) {
		trc.dprint("vspline3 with x=",x,", y=",y,", z=",z);
	 	ip = new vspline3(x,y,z,data);
	}
	if (ip == nullptr)
		throw runtime_error(vastr("cannot iterpolate ",theparams.size()," parameters"));
	spline.add(ip);
}

// register class IntPz for Fio::get/put serialization
bool IntPz::regd = Fio::register_class(IntPz::id(), IntPz::get);

IntPz::
IntPz (Receiver& s) : Pz(s) {
	Trace trc(1,"IntPz Receiver constructor");
	s.serialize(extrap);
	Fio* op = Interp::get(s);
	Interp* bp = dynamic_cast<Interp*>(op);
	assert(bp != nullptr);
	spline = *bp;
	delete bp;
}

void
IntPz::
put (Sender& s) const {
	Trace trc(1,"IntPz::put");
	Pz::put(s);
	s.serialize(extrap);
	spline.put(s);
}

vector<Ad>
IntPz::
get_par_values (pset& plt) {
//
// Helper function for IntPz eval() function.
//
// Put all the parameter values in an array, checking to see
// that they are in range; if out of range the value is set
// to the limit to avoid extrapolation.
	Trace trc(1,"get_par_values");
	vector<Ad> rval;
	vector<string> dep = this->spline.dependson();

	trc.dprint("dep parameter(s): ",dep);

	// Get the parameter values by calling parval with "inrange" set
	bool inrange{true};
	for (auto& di : dep) {
		Ad pv = plt.parval(di, inrange);
		rval.push_back(pv);
	}
	trc.dprint("returning ",rval);
	return rval;
}

bool
IntPz::
eval(pset& plt, vector<complex<Ad>>& result, size_t nr, size_t nc) {
// evaluate a B-spline at the current values of the parameters
	Trace trc(2,"IntPz::eval");

	// Put all the parameter values in an array,
	// with the ordering determined by my parameters (dep),
	// checking to see that they are in range.
	vector<Ad> x = get_par_values(plt);

	// evaluate the Interp
	spline.eval (x, result);

	return true;
}

void
IntPz::
osprint(std::ostream& s) const {
	s << "interpolation wrt ";
	vector<string> dep = this->spline.dependson();
	string sep;
	for (auto& di : dep) {
		s << sep << di;
		sep = ", ";
	}
}

void
IntPz::
plot(vector<int> elements, string const& var, int nstep) {
// not implemented
}

//------------------------------------------------------------------
// Matrix eLeMent: LMPz
//------------------------------------------------------------------

// register class LMPz for Fio::get/put serialization
bool LMPz::regd = Fio::register_class(LMPz::id(), LMPz::get);

LMPz::
~LMPz() { /* virtual desctructor */ }

LMPz::
LMPz(size_t nr, size_t nc, Elem& el, std::string const& equation, char oper) {
	Trace trc(1,"LMPz constructor");

	trc.dprint(el," = ",equation," op ",oper);
	elem = el;
	eqn = equation;
	this->nrow = nr;
	this->ncol = nc;
	op = oper;
	bool is_complex;
	// get the dependent parameter names
	dep_par = adeqn::dependson(gpset::get(), eqn, is_complex);
	trc.dprint("created: ",*this);
}

LMPz::
LMPz (Receiver& s) : Pz(s) {
	Trace trc(2,"LMPz Receiver constructor");
	s.serialize(eqn);
	s.serialize(elem.row);
	s.serialize(elem.col);
	int iop;
	s.serialize(iop);
	op = char(iop);
	s.serialize(dep_par);
}

void
LMPz::
put (Sender& s) const {
	Trace trc(2,"LMPz::put");
	Pz::put(s);
	s.serialize(eqn);
	s.serialize(elem.row);
	s.serialize(elem.col);
	int iop{op};
	s.serialize(iop);
	s.serialize(dep_par);
}

//------------------------------------------------------------------
// Evaluating an LMPz means evaluating the equation and setting
// the matrix element to the result or modifying the element with
// op=, where op is +-*/
//------------------------------------------------------------------
bool
LMPz::
eval(pset& plt, vector<complex<Ad>>& result, size_t nr, size_t nc) {
	Trace trc(1,"LMPz::eval");
	std::ostringstream os;

	int i = elem.row;  // i is 0b
	int j = elem.col;  // j is 0b
	int k = i + nr*j;	// vector element #

	Ad val = adeqn::eval(plt, this->eqn);

	trc.dprint(elem," = ",this->eqn," = ",val,", op: ",op);
	trc.dprint("result before op: ",result[k]);

	if (op == '*')
		result[k] *= val;
	else if (op == '/')
		result[k] /= val;
	else if (op == '+')
		result[k] += val;
	else if (op == '-')
		result[k] -= val;
	else
		result[k] = val;
	trc.dprint("result after op: ",result[k]);

	return true;
}

void
LMPz::
osprint(std::ostream& s) const {
	s << elem;
	if (op != 0)
		s << " " << op << "= " << eqn;
	else
		s << " = " << eqn;
	s << '(';
	string comma{""};
	for (auto& di : dep_par) {
		s << comma << di;
		comma = ",";
	}
	s << ')';
}

vector<string>
LMPz::
dependson(pset& plt) {
	dep_par = get_dependson(plt);
	return dep_par;
}

//------------------------------------------------------------------
// RFAPz: rational-function approximation
//------------------------------------------------------------------

// register class RFAPz for Fio::get/put serialization
bool RFAPz::regd = Fio::register_class(RFAPz::id(), RFAPz::get);

// double (beta & R) version
void
RFAApprox (std::vector<std::complex<double> > const& pval,
		std::vector<double> const& beta,
		size_t n, std::vector<std::complex<double>*> const& matrices,
	   std::vector<std::vector<double>>& R, bool kscale=false, bool forceR0=false);

RFAPz::
RFAPz (std::vector<Matrix*> const& matrices, std::vector<double> betas,
			 bool kscale, bool forceR0, bool optbeta) {
	Trace trc(1,"RFAPz constructor");
	vector<complex<double> > pval;
	vector<complex<double>*> cmlist;
	size_t n = matrices[0]->rsize();
	double rfmin = 1.0e+20;
	double rfmax = 0.0;

	trc.dprint(matrices.size()," matrices", "betas:\n",betas);

	// check that all matrices are complex - cast if not
	for (auto mp : matrices) {
		if (!mp->is_complex()) {
			throw runtime_error(vastr("RFA gaf matrices must be complex: ",
						mp->summary()));
		} else {
			cmlist.push_back(mp->celem());
		}
		// get the parameter value: must be a Pz_const with rf
		// and possibly rsf
		if (mp->pz.size() != 1)
			throw runtime_error(vastr("matrices must have exactly 1 Pz_const: ",
					mp->summary()));

		Pz_const* np = dynamic_cast<Pz_const*>(mp->pz[0]);
		// vector<string> dep = np->dependson();
		double rsf = 0.0;
		double rf = 0.0;
		bool hasrf = false;
		for (auto pp : np->params) {
			if (pp->name == "rf") {
				rf = pp->value();
				hasrf = true;
			}
			if (pp->name == "rsf")
				rsf = pp->value();
		}
		if (!hasrf) {
			throw runtime_error(vastr("matrices must be functions of reduced "
						" frequency for RFA approximation: ", mp->summary()));
		}
		pval.push_back(complex<double>(rsf, rf));
		// gather the limits on rf
		rfmin = std::min(rf, rfmin);
		rfmax = std::max(rf, rfmax);
	}
	// create a RFAPz...
	*this = RFAPz(n, cmlist, pval, betas, kscale, forceR0, optbeta);

	// allow for extrapolation in rf but make sure it has limits
	// for plotting
	Par* rfp = gpset::find("rf");
	if (!rfp->has_min()) {
		rfp->min(rfmin);
	}
	if (!rfp->has_max()) {
		rfp->max(rfmax);
	}
}

RFAPz::
RFAPz (size_t n, std::vector<std::complex<double>*> const& matrices,
			 std::vector<std::complex<double> > pval, std::vector<double> betas,
			 bool kscale, bool forceR0, bool optbeta) {
/*------------------------------------------------------------------
 * RFAPz constructor
 * Method (latex) is in "rfa.tex"
 *------------------------------------------------------------------*/
	Trace trc(1,"RFAPz constructor");
	size_t i;

	trc.dprint(matrices.size()," matrices, ",betas.size()," betas");
	trc.dprint("pvals:\n",pval);
	trc.dprint("betas:\n",betas);

	if (matrices.size() == 0)
		throw runtime_error("attempt to construct a RFA with no matrices");

	if (matrices.size() != pval.size())
		throw runtime_error(vastr("attempt to construct a RFA with ", pval.size(),
			" p-values, but ", matrices.size(), " associated matrices"));

	this->nrow = n;
	this->ncol = n;

	// The one and only parameter is p = s/vtas = complex(rsf,rf)/vtas
	double rsfmin = 1.0e+10;
	double rsfmax = -1.0e+10;
	double rfmin = 1.0e+10;
	double rfmax = -1.0e+10;
	for (i=0; i<pval.size(); i++) {
		double rsf = pval[i].real();
		double rf = pval[i].imag();
		rsfmin = std::min(rsf, rsfmin);
		rsfmax = std::max(rsf, rsfmax);
		rfmin = std::min(rf, rfmin);
		rfmax = std::max(rf, rfmax);
	}

	if (optbeta) {
		optimize (pval, betas, n, matrices, R);
		beta = betas;
	} else {
		beta = betas;
		RFAApprox (pval, beta, n, matrices, R, kscale, forceR0);
	}
}

double
RFAerror(vector<double> const& beta, vector<vector<double>>& R,
	vector<complex<double>> const& ps, vector<complex<double>*> matrices) {
// compute the error in an RFA at "pval" relative to (n,n) complex "matrix"
	int np = ps.size();

	double err{0.0};
	for (int j=0; j<np; j++) {
		vector<complex<double>> rfa;
		RFAPz::eval(ps[j], beta, R, rfa);
		int nsq = rfa.size();
		for (int i=0; i<nsq; i++) {
			complex<double> t = rfa[i] - matrices[j][i];
			err += t.real()*t.real() + t.imag()*t.imag();
		}
	}
	err = sqrt(err);
	return err;
}

double*
open_octlab(string const& file, size_t n) {
// open shared memory file for communicating with octlab, with n doubles;
// returns a pointer to the shared memory with n+1 doubles - the first is
// a semaphore: if it is zero the data has been read and can be overwritten,
// otherwise it is the number of doubles following.
	Trace trc(1,"open_octlab");

	int fd = open(file.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRWXU);
	if (fd == -1)
		throw std::system_error(errno,std::generic_category(),
			vastr("cannot open ",file,": ",strerror(errno)));

	// create the file with n+1 doubles: the first double is the length
	size_t nbytes = (n+1)*sizeof(double);
	vector<double> buf(n+1, 0.0);
	write(fd, buf.data(), nbytes);

	// create a memory-map linked to the file
	double* rval = (double*)mmap(0, nbytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	if (rval == MAP_FAILED) {
		string exc{vastr("cannot map ",file,": ",strerror(errno))};
		trc.dprint("throwing exception: ",exc);
		throw std::system_error(errno,std::generic_category(),
			vastr("cannot map ",file,": ",strerror(errno)));
	}
	return rval;
}
	
string
script(const vector<double>& betas) {
// create a script to run in matlab/octave
// XXX run it in a different thread?
	int maxfunevals{60};

	string cmd("matlab");
	ofstream ofs("optbeta.m", std::ios::trunc);
	ofs << "cd " << getftmp() << endl;
	ofs << "x = [";
	string sep;
	for (auto bi : betas) {
		ofs << sep << bi;
		sep = "; ";
	}
	ofs << "];\n";
	ofs << "opt = optimset(\'MaxFunEvals\'," << maxfunevals << ",\'TolFun\',0.1);\n";
	ofs << "fminsearch(@optfcn,x,opt)\n";
	ofs << "quit\n";
	// create the command line to run: only matlab works - octave doesn't have memmap
	string cmdline{vastr(cmd," -batch \"optbeta\"&")};
	//!! if (matlab_cmd == "octave")
		//!! cmdline = vastr(matlab_cmd, " optbeta.m");

 // unset LD_PRELOAD while running matlab (libstdc++ conflict)
   string ld_preload = getEnv("LD_PRELOAD");
   if (!ld_preload.empty())
      unsetenv("LD_PRELOAD");

	// run matlab/octave in the background
	errno = 0;
	int status = system(cmdline.c_str());
	if (status != 0)
		throw std::system_error(errno,std::generic_category(),
			vastr("could not run matlab/octave: ",strerror(errno)));

	return cmdline;
}

bool
RFAPz::
optimize (vector<complex<double>> const& pval, vector<double> const& beta,
		size_t n, vector<complex<double>*> const& matrices, vector<vector<double>>& R) {
	// compute an optimal set of beta to minimize the error between the RFA
	// and the input matrices
	Trace trc(1,"optimize");

	// change to temp directory using RAII class Chdir
	Chdir tmp(getftmp());

	int nbeta = beta.size();

	vector<double> x = beta;
	trc.dprint("initial x:",x);

	// open memory-mapped files to send/receive data to matlab/octave
	double* enviar = open_octlab("function", 1);
	double* recibir = open_octlab("betas", nbeta);

	// start matlab/octave with the initial betas
	string file = script(beta);

	// loop evaluating betas in new RFA and it's error
	int maxiter{100};
	for (int iter=0; iter<maxiter; iter++) {
		// wait for a recibir message
		while (recibir[0] == 0.0)
			sleep(1);
		int len = recibir[0];
		cout << "got length = " << len << endl;
		if (len != nbeta)
			throw runtime_error(vastr("expected ",nbeta," betas, got ",len));
		for (int i=0; i<len; i++)
			x[i] = recibir[i+1];
		cout << "got x = " << x << endl;
		// signal that we got the message, then optfcn will wait for a reply
		recibir[0] = 0.0;
		trc.dprint("set recibir[0] to 0");
		// compute the new RFA ...
		R.clear();
		RFAApprox (pval, x, n, matrices, R);
		// ... and the error
		double err = RFAerror (x, R, pval, matrices);

		trc.dprint("sending err = ",err);
		// return the function value (err)
		enviar[1] = err;
		enviar[0] = 1.0;	// only one return value
		// wait until optfcn reads the reply
		while(enviar[0] == 1.0)
			sleep(1);
	}

	return true;
}

bool
RFAPz::
eval(pset& plt, vector<complex<Ad>>& result, size_t nr, size_t nc) {
// Add result = cbeta*result + alpha*Q
// where Q = R_0 + pR_1 + p^2 R_2 + (p/(p + beta_1)) R_3 + ...
// p = s/v = complex<double>(rsf, rf)
// result is dimensioned (nr,nc) where n = A->rsize() and nr >= n
	Trace trc(2,"RFAPz::eval");
	std::ostringstream os;

	trc.dprint(R.size()," R matrices");

	if (R.size() < 3) {
		throw runtime_error(vastr("invalid RFA representation: only ",
				R.size()," R matrices"));
	}

	// the R matrices may be dimensioned smaller than "result" (nr,nc)
	size_t n = this->nrow;

	assert(n > 0 && n <= nr);

	// Get the current values of p = (rsf,rf) from gpset
	// If p is out of range, take closest value in range, assuming
	// the correct limits have been set in gpset
	bool inrange{true};
	Ad rsf = plt.parval("rsf", inrange);
	Ad rf = plt.parval("rf", inrange);
	complex<Ad> p(rsf, rf);

	trc.dprint("p = ",p);

	// result = R0 + p*R1 + p*p*R2 = r^t * pv
	// XXX 3 saxpy's?
	vector<complex<Ad>> pv{{Ad(1.0)}, {p}, {p*p}};
	vector<double> r(3);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<n; j++) {
			for (int ir=0; ir<3; ir++)
				r[ir] = R[ir][i+n*j];
			// specialized (double,complex<Ad>) dot
			blas::dot(3, &r[0], 1, &pv[0], 1, result[i+nr*j]);
		}
	}

	// result += R3*(p/(p+beta0)) + R4*(p/(p+beta1)) + ...
	//           = r^t * b
	vector<complex<Ad>> b(beta.size());
	complex<Ad> z;
	r.resize(beta.size());
	for (size_t k=0; k<beta.size(); k++)
		b[k] = p/(p + beta[k]);
	for (size_t i=0; i<n; i++) {
		for (size_t j=0; j<n; j++) {
			for (size_t k=0; k<beta.size(); k++)
				r[k] = R[k+3][i+j*n];
			// use specialized (double,Ad) version of dot
			blas::dot(beta.size(), &r[0], 1, &b[0], 1, z);
			result[i+nr*j] += z;
		}
	}
	return true;
}

bool
RFAPz::
eval(complex<double> const& p, vector<double> const& beta,
	vector<vector<double>> const& R, vector<complex<double>>& result) {
// Evaluate an RFA as a complex (nr,nc) matrix at "p"
// where Q = R_0 + pR_1 + p^2 R_2 + (p/(p + beta_1)) R_3 + ...
// p = s/v = complex<double>(rsf, rf)
	Trace trc(2,"RFAPz::eval(complex<double>)");
	std::ostringstream os;

	trc.dprint(R.size()," R matrices");

	int nR = R.size();
	int nbeta = beta.size();
	if (nR != 3+nbeta)
		throw runtime_error(vastr("invalid RFA representation: ",nbeta," betas, ",
				nR, " R matrices"));

	// result = R0 + p*R1 + p*p*R2
	int nsq = R[0].size();
	result = vector<complex<double>>(nsq);
	blas::axpy(nsq, p, R[1].data(), 1, result.data(), 1);
	complex<double> pp = p*p;
	blas::axpy(nsq, pp, R[2].data(), 1, result.data(), 1);

	// result += R3*(p/(p+beta0)) + R4*(p/(p+beta1)) + ...
	//           = r^t * b
	for (int i=0; i<nbeta; i++) {
		complex<double> t = p/(p+beta[i]);
		blas::axpy(nsq, t, R[3+i].data(), 1, result.data(), 1);
	}
	return true;
}

RFAPz::
RFAPz (Receiver& s) : Pz(s) {
	Trace trc(1,"RFAPz Receiver constructor");
	size_t i;
	size_t nR;
	s.serialize(beta);

	s.serialize(nR);
	trc.dprint("deserializing ",nR," R matrices");
	for (i=0; i<nR; i++) {
		vector<double> ri;
		s.serialize(ri);
		R.push_back(ri);
	}
}

void
RFAPz::
put (Sender& s) const {
	Trace trc(2,"RFAPz::put");
	size_t i;
	size_t nR = R.size();
	Pz::put(s);
	s.serialize(beta);

	s.serialize(nR);
	for (i=0; i<nR; i++)
		s.serialize(R[i]);
}


extern "C" void dgelss_(int *m, int *n, int *nrhs, double *a, int *lda,
				double *b, int *ldb, double *s, double *rcond, int *rank,
				double *work, int *lwork, int *info);

int
dgelss (int m, int n, int nrhs, double* a, int lda,
	double* b, int ldb, double* s, double rcond, int* rank) {

	int info;

	// workspace query
	int lwork{-1};
	double query{0.0};
	dgelss_ (&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond,
		rank, &query, &lwork, &info);
	lwork = query;
	assert(lwork > 0);
	assert(info == 0);

	vector<double> work(lwork, 0.0);

	try {
		dgelss_ (&m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond,
			rank, &work[0], &lwork, &info);
	} catch (runtime_error& s) {
		MM::exporter("dgelss_a.mm",s.what(),a,m,n);
		MM::exporter("dgelss_b.mm",s.what(),b,m,nrhs);
		throw runtime_error(vastr("fp error: see dgelss_a/b.mm: ",s.what()));
	}
	return info;
}

void
RFAApprox (vector<complex<double> > const& pval,
		vector<double> const& beta,
		size_t n, vector<complex<double>*> const& matrices,
	   vector<vector<double>>& R, bool kscale, bool forceR0) {
	Trace trc(1,"RFAApprox(double)");
	size_t nsq = n*n;
	std::ostringstream os;

	trc.dprint(pval.size()," parameter values:",pval,", forceR0? ",forceR0);
	trc.dprint(beta.size(), " betas:",beta,", kscale? ",kscale);

	// allocate the (3+l) R matrices
	int nbeta = beta.size();
	for (int i=0; i<nbeta+3; i++)
		R.push_back(vector<double>(nsq,0.0));

	size_t np = matrices.size();
	int mca = np;               // number of rows in complex a
	int ma = 2*mca;             // number of real rows in a
	int na = 3 + beta.size();  // number of R matrices
	int nrhs = n*n;
	double k0scaleval{1.0};

	// A is (mca,na) complex
	vector<complex<double> > a(mca*na, complex<double>(0.0));
	// B is (mca,n*n) complex
	vector<complex<double> > b(mca*nrhs, complex<double>(0.0));

	//index for pval and matrices
	int k0index = -1;

	if (forceR0 || kscale) {
	    for (size_t i=0; i < pval.size(); ++i){
			if (k0index < 0 && pval[i].imag() == 0) {
				k0index = i; 
			} else if (k0index >= 0 && pval[i].imag() == 0 ){
				flaps::warning("More than one k=0 gaf matrices exist. Using matrix ",
						k0index, "(0b) as R0 matrix.");
			}		    
	    }

	    // if still negative that means not found. so set to something
	    // out of bounds
	    if( k0index < 0 )
			k0index = pval.size();
	}

	// error check for forceR0 and K=0 matrix
	// matrices.size() and pval.size() should be equal
	if (k0index > (int)matrices.size() && forceR0){
		flaps::warning("k=0 matrix is not found. Cannot force R0 to be equal to k=0.",
			 "Continuing with RFA fit without this option.");
	    forceR0 = false;
	}

	// if you have a k=0 matrix and you would like to scale by 1/k find a suitable
	// k value.  In this case it is the min non zero k value over 10
	if (k0index < (int)matrices.size() && kscale){
		std::complex<double> minNonzeroK ( 0, std::numeric_limits<double>::max());
	    for (size_t i=0; i < pval.size(); i++) {
			if (pval[i].imag() < minNonzeroK.imag() && pval[i].imag() != 0 )
				minNonzeroK = pval[i];
	    }
	    k0scaleval = minNonzeroK.imag() / 10.0;
	}

	// set up the A & B matrices
	for (size_t i=0; i<np; i++) {
		std::complex<double> p = pval[i];
		a[IJ(i,0,mca)] = 1.0; 
		
		a[IJ(i,1,mca)] = p;
		a[IJ(i,2,mca)] = p*p;
		for (size_t j=0; j<beta.size(); j++) {
			a[IJ(i,3+j,mca)] = p/(p + beta[j]);
		}
		for (int j=0; j<nrhs; j++) {
		    b[IJ(i,j,mca)] = matrices[i][j];

		    if(forceR0)
				b[IJ(i,j,mca)] -= matrices[k0index][j];
		}

		if(kscale) {
		    double scaleval;
		    if(p.imag() == 0)
				scaleval = k0scaleval;
		    else
				scaleval = p.imag();

			a[IJ(i,0,mca)] /= scaleval;
		    a[IJ(i,1,mca)] /= scaleval;
		    a[IJ(i,2,mca)] /= scaleval;
		    for (size_t j=0; j<beta.size(); j++) {
				a[IJ(i,3+j,mca)] /= scaleval;
		    }
		    for (int j=0; j<nrhs; j++) {
				b[IJ(i,j,mca)] /= scaleval;
		    }
		}
	}

	// Remove the rows and columns here. Update other parameters prior to *gelss.
	// Be careful when removing rows as k=0 matrix may not be in first row
	if(forceR0){
		std::vector<int> rows;
		std::vector<int> cols;
		std::vector<std::complex<double> > newmat;
		int newnr, newnc;

	    for (int i=0; i < mca; i++) {
			if ((int)i == k0index)
				continue;
			else
				rows.push_back(i+1);
	    }
	    for (int i=1; i < na; i++) cols.push_back(i+1);

	    newnr = rows.size();
	    newnc = cols.size();

	    newmat.resize(newnr*newnc, std::complex<double>(0,0));
	    flaps::extract<std::complex<double> >(mca, na, &a[0], newnr, &rows[0],
							 newnc, &cols[0], newnr, &newmat[0]);
	    a = newmat;
	    	    
	    cols.clear();
	    newmat.resize(newnr*nrhs, std::complex<double>(0,0));
	    flaps::extract<std::complex<double> >(mca, nrhs, &b[0], newnr, &rows[0],
							 0, NULL, newnr, &newmat[0]);
	    b = newmat;

	    mca--;
	    ma = 2*mca;
	    na--;
	    np--;
	}

	trc.dprintm(mca,na,ma,a,"least-squares A");
	trc.dprintm(mca,nrhs,ma,b,"least-squares B");
	// solve the linear system a*x = b
	// double rcond = sqrt(std::numeric_limits<double>::epsilon());
	double rcond = -1.0;  // use machine eps
	int rank;
	std::vector<double> s(ma, 0.0);
	int info = dgelss(ma, na, nrhs, reinterpret_cast<double*>(&a[0]), ma,
		reinterpret_cast<double*>(&b[0]), ma, &s[0], rcond, &rank);

	trc.dprint("dgelss returned ",info,", rank ",rank, " (",ma,',',na,"), rcond ",rcond);
	if (rank < (int)na) {
		flaps::warning( "ill-conditioned least-squares problem: ", na-rank,
				" rank-deficiencies");
	}

	// ptr to a real (na,nrhs) matrix with ldb=ma
	double* bp = reinterpret_cast<double*>(&b[0]);

	if (forceR0) {
		//assign R0 matrix
		for (int j=0; j<nrhs; j++) {
			R[0][j] = matrices[k0index][j].real();
		}

		// assign R1-Rn matrices. na is equal to 2+nlags since one column was 
		// removed
		for (int i=0; i<na; i++) {
			for (int j=0; j<nrhs; j++) {
				R[i+1][j] = bp[IJ(i,j,ma)];
			}
			trc.dprintm(n,n,n,R[i],vastr("R",i));
		}
	} else {
	    for (int i=0; i<na; i++) {
			for (int j=0; j<nrhs; j++) {
				R[i][j] = bp[IJ(i,j,ma)];
			}
			trc.dprintm(n,n,n,R[i],vastr("R",i));
	    }
	}
}

void
RFAPz::
osprint(std::ostream& s) const {
	s << "Rational-function approximation with " << beta.size()
		<< " betas";
}

//------------------------------------------------------------------
// CustomPz: custom matrix evaluation code
//------------------------------------------------------------------

// register class CustomPz for Fio::get/put serialization
bool CustomPz::regd = Fio::register_class(CustomPz::id(), CustomPz::get);

CustomPz::
CustomPz (string const& name, size_t nr, size_t nc, int extra) : Pz(nr,nc) {
/*------------------------------------------------------------------
 * CustomPz main constructor. CustomPz code should already be compiled
 * into custom.so on FTMP and made available here with LD_PRELOAD.
 * Input:
 *   name  of the custom function entry point
 *   nr, nc  dimensions of the matrix
 *   extra   number of rows/columns above the other matrices this matrix
 *           will be used with
 *------------------------------------------------------------------*/
	Trace trc(2,"CustomPz constructor");

	entrypt = name;
	// sanity check, arbitrary upper limit:
	assert(extra >= 0 && extra < 1000);
	nExtra = extra;
}

CustomPz::
CustomPz (CustomEval fcn, size_t nr, size_t nc, int extra) : Pz(nr,nc) {
// alternative constructor for the case where the function is not
// "user-written" and needs to be available in the current process,
// e.g. in fem or flut (default_dmatrix)
	Trace trc(2,"CustomPz(function) constructor");
	// ok to leave entrypt empty?
	function = fcn;
	assert(extra >= 0 && extra < 1000);
	nExtra = extra;
}

CustomPz::
CustomPz (Receiver& s) : Pz(s) {
	Trace trc(2,"CustomPz Receiver constructor");
	s.serialize(entrypt);
	s.serialize(nExtra);
	function = nullptr;
	s.serialize(dep_par);
}

void
CustomPz::
put (Sender& s) const {
	Trace trc(2,"CustomPz::put");
	// invoke (call) parent put
	Pz::put(s);
	s.serialize(entrypt);
	s.serialize(nExtra);
	s.serialize(dep_par);
}

// replace these with default?
CustomPz::
CustomPz (const CustomPz& pz) : Pz(pz) {
	Trace trc(2,"CustomPz copy constructor");
	dep_par = pz.dep_par;
	entrypt = pz.entrypt;
	nExtra = pz.nExtra;
	function = pz.function;
}

CustomPz&
CustomPz::
operator=(const CustomPz& rhs) {
// assignment operator
	if (this != &rhs) {
		// invoke base class operator=
		(static_cast<Pz*>(this))->operator=(rhs);
		dep_par = rhs.dep_par;
		entrypt = rhs.entrypt;
		nExtra = rhs.nExtra;
		function = rhs.function;
	}
	return *this;
}



bool
CustomPz::
eval(pset& plt, vector<complex<Ad>>& result, size_t nr, size_t nc) {
// Evaluate a matrix by calling a custom function
	Trace trc(2,"CustomPz::eval");
	if (nr != Pz::nrow) {
		throw runtime_error(vastr("CustomPz::eval: nr(", nr,
					") != Pz::nrow(", Pz::nrow, ')'));
	}

	// if member "function" has not been set get it from CustomPz::find()
	if (function == nullptr) {
		Custom::get(entrypt, function);
		if (function == nullptr) {
			ostringstream os;
			os << "user subroutine " << entrypt << " has not been added. ";
			char* cp = getenv("LD_PRELOAD");
			if (cp == nullptr) {
				os << "LD_PRELOAD has not been set";
			} else {
				os << "LD_PRELOAD = " << cp << endl;
			}
			throw runtime_error(os.str());
		}
	}
	function (plt, (int)nr, (int)nc, result);

	return true;
}

void
CustomPz::
osprint(ostream& s) const {
	s << "custom function \"" << entrypt << "\"";
	if (nExtra > 0)
		s << ", " << nExtra << " extra rows/columns";
	if (!dep_par.empty()) {
		ostringstream os;
		os << "function of " << dep_par.size() << " parameters: ";
		string sep;
		for (auto& di : dep_par) {
			os << sep << di;
			sep = ", ";
		}
		s << endl << roff(os.str(), 4, 80);
	}
}

vector<string>
CustomPz::
dependson(pset& plt) {
// returns the list of dependent parameter names and sets
// member "dep_par"
	dep_par = get_dependson(plt);
	return dep_par;
}

void
whirl_deriv(const Ad& beta, Ad& cmq, Ad& czr,
   Ad& czpsi, Ad& cztheta, Ad& cmpsi) {
// Computes the unsteady aero derivatives as a function of
// blade angle at .75 span. Data is from [1], figure 5.
// XXX move to gaffcn?
   vector<plinear> interp;
   vector<double> xs{34.0, 58.0};

#ifdef NEVER // plinear extrap
	// check for beta in range
	double eps = std::numeric_limits<double>::epsilon();
	double bv = beta.value();
	double x0 = xs[0]*(1.0-eps);
	double xn = xs[xs.size()-1]*(1.0+eps);
	if (bv < x0 || bv > xn)
		throw runtime_error(vastr("beta (",beta,") is out of the range [",x0, ":",xn,"]"));
#endif // NEVER // plinear extrap

   if (interp.empty()) {
      interp.push_back(plinear(xs,{-0.11,-0.03})); // mq
      interp.push_back(plinear(xs,{-0.23,-0.15})); // zr
      interp.push_back(plinear(xs,{0.08,0.09}));   // zpsi
      interp.push_back(plinear(xs,{-0.38,-0.55})); // ztheta
      interp.push_back(plinear(xs,{0.12,0.08}));   // mpsi
   }
   cmq = interp[0].eval(beta);
   czr = interp[1].eval(beta);
   czpsi = interp[2].eval(beta);
   cztheta = interp[3].eval(beta);
   cmpsi = interp[4].eval(beta);
   static int visit{0};
   if (visit == 0) {
      cout << "cmq = " << cmq.value() << endl;
      cout << "czr = " << czr.value() << endl;
      cout << "czpsi = " << czpsi.value() << endl;
      cout << "cztheta = " << cztheta.value() << endl;
      cout << "cmpsi = " << cmpsi.value() << endl;
      visit++;
   }
}

