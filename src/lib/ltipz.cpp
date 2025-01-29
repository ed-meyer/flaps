//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "Pz.h"

using namespace std;
/*------------------------------------------------------------------
 * LTIPz implementation
 *------------------------------------------------------------------*/

#ifdef NEVER // do this in exim
// smoothSplines(int nb, int nix, Real* xb, Real* coef) {
void
LTIPz::
smoothSplines() {
/*
 * Replace the spline coeff that come in the file from
 * matlab/simulink with splines that are smoother.
 * Arguments to this function are a major kludge - passed
 * via environment variables:
 *   SMOOTH_RHO     - smoothness factor: 0 gives orig splines,
 *                    5-10 gives reasonable smoothness, >1000
 *                    gives closer to a straight line
 *   SMOOTH_NSTEP   - make plots of before and after the
 *                    smoothing process if nstep>0
 *
 * coef    a (nb-1, 4, nix) array
 * ABCD:
 *    a(t) = ((c(i,0)*d + c(i,1))*d + c(i,2))*d + c(i,3)
 */
	double eps = sqrt(std::numeric_limits<double>::epsilon());
	int nbm = nb - 1;
	/*
	 * need at least three breakpoints for smoothing
	 */
	if (nb < 3) {
		warningMsg("at least three breakpoints required for smoothing");
		return;
	}

	int nstep = 0;
	char* nstep_ev = getenv("SMOOTH_NSTEP");
	if (nstep_ev)
		nstep = atoi(nstep_ev);

	int i, k;
	vector<double> f(nb, 0.0);
	vector<double> absdiff(nb, 0.0);
	vector<double> reldiff(nb, 0.0);

	double rho = 0.0;    // smoothness parameter
	char* rho_ev = getenv("SMOOTH_RHO");
	if (rho_ev) {
		rho = atof(rho_ev);
	}
	D1(T("using %g for smoothing parameter\n",rho);)

	double maxreldiff = 0.0;
	double maxxb = 0.0;
	size_t maxk = 0;
	double oldy = 0.0, newy = 0.0;

	/*
	 * loop over all variable elements of A
	 */
	for (k=0; k<nix; k++) {
		if (nstep > 0)
			plotElem(k, nstep, "orig", false);
		for (i=0; i<nb; i++)
			f[i] = interpA(xb[i], k, 0);
		int ier = cubgcv(xb, &f[0], nb, rho, &coef[IJK(0,0,k,nbm,4)]);
		if (ier != 0) {
			ostringstream os;
			os << "cubgcv failed, ier = " << ier;
			throw Err(os.str());
		}

		// compute largest difference between breakpoints before
		// and after smoothing relative to the largest breakpoint
		double dmag = Csnrminf(nbm, &f[0], 1);
		for (i=0; i<nbm; i++) {
			absdiff[i] = fabs(coef[IJK(i,3,k,nbm,4)] - f[i]);
			if (dmag > eps)
				reldiff[i] = absdiff[i]/dmag;
			else
				reldiff[i] = 0.0;
			if (reldiff[i] > maxreldiff) {
				maxreldiff = reldiff[i];
				maxxb = xb[i];
				maxk = k;
				oldy = f[i];
				newy = coef[IJK(i,3,k,nbm,4)];
			}
		}
		if (nstep > 0)
			plotElem(k, nstep, "smoothed", true);
	}
	if (nstep > 0)
		printExecInfo("the interpolation for A has been smoothed - "
			"you should compare the original with the smoothed");
	printExecInfo("the max relative change is %g at %s = %g in A[%d,%d]: %g -> %g\n",
			maxreldiff, paramA.c_str(), maxxb, elem[maxk].row+1, elem[maxk].col+1,
			oldy, newy);
}
#endif // NEVER // do this in exim

	// main constructor:
LTIPz::
LTIPz(string const& path) {
#ifdef NEVER //  old args
	string const& psiName, string const& EName,
		string const& KEName, vector<Egc>& egcs, vector<Egc>& sgcs,
		vector<string> ig, vector<string> ip, vector<string> og, vector<string> op,
		string const& nameITD) {
#endif // NEVER //  old args
	//char const* here = "LTIPz constructor";
	vector<string> tok;

 	// import the LTI file from simulink
	LTI::importer(path, A, B, C, D, Atd, itd, otd, Sidx);

	// Egc is an alternative to supplying a KE matrix
	egc = egcs;

	// sgc are scale factors
	sgc = sgcs;

	// Parameters: gain, phase for each input and output
	igain = ig;
	iphase = ip;
	ogain = og;
	ophase = op;

	 // First non-comment line: 6 integers
	 //   ns:   order of A
	 //   ni:   number of input channels
	 //   no:   number of ouput channels
	 //   nb:   number of interpolation breakpoints
	 //         for the elements of A that vary
	 //   nix:  number of elements of A that vary
	 //   nin:  number of elements of A that are constant (this is redundant)
	string line;
	in(line,"#");
	tok = string2tok(line, " ");
	if (tok.size() != 6)
		throw runtime_error(vastr("expecting 6 integers on line ",in.line_number()));

	ns = str2int(tok[0]);
	ni = str2int(tok[1]);
	no = str2int(tok[2]);
	nb = str2int(tok[3]);
	nix = str2int(tok[4]);
	nin = str2int(tok[5]);

	// Second (non-comment) line: interpolation variable name
	in(line, "#");
	tok = string2tok(line, " ");
	if (tok.size() != 1)
		throw runtime_error(vastr("expecting interpolation variable name on line ",
			in.line_number()));

	paramA = tok[0];
	// check that the parameter is in the stdpl list
	// Replace paramA with the real name from stdpl: allow different
	// case and longer name, e.g. "Altitude" instead of "alt"
	// Also put this parameter on the Pz parameter list (par)
	Par* parA{nullptr};
#ifdef NEVER // new
	for (i=0; i<stdpl.size(); i++) {
		if (stringCmpOpt(paramA, stdpl[i]->name)) {
			parA = stdpl[i];
			paramA = stdpl[i]->name;
			P14n::p14npar.push_back(parA->clone());
			break;
		}
	}
#else // NEVER // new
	for (auto& elem : gpset::get()) {
		if (compare(paramA, elem.first)) {
			paramA = elem.first;
			parA = elem.second;
			break;
		}
	}
#endif // NEVER // new
		
	if (parA == nullptr)
		throw runtime_error(vastr("the A-matrix interpolation parameter (",
			paramA,") from the LTI input file (",filename,") has not been defined"));

	// nb double: interpolation breakpoints (alt)
	// If parA has a conversion factor, don't apply it
	xb = new double[nb];
	for (i=0; i<nb; i++) {
		getLine(in, line);
		xb[i] = string2Real(line);
	}
	
	/*
	 * nix pairs of integers: row,col of elements of A that vary
	 */
	int row, col;
	elem = new Elem[nix];
	for (i=0; i<nix; i++) {
		getLine(in, line);
		tok = string2tok(line, " ");
		if (tok.size() != 2) {
			os << "expected 2 integers, line " << getLineNumber()
				<< ", got \"" << line << '\"';
			throw Err(os.str());
		}
		row = string2Int(tok[0]);
		col = string2Int(tok[1]);
		elem[i].row = row - 1;  // zero-based
		elem[i].col = col - 1;  // zero-based
	}
	/*
	 * indices of constant terms of A - ignore
	 */
	for (i=0; i<nin; i++)
		getLine(in, line);
	/*
	 * Input, output time delays
	 */
	if (ni > 0) {
		deltati = new double[ni];
		for (i=0; i<ni; i++) {
			getLine(in, line);
			deltati[i] = string2Real(line);
		}
		D1(printReal(0,"deltati",deltati,ni,1);)
	}
	if (no > 0) {
		deltato = new Real[no];
		for (i=0; i<no; i++) {
			getLine(in, line);
			deltato[i] = string2Real(line);
		}
		D1(printReal(0,"deltato",deltato,no,1);)
	}
	/*
	 * exponents for S: come as floats, convert to ints
	 */
	if (ni > 0) {
		sinput = new int[ni];
		Real si;
		for (i=0; i<ni; i++) {
			getLine(in, line);
			si = string2Real(line);
			if (isEqual(si, 0.0, 5))
				sinput[i] = 0;
			else if (isEqual(si, 1.0, 5))
				sinput[i] = 1;
			else if (isEqual(si, 2.0, 5))
				sinput[i] = 2;
			else {
				os << "illegal value for an S-matrix exponent: " << si;
				throw Err(os.str());
			}
		}
		D1(printInt(0,"sinput",sinput,ni,1);)
	}
	/*
	 * A matrix (constant terms)
	 */
	if (ns > 0) {
		A = new Real[ns*ns];
		for (i=0; i<ns; i++) {
			for (j=0; j<ns; j++) {
				getLine(in, line);
				A[IJ(i,j,ns)] = string2Real(line);
			}
		}
		D1(matrixMarketReal(0,"A",ns,ns,A);)
	}
	// B matrix: (ns,ni)
	if (ni > 0 && ns > 0) {
		B = new Real[ns*ni];
		for (i=0; i<ns; i++) {
			for (j=0; j<ni; j++) {
				getLine(in, line);
				B[IJ(i,j,ns)] = string2Real(line);
			}
		}
		D1(matrixMarketReal(0,"B",ns,ni,B);)
	}
	// C matrix: (no,ns)
	if (no > 0 && ns > 0) {
		C = new Real[no*ns];
		for (i=0; i<no; i++) {
			for (j=0; j<ns; j++) {
				getLine(in, line);
				C[IJ(i,j,no)] = string2Real(line);
			}
		}
		D1(matrixMarketReal(0,"C",no,ns,C);)
	}
	// D matrix: (no,ni)
	if (no > 0 && ni > 0) {
		D = new Real[no*ni];
		for (i=0; i<no; i++) {
			for (j=0; j<ni; j++) {
				getLine(in, line);
				D[IJ(i,j,no)] = string2Real(line);
			}
		}
		D1(matrixMarketReal(0,"D",no,ni,D);)
		int imax = ftn_isamax(no*ni, D, 1);
		if (D[imax-1] == 0.0) {
			delete[] D;
			D = 0;
		}
	}
	/*
	 * cubic-spline coeff for the variable terms of A:
	 * a (nb-1, 4, nix) array
	 */
	if (nb > 1 && nix > 0) {
		size_t nbm = nb - 1;
		coef = new Real[nbm*4*nix];
		for (k=0; k<nix; k++) {
			for (i=0; i<nbm; i++) {
				for (j=0; j<4; j++) {
					getLine(in, line);
					coef[IJK(i,j,k,nbm,4)] = string2Real(line);
				}
			}
		}
		D1(printReal(0,"coef",coef,nbm,4*nix);)
	}
	/*
	 * smooth the splines
	 */
	if (getenv("SMOOTHSPLINES"))
		smoothSplines();
	/*
	 * Optional internal time delays
	 */
	size_t nia = 0;
	try {
		getLine(in, line);
		nia = string2Int(line);
	} catch (Err& s) {
		os.str("");
		os << "no internal time delays in " << path;
		printExecInfo(os.str().c_str());
	}

	for(i=0; i<nia; i++) {
		Itd itd;
		size_t ninda;
		getLine(in, line);
		ninda = string2Int(line);
		for (j=0; j<ninda; j++) {
			getLine(in, line);
			tok = string2tok(line, " ");
			if (tok.size() != 2) {
				os << "expected 2 integers, line " << getLineNumber()
					<< ", got \"" << line << '\"';
				throw Err(os.str());
			}
			row = string2Int(tok[0]);
			col = string2Int(tok[1]);
			itd.elem.push_back(Elem(row-1,col-1));  // zero-based
		}
		internalTimeDelays.push_back(itd);
	}
	for(i=0; i<nia; i++) {
		getLine(in, line);
		internalTimeDelays[i].deltat = string2Real(line);
	}
	/*
	 * Read KE, convert to complex if necessary (the cast takes
	 * care of this)
	 */
	KE = 0;
	if (!KEName.empty()) {
		Array* ka = axdm::fetchMatrix(KEName);
		KE = dynamic_cast<ComplexMatrix*>(ka);
		RealMatrix* RKE = dynamic_cast<RealMatrix*>(ka);
		if (RKE) {
			KE = new ComplexMatrix(RKE);
		} else if (!KE) {
			os.str("");
			os << "unrecognized datatype (" << ka->id()
				<< ") for KE matrix";
			throw Err(os.str());
		}
		D1(printComplex(0,"LTI KE",KE->gm,KE->rsize(),KE->csize());)
	}
	/*
	 * Read the E matrix if name was given, cast to complex...
	 */
	E = 0;
	if (!EName.empty()) {
		Array* e = axdm::fetchMatrix(EName);
		//E = NumericMatrix<Complex>::cast(ke);
		E = dynamic_cast<NumMat*>(e);
		D1(E->display("cerr","LTI E");)
	}
	/*
	 * Read psi
	 */
	psi = 0;
	if (!psiName.empty()) {
		Array* pa = axdm::fetchMatrix(psiName);
		psi = dynamic_cast<RealMatrix*>(pa);
		D1(printReal(0,"LTI psi",psi->gm,psi->rsize(),psi->csize());)
	}

	// nameITD is the name of a parameter which multiplies all internal
	// time delays
	if (!nameITD.empty()) {
		paramITD = nameITD;
		Par* parITD = 0;
		for (i=0; i<stdpl.size(); i++) {
			if (stringCmpOpt(paramITD, stdpl[i]->name)) {
				parITD = stdpl[i];
				paramITD = stdpl[i]->name;
				break;
			}
		}
			
		if (!parITD) {
			os << "the internal time-delay parameter (" << paramITD
				<< ") has not been defined";
			throw Err(os.str());
		}
	}

	check();

	D1(vector<Elem> elements;\
		plot(elements,1001);)
	D1(plotEigen();)

	//D1(X("%s",here);)
}	// LTIPz constructor

void
throwWrongSize(size_t ne, size_t ng) {
	ostringstream os;
	os << "reading LTI: expected " << ne << ", got " << ng;
	throw Err(os.str());
}

LTIPz::
LTIPz(Receiver& s) : P14n(s) {
	size_t nel;
	size_t i, j;
	ostringstream os;

	s.serialize(ns);
	s.serialize(ni);
	s.serialize(no);
	s.serialize(nb);
	s.serialize(nix);
	xb = 0;
	s.serialize(nel, &xb);
	if (nel != nb)
		throwWrongSize(nb, nel);
	int* row = 0;
	s.serialize(nel, &row);
	if (nel != nix)
		throwWrongSize(nix, nel);
	int* col = 0;
	s.serialize(nel, &col);
	if (nel != nix)
		throwWrongSize(nix, nel);
	elem = new Elem[nix];
	for (i=0; i<nix; i++) {
		elem[i].row = row[i];
		elem[i].col = col[i];
	}
	delete[] row;
	delete[] col;

	// deltai: input time-delays
	OIO* op;
	Par* pp;
	deltati = 0;
	s.serialize(nel, &deltati);
	if (nel != ni)
		throwWrongSize(ni, nel);
	// input gains
	s.serialize(nel);
	for (i=0; i<nel; i++) {
		op = getOIO(s);
		pp = dynamic_cast<Par*>(op);
		if (!pp) {
			os << "deserializing LTIPz: unrecognized obj: " << op->vid();
			throw Err(os.str());
		}
		igain.push_back(pp);
	}
	// input phases
	s.serialize(nel);
	for (i=0; i<nel; i++) {
		op = getOIO(s);
		pp = dynamic_cast<Par*>(op);
		if (!pp) {
			os << "deserializing LTIPz: unrecognized obj: " << op->vid();
			throw Err(os.str());
		}
		iphase.push_back(pp);
	}
	// deltao: output time-delays
	deltato = 0;
	s.serialize(nel, &deltato);
	if (nel != no)
		throwWrongSize(no, nel);
	// output gains
	s.serialize(nel);
	for (i=0; i<nel; i++) {
		op = getOIO(s);
		pp = dynamic_cast<Par*>(op);
		if (!pp) {
			os << "deserializing LTIPz: unrecognized obj: " << op->vid();
			throw Err(os.str());
		}
		ogain.push_back(pp);
	}
	// output phases
	s.serialize(nel);
	for (i=0; i<nel; i++) {
		op = getOIO(s);
		pp = dynamic_cast<Par*>(op);
		if (!pp) {
			os << "deserializing LTIPz: unrecognized obj: " << op->vid();
			throw Err(os.str());
		}
		ophase.push_back(pp);
	}
	// S matrix diagonal codes:
	sinput = 0;
	s.serialize(nel, &sinput);
	if (nel != ni)
		throwWrongSize(ni, nel);
	A = 0;
	s.serialize(nel, &A);
	if (nel != ns*ns)
		throwWrongSize(ns*ns, nel);
	B = 0;
	s.serialize(nel, &B);
	if (nel != ns*ni)
		throwWrongSize(ns*ni, nel);
	C = 0;
	s.serialize(nel, &C);
	if (nel != no*ns)
		throwWrongSize(no*ns, nel);
	D = 0;
	s.serialize(nel, &D);
	if (nel > 0 && nel != no*ni)
		throwWrongSize(no*ni, nel);
	coef = 0;
	s.serialize(nel, &coef);
	if (nel != (nb-1)*4*nix)
		throwWrongSize((nb-1)*4*nix, nel);
	// psi matrix: t is guaranteed to be a MatrixList*
	OIO* t = MatrixList::get(s);
	MatrixList* ml = dynamic_cast<MatrixList*>(t);
	if (ml->size() > 0)
		psi = dynamic_cast<RealMatrix*>((*ml)[0]);
	else
		psi = 0;
	// E matrix: t is guaranteed to be a MatrixList*
	t = MatrixList::get(s);
	ml = dynamic_cast<MatrixList*>(t);
	if (ml->size() > 0)
		E = dynamic_cast<NumMat*>((*ml)[0]);
	else
		E = 0;
	// KE matrix: t is guaranteed to be a MatrixList*
	t = MatrixList::get(s);
	ml = dynamic_cast<MatrixList*>(t);
	if (ml->size() > 0)
		KE = dynamic_cast<ComplexMatrix*>((*ml)[0]);
	else
		KE = 0;
	// egc
	s.serialize(nel);
	for (i=0; i<nel; i++) {
		size_t gc;
		Real fac;
		s.serialize(gc);
		s.serialize(fac);
		egc.push_back(Egc(gc,fac));
	}
	// sgc
	s.serialize(nel);
	for (i=0; i<nel; i++) {
		size_t gc;
		Real fac;
		s.serialize(gc);
		s.serialize(fac);
		sgc.push_back(Egc(gc,fac));
	}
	// internal time delays
	size_t nitd;
	s.serialize(nitd);
	for (i=0; i<nitd; i++) {
		Itd itd;
		s.serialize(nel);
		for (j=0; j<nel; j++) {
			size_t ri, cj;
			s.serialize(ri);
			s.serialize(cj);
			itd.elem.push_back(Elem(ri,cj));
		}
		s.serialize(itd.deltat);
		internalTimeDelays.push_back(itd);
	}
	// name of the parameter A is a function of
	s.serialize(paramA);
	// name of the parameter multiplying internal time delays
	s.serialize(paramITD);
}

void
LTIPz::
put (Sender& s) {
	size_t nel;
	size_t i, j;

	P14n::put(s);

	s.serialize(ns);
	s.serialize(ni);
	s.serialize(no);
	s.serialize(nb);
	s.serialize(nix);
	s.serialize(nb, &xb);
	int* row = new int[nix];
	int* col = new int[nix];
	for (i=0; i<nix; i++) {
		row[i] = elem[i].row;
		col[i] = elem[i].col;
	}
	s.serialize(nix, &row);
	s.serialize(nix, &col);
	delete[] row;
	delete[] col;
	// deltati: input time-delays
	s.serialize(ni, &deltati);
	// igain: input gains
	nel = igain.size();
	s.serialize(nel);
	for (i=0; i<igain.size(); i++)
		putOIO(s, *igain[i]);
	// iphase: input phases
	nel = iphase.size();
	s.serialize(nel);
	for (i=0; i<iphase.size(); i++)
		putOIO(s, *iphase[i]);
	// deltato: output time-delays
	s.serialize(no, &deltato);
	// ogain: output gains
	nel = ogain.size();
	s.serialize(nel);
	for (i=0; i<ogain.size(); i++)
		putOIO(s, *ogain[i]);
	// ophase: output phases
	nel = ophase.size();
	s.serialize(nel);
	for (i=0; i<ophase.size(); i++)
		putOIO(s, *ophase[i]);
	// S matrix diagonal codes:
	s.serialize(ni, &sinput);
	// A matrix: (ns,ns) Real
	nel = ns*ns;
	s.serialize(nel, &A);
	// B matrix: (ns,ni) Real
	nel = ns*ni;
	s.serialize(nel, &B);
	// C matrix: (no,ns) Real
	nel = no*ns;
	s.serialize(nel, &C);
	// D matrix: (no,ni) Real - might be null!
	if (D)
		nel = no*ni;
	else
		nel = 0;
	s.serialize(nel, &D);
	// A matrix interpolation coeff: (nb-1,4,nix) Real
	nel = (nb - 1)*4*nix;
	s.serialize(nel, &coef);
	// psi matrix
	MatrixList ml;
	if (psi)
		ml.push_back(psi);
	ml.put(s);
	// E matrix: an Apex matrix - use MatrixList::put
	ml.clear();
	if (E)
		ml.push_back(E);
	ml.put(s);
	// KE matrix: an Apex matrix - use MatrixList::put
	ml.clear();
	if (KE)
		ml.push_back(KE);
	ml.put(s);
	// egc
	nel = egc.size();
	s.serialize(nel);
	for (i=0; i<egc.size(); i++) {
		s.serialize(egc[i].gc);
		s.serialize(egc[i].factor);
	}
	// sgc
	nel = sgc.size();
	s.serialize(nel);
	for (i=0; i<sgc.size(); i++) {
		s.serialize(sgc[i].gc);
		s.serialize(sgc[i].factor);
	}
	// internal time delays
	size_t nitd = internalTimeDelays.size();
	s.serialize(nitd);
	for (i=0; i<nitd; i++) {
		nel = internalTimeDelays[i].elem.size();
		s.serialize(nel);
		for (j=0; j<nel; j++) {
			s.serialize(internalTimeDelays[i].elem[j].row);
			s.serialize(internalTimeDelays[i].elem[j].col);
		}
		s.serialize(internalTimeDelays[i].deltat);
	}
	// name of the parameter A is a function of
	s.serialize(paramA);
	// name of the parameter multiplying internal time delays
	s.serialize(paramITD);
}

Real*
cloneReal(size_t n, Real* from) {
	if (n == 0 || !from) return 0;
	Real* rval = new Real[n];
	for (size_t i=0; i<n; i++)
		rval[i] = from[i];
	return rval;
}

int*
cloneInt(size_t n, int* from) {
	if (n == 0 || !from) return 0;
	int* rval = new int[n];
	for (size_t i=0; i<n; i++)
		rval[i] = from[i];
	return rval;
}

LTIPz::
LTIPz(const LTIPz& rhs) : P14n(rhs) {
// copy constructor
	size_t i;
	ns = rhs.ns;
	ni = rhs.ni;
	no = rhs.no;
	nb = rhs.nb;
	nix = rhs.nix;
	xb = cloneReal(nb, rhs.xb);
	if (nix > 0) {
		elem = new Elem[nix];
		for (i=0; i<nix; i++)
			elem[i] = rhs.elem[i];
	}
	deltati = cloneReal(ni, rhs.deltati);
	for (i=0; i<rhs.igain.size(); i++)
		igain.push_back(rhs.igain[i]->clone());
	for (i=0; i<rhs.iphase.size(); i++)
		iphase.push_back(rhs.iphase[i]->clone());
	
	deltato = cloneReal(no, rhs.deltato);
	for (i=0; i<rhs.ogain.size(); i++)
		ogain.push_back(rhs.ogain[i]->clone());
	for (i=0; i<rhs.ophase.size(); i++)
		ophase.push_back(rhs.ophase[i]->clone());

	sinput = cloneInt(ni, rhs.sinput);
	A = cloneReal(ns*ns, rhs.A);
	B = cloneReal(ns*ni, rhs.B);
	C = cloneReal(no*ns, rhs.C);
	D = cloneReal(no*ni, rhs.D);
	coef = cloneReal((nb-1)*4*nix, rhs.coef);
	if (rhs.psi)
		psi = new RealMatrix(*rhs.psi);
	else
		psi = 0;
	if (rhs.E)
		E = rhs.E->clone();
	else
		E = 0;
	if (rhs.KE)
		KE = new ComplexMatrix(*rhs.KE);
	else
		KE = 0;
	egc = rhs.egc;
	sgc = rhs.sgc;
	internalTimeDelays = rhs.internalTimeDelays;
	// name of the parameter A is a function of
	paramA = rhs.paramA;
	// name of the parameter multiplying internal time delays
	paramITD = rhs.paramITD;
}

LTIPz::
~LTIPz() {
	if (xb)
		delete[] xb;
	if (elem)
		delete[] elem;
	if (deltati)
		delete[] deltati;
	if (deltato)
		delete[] deltato;
	if (sinput)
		delete[] sinput;
	if (A)
		delete[] A;
	if (B)
		delete[] B;
	if (C)
		delete[] C;
	if (D)
		delete[] D;
	if (coef)
		delete[] coef;
	if (psi)
		delete psi;
	if (E)
		delete E;
	if (KE)
		delete KE;
}

LTIPz&
LTIPz::
operator=(const LTIPz& rhs) {
	size_t i;
// assignment operator
	if (this != &rhs) {
		// invoke base class operator=
		(static_cast<P14n*>(this))->operator=(rhs);
		ns = rhs.ns;
		ni = rhs.ni;
		no = rhs.no;
		nb = rhs.nb;
		nix = rhs.nix;
		xb = cloneReal(nb, rhs.xb);
		if (nix > 0) {
			elem = new Elem[nix];
			for (i=0; i<nix; i++)
				elem[i] = rhs.elem[i];
		}
		deltati = cloneReal(ni, rhs.deltati);
		for (i=0; i<rhs.igain.size(); i++)
			igain.push_back(rhs.igain[i]->clone());
		for (i=0; i<rhs.iphase.size(); i++)
			iphase.push_back(rhs.iphase[i]->clone());

		deltato = cloneReal(no, rhs.deltato);
		for (i=0; i<rhs.ogain.size(); i++)
			ogain.push_back(rhs.ogain[i]->clone());
		for (i=0; i<rhs.ophase.size(); i++)
			ophase.push_back(rhs.ophase[i]->clone());

		sinput = cloneInt(ni, rhs.sinput);
		A = cloneReal(ns*ns, rhs.A);
		B = cloneReal(ns*ni, rhs.B);
		C = cloneReal(no*ns, rhs.C);
		D = cloneReal(no*ni, rhs.D);
		coef = cloneReal((nb-1)*4*nix, rhs.coef);
		if (rhs.psi)
			psi = new RealMatrix(*rhs.psi);
		else
			psi = 0;
		if (rhs.E)
			E = rhs.E->clone();
		else
			E = 0;
		if (rhs.KE)
			KE = new ComplexMatrix(*rhs.KE);
		else
			KE = 0;
		egc = rhs.egc;
		sgc = rhs.sgc;
		internalTimeDelays = rhs.internalTimeDelays;
		// name of the parameter A is a function of
		paramA = rhs.paramA;
		// name of the parameter multiplying internal time delays
		paramITD = rhs.paramITD;
	}
	return *this;
}


	// Evaluate an LTI parameterization...
#define REALX Realxp
Real
LTIPz::
interpA(Real xh, size_t eleno, int ider) const {
/*
 * Compute element "eleno" (index into LTIPz::elem) of A
 * or a derivative of A wrt the interpolation parameter at xh.
 * If xh is outside the limits of the interpolant
 * return the value at the limit (ider=0), or zero (ider=1)
 */
	size_t index = 0;
	REALX dx = 0.0;
	Real x = xh;
	size_t i;
	size_t nbm = nb - 1;

	if (x <= xb[0]) {
		index = 0;
		dx = 0.0;
	} else if (x >= xb[nbm]) {
		index = nb - 2;
		dx = xb[nbm] - xb[index];
	} else {
		for (i=0; i<nbm; i++) {
			if (x >= xb[i] && x <= xb[i+1]) {
				dx = x - xb[i];
				index = i;
				break;
			}
		}
	}

	REALX rval = 0.0;
	REALX two(2.0);
	REALX three(3.0);
	if (ider == 0) {
		rval = coef[IJK(index,0,eleno,nbm,4)];
		for (i = 1; i<4; i++) {
			rval = dx*rval + REALX(coef[IJK(index,i,eleno,nbm,4)]);
		}
	} else if (ider == 1) {
		if (xb[0] <= x && x <= xb[nbm])
			rval = REALX(coef[IJK(index,2,eleno,nbm,4)]) +
				two*dx*REALX(coef[IJK(index,1,eleno,nbm,4)]) +
				three*dx*dx*REALX(coef[IJK(index,0,eleno,nbm,4)]);
	} else if (ider == 2) {
		if (xb[0] <= x && x <= xb[nbm])
			rval = 2.0*coef[IJK(index,1,eleno,nbm,4)] +
				6.0*dx*coef[IJK(index,0,eleno,nbm,4)];
		else
			rval = 0.0;
	} else if (ider == 3) {
		if (xb[0] <= x && x <= xb[nbm])
			rval = 6.0*coef[IJK(index,0,eleno,nbm,4)];
		else
			rval = 0.0;
	} else {
		ostringstream os;
		os << "interpA: expected ider=0, 1 or 2, got " << ider;
		throw Err(os.str());
	}
	return rval;
}

bool
limitAlt (Real& alt) {
/*
 * Limit alt to LTI_ALTMAX feet if set in environment
 */
	bool rval = false;
	static bool first = true;
	static Real altmax = 250000.0;
	if (first) {
		char* cp = getenv("LTI_ALTMAX");
		ParList& stdpl = stdParamListRef();
		Param* altp = dynamic_cast<Param*>(stdpl.find("alt"));
		if (cp) {
			altmax = atof(cp);
			printExecInfo("limiting LTI altitude to %g ft",altmax);
		} else if (altp->hasMax()) {
			altmax = altp->const_max();
		} else {
			altmax = 250000.0;
		}
		first = false;
	}
	if (alt > altmax) {
		rval = true;
		alt = altmax;
	}
	return rval;
}

void
LTIPz::
evaluateA(size_t na, Real* a, int ider) const {
	Timer tim("evaluateA", 2);
	ParList& stdpl = stdParamListRef();
	Param* parA = dynamic_cast<Param*>(stdpl.find(paramA));
	if (!parA) {
		ostringstream os;
		os << paramA << " is not a parameter in this analysis";
		throw Err(os.str());
	}
	/*
	 * get the value in external units - that is what the
	 * spline.dat file is based on
	 */
	Real val = parA->prevalue();

	limitAlt(val);
	/*
	 * interpolate each variable element in the A matrix
	 * XXX should I mult by internal time delays here? that would
	 * mean "a" would be complex NO! cant because evalAD does this
	 * using ADs
	 */
	for (size_t i=0; i<nix; i++) {
		if (elem[i].row >= na || elem[i].col >= na) {
			ostringstream os;
			os << elem[i].toString() << " is out of range";
			throw Err(os.str());
		}
		a[IJ(elem[i].row, elem[i].col, na)] = interpA(val, i, ider);
	}
}

void
multRC(size_t nr, size_t nc, size_t nca, Real* a, Complex* b, Complex* c) {
/*
 * Mult C = A*B
 * where A is (nr,nca) real, B is (nca,nc) complex, and C is (nr,nc) complex
 */
	for (size_t i = 0; i<nr; i++) {
		for (size_t j=0; j<nc; j++) {
			Complex t(0.0);
			for (size_t k=0; k<nca; k++) {
				t += a[IJ(i,k,nr)]*b[IJ(k,j,nca)];
			}
			c[IJ(i,j,nr)] = t;
		}
	}
}

void
multCR(size_t nr, size_t nc, size_t nca, Complex* a, Real* b, Complex* c) {
/*
 * Mult C = A*B
 * where A is (nr,nca) complex, B is (nca,nc) real, and C is (nr,nc) complex
 */
	for (size_t i = 0; i<nr; i++) {
		for (size_t j=0; j<nc; j++) {
			Complex t(0.0);
			for (size_t k=0; k<nca; k++) {
				t += a[IJ(i,k,nr)]*b[IJ(k,j,nca)];
			}
			c[IJ(i,j,nr)] = t;
		}
	}
}

void
LTIPz::
gainValues(pset& plt, vector<complex<Ad>>& igainval, vector<complex<Ad>>& ogainval) {
// Evaluate the (complex) gains: find the gain and phase parameters in
// the std param list, evaluate, save in igainval, ogainval.
// Missing gains are set to 1.0
// c_k = g_k e^{i\phi_k}
// where $c_k$ is the complex gain, $g_k$ and $\phi_k$ are the
// gain and phases for the $k^{th}$ input or output, respectively.
	Par* gp;
	Par* pp;
	CParam* cgp;
	Param* rgp;
	Param* rpp;
	Complex gain;
	Real phase;
	size_t i;

	for (i=0; i<ni; i++) {
		complex<Ad> gain(Ad(1.0));
		Ad phase{0.0};
		// if i >= igain.size() use last one...
		if (i < igain.size())
			real(gain) = plt.parval(igain[i]);	// Ad::real
		if (i < iphase.size())
			phase = plt.parval(iphase[i]);
		gain *= exp(Complex(0.0, phase));
		igainval.push_back(gain);
	}
	for (i=0; i<no; i++) {
		complex<Ad> gain(Ad(1.0));
		Ad phase{0.0};
		// if i >= ogain.size() use last one...
		if (i < ogain.size())
			real(gain) = plt.parval(ogain[i]);	// Ad::real
		if (i < ophase.size())
			phase = plt.parval(ophase[i]);
		gain *= exp(Complex(0.0, phase));
		ogainval.push_back(gain);
	}
}

#ifdef NEVER // XXX make KE a parameterized matrix
ComplexMatrix*
LTIPz::
makeKE(size_t ne) {
// Create the KE matrix from the stiffness matrix and either
// 1) a set of g.c. and factors (egc)
// 2) an E matrix
// Or just return the KE pointer if there is one
	static ComplexMatrix* rval = 0;

	// quick return if no outputs or structure, or if
	// no E matrix and no egc...
	// or if there is a KE matrix
	if (KE)
		return KE;
	if (no == 0 || ne == 0 || (egc.empty() && !E)) {
		return 0;
	}

	ostringstream os;
	Complex zero(0.0);
	size_t j;
	int k;

	NumMatList& ml = stdNumMatList();
	NumMat* stif = ml.findDesc("stif");
	if (!stif)
		return 0;

	RealMatrix* rstif = dynamic_cast<RealMatrix*>(stif);
	ComplexMatrix* cstif = dynamic_cast<ComplexMatrix*>(stif);
	/*
	 * Evaluate the stiffness matrix - this allows stiffness parameters
	 * on ouput generalized coordinates
	 */
	stif->eval();

	// create a KE matrix if not already there
	if (!rval)
		rval = new ComplexMatrix("KE","KE",ne,no);

	/*
	 * If an E matrix was given instead of egc's, just
	 * do the multiplication and return; but watch out for
	 * a stiffness matrix that has the extra rows/columns...
	 */
	if (E) {
		NumMat* ke = stif->mult("n", "n", E);
		RealMatrix* rke = dynamic_cast<RealMatrix*>(ke);
		ComplexMatrix* cke = dynamic_cast<ComplexMatrix*>(ke);
		if (rke) {
			ftn_scopy(rke->rsize()*rke->csize(), rke->gm, 1, (Real*)rval->gm, 2);
		} else if (cke) {
			ftn_ccopy(cke->rsize()*cke->csize(), cke->gm, 1, rval->gm, 1);
		} else {
			os << "expecting real or complex from mult K*E, got " << ke->summary();
			throw Err(os.str());
		}
		delete ke;
		return rval;
	}
	/*
	 * generalized-coordinates were given: extract from the
	 * stiffness matrix
	 */
	// check stiffness matrix size against input ne; it may
	// have ns rows/cols added to it...
	if (stif->rsize() != ne && stif->rsize() != ne+ns) {
		os << "expecting a stiffness matrix with " << ne
			<< " rows, got " << stif->summary();
		throw Err(os.str());
	}
	// check that all specified gc are in range
	for (j=0; j<egc.size(); j++) {
		k = egc[j].gc - 1;  // zero-based
		if (k < 0 || k >= (int)ne) {
			os << "a g.c. number (" << egc[j].gc << ") was "
				<< "specified for the E matrix that is out "
				<< "of the range of the stiffness matrix: "
				<< stif->summary();
			throw Err(os.str());
		}
	}

	// real stiffness matrix:
	if (rstif) {
		for (j=0; j<egc.size(); j++) {
			k = egc[j].gc - 1;  // zero-based
			ftn_scopy(ne, &rstif->gm[IJ(0,k,rstif->rsize())], 1, (Real*)&rval->gm[IJ(0,j,ne)], 2);
			if (egc[j].factor != 1.0)
				ftn_sscal(ne, egc[j].factor, (Real*)&rval->gm[IJ(0,j,ne)], 2);
		}
	}

	// Complex stiffness matrix:
	if (cstif) {
		for (j=0; j<egc.size(); j++) {
			k = egc[j].gc - 1;  // zero-based
			ftn_ccopy(ne, &cstif->gm[IJ(0,k,cstif->rsize())], 1, &rval->gm[IJ(0,j,ne)], 1);
			if (egc[j].factor != 1.0) {
				Complex f(egc[j].factor, 0.0);
				ftn_cscal(ne, f, &rval->gm[IJ(0,j,ne)], 1);
			}
		}
	}
	return rval;
}
#endif // NEVER // XXX make KE a parameterized matrix

void
LTIPz::
scaleX(int ne, int ns, Complex* a) {
	int n = ne + ns;

	if (sgc.empty())
		return;

	for (size_t j=0; j<sgc.size(); j++) {
		Complex sc(sgc[j].factor, 0.0);
		int k = sgc[j].gc - 1;  // zero-based
		if (k < ne || k >= n) {
			ostringstream os;
			os << "a g.c. number (" << sgc[j].gc << ") was "
				<< "specified for a scale factor that is out "
				<< "of the range " << ne+1 << " to " << n;
			throw Err(os.str());
		}
		ftn_cscal(n, sc, &a[IJ(0,k,n)], 1);
	}
}

void
LTIPz::
approx (vector<Complex> const& sval, vector<Complex>& beta, vector<ComplexMatrix*>& R) {
/*------------------------------------------------------------------
 * Compute a least-squares RFA fit to the T matrix evaluated at a number
 * of s values (sval), using a set of given beta values.
 *------------------------------------------------------------------*/
	vector<ComplexMatrix*> T;
	Complex s;
	ParList const& stdpl =  stdParamListRef();
	Par* sp = stdpl.find("s");
	Param* sigmap = dynamic_cast<Param*>(stdpl.find("sigma"));
	Param* freqp = dynamic_cast<Param*>(stdpl.find("freq"));
	ostringstream os;
	size_t n = ns;

	if (psi)
		n += psi->csize();

	for (size_t i=0; i<sval.size(); i++) {
		s = sval[i];
		if (sp) {
			CParam* csp = dynamic_cast<CParam*>(sp);
			csp->value() = s;
		} else {
			sigmap->value() = s.real();
			freqp->value() = s.imag();
		}
		os << "T" << i;
		string mid = os.str();
		ComplexMatrix* tp = new ComplexMatrix(mid,mid,n,n);
		eval (Complex(1.0), Complex(0.0), tp->gm, n, n);
		T.push_back(tp);
	}

	RFAApprox (sval, beta, T, R);
}

vector<Complex>
approxExp(Real deltat) {
/*
 * Approximate e^{-s\delta t} between s_{min}
 * and s_{max} by fitting a polynomial minimizing
 * the least-squares error
 */
// exp(-s*deltat) = a0 + a1*s + a2*s^2
//    Ax = b
//    A:   | 1    s_i    s_1^2 |
//    b:   | exp(-s*deltat) |
	ParList const& stdpl =  stdParamListRef();
	Par* sp = stdpl.find("s");
	Param* sigmap = dynamic_cast<Param*>(stdpl.find("sigma"));
	Param* freqp = dynamic_cast<Param*>(stdpl.find("freq"));
	int nx = 3;
	int ny = 5;
	int i, j;
	int nrhs = 1;
	int nra = nx*ny;
	int nca = 3;
	vector<Complex> A(nra*3, Complex(0.0));
	vector<Complex> b(nra, Complex(0.0));
	Complex smin;
	Complex smax;
	ostringstream os;

	if (sp) {
		CParam* cp = dynamic_cast<CParam*>(sp);
		smin = cp->const_min();
		smax = cp->const_max();
	} else {
		if (!(sigmap && freqp)) {
			os << "both sigma and freq must be in the analysis";
			throw Err(os.str());
		}
		Real sigmaMin = -100.0;
		Real sigmaMax = 100.0;
		if (sigmap->hasMin())
			sigmaMin = sigmap->const_min();
		if (sigmap->hasMax())
			sigmaMax = sigmap->const_max();
		smin = Complex(sigmaMin, freqp->const_min());
		smax = Complex(sigmaMax, freqp->const_max());
	}
	Real xmin = smin.real();
	Real xmax = smax.real();
	Real ymin = smin.imag();
	Real ymax = smax.imag();
	Real dx = (xmax - xmin)/(nx-1);
	Real dy = (ymax - ymin)/(ny-1);

	for (i=0; i<nx; i++) {
		Real x = xmin + i*dx;
		for (j=0; j<ny; j++) {
			Real y = ymin + j*dy;
			Complex s = Complex(x,y);
			A[i+j*nx] = 1.0;
			A[i+j*nx+nra] = s;
			A[i+j*nx+2*nra] = s*s;
			b[i+j*nx] = exp(-s*deltat);
		}
	}
	D1(printComplex(0, "ls A", &A[0], nra, nca);)
	D1(printComplex(0, "ls b", &b[0], nra, nrhs);)

	vector<Real> s(nx*ny, 0.0);
	Real rcond = -1.0;  // use machine eps
	int rank;
	int info = ftn_cgelss(nra, nca, nrhs, &A[0], nra, &b[0], nra,
			&s[0], rcond, &rank);
	if (info != 0) {
		os << "least-squares fit to exponential failed: info = " << info;
		throw Err(os.str());
	}
	printReal(0, "ls s", &s[0], 3, 1);
	cout << "a0 = " << b[0] << endl;
	cout << "a1 = " << b[1] << endl;
	cout << "a2 = " << b[2] << endl;
	vector<Complex> rval(3, Complex(0.0));
	rval[0] = b[0];
	rval[1] = b[1];
	rval[2] = b[2];
	return rval;
}

class ExpApprox {
public:
	Real deltat;
	vector<Complex> a;

	ExpApprox(Real dt) {
		deltat = dt;
		a = approxExp(deltat);
	}
	Complex eval(Complex s, int ider) {
		size_t i;
		Complex ct(1.0);
		Complex rval(0.0);
		if (ider == 0) {
			for (i=0; i<a.size(); i++) {
				rval += ct*a[i];
				ct *= s;
			}
		} else if (ider == 1) {
			for (i=1; i<a.size(); i++) {
				rval += double(i)*ct*a[i];
				ct *= s;
			}
		} else if (ider == 2) {
			for (i=2; i<a.size(); i++) {
				rval += double(i)*double(i-1)*ct*a[i];
				ct *= s;
			}
		} else {
			ostringstream os;
			os << "attempt to evaluate derivative " << ider << " of exp(-s*deltat)";
			throw Err(os.str());
		}
		return rval;
	}
}; // ExpApprox


bool UseApproxExp = false;

bool
LTIPz::
useApproxExp(bool yn) {
	bool rval = UseApproxExp;
	D0(T("changing UseApproxExp from %s to %s\n",\
		rval?"true":"false", yn?"true":"false");)
	UseApproxExp = yn;
	return rval;
}

static
Complex
exponential (Complex s, Real deltat, int ider) {
	static vector<ExpApprox*> expapprox;

	if (UseApproxExp) {
		ExpApprox* ep = 0;
		for (size_t i=0; i<expapprox.size(); i++) {
			if (isEqual(deltat, expapprox[i]->deltat, 5)) {
				ep = expapprox[i];
				break;
			}
		}
		if (!ep) {
			ep = new ExpApprox(deltat);
			D0(T("created new exp approx for deltat = %g\n",deltat);)
			expapprox.push_back(ep);
		}
		return ep->eval(s, ider);
	} else {
		if (ider == 0) {
			return exp(-deltat*s);
		} else if (ider == 1) {
			return -deltat*exp(-deltat*s);
		} else if (ider == 2) {
			return deltat*deltat*exp(-deltat*s);
		} else {
			ostringstream os;
			os << "attempt to evaluate derivative " << ider << " of exp(-s*deltat)";
			throw Err(os.str());
		}
	}
}

bool
LTIPz::
eval (pset& plt, vector<complex<Ad>>& result, int nr, int nc) {
// Evaluates an ABCD p14n
	ostringstream os;
	size_t i, j, k;
	vector<complex<Ad>> igainval;
	vector<complex<Ad>> ogainval;
	complex<Ad> zero(Ad(0.0),Ad(0.0));
	complex<Ad> t;
	Timer tim("evalAD", 2);
	
	int ne{0};
	if (psi)
		ne = psi->csize();

	Ad sigma = plt.parval("sigma");
	Ad freq = plt.parval("freq");
	complex<Ad> s(sigma,freq);

#ifdef NEVER // pz'd KE
	// If we have an E matrix (defined by egc) multiply
	// it by the stiffness matrix (alternative to KE matrix)
	Matrix* ke = makeKE(ne);
#endif // NEVER // pz'd KE

	// Evaluate the gains: find the parameters in the std param list,
	// evaluate, save in igainval, ogainval. Missing gains are set to 1.0
	gainValues(plt, igainval, ogainval);

	if (nr != (int)(ns + ne))
		throw runtime_error(vastr("Ltipz eval got nr ",nr,", expected ",ns+ne));

	//!! check();
	
	// (2,2) block: \prod_{k=1}^{n_t} A - sI
	// Do this first since if psi and KE were not specified we are done
	// XXX what does this mean?
	vector<complex<Ad>> work(ns*ns);
	A->eval(plt, work);
	for (j=0; j<ns; j++) {
		for (i=0; i<ns; i++) {
			result[IJ(ne+i,ne+j,nr)] = A[IJ(i,j,ns)];
		}
	}
	// first mult by internal time delays...
	Ad parITDval{1.0};
	if (!paramITD.empty())
		parITDval = parval(paramITD);
	if (parITDval.value() != 0.0) {
		for (k=0; k<internalTimeDelays.size(); k++) {
			Ad ct = exp(-s*internalTimeDelays[k].deltat*parITDval);
			for (size_t l=0; l<internalTimeDelays[k].elem.size(); l++) {
				i = internalTimeDelays[k].elem[l].row;
				j = internalTimeDelays[k].elem[l].col;
				result[IJ(ne+i,ne+j,nr)] *= ct;
			}
		}
	}
	// ... then subtract sI
	for (j=0; j<ns; j++)
		result[IJ(ne+j,ne+j,nr)] -= s;

	// if neither psi nor KE were specified we are done XXX???
	if (psi.empty() && KE == nullptr)
		return true;

	// S:
	vector<complex<Ad>> S(ni);
	for (i=0; i<ni; i++) {
		if (sinput[i] == 0)
			S[i] = complex<Ad>(1.0);
		else if (sinput[i] == 1)
			S[i] = s;
		else
			S[i] = s*s;
	}
	// G, H: input and output gain and time delays
	vector<complex<Ad>> G(ni);
	vector<complex<Ad>> H(no);
	for (i=0; i<ni; i++) {
		G[i] = igainval[i]*exp(-s*deltati[i]);
	}
	// change sign on H
	for (i=0; i<no; i++) {
		H[i] = -ogainval[i]*exp(-s*deltato[i]);
	}

	// evaluate KE -> ke
	vector<complex<Ad>> ke(ne*no);
	KE->eval(plt, ke);

	if (!D.empty()) {
		// (1,1) block: -KE*H*D*G*S*Psi
		vector<complex<Ad>> HDGS(no*ni);
		for (i=0; i<no; i++) {
			for (j=0; j<ni; j++) {
				HDGS[IJ(i,j,no)] = H[i]*D[IJ(i,j,no)]*G[j]*S[j];
			}
		}
		vector<complex<Ad>> KEHDGS(ne*ni);
		// ftn_cgemm("n", "n", ne, ni, no, Complex(1.0), ke->gm, ne,
		// 		&HDGS[0], no, Complex(0.0), &KEHDGS[0], ne);
		for (i=0; i<ne; i++) {
			for (j=0; j<ni; j++) {
				t = zero;
				for (k=0; k<no; k++) {
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
		for (i=0; i<ne; i++) {
			for (j=0; j<ne; j++) {
				t = zero;
				for (k=0; k<ni; k++) {
					t += KEHDGS[IJ(i,k,ne)]*psi[IJ(k,j,ni)];
				}
				result[IJ(i,j,nr)] = t;
			}
		}
	}

	// (1,2) block: -KE*H*C
	vector<complex<Ad>> HC(no*ns, zero);
	for (i=0; i<no; i++) {
		for (j=0; j<ns; j++) {
			HC[IJ(i,j,no)] = H[i]*C[IJ(i,j,no)];
		}
	}
	// vector<ADs> KEHC(ne*ns, zero);
	//ftn_cgemm("n", "n", ne, ns, no, Complex(1.0), ke->gm, ne,
	//		&HC[0], no, Complex(0.0), &KEHC[0], ne);
	//for (j=0; j<ns; j++) {
	//	ftn_ccopy (ne, &KEHC[IJ(0,j,ne)], 1, &result[IJ(0,ne+j,nr)], 1);
	//}
	for (i=0; i<ne; i++) {
		for (j=0; j<ns; j++) {
			t = zero;
			for (k=0; k<no; k++) {
				t += ke[IJ(i,k,ne)]*HC[IJ(k,j,no)];
			}
			result[IJ(i,j+ne,nr)] = t;
		}
	}

	// (2,1) block: B*G*S*Psi
	vector<complex<Ad>> GSPsi(ni*ne, zero);
	for (i=0; i<ni; i++) {
		for (j=0; j<ne; j++) {
			GSPsi[IJ(i,j,ni)] = G[i]*S[i]*psi[IJ(i,j,ni)];
		}
	}
	//vector<ADs> BGSPsi(ns*ne, zero);
	// multRC(ns, ne, ni, &B[0], &GSPsi[0], &BGSPsi[0]);
	// for (j=0; j<ne; j++) {
	// 	ftn_ccopy (ns, &BGSPsi[IJ(0,j,ns)], 1, &result[IJ(ne,j,nr)], 1);
	// }
	for (i=0; i<ns; i++) {
		for (j=0; j<ne; j++) {
			t = zero;
			for (k=0; k<ni; k++) {
				t += B[IJ(i,k,ns)]*GSPsi[IJ(k,j,ni)];
			}
			result[IJ(i+ne,j,nr)] = t;
		}
	}
	/*
	 * Scale the last ns columns
	 */
	// scaleX(ne, ns, result);

	return true;
}

#ifdef NEVER // unused
	// printing stuff...
string
LTIPz::
desc() const {
	return string("ABCD Control Law");
}

void
LTIPz::
print(FILE* stream, const char* title) const {
}

string
LTIPz::
display (string const& path, string const& title) const {
	// XXX when we change to NumericMatrix<Real> for A, B, C, and D
	// we can change these to display
#endif // NEVER // unused

#ifdef NEVER // use this as osprint?
void
LTIPz::
osprint(ostream& s) const {
	size_t i, j, k;
	string offset("    ");
	ostringstream os;

	s << endl;
	if (ni > 0) {
		s << ni << " input time delays:\n";
		for (i=0; i<ni; i++)
			s << offset << "input time delay " << i+1
				<< " = " << deltati[i]*1000.0 << " ms\n";
	} else
		s << "no input time delays";
	s << endl;
	if (no > 0) {
		s << no << " output time delays:\n";
		for (i=0; i<no; i++)
			s << offset << "output time delay " << i+1
				<< " = " << deltato[i]*1000.0 << " ms\n";
	} else
		s << "no output time delays";
	s << endl;

	if (internalTimeDelays.empty()) {
		s << "no internal time delays\n";
	} else {
		s << internalTimeDelays.size() << " internal time delays:\n";
		for (i=0; i<internalTimeDelays.size(); i++) {
			s << offset << "internal time delay " << i+1
				<< ": " << internalTimeDelays[i].deltat*1000.0 << " ms on A:";
			for (j=0; j<internalTimeDelays[i].elem.size(); j++) {
				s << " (" << internalTimeDelays[i].elem[j].row+1
					<< ',' << internalTimeDelays[i].elem[j].col+1 << ')';
				if ((j+1)%10 == 0)
					s << endl << offset << "               ";
			}
			s << endl;
		}
	}

	vector<vector<int> > multitd(ns*ns);
//	for (i=0; i<internalTimeDelays.size()-1; i++) {
//		for (j=0; j<internalTimeDelays[i].elem.size(); j++) {
//			Elem elej = internalTimeDelays[i].elem[j];
//			for (k=i+1; k<internalTimeDelays.size(); k++) {
//				for (l=0; l<internalTimeDelays[k].elem.size(); l++) {
//					if (elej == internalTimeDelays[k].elem[l]) {
//						multitd[IJ(elej.row,elej.col,ns)].push_back(i);
//						multitd[IJ(elej.row,elej.col,ns)].push_back(k);
//					}
//				}
//			}
//		}
//	}
	// for each element of each td add an entry in multitd...
	for (i=0; i<internalTimeDelays.size(); i++) {
		for (j=0; j<internalTimeDelays[i].elem.size(); j++) {
			Elem elej = internalTimeDelays[i].elem[j];
			multitd[IJ(elej.row,elej.col,ns)].push_back(i);
		}
	}
	// count the non-zero td entries for matrix market file...
	int ntd = 0;
	for(i=0; i<ns*ns; i++)
		if (!multitd[i].empty())
				ntd++;

	// how wide is the largest set of td numbers?
	size_t tdfw = 0;
	for (i=0; i<ns; i++) {
		for (j=0; j<ns; j++) {
			if (!multitd[IJ(i,j,ns)].empty()) {
				os.str("");
				for (k=0; k<multitd[IJ(i,j,ns)].size(); k++) {
					os << "  " << setw(3) << multitd[IJ(i,j,ns)][k]+1;
				}
				tdfw = std::max(tdfw, os.str().size());
			}
		}
	}

	string tdmmFile("internalTimeDelays.mm");
	ofstream tdmm(tdmmFile.c_str());
	tdmm << "%%MatrixMarket matrix coordinate real general\n";
	tdmm << ns << " " << ns << " " << ntd << endl;
	s << endl;
	s << ntd << " elements of A have time delays"
		<< " (also in matview file " << tdmmFile << "):\n";
	for (i=0; i<ns; i++) {
		for (j=0; j<ns; j++) {
			if (!multitd[IJ(i,j,ns)].empty()) {
				s << offset << "A[" << setw(3) << i+1
					<< ',' << setw(3) << j+1 << "] has "
					<< multitd[IJ(i,j,ns)].size() << " internal time delays:";
				Real totaltd = 0.0;
				os.str("");
				for (k=0; k<multitd[IJ(i,j,ns)].size(); k++) {
					os << "  " << setw(3) << multitd[IJ(i,j,ns)][k]+1;
					totaltd += internalTimeDelays[multitd[IJ(i,j,ns)][k]].deltat;
				}
				s << setw(tdfw) << left << os.str()
					<< " (total: " << totaltd*1000.0 << " ms)\n";
				tdmm << i+1 << " " << j+1 << " " << totaltd*1000.0 << endl;
			}
		}
	}
	// variable elements of A
	vector<Real> Av(ns*ns, 0.0);
	evaluateA(ns, &Av[0], 0);
	matrixMarketReal(0, "A variable", ns, ns, &Av[0]);

	if (A)
		matrixMarketReal(0, "A", ns, ns, A);
	if (B)
		matrixMarketReal(0, "B", ns, ni, B);
	if (C)
		matrixMarketReal(0, "C", no, ns, C);
	if (D)
		matrixMarketReal(0, "D", no, ni, D);
	return "";
}	// osprint
#endif // NEVER // use this as osprint?

void
LTIPz::
osprint(ostream& s) const {
	size_t ne = 0;
	size_t i;
	if (psi)
		ne = psi->csize();
	s << desc() << " A(" << paramA << '[' << xb[0] << ':'
		<< xb[nb-1] << "]) ns=" << ns << ", ni=" << ni << ", no=" << no
		<< ", nb = " << nb << ", nix=" << nix << ", struct dof = " << ne << endl;
	s << P14n::summary();

	if (igain.empty()) {
		s << "  no input gains\n";
	} else {
		for (i=0; i<igain.size(); i++)
			s << "   input gain: " << igain[i]->summary() << endl;
	}
	if (iphase.empty()) {
		s << "  no input phases\n";
	} else {
		for (i=0; i<iphase.size(); i++)
			s << "   input phase: " << iphase[i]->summary() << endl;
	}

	if (ogain.empty()) {
		s << "  no output gains\n";
	} else {
		for (i=0; i<ogain.size(); i++)
			s << "   output gain: " << ogain[i]->summary() << endl;
	}
	if (ophase.empty()) {
		s << "  no output phases\n";
	} else {
		for (i=0; i<ophase.size(); i++)
			s << "   output phase: " << ophase[i]->summary() << endl;
	}

	if (sgc.empty())
		s << "   no scaling";
	else {
		s << ", sgc: ";
		for (size_t i=0; i<sgc.size(); i++) {
			s << sgc[i].gc << " (" << sgc[i].factor << ')';
		}
	}
}

#ifdef NEVER // unused
std::string
LTIPz::
table() const {
	ostringstream os;
	size_t ne = 0;
	size_t i;
	if (psi)
		ne = psi->csize();
	os << desc() << " A(" << paramA << '[' << xb[0] << ':'
		<< xb[nb-1] << "]) ns=" << ns << ", ni=" << ni << ", no=" << no
		<< ", nb = " << nb << ", nix=" << nix << ", struct dof = " << ne << endl;
	os << P14n::table();

	if (igain.empty()) {
		os << "  no input gains\n";
	} else {
		for (i=0; i<igain.size(); i++)
			os << "   input gain: " << igain[i]->summary() << endl;
	}
	if (iphase.empty()) {
		os << "  no input phases\n";
	} else {
		for (i=0; i<iphase.size(); i++)
			os << "   input phase: " << iphase[i]->summary() << endl;
	}

	if (ogain.empty()) {
		os << "  no output gains\n";
	} else {
		for (i=0; i<ogain.size(); i++)
			os << "   output gain: " << ogain[i]->summary() << endl;
	}
	if (ophase.empty()) {
		os << "  no output phases\n";
	} else {
		for (i=0; i<ophase.size(); i++)
			os << "   output phase: " << ophase[i]->summary() << endl;
	}

	if (sgc.empty())
		os << "   no scaling";
	else {
		os << ", sgc: ";
		for (size_t i=0; i<sgc.size(); i++) {
			os << sgc[i].gc << " (" << sgc[i].factor << ')';
		}
	}

	return os.str();
}

std::vector<Real>
LTIPz::
origValues(std::string const& name) const {
	std::vector<Real> rval;
	for (size_t i=0; i<nb;  i++)
		rval.push_back(xb[i]);
	return rval;
}
#endif // NEVER // unused

void
LTIPz::
plot(std::vector<Elem> const& elements, int nstep) {
//plot(std::vector<int> elements, std::string const& var) {
// for now ignore element, plot all varying A terms...
// XXX do not do this if we smoothed - plotting done there
	char* nstep_ev = getenv("SMOOTH_NSTEP");
	if (nstep_ev)
		return;
	for (size_t i=0; i<nix; i++) {
		plotElem(i, nstep, "unsmoothed", false);
	}
}

void
LTIPz::
plotElem(size_t eleno, int nstep, string const& runid, bool append) {
	ostringstream os;
	size_t i;

	// os << fn << elem[eleno].toString() << ".esa";
	os << "A" << elem[eleno].toString() << ".esa";
	string path(os.str());

	std::ios_base::openmode om;
	if (append)
		om = std::ios::app;
	else
		om = std::ios::trunc;

	ofstream file(path.c_str(), om);
	if (!file) {
		os.str("");
		os << "cannot open " << path;
		throw Err(os.str());
	}

	file << "$ A" << elem[eleno].toString() << endl;
	file << "+1 " << paramA << endl;
	file << "+2 A\n";
	file << "+3 first deriv of A\n";
	file << "+4 second deriv of A\n";
	file << "+5 third deriv of A\n";
	/*
	 * original values: we have nb xb's but only nb-1 coeff
	 */
	file << runid << "-breakpts\n";
	file << "x y dy d2y d3y\n";
	size_t nbm = nb - 1;
	for (i=0; i<nbm; i++) {
		file << xb[i] << " " << coef[IJK(i,3,eleno,nbm,4)]
			<< " " << coef[IJK(i,2,eleno,nbm,4)]
			<< " " << coef[IJK(i,1,eleno,nbm,4)]*2.0
			<< " " << coef[IJK(i,0,eleno,nbm,4)]*6.0 << endl;
	}
	// last breakpoint must use coef
	file << std::setprecision(8) << xb[nbm]
		<< " " << interpA(xb[nbm],eleno,0)
		<< " " << interpA(xb[nbm],eleno,1)
		<< " " << interpA(xb[nbm],eleno,2)
		<< " " << interpA(xb[nbm],eleno,3) << endl;
	file << "*EOF\n";

	/*
	 * Interpolated values
	 */
	file << runid << "-interp\n";
	file << "x y dy d2y d3y\n";
	Complex alpha(1.0);
	Complex beta(0.0);
	int npoint = 4001;
	if (nstep > 0)
		npoint = nstep;
	Real pmin = xb[0];
	Real pmax = xb[nb-1];
	ParList& stdpl = stdParamListRef();
	Param* parA = dynamic_cast<Param*>(stdpl.find(paramA));
	if (parA->hasMin())
		pmin = parA->premin();
	if (parA->hasMax())
		pmax = parA->premax();

	Real factor = 1.0;
	if (parA->conv)
		factor = parA->conv->factor;
#define UNIFORM 1
#ifdef UNIFORM
	Real dp = (pmax - pmin)/(npoint-1);
	for (i=0; i<(size_t)npoint; i++) {
		Real val = pmin + i*dp;
		parA->value() = val/factor;
		file << std::setprecision(8) << val
			<< " " << interpA(val,eleno,0)
			<< " " << interpA(val,eleno,1)
			<< " " << interpA(val,eleno,2)
			<< " " << interpA(val,eleno,3) << endl;
	}
#else
	for (i=0; i<nb-1; i++) {
		Real val = (xb[i] + xb[i+1])/2.0;
		parA->changeFixedValue(val/factor);
		stdpl.touch();
		file << std::setprecision(8) << val
			<< " " << std::setprecision(8) << interpA(val,eleno,0)
			<< " " << std::setprecision(8) << interpA(val,eleno,1)
			<< " " << std::setprecision(8) << interpA(val,eleno,2)
			<< " " << std::setprecision(8) << interpA(val,eleno,3) << endl;
	}
#endif
	file << "*EOF\n";
}

class CompareReals {
public:
	bool operator() (CEigenPtr a, CEigenPtr b) {
		return a->s.real() > b->s.real();
	}
};

void
LTIPz::
plotEigen() {
	ostringstream os;
	size_t i, j;
	string path("ABCDeigen.esa");
	ofstream file(path.c_str(), std::ios::trunc);
	if (!file) {
		os.str("");
		os << "cannot open " << path;
		throw Err(os.str());
	}

	ParList& stdpl = stdParamListRef();
	Param* parA = dynamic_cast<Param*>(stdpl.find(paramA));
	Real factor = 1.0;
	if (parA->conv)
		factor = parA->conv->factor;

	size_t npoint = 101;
	Real xmin = xb[0];
	Real xmax = xb[nb-1];
	Real dx = (xmax - xmin)/(npoint-1);
	vector<vector<Complex> > eigenvalues;
	vector<Complex> eigi;

	for (i=0; i<npoint; i++) {
		Real val = xmin + i*dx;
		parA->changeFixedValue(val/factor);
		stdpl.touch();
		evaluateA(ns, A, 0);
		CEigenList eig = CEigenList::ceigenRealSoln(ns, A);
		eig.sort();
		std::sort (eig.begin(), eig.end(), CompareReals());
		eigi.clear();
		for (j=0; j<eig.size(); j++) {
			eigi.push_back(eig[j]->s);
		}
		for (j=eig.size(); j<ns; j++)
			eigi.push_back(Complex(0.0));
		eigenvalues.push_back(eigi);
	}

	cout << "    A plot file of the eigenvalues of A has been written on "
		<< path << endl;

	file << "$ Eigenvalues of A\n";
	file << "+1 " << paramA << endl;

	for (j=0; j<ns; j++) {
		file << "eig" << j+1 << endl;
		file << paramA << " real imag\n";
		for (i=0; i<npoint; i++) {
			Real val = xmin + i*dx;
			file << val << " " << eigenvalues[i][j].real()
				<< " " << eigenvalues[i][j].imag() << endl;
		}
		file << "*EOF\n";
	}
}
