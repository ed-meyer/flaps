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
#include <sstream>
#include <vector>

#include "config.h"
#include "adeqn.h"
#include "exim.h"
#include "lexer.h"
#include "matrix.h"
#include "settings.h"
#include "specs.h"


using namespace std;

// the Specs are kept in "thespecs" accessed with:
Specs& specs() {
	static Specs thespecs;
	return thespecs;
}

static size_t get_size(string const& str, int dim);
static void change_units(vector<Matrix*>& matrices, vector<double>& betas);

vector<Matrix*>
get_input(vector<tuple<string,string>>& input);
static void
get_params (Matrix* mp, string const& options);
static Matrix*
do_RFA(vector<Matrix*> const& matrices, vector<double> const& betas,
	bool forcer0, bool kscale, bool optbeta, string const& outputname);
static Matrix*
do_interp(vector<Matrix*> const& matrices, string const& outputname,
	int extrap, double smooth);


// Handler functions: *_f
#ifdef NEVER // move to flaps
static bool lti (const Tok& opt);
#endif // NEVER // move to flaps
static bool code_f (const Tok& opt, Matrix*& output_matrix,
		const string& outputname, vector<Matrix*> input);
static bool param_f (const Tok& opt);
static bool plot_f (const Tok& opt, Matrix* output_matrix);

Matrix*
parser () {
//------------------------------------------------------------------
// parse the user's preferences from stdin
//
// parameters may be associated with an input matrix by
// enclosing them in curly braces ({}) after the matrix name, e.g.
//
//  gaf2 {rsf=0.01, rf=0.001, mach=0.4}
//
// CustomPz parameterized matrices are declared with the code option:
//    ..., code=file, ...
// where "file" is either the path of a file or a submit-file
// block containing the user-written function(s)

// Options are parsed in stages by calls to flaps::lexer():
//    1) get the units=uscs option first, reading all prefs from cin
//    2) input, output name and betas
//    3) interpolation & RFA options if multiple input matrices
//    4) custom-code parameterization: can be either generative or modifier
//    5) matrix-element parameterizations
//    6) plotting
//    7) parameter definitions
//------------------------------------------------------------------
	Trace trc(1,"parser");
	Specs& sp = specs();
	
	// 1) get the units=uscs option first, reading all prefs from cin
	//    all other options will be returned as a vector<Tok*>
	vector<Tok*> unrec = flaps::lexer("", {
		{"units", [&](const Tok& p) {return  sp.units = true; }}});

	// 2) now get input, output name, betas - lambda handlers
	//    inputs may have options attach - put the name & options in tuples
	unrec = flaps::lexer(unrec, {
		{"i(nput)?", [&](const Tok& p){
			for (size_t i=0; i<p.svec.size(); i++)
				sp.input.push_back({p.svec[i], p.roptvec[i]});
			return true;
		}},
		{"o(utput)?", [&](const Tok& p) { sp.outputname = p.srhs; return true;}},
		{"out", [&](const Tok& p) { sp.outputname = p.srhs; return true;}},
		{"beta", [&](const Tok& p) {sp.betas = p.rvec; return true; }}
	});

	// fetch the input matrices, change units (gaf and betas only)
	sp.matrices = get_input(sp.input);
	if (sp.units)
		change_units(sp.matrices, sp.betas);

	// if output name was not set, default it to input
	if (sp.outputname.empty()) {
		if (sp.matrices.size() != 1)
			throw runtime_error(vastr("no output matrix name given and ",
				sp.matrices.size()," inputs"));
		sp.outputname = sp.matrices[0]->mid();
	}

	// 3) then look for interp & RFA options if multiple input matrices
	if (sp.matrices.size() > 1) {
		unrec = flaps::lexer(unrec, {
			{"forceR0", [&](const Tok& p) {sp.forcer0 = true; return true; }},
			{"kscale", [&](const Tok& p) {sp.kscale = true; return true; }},
			{"optbeta", [&](const Tok& p) {sp.optbeta = true; return true; }},
			{"^extrap(olate)?", [&](const Tok& p) {return sp.extrap = true;}},
			{"^smooth", [&](const Tok& p) { sp.smooth = p.rvec[0]; return true; }}
		});

		// XXX change these args to Spec&
		if (sp.betas.empty())
			sp.output_matrix = do_interp(sp.matrices, sp.outputname, sp.extrap, sp.smooth);
		else
			sp.output_matrix = do_RFA(sp.matrices, sp.betas, sp.forcer0,
				sp.kscale, sp.optbeta, sp.outputname);
	} // multiple input matrices

#ifdef NEVER // move to flaps
	// 4) LTI parameterization: generative
	unrec = flaps::lexer(unrec, {
		{"lti", [&](const Tok& p) { return lti (p); }}});
#endif // NEVER // move to flaps

	// 4) custom-code parameterization: can be either generative or modifier
	unrec = flaps::lexer(unrec, {
		{"code", [&](const Tok& p) { return code_f (p, sp.output_matrix,
			sp.outputname, sp.matrices); }}});

	// no output_matrix defined yet? must be a single input matrix
	if (sp.output_matrix == nullptr) {
		if (sp.matrices.size() != 1)
			throw runtime_error("at least one input matrix must be specified");
		sp.output_matrix = sp.matrices[0];
		sp.output_matrix->mid(sp.outputname);
	}

	// 5) matrix-element parameterizations
	size_t nr = sp.output_matrix->rsize();
	size_t nc = sp.output_matrix->csize();
	unrec = flaps::lexer(unrec, {
		{"^\\[[0-9]*,[ \t]*[0-9]*\\]$", [&](const Tok& p) {
			Elem el(p.lhs);
			sp.output_matrix->pz.push_back(new LMPz(nr, nc, el, p.srhs, p.op));
			return true;
		}},
	});

	Matrix* rval = sp.output_matrix;

	// 6) look for plot
	unrec = flaps::lexer(unrec, {
		{"plot", [&](const Tok& p) { return plot_f(p, sp.output_matrix); } }});

	// 7) no options left - maybe a parameter defn?
	if (!unrec.empty())
		unrec = flaps::lexer(unrec, {
			{"^.*", [&](const Tok& p) { return param_f(p); }}});

	if (!unrec.empty()) {
		string exc{vastr(unrec.size()," options were not recognized: ",unrec)};
		throw runtime_error(exc);
	}

	trc.dprint("returning ",rval->summary());
	return rval;
}

Matrix*
do_RFA(vector<Matrix*> const& matrices, vector<double> const& betas,
	bool forcer0, bool kscale, bool optbeta, string const& outputname) {
// Create a new matrix with an RFA parameterization

	// Create a new matrix with this RFA
	size_t nr = matrices[0]->rsize();
	size_t nc = matrices[0]->csize();
	Matrix* rval = new Matrix(outputname, "RFA", nr, nc, true);
	rval->pz.push_back(new RFAPz(matrices, betas, kscale, forcer0, optbeta));
	return rval;
}

#ifdef NEVER // move to flaps
static bool
lti (const Tok& opt) {
// Parse options for the "lti" option and create the output matrix with an
// LTIPz; options:
//   lti{file=path, psi=name, e=(list{scale})|name, igain=(list)|parname, ogain=(list)|parname}
	Specs& sp = specs();
#ifdef NEVER // use pz
	string path;
	string psiname;
	string stifname;
	vector<pair<int,double>> ecolk;
	vector<string> igains;
	vector<string> iphases;
	vector<string> ogains;
	vector<string> ophases;
#endif // NEVER // use pz
	LTIPz* pz = new LTIPz();

	// parse the lopt part of opt
	vector<Tok*> unrec = flaps::lexer(opt.lopt, {
		{"file", [&](const Tok& p) {
			LTI::importer (p.srhs, pz->A, pz->B, pz->C, pz->D, pz->Atd,
				pz->itd, pz->otd, pz->Sidx);
			return true; }},
		{"psi", [&](const Tok& p) { pz->psiname = p.srhs; return true; }},
		{"stif", [&](const Tok& p) { pz->Kname = p.srhs; return true; }},
		{"e", [&](const Tok& p) {
			if (p.srhs.empty()) return false;
			if (isdigit(p.srhs[0])) {	// pairs
				for (size_t i=0; i<p.ivec.size(); i++) {
					double sf{1.0};
					if (!p.roptvec.empty())
						sf = stod(p.roptvec[i]);
					pz->ecolk.push_back({p.ivec[i],sf});
				}
			} else {
				return false;	// e=name disabled for now
			}
			return true;
		}},
		{"igain", [&](const Tok& p) {		// parameter name or value
			pz->igains = p.svec; return true;
		}},
		{"iphase", [&](const Tok& p) {		// parameter name or value
			pz->iphases = p.svec; return true;
		}},
		{"ogain", [&](const Tok& p) {		// parameter name or value
			pz->ogains = p.svec; return true;
		}},
		{"ophase", [&](const Tok& p) {		// parameter name or value
			pz->ophases = p.svec; return true;
		}},
	});

	// do some checking
	if (pz->Kname.empty())
		throw runtime_error("a stiffness matrix was not specified");
	
	// dependent parameters: sigma, freq + K->dependson + A->dependson
	pz->dep_par.push_back("sigma");
	pz->dep_par.push_back("freq");
	auto pbegin = pz->dep_par.begin();
	auto pend = pz->dep_par.end();
	vector<string> kdep = K->dependson();
	for (d : kdep) {
		if (find(pbegin,pend,d) == pend)
			pz->dep_par.push_back(d);
	}
	vector<string> adep = A->dependson();
	for (d : adep) {
		if (find(pbegin,pend,d) == pend)
			pz->dep_par.push_back(d);
	}

	// read the stiffness matrix (fetch constructor) and get it's size
	Matrix K(pz->Kname);
	int ne = K.rsize();
	// create the output matrix: (ne+ns,ne+ns) complex...
	int ns = pz->A->rsize();
	int m = ne + ns;
	sp.output_matrix = new Matrix(sp.outputname, "LTI T Matrix", m, m, true);
	// ... and add the pz
	sp.output_matrix->pz.push_back(pz);

	return true;
}
#endif // NEVER // move to flaps

static bool
code_f (const Tok& opt, Matrix*& output_matrix,
	const string& outputname, vector<Matrix*> input) {
/*
 * The "code" option triggers the creation of a CustomPz.
 * Syntax (allows for multiple pz):
 *    code=filename{size=string, extra=int}
 *       filename:  name of the file or submit-file block containing
 *                  C++ code. This name must match the entry-point
 *                  name
 *       sizestr:   (string) name of a matrix which is the same size
 *                  as what the code expects
 */
	Trace trc(1,"code_f");
	ostringstream os;

	if (opt.svec.empty())
		throw runtime_error(vastr("option (",opt,") has no rhs"));

	string filename = opt.svec[0];
	trc.dprint("CustomPz code is in ",filename);

	// check the rhs options for "extra" and "size"
	// The actual size of the output matrix is size+extra
	//  size = n            matrix is (n,n)
	//  size = (n,m)        matrix is (n,m)
	//  size = mid          matrix has the same dimensions as "mid"
	//  size = (nr, mid):   matrix has nr rows and same number of
	//                      columns as "mid"
	//  size = (mid, nc):   matrix has nc columns and same number of
	//                      rows as "mid"
	int nextra{0};
	size_t rsize{0};
	size_t csize{0};
	if (!opt.ropt.empty()) {
		vector<Tok*> prefs = flaps::lexer(opt.ropt);
		for (auto pp : prefs) {
			if (pp->lhs == "extra")
				nextra = pp->ivec[0];
			else if (pp->lhs == "size") {
				// there should be either one or two rhs...
				if (pp->svec.size() == 1) {
					rsize = get_size(pp->svec[0], 1);
					csize = get_size(pp->svec[0], 2);
				} else if (pp->svec.size() == 2) {
					rsize = get_size(pp->svec[0], 1);
					csize = get_size(pp->svec[1], 2);
				} else {
					throw runtime_error("illegal size specification: only 2 values allowed");
				}
			}
		}
	}
	if (csize == 0)
		csize = rsize;
	// if the output dimensions were not specified take them
	// from the input matrix if there is one
	// if the output matrix hasn't been assigned, either take the
	// input (modifier) or create a new one (generative) with outputname
	if (output_matrix == nullptr) {
		if (!input.empty()) {
			if (input.size() > 1)
				throw runtime_error("multiple matrices input with the \"code\" option");
			if (rsize != 0 || csize != 0)
				throw runtime_error("cannot specify \"size\" with an input matrix");
			output_matrix = input[0];
			output_matrix->mid(outputname);
			rsize = output_matrix->rsize();
			csize = output_matrix->csize();
		} else {
			if (rsize == 0)
				throw runtime_error("row size not specified");
			output_matrix = new Matrix(outputname, filename, rsize+nextra, csize+nextra, true);
		}
	}

	size_t nr{rsize+nextra};
	size_t nc{csize+nextra};
	// the file is assumed to already be in custom.so, the entry point is the base name ...
	string entrypt = stringBasename(filename);
	// ... create a CustomPz ...
	output_matrix->pz.push_back(new CustomPz(entrypt, nr, nc, nextra));

	return true;
}

vector<Matrix*>
get_input(vector<tuple<string,string>>& input) {
// Fetch input matrices and add Pz_const if the second of the tuple
// is a set of param values
	Trace trc(1,"get_input");
	vector<Matrix*> rval;
	for (auto& in : input) {
		vector<Matrix*> mps = Matrix::fetch(get<0>(in));
		//!! Matrix* mp = new Matrix(get<0>(in));
		//!! if (mp == nullptr)
			//!! throw runtime_error(vastr(get<0>(in), " is not available"));
		for (auto mp : mps) {
			get_params (mp, get<1>(in));
			rval.push_back(mp);
		}
	}
	// put on the list of input matrices (Matrix::collection())
	for (auto mp : rval)
		Matrix::insert(mp);

	flaps::info ("Input matrix(s):\n", Matrix::inventory());
	trc.dprint("now have input matrix list:\n",Matrix::inventory());

	return rval;
}

static void
change_units(vector<Matrix*>& matrices, vector<double>& betas) {
// change units from USCS to SI for each matrix in "matrices"
// but only if they are aero matrices, functions of rf and
// possibly Mach and rsf.
// Rationale:
// according to table 2.3 in the manual, to convert the flutter equation
// from USCS to SI, all matrices except aero are multiplied by
// 0.112984829 N-m/lb_f-in, and the aero matrix is multiplied by
// 1.6387064(-5) m^3/in^3, so rather than scaling all matrices it is easier
// and safer to divide the flutter equation by 0.112984829 and scale
// the aero matrix by (1.6387064e-5 m^2/in^3)/0.112984829 = 1.45037737766e-4,
// or divided by 6894.757291467.
// Reduced frequencies are scaled by 1/0.0254 m/in
// Betas are also scaled by 1/0.0254 if not empty
	Trace trc(1,"change_units");
	double scale{1.45037737766e-4};
	double rfscale = 1.0/0.0254;

	for (auto mp : matrices) {
		assert(mp->pz.size() == 1);
		Pz_const* pzi = dynamic_cast<Pz_const*>(mp->pz[0]);
		assert(pzi != nullptr);
		bool is_aero{false};
		// to convert rf = \omega b/V and rsf we just need to
		// convert the ref length "b" which is assumed to be 1", so 
		// scale them so b=1m: 1/0.0254 in/m
		for (auto& pi : pzi->params) {
			if (pi->name == "rf" || pi->name == "rsf") {
				double t = pi->value()*rfscale;
				pi->value(t);
				is_aero = true;
			}
		}
		if (!is_aero)
			throw runtime_error(vastr("input matrix \"",mp->mid(),
					" is not a function of rf"));
		// scale the matrix by 1.45037737766(-4)
		vector<double>& mdata = mp->data();
		blas_scal(mdata.size(), scale, &mdata[0], 1);
	}
	// scale betas
	if (!betas.empty())
		blas_scal(betas.size(), 1.0/0.0254, betas.data(), 1);

	flaps::info("units have been changed from USCS to SI: scaled by ",scale,
			" (divided by 6894.757291467),");
	if (!betas.empty())
		flaps::info("reduced-frequencies and reduced stability factors have "
			"been multiplied by ",rfscale);
	else
		flaps::info("reduced-frequencies, reduced stability factors, "
			"and betas have been multiplied by ", rfscale);
}

static Matrix*
do_interp(vector<Matrix*> const& matrices, string const& outputname,
	int extrap, double smooth) {
// Interpolate the input matrices wrt 1, 2, or 3 parameters

	// are any of the matrices complex?
	bool cmplx{false};
	for (auto mp : matrices) {
		if (mp->is_complex()) {
			cmplx = true;
			break;
		}
	}
	if (cmplx) {
		for (auto mp : matrices) {
			if (!mp->is_complex())
				mp->cast_complex();
		}
	}

	size_t nr = matrices[0]->rsize();
	size_t nc = matrices[0]->csize();
	Matrix* rval = new Matrix(outputname, "Interp", nr, nc, cmplx);
	// do the interpolation
	rval->pz.push_back(new IntPz(matrices, extrap, smooth));
	return rval;
}

static bool
plot_f (const Tok& opt, Matrix* output_matrix) {
// Plot matrix elements; all preferences are in curly
// braces following "plot":
// Certain diagonals:
//   plot{diag=(1,3,5 to 10), ...}
// A range of rows, columns:
//   plot{rows=(1 to 10), col=(14,15)}
// individual elements:
//   plot{[2,3], [4,6], ...}
// Take ns steps in plot
//   plot{nstep=ns, ... }
// Use "parname" as independent variable:
//   plot{..., indep=parname}
// hold "parname" constant:
//   plot{..., fixed=parname}
	Trace trc(1,"plot_f");
	size_t i, j;
	vector<int> rows, cols;
	vector<Elem> elements;
	vector<string> indep_par;
	string fixed_par;
	size_t nstep{101};
	Specs& sp = specs();

	// Create a copy of the gpset so we can change things and
	// not disturb it
	pset altpl(gpset::get());
	altpl.desc("Matrix::plot");

	// all options must be in curly braces: plot{ ... }
	if (opt.lopt.empty()) {
		flaps::warning("no options given with the \"plot\" option");
		return true;
	}

	// parse the lopt part of opt
	vector<Tok*> unrec = flaps::lexer(opt.lopt, {
		{"diag(onal)?", [&](const Tok& p) {
			if (!p.ivec.empty()) {
				for (j=0; j<p.ivec.size(); j++)
					elements.push_back(Elem(p.ivec[j]-1, p.ivec[j]-1));
			// just the word "diagonal": plot all diags
			} else {
				size_t nr = sp.output_matrix->rsize();
				for(size_t k=0; k<nr; k++)
					elements.push_back(Elem(k,k));
			}
			return true;
		}},
		{"nstep", [&](const Tok& p) { nstep = p.ivec[0]; return true; }},
		{"^\\[[0-9]*,[ \t]*[0-9]*\\]$", [&](const Tok& p) {
			elements.push_back(Elem(p.lhs));  return true;}},
		{"row", [&](const Tok& p) {
			for (size_t j=0; j<p.ivec.size(); j++)
				rows.push_back(p.ivec[j]-1);
			return true;
		}},
		{"col", [&](const Tok& p) {
			for (size_t j=0; j<p.ivec.size(); j++)
				cols.push_back(p.ivec[j]-1);
			return true;
		}},
		{"indep", [&](const Tok& p) {
			for (auto& par : p.svec) {
				if (std::find(indep_par.begin(), indep_par.end(), par)
						== indep_par.end())
					indep_par.push_back(par);
			}
			return true;
		}},
		{".*", [&](const Tok& p) {
			// Must be a fixed parameter; if one or more rhs included,
			// hold this parameter const, possibly at multiple values,
			// otherwise it is an indep_par
			Par* ep = altpl.findp(p.lhs);
			if (ep == nullptr)
				return false;
			if (p.rvec.empty()) {
				indep_par.push_back(p.lhs);
			} else {
				ep->valuef(p.rvec[0]);
				ep->altval.clear();
				if (p.rvec.size() > 1) {
					for (auto x : p.rvec)
						ep->altval.push_back(x);
				}
				ep->set_fixed();
				fixed_par = p.lhs;
			}
			return true;
		}} });

	if (!rows.empty()) {
		if (!cols.empty()) {
			for (i=0; i<rows.size(); i++) {
				for (j=0; j<cols.size(); j++) {
					elements.push_back(Elem(rows[i],cols[j]));
				}
			}
		}
	}

	trc.dprint("got ",elements.size()," plot elements");

	// do the plotting
	sp.output_matrix->plot(altpl, elements, indep_par, fixed_par, nstep);

	return true;
}

static size_t
get_size(string const& str, int dim) {
//------------------------------------------------------------------
// Get a matrix dimension from a string which is either
// an integer or a matrix name. If it is a matrix name, read
// the matrix and return either the first or second dimension
// according to "dim"
//------------------------------------------------------------------
	size_t rval{0};
	int k;

	if (str2int(str, k)) {
		rval = (size_t)k;
	} else {
		vector<Matrix*> mps = Matrix::fetch(str);
		if (mps.empty()) {
			string exc = vastr("matrix \"",str,"\" is not available");
			throw runtime_error(exc);
		}
		Matrix* mp = mps[0];
		if (dim == 1)
			rval = mp->rsize();
		else
			rval = mp->csize();
		delete mp;
	}
	return rval;
}

static bool
param_f (const Tok& opt) {
//------------------------------------------------------------------
// Handle any options not yet recognized: the only legal possibility
// is a parameter definition
//------------------------------------------------------------------
	Trace trc(1,"param_f");
	ostringstream os;

	Par* param{nullptr};
	try {
		double x;
		if (str2double(opt.srhs, x)) {
			// rhs is a double so it is fixed
			param = new Par(opt.lhs, opt.srhs);
			param->set_fixed();
		} else {
			// rhs is an equation
			param = new Par(opt.lhs, {opt.srhs});
		}
	} catch (runtime_error& s) {
		trc.dprint("not a parameter: return false");
		return false;
	}
	param->pref = true;
	// Add this new parameter to the gpset
	Par* existing = gpset::find(param->name);
	if (existing == nullptr) {
		param = gpset::get().add(param);
		flaps::info("adding new parameter (", param->summary(),") to standard list");
	}

	return true;
}

/*------------------------------------------------------------------
 * end of handler functions
 *------------------------------------------------------------------*/


static
void
get_params (Matrix* mp, string const& options) {
/*------------------------------------------------------------------
 * Parse a set of options associated with a matrix, e.g.
 * gm2{cg=4}. Add a Pz_const to the matrix with the parameter(s)
 *
 * The legal options for a matrix are:
 *  par=val    value of a parameter to be associated with this matrix
 *------------------------------------------------------------------*/
	Trace trc(1,"get_params ", options);
	ostringstream os;

	if (options.empty()) {
		return;
	}

	// parse with no handlers - just get Toks
	vector<Tok*> toks = flaps::lexer(options);

	// the only options allowed are parameter definitions
	//! vector<string> pnames;
	//! vector<double> pvals;
	vector<Par*> params;
	for (auto ti : toks) {
		Par* pp = new Par(ti->lhs, ti->svec[0]);
		params.push_back(pp);
		// first add it to the gpset in case it's not already there...
		gpset::get().add(pp);
		// ... save the name & value
		//! pnames.push_back(pp.name);
		//! pvals.push_back(pp.value());
	}

	// ... then create a Pz_const
	size_t nr = mp->rsize();
	size_t nc = mp->csize();
	//! Pz_const* op = new Pz_const(pnames, pvals, nr, nc);
	Pz_const* op = new Pz_const(params, nr, nc);
	mp->pz.push_back(op);
	trc.dprint("added Pz_const to ",mp);
	return;
}

