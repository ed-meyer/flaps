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

#include "conv.h"
#include "lapack.h"
#include "lexer.h"
#include "matrix.h"
#include "specs.h"
#include "settings.h"

using namespace std;

bool
is_flutmat (string const& desc) {
/*------------------------------------------------------------------
 * Is "desc" one of the allowable matrices for flutter runs?
 * e.g. "mass", "stif", etc.
 * Returns the matched description if found, empty string otherwise
 *------------------------------------------------------------------*/
	Trace trc(1,"is_flutmat");

	trc.dprint("is ",desc," a flutter matrix?");

	// These are the (only) options for flutter matrices
	static vector<string> valid;
	if (valid.empty()) {
		valid.push_back("mass");
		valid.push_back("stif");
		valid.push_back("gyro");
		valid.push_back("vdamp");
		valid.push_back("cstif");
		valid.push_back("gaf");
		valid.push_back("controls");
	}

	// check for an exact match with one of the allowable names
	if (find(valid.begin(), valid.end(), desc) != valid.end()) {
		return true;
	}
	trc.dprint("returning false: ",desc," is not a flut matrix");
	return false;
}

/*------------------------------------------------------------------
 * specs handler functions: names end in _f
 *------------------------------------------------------------------*/

static bool
checkjac_f(const Tok& opt) {
// handle the "checkjacobian" option:
//    checkjacobian{p1,p2,...} = nstep: check the jacobian every
// "nstep" steps, only columns corresponding to parameter names
//     "p1", "p2", ...
	int stepsper{5};  // default steps between checks

	if (!opt.ivec.empty())
		stepsper = opt.ivec[0];

	// check for parameter names, e.g.
	//    checkjac{freq,sigma,vtas} = 10
	vector<string> pnames;
	if (!opt.lopt.empty()) {
		vector<Tok*> prefopt = flaps::lexer(opt.lopt);
		for (size_t i=0; i<prefopt.size(); i++) {
			pnames.push_back(prefopt[i]->lhs);
		}
	}
	// set checkjacobian
	Settings::defaults.set("checkjacobian",stepsper);
	if (!pnames.empty()) {
		Settings::defaults.set("checkjacobian_par",pnames);
	}

	ostringstream os;
	os << "The Jacobian ";
	if (!pnames.empty()) {
		os << "(only " << pnames[0];
		for (size_t i=1; i<pnames.size(); i++)
			os << ", " << pnames[i];
		os << ") ";
	}
	os << "will be checked by differencing";
	os << " every " << stepsper << " steps";
	flaps::info (os.str());

	return true;
}

static bool
debug_f(const Tok& opt) {
/*
 * debugging options:
 * options may have either "check" or "debug" or "d" on the lhs
 * the rhs may be followed by the debug level in braces, e.g.
 *    debug = start{3}
 * which sets the debug level to 3 for the start process.
 * The rhs may be any of the names that are used in a statement
 * like
 *     Debug dbg("name");
 * The current list of these (found by "grep Debug *.c"):
 * 	Debug dbg("bifurcation");
 *  	Debug dbg("start");
 * 	Debug dbg("stepsize");
 * 	Debug dbg("constrain");
 * 	Debug dbg("checkjacobian");
 * 	Debug dbg("jacobian");
 * 	Debug dbg("homotopy");
 * 	Debug dbg("make_equations");
 * and one that works a little differently (see main.cpp):
 *    debug=fpexceptions
 * So for example the option
 *    d=start{3}
 * will set the debug level to 3 for the start process
 *
 */
	Trace trc(1,"debug");

	// the rhs may be followed by {n}: the debug level (default 2)
	int level{2};
	if (!opt.ropt.empty()) {
		str2int(opt.ropt, level);
	}

	// the rhs may be any of the routines that include a statement
	// like
	//    Debug dbg("fcn");
	// the rhs might be an integer: the global debug level
	if (!opt.ivec.empty()) {
		debug(opt.ivec[0]);
		flaps::info("Debug level ",debug());
	} else {
		// allow for multiple functions, e.g. debug=(start,stepsize)
		for (auto& fcn : opt.svec) {
			debug(fcn, level);
			if (fcn == "start")
				Settings::defaults.set("plot_homotopy", true);
		}
	}
	return true;
}

static bool
indep_f(const Tok& opt) {
	Trace<> trc("indep ",opt.svec.size()," independents");
	bool rval = true;

	if (opt.svec.empty()) {
		string exc = vastr("no independent parameters specified in ",opt);
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	// parse each rhs
	vector<Par*> indep;
	for (auto srhs : opt.svec) {
		Par* pp;
		try {
			// the rhs might be more than just a parameter name
			// so call the Par parsing constructor:
			pp = new Par(srhs, "");
			// change growth to sigma
			if (pp->name == "growth") {
				pp->name = "sigma";
				flaps::info ("changing indep \"growth\" to \"sigma\"");
			}
		} catch (runtime_error& s) {
			throw runtime_error(vastr("bad indep parameter definition (",srhs,"): ",s.what()));
		}

		// keep track of the independents to determine if this is VOE
		indep.push_back(pp);
		// ... and note that these were set according to the user's specs
		pp->pref = true;
		// adding it to gpset will change its state if it already exists
		pp = gpset::get().add(pp);
		// make the parameter Indep...
		pp->set_indep();

		// set the value to it's min
		if (pp->has_min())
			pp->value(pp->min());
	}
	// is this a VOE? used to automatically turn on lcostab
	if (indep.size() == 3) {
		bool v{false};
		bool f{false};
		bool e{false};
		for (auto pp : indep) {
			string nm = pp->name;
			if (nm == "vtas" || nm == "alt" || nm == "veas" || nm == "vcas"
				|| nm == "dpress" || nm == "cdpress" || nm == "rf")
				v = true;
			if (nm == "freq")
				f = true;
			if (nm == "gcnorm")
				e = true;
		}
		// if this is a VOE analysis turn on lcostab
		if (v && f && e)
			Settings::defaults.set("lcostab",true);
	}
	return rval;
}

static bool
matrix_f (const Tok& opt) {
// treat matrix input options, e.g. "mass=KHH", add the matrix to Matrix::collection,
// so it is accessible with, e.g.
//   Matrix* mass = Matrix::get_desc("mass");
	Trace trc(1,"matrix ",opt);
	ostringstream os;
	
	// Check for a standard matrix name (e.g. mass= )
	string descrip = opt.lhs;
	Matrix* mp{nullptr};
	if (is_flutmat(descrip) && !opt.svec.empty()) {
		// fetch might throw an exception - let someone else catch it
		mp = new Matrix(opt.svec[0]);
		if (mp == nullptr)
			throw runtime_error(vastr("matrix \"",opt.svec[0],"\" is not available"));
		mp->desc(descrip);
		Matrix::insert(mp);
		trc.dprint("returning true: got ",mp->summary());
		return true;
	}
	trc.dprint("returning false: not a flutter matrix");
	return false;
}

static bool
plot_f(const Tok& opt) {
/*
 * plot=(list)
 * Special items in "list":
 *   stepsize     include a number of parameters associated with tracking
 *                quality, e.g. stepsize, rcond, quality, etc.
 *   iterations   if the correction phase at a point fails, plot values
 *                at each corrector iteration
 *   homotopy     plot any homotopy curve (usually from normalModes)
 */
	Trace trc(1,"plot");
	ostringstream os;

	trc.dprint(opt.svec.size()," rhs");

	string rhs;
	vector<string> toplot;
	for (auto& rhs : opt.svec) {
		trc.dprint("got plot \"",rhs,"\"");


		// All special rhs's tried - must be a parameter. Create
		// a new parameter, add it to the gpset
		try {
			Par newpar(rhs, "0.0");
			gpset::get().add(&newpar);
			toplot.push_back(rhs);
		} catch (runtime_error& s)  {
			os.str("");
			os << "bad plot parameter definition (" << rhs
				<< "): " << s;
			trc.dprint("throwing exception: ",os.str());
			throw runtime_error(os.str());
		}
	}

	// if the user specified parameters to plot make them the default
	if (!toplot.empty()) {
		Settings::defaults.set("toplot", toplot);
		trc.dprint("now have ",toplot.size()," plot parameters");
	}

	return true;
}

static bool
plotfile_f(const Tok& opt) {
	if (opt.svec.empty())
		throw runtime_error(vastr("the ",opt.lhs," option has no rhs"));
	Settings::defaults.set("plotfile", opt.svec[0]);
	return true;
}

static bool
toprint_f(const Tok& opt) {
/*------------------------------------------------------------------
 * Parse the print option. The general form is
 *    print{options} = (rhs_keywords, parameter_names)
 * Print options:
 *    width    set output file width
 * rhs_keywords:
 *    matrices    print all matrices once
 *    gc          print summary of generalized coordinates at each soln
 *    vec           "
 *    tangent     print summary of tangent at each soln
 *    summary     only print summary of soln
 *    nosummary   print all solutions
 *    full
 *    all
 * Any other rhs is assumed to be a parameter name which will
 * be added to the list of parameters to print in the solution
 * summary
 *------------------------------------------------------------------*/
	Trace trc(1,"toprint");
	size_t i;
	ostringstream os;

	//  lhs options (print{prefopt} = ...): none currently
	//  XXX add width option
	if (!opt.lopt.empty()) {
		vector<Tok*> prefopt = flaps::lexer(opt.lopt);
		for (i=0; i<prefopt.size(); i++) {
			flaps::warning("unrecognized print option: ", prefopt[i]);
		}
	}

	// the rhs consists of a list of strings: keywords and parameter names
	string rhs;
	vector<string> toprint;
	// for (i=0; i<std::max((size_t)1, opt.svec.size()); i++) 
	for (auto& rhs : opt.svec) {
		if (rhs.substr(0,4) == "summ") {
			Settings::defaults.set("printsummary",true);
		} else if ((rhs.substr(0,5) == "nosum") || (rhs.substr(0,4) == "full") ||
				(rhs.substr(0,3) == "all")) {
			Settings::defaults.set("printsummary",false);
		} else if (rhs.substr(0,3) == "mat") {
			Settings::defaults.set("printmatrices",true);
		} else {
			// throw an exception if the parameter is undefined unless it
			// is an eigenvector component
			if (gpset::find(rhs) == nullptr && evparse(rhs) == -1) {
				string exc = vastr("parameter (",rhs,") has not been defined");
				throw runtime_error(exc);
			}

			toprint.push_back(rhs);

			trc.dprint(rhs," will be printed");
		}
	}

	// if the user specified parameters to print set them as the default
	if (!toprint.empty()) {
		Settings::defaults.set("toprint", toprint);
	}
	return true;
}

static bool
param_f(const Tok& opt) {
// all other keywords have been filtered.
// Check for a parameter declaration in the form
//		name(desc)[ min : max ]<conv PU/CU> = value
// where "value" is in PU, or
//		name(desc)[ min : max ]<conv PU/CU> = equation
// where "equation" is a string defining the equation to be used
// to evaluate the parameter. XXX the equation must be in
// PU, so this must be accounted for in adeqn
//
// If the user specifies a parameter value => the parameter is Fixed.
// If the user does not specify units they default to presentation
// (aka "external")
//
// XXX Allow the user to say something like "alt=Fixed" instead
// of specifying a value - this means to take the parameter value
// from a previous analysis
	Trace trc(1,"param");
	ostringstream os;

	trc.dprint("opt<",opt,"> ",opt.svec.size()," rhs");

	Par* par(nullptr);
	try {
		if (opt.svec.size() == 1 && opt.svec[0] == "fixed") {
			par = new Par(opt.lhs, "0.0");
			par->set_fixed();
			par->pref = true;
		} else {
			par = new Par(opt.lhs, opt.svec);
			par->pref = 1;
			if (opt.is_int() || opt.is_real() || opt.is_complex()) {
				par->set_fixed();
			}
			if (!par->equation.empty()) {
				par->set_derived();
			}
		}
		// some parameters are "known" and Aux
		if (par->name == "milepost" || par->name == "coord")
			par->set_aux();
	} catch (runtime_error& s) {
		string exc = vastr("illegal parameter definition: ",s.what());
		trc.dprint("throwing exception: ",exc);
		throw runtime_error(exc);
	}

	// add this parameter to gpset or upgrade the existing parameter
	// XXX should gpset::fetch set all states to nostate?
	Par* existing = gpset::find(par->name);
	if (existing == nullptr)
		throw runtime_error(vastr("parameter \"",par->name,"\" is undefined"));
	if (existing != nullptr && existing->is_indep())
		throw runtime_error(vastr(opt,": ",par->name,
			" was already declared ",existing->state));

	par = gpset::get().add(par);

	trc.dprint("got parameter ",par->longsummary());
	return true;
}
/*------------------------------------------------------------------
 * End of specs handler functions: names end in _f
 *------------------------------------------------------------------*/

static void
invert(vector<double>& P) {
// invert an (n,n) matrix, return it in the input P
	Trace trc(1,"invert");

	int n = sqrt(P.size());
	trc.dprintm(n,n,n,P,"inverting");
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
	trc.dprint("rcond ",rcond);
	trc.dprintm(n,n,n,x,"inverted");
	// check by multiplying x*P
	blas_sgemm("n", "n", n, n, n, 1.0, psave.data(), n,
			x.data(), n, 0.0, pf.data(), n);
	trc.dprintm(n,n,n,pf,"P*x");
	P = x;
}

void
parser (string const& prog, Specs& specs) {
	Trace<> trc("parser");

	// gcnorm defaults to Fixed zero
	Par* gcnorm = gpset::find("gcnorm");
	assert(gcnorm != nullptr);
	gcnorm->set_nostate();
	gcnorm->value(0.0);
	gcnorm->set_fixed();
	// sdamp defaults to Fixed 0.0
	Par* sdamp = gpset::find("sdamp");
	assert(sdamp != nullptr);
	sdamp->set_nostate();
	sdamp->value(0.0);
	sdamp->set_fixed();

	// first get the "aid" pref
	// parse these two options taking the input from cin (empty string)
	vector<Tok*> unrec = flaps::lexer("", {
		{"^(a)?id", [&](const Tok& p) {
			Settings::defaults.set("aid",p.svec[0]);
			// write the id on cerr
			cerr << "analysis id " << p.svec[0];
			return true;
		}}}
	);

	// ... then change the states of any indep par to derived
	for (auto pj : gpset::get().pmap()) {
		Par* pp = pj.second;
		if (pp->is_indep())
			pp->set_derived();
	}

	// parse the options taking the input from "unrec": the options
	// not recognized by previous calls to flaps::lexer
	unrec = flaps::lexer(unrec, {
			{"^checkjac(obian)?", checkjac_f},
			{"^d(ebug)?$", debug_f},
			{"^indep(endent)?", indep_f},
			{".*", matrix_f},
			{"maxstep", [&](const Tok& p) { specs.maxstep = p.ivec[0]; return true; }},
			{"^plotf(ile)?", plotfile_f},
			{"^plot$", plot_f},
			{"^print", toprint_f}
		});

	// ... and lastly check for parameter defn
	unrec = flaps::lexer(unrec, {{".*", param_f}});

	if (!unrec.empty()) {
		string exc{vastr(unrec.size()," options were not recognized: ",unrec)};
		flaps::warning(exc);
	}

	// get n, the order of mass, stif, etc...
	Matrix* mass = Matrix::find_desc("mass");
	if (mass == nullptr)
		throw runtime_error("no mass matrix specified");
	int n = mass->rsize();
	specs.mass = mass->data();

	specs.stif = Matrix::find_desc("stif");
	if (specs.stif == nullptr)
		throw runtime_error("no stif matrix specified");

	// ... and nbeta
	Matrix* gaf = Matrix::find_desc("gaf");
	if (gaf == nullptr)
		throw runtime_error("no gaf matrix specified");
	int na = gaf->rsize();
	if (na != n)
		throw runtime_error(vastr(gaf->mid()," is order ",na,", but ",
				mass->mid()," is ",n));
	if (gaf->pz.size() != 1)
		throw runtime_error(vastr(gaf->mid()," has ",gaf->pz.size()," parameterizations"));
	RFAPz* rfa = dynamic_cast<RFAPz*>(gaf->pz[0]);
	if (!rfa)
		throw runtime_error("gaf matrix is not a rational-function approximation (RFA)");
	specs.beta = rfa->beta;
	specs.R = rfa->R;
	// compute P = (M - \rho/2 R_2)^{-1} assuming rho is constant
	specs.P = mass->data();
	Par* rho = gpset::find("rho");
	double t = rho->value()/2.0;
	blas_axpy(n*n, t, specs.R[2].data(), 1, specs.P.data(), 1);
	invert(specs.P);

	return;
}
