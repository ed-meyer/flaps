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
#include "concurrent.h"
#include "conv.h"
#include "flutcurve.h"
#include "lexer.h"
#include "matrix.h"
#include "process.h"
#include "qr.h"
#include "settings.h"
#include "specs.h"
#include "startpts.h"
#include "trace.h"

using namespace std;

void
Specs::
update(const string& nm, double val) {
// update Cu2eu "nm" with val if it exists, otherwise add
// it to my vector<Cu2eu>
	T_(Trace trc(2,"Specs::update");)
	for (auto& ci : cu2eu) {
		if (ci.name == nm) {
			ci.value = val;
			T_(trc.dprint("updated cu2eu: ",ci);)
			return;
		}
	}
	cu2eu.push_back({nm, val});
	T_(trc.dprint("added new cu2eu: ",cu2eu[cu2eu.size()-1]);)
}

Specs&
flutspecs() {
	static Specs thespecs;
	return thespecs;
}

bool
get_source(string const& sid);

static void
e9n_create() {
// create parameters corresponding to the Step members & add them to gpset
	const vector<pair<string,string>>& par = Step::parnames;
	for (auto& pi : par) {
		string name = std::get<0>(pi);
		Par* pp = gpset::find(name);
		if (pp == nullptr) {
			Par newpar(vastr(name,"(",get<1>(pi),")"), "0");
			newpar.set_aux();
			pp = gpset::get().add(&newpar);
		}
	}
	// also add Pac::stepsize...
	Par* pp = gpset::find("stepsize");
	if (pp == nullptr) {
		Par newpar("stepsize","0");
		newpar.set_aux();
		pp = gpset::get().add(&newpar);
	}
	// ... projnorm
	pp = gpset::find("projnorm");
	if (pp == nullptr) {
		Par newpar("projnorm","0");
		newpar.set_aux();
		pp = gpset::get().add(&newpar);
	}
	// ... and extended convergence criterum
	pp = gpset::find("ecc");
	if (pp == nullptr) {
		Par newpar("ecc","0");
		newpar.set_aux();
		pp = gpset::get().add(&newpar);
	}
}
	
string
incr_mvf(int init) {
// Set multiple-value-fixed (mvf) parameter values in the gpset to
// the values for the "visit" set of its mvf parameter values.
// To initialize the incrementation, call with a starting argument,
// usually 0; with no argument init=-1
// The first time this function is called (with an argument) the mvf
// parameter values are set to their first value; subsequent calls
// increment the first parameter value through all multi-values,
// then set the second parameter to its second value and
// go through all of the first parameter values, etc.
// That is, if there are m multi-valued parameters, the values
// are treated as an m-dimensional array and the values
// are set "fortran-style": first index increasing most rapidly.
//
// Returns a string like cbba where the letters are the
// ordinal numbers of the mvf parameters
// (in this example there are four mvf parameters, the first
// is on it's 3rd mv value, the 2nd & 3rd on their second, and
// the 4th on its first.
//
// Note: this function is not threadsafe: contains statics
	T_(Trace trc(1,"incr_mvf");)
#ifdef NEVER // use get_mvf()
	vector<Par*> fixed = gpset::get().get_fixed();
#else // NEVER // use get_mvf()
	vector<Par*> mvfpar = gpset::get().get_mvf();
#endif // NEVER // use get_mvf()
	size_t i;
	ostringstream os;
	static int visit{0};
	string rval;

	T_(trc.dprint("init ",init,", visit ",visit);)

#ifdef NEVER // use get_mvf()
	if (fixed.empty())
#else // NEVER // use get_mvf()
	if (mvfpar.empty()) {
#endif // NEVER // use get_mvf()
		T_(trc.dprint("returning empty string: no fixed parameters");)
		return "";
	}

	// if init >= 0 => initialize
	if (init >= 0) {
		visit = init;
	}

	// Create a vector of sizes of each "altval" member of
	// the fixed parameters...
	vector<size_t> ld;
	size_t nmv{1};  // number of mv combinations
#ifdef NEVER // use get_mvf()
	vector<Par*> mvfpar;
	for (auto pp : fixed) {
		size_t sz = pp->altval.size();
		if (sz > 0) {
			mvfpar.push_back(pp);
			ld.push_back(sz);
			T_(trc.dprint("mvf parameter (",pp->name,") has ",sz," altval");)
			nmv *= sz;
		}
	}
#else // NEVER // use get_mvf()
	for (auto pp : mvfpar) {
		size_t sz = pp->altval.size();
		ld.push_back(sz);
		T_(trc.dprint("mvf parameter (",pp->name,") has ",sz," altval");)
		nmv *= sz;
	}
#endif // NEVER // use get_mvf()

	if (visit >= (int)nmv) {
		T_(trc.dprint("returning empty string: all mvf par used");)
		return "";
	}

	// Create a vector of indices into the "altval" member of
	// each mvf parameter...
	static vector<size_t> mvfindex;	// not threadsafe
	mvfindex = vec2mdim(visit, ld);
	T_(trc.dprintv(mvfindex,"mvfindex");)
	// ... then set each mvf parameter to that value
	for (i=0; i<mvfpar.size(); i++) {
		Par* pp = mvfpar[i];
		double pi = pp->altval[mvfindex[i]];
		T_(trc.dprint("set ",pp->name," to ",pi," @",pp);)
		// use valuef to force the change; note this does not change
		// the "constant" member
		pp->valuef(pi);
	}

	// create the return string
	os.str("");
	for (auto i : mvfindex)
		os << char('a'+i);
	rval = os.str();

	// mark all parameters and matrices out-of-date
	gpset::get().touch();

	// increment visit for the next call
	visit++;
	T_(trc.dprint("returning ",rval);)
	return rval;
}


bool
is_flutmat (string const& desc) {
/*------------------------------------------------------------------
 * Is "desc" one of the allowable matrices for flutter runs?
 * e.g. "mass", "stif", etc.
 * Returns the matched description if found, empty string otherwise
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"is_flutmat");)

	T_(trc.dprint("is ",desc," a flutter matrix?");)

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
	T_(trc.dprint("returning false: ",desc," is not a flut matrix");)
	return false;
}

ostream&
operator<<(ostream& s, const Cu2eu& t) {
	s << t.name << ": " << t.value;
	return s;
}

ostream&
operator<<(ostream& s, const Specs& t) {
	s << "aid " << t.aid << endl;
	s << "sid " << t.sid << endl;
	s << "checkjacobian " << t.checkjacobian << endl;
	s << "checkjacobian_par " << t.checkjacobian_par << endl;
	s << "checklooping " << t.checklooping << endl;
	s << "compare_det " << t.compare_det << endl;
	s << "constrained " << t.constrained << endl;
	if (!t.cu2eu.empty())
		s << "cu2eu: " << t.cu2eu << endl;
	s << "dmatrix " << t.dmatrix << endl;
	s << "free_vibration " << t.free_vibration << endl;
	s << "gmethod " << t.gmethod << endl;
	s << "lcostab " << t.lcostab << endl;
	s << "linearize " << t.linearize << endl;
	s << "nlev " << t.nlev << endl;
	s << "plotfile " << t.plotfile << endl;
	s << "plot_homotopy " << t.plot_homotopy << endl;
	s << "printmatrices " << t.printmatrices << endl;
	s << "printsummary " << t.printsummary << endl;
	s << "process_targets " << t.process_targets << endl;		// name of a Custom function
	s << "project " << t.project << endl;
	s << "projdir " << t.projdir << endl;
	s << "reflen " <<  t.reflen << endl;
	s << "rflimits " <<  t.rflimits << endl;
	s << "sequential " << t.sequential << endl;
	s << "stopatdsc " << t.stopatdsc << endl;
	s << "svd " << t.svd << endl;
	s << "e9n " << t.e9n << endl;
	s << "targets " << t.targets << endl;
	s << "title " << t.title << endl;
	s << "toprint " << t.toprint << endl;
	s << "vzid " << t.vzid << endl;
	return s;
}

/*------------------------------------------------------------------
 * Parser functions: names end in _f
 *------------------------------------------------------------------*/
static bool
checkjac_f(const Tok& opt) {
// handle the "checkjacobian" option:
//    checkjacobian{p1,p2,...} = nstep: check the jacobian every
// "nstep" steps, only columns corresponding to parameter names
//     "p1", "p2", ...
	int stepsper{5};  // default steps between checks
	Specs& sp = flutspecs();

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
	sp.checkjacobian = stepsper;
	sp.checkjacobian_par = pnames;

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
vacm_f(double mach_value) {
// treat a shorthand, either
//   vacm=mach_number
// or
//   cmvd=mach_number     (old syntax)
// by setting up to do a "variable alt-constant Mach" solution
	Par* pp = gpset::find("mach");
	assert(pp != nullptr);
	pp->value(mach_value);
	pp->set_fixed();
	// independents are alt, freq, sigma
	for (auto pi : {"alt", "freq", "sigma"}) {
		pp = gpset::find(pi);
		assert(pp != nullptr);
		pp->set_indep();
	}
	// add veas to print and plot
	Specs& sp = flutspecs();
	add2vector(string("veas"), sp.toprint);
	return true;
}

static bool
vvca_f(double alt_value) {
// treat a shorthand
//   vvca=alt
// by setting up to do a "variable vtas-constant alt" solution within
// the limits of rf
	Par* pp = gpset::find("alt");
	assert(pp != nullptr);
	pp->value(alt_value);
	pp->set_fixed();
	// independents are rf, freq, sigma
	for (auto pi : {"rf", "freq", "sigma"}) {
		pp = gpset::find(pi);
		assert(pp != nullptr);
		pp->set_indep();
	}
	// add vtas to print and plot
	Specs& sp = flutspecs();
	add2vector(string("vtas"), sp.toprint);
	return true;
}

static bool
dmatrix_f (const Tok& opt) {
// assign a custom function to be called after the
// dynamic matrix has been evaluated to modify the dynamic
// matrix
	flutspecs().dmatrix = opt.srhs;
	flaps::info("custom function \"",opt.srhs,"\" will be called"
		" to modify the dynamic matrix");
	return true;
}

static bool
indep_f(const Tok& opt) {
	T_(Trace trc(1,"indep ",opt.svec.size()," independents");)
	bool rval = true;

	if (opt.svec.empty()) {
		string exc = vastr("no independent parameters specified in ",opt);
		T_(trc.dprint("throwing exception: ",exc);)
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
		} catch (runtime_error& s) {
			throw runtime_error(vastr("bad indep parameter definition (",srhs,"): ",s.what()));
		}

		// keep track of the independents to determine if this is VOE
		indep.push_back(pp);
		// ... and note that these were set according to the user's preference
		pp->pref = true;
		// adding it to gpset will change its state if it already exists
		pp = gpset::get().add(pp);
		// make the parameter independent...
		pp->set_indep();

		// set the value to it's min...
		if (pp->has_min())
			pp->value(pp->min());
		// ... unless it is alt or rf: then max
		if (pp->name == "alt" || pp->name == "rf") {
			if (pp->has_max())
				pp->value(pp->max());
		}
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
			flutspecs().lcostab = true;
	}
	return rval;
}

static bool
matrix_f (const Tok& opt) {
// treat matrix input options, e.g. "mass=KHH"
	T_(Trace trc(1,"matrix ",opt);)
	ostringstream os;
	Specs& sp = flutspecs();
	
	// Check for a standard matrix name (e.g. mass= )
	string descrip = opt.lhs;
	Matrix* mp{nullptr};
	if (is_flutmat(descrip) && !opt.svec.empty()) {
		// fetch might throw an exception - let someone else catch it
		mp = new Matrix(opt.svec[0]);
		if (mp == nullptr)
			throw runtime_error(vastr("matrix \"",opt.svec[0],"\" is not available"));
		mp->desc(descrip);
		// add it to the collection if source was not included
		if (sp.sid.empty())
			Matrix::insert(mp);
		T_(trc.dprint("returning true: got ",mp->summary());)
		return true;
	}
	T_(trc.dprint("returning false: not a flutter matrix");)
	return false;
}

static bool
optimize_f (const Tok& opt) {
// optimize|decrease|increase = parameter_name
// or
// optimize{decrease|increase} = parameter_name
// or
// -------------------------------
// optimize{vtas=0.5, gcnorm=-0.8}
// -------------------------------
// optimize{increase=vtas{0.5}, decrease=gcnorm{0.8}}
// or
// increase{vtas=0.4, gcnorm=0.1, sigma=-0.4}
	T_(Trace trc(1,"optimize");)
	Specs& sp = flutspecs();

	// parse the left-hand side options: split into Toks by
	// calling parse with no Parsers
	// check for optimize{par[=dir]}
	vector<string> parnames;
	vector<double> pardir;
	if (!opt.lopt.empty()) {
		vector<Tok*> prefopt = flaps::lexer(opt.lopt);
		for (auto op : prefopt) {
			parnames.push_back(op->lhs);
			if (!op->rvec.empty())
				pardir.push_back(op->rvec[0]);
			else
				pardir.push_back(1.0);
		}
	} else {
		// treat an option like
		//    optimize=vtas
		parnames.push_back(opt.srhs);
		// both directions
		pardir.push_back(1.0);
	}

	// add new parameters & directions
	for (size_t i=0; i<parnames.size(); i++) {
		sp.project.push_back(parnames[i]);
		sp.projdir.push_back(pardir[i]);
	}
	T_(trc.dprint("projection parameters: ",sp.project);)
	T_(trc.dprint("projection directions: ",sp.projdir);)

	return true;
}

static bool
start_f(const Tok& opt) {
/*------------------------------------------------------------------
 * Start region:
 *   start{freq[0:4]}
 * or
 *   start{freq[0:4]=0.1, sigma[-1:1]=0.1}, ...
 * or
 *   start{freq[0:4]=4, sigma[-1:1]=5 }, ...
 * or
 *   start{vtas=100}
 * or
 *   start{vtas=minimum}
 * or
 *   start{modes=(implicit_list)}
 * or
 *   start{mode_1[0-9].a, ...}  regular-expression curveid
 * or
 *   start{coord=122, ... }     coord number
 * or
 *   start{milepost=3.1415,...} distance along the curve
 * or
 *   start{ordinal=3}           3rd place on the curve where a fixed par
 *                              occurs; e.g. the 3rd flutter crossing in
 *                              a particular mode
 * or
 *   start{sigma=closest{0}}
 * See region.cpp:Region(Tok) constructor for precise parsing options
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"start_f");)

	T_(trc.dprint("Tok: ",opt);)

	// all options must be in the lopt
	if (opt.lopt.empty()) {
		string exc = vastr("illegal startregion: ",opt);
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}
	// add these options to the startregions
	// Let Region parse the option
	startregions(opt.lopt);

	return true;
}

static bool
target_f(const Tok& p) {
// Create one or more Targets from a set of options, usually
// the sub-options (enclosed by braces) to a target option, e.g.
//    target{growth=(0,.03), sigma=1, vtas[10:100]}
// which creates three targets (growth=0, growth=0.03,
// and sigma=1) each limited to the interval vtas[10:100].
// This function can be called multiple times - Targets are
// accumulated in flutspecs().
	T_(Trace trc(1,"target_f");)
	vector<Window> windows;
	string sortpar;
	Specs& sp = flutspecs();

	if (p.lopt.empty())
		throw runtime_error(vastr("the \"target\" option must have the form ",
			"(for example) \"target{freq=(0,3), vcas=233, vtas[20:]}\""));

	vector<Target*> newtargets;	// all the targets created here

	// break the lhs option list into Tok*
	vector<Tok*> opt = flaps::lexer(p.lopt, {
		{"^sort$", [&](const Tok& p) { sortpar = p.srhs; return true; }},
		{"^process$", [&](const Tok& p) { sp.process_targets = p.srhs; return true; }},
		{".*", [&](const Tok& p) {
			if (!p.svec.empty()) {
				// has rhs: must be a target parameter and value(s); create a
				// Target for each value on the rhs. Watch for par=closest{v}
				for (auto& rhs : p.svec) {
					if (rhs == "closest")
						newtargets.push_back(new Target(p.lhs, p.ropt, true));
					else
						newtargets.push_back(new Target(p.lhs, rhs, false));
				}
			} else {
				// no rhs: must be a Window, e.g. vtas[20:50]
				Par pp(p.lhs, "0");
				Par* stdp = gpset::find(pp.name);
				if (stdp == nullptr)
					throw runtime_error(vastr("target parameter \"",pp.name,
							"\" has not been defined"));
				// convert limits from PU to EU
				double cf = stdp->conv.factor;
				windows.push_back({pp.name,pp.min(),pp.max(),cf});
			}
			return true;
		}}
	});

	// Add windows to each new target...
	if (!windows.empty()) {
		for (auto t : newtargets)
			t->windows = windows;
	}
	
	// ... and sort parameter
	if (!sortpar.empty()){
		for (auto t : newtargets)
			t->sortpar = sortpar;
	}
	// now add these new targets to Specs
	for (auto t : newtargets)
		sp.targets.push_back(t);

	return true;
}

static bool
toprint_f(const Tok& opt) {
/*------------------------------------------------------------------
 * Parse the print option. The general form is
 *    print{options} = (rhs_keywords, parameter_names)
 * Print options:
 *    width    set output file width
 *    summary     only print summary of soln		XXX move to here
 *    nosummary   print all solutions				XXX move to here
 * rhs_keywords:
 *    matrices    print all matrices once
 *    gc          print summary of generalized coordinates at each soln
 *    vec           "
 *    tangent     print summary of tangent at each soln
 *    summary     only print summary of soln		XXX from here
 *    nosummary   print all solutions				XXX from here
 *    full
 *    all
 * Any other rhs is assumed to be a parameter name which will
 * be added to the list of parameters to print in the solution
 * summary
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"toprint");)
	size_t i;
	ostringstream os;
	Specs& sp = flutspecs();

	//  lhs options (print{prefopt} = ...): none currently
	//  XXX add width option
	if (!opt.lopt.empty()) {
		vector<Tok*> toks = flaps::lexer(opt.lopt);
		for (i=0; i<toks.size(); i++) {
			flaps::warning("unrecognized print option: ", toks[i]);
		}
	}

	// the rhs consists of a list of strings: keywords and parameter names
	string rhs;
	// for (i=0; i<std::max((size_t)1, opt.svec.size()); i++) 
	for (auto& rhs : opt.svec) {
		if (rhs.substr(0,4) == "summ") {
			sp.printsummary = true;
		} else if ((rhs.substr(0,5) == "nosum") || (rhs.substr(0,4) == "full") ||
				(rhs.substr(0,3) == "all")) {
			sp.printsummary = false;
		} else if (rhs.substr(0,3) == "mat") {
			sp.printmatrices = true;
		} else if (gpset::find(rhs) != nullptr) {
			add2vector(rhs, sp.toprint);
		} else {
			// throw an exception if the parameter is undefined unless it
			// is an eigenvector component: name starting with "ev", "gc", or "eig.."
			if (rhs != "ev" && rhs != "gc" && rhs != "eigenvector") {
				if (evparse(rhs) == -1) {		// not an eigenvector component?
					return false;
					//!! throw runtime_error(vastr("parameter (",rhs,") is undefined"));
				} else
					add2vector(rhs, sp.toprint);
					//!! sp.toprint.push_back(rhs);
			} else {
				add2vector(rhs, sp.toprint);
				//!! sp.toprint.push_back(rhs);
				T_(trc.dprint(rhs," will be printed");)
			}
		}
	}
	return true;
}

static bool
pac_f(const string& opt) {
// set some Pac parameters in Flutspecs; the parameters in Pacspecs
// with default values (see class Pacspecs in pac.h):
//	double curvaturefac{4.0};	// increase curvature, decrease stepsize
//	int ecc{0};						// -1/0/1: never/default/always use minf
//	double ecc_eps{1.0e-13};	// used in ecc()
//	double epsabs{0.0};				// error tolerance for fxnorm
//	double growthfac{3.0};			// increase stepsize time secant if distance==0
//	double initial_stepsize{0.001};	// if hk=0 in trace
//	double maxangle{0.5235};		// angle between tans: 30 deg
//	double maxcurvature{100.0};		// max curvature
//	int maxiter{10};				// Newton - use 20 for Modified Newton
//	int maxsteps{1000};
//	double maxstepsize{0.2};
//	std::vector<double> minf;		// smallest possible value for f_j j=1,nf
//	double minstepsize{1.0e-12};
//	double reductionfac{3.0};

	Specs& sp = flutspecs();
	Pacspecs pp;

	vector<Tok*> unrec = flaps::lexer(opt, {
		{"^curvaturefac", [&](const Tok& p) { pp.curvaturefac = p.rvec[0]; return true; }},
		{"^ecc$", [&](const Tok& p) { pp.ecc = p.ivec[0]; return true; }},
		{"^ecc_eps$", [&](const Tok& p) { pp.ecc_eps = p.rvec[0]; return true; }},
		{"^epsabs$", [&](const Tok& p) { pp.epsabs = p.rvec[0]; return true; }},
		{"^growthfac$", [&](const Tok& p) { pp.growthfac = p.rvec[0]; return true; }},
		{"^initial_stepsize$", [&](const Tok& p) { pp.initial_stepsize = p.rvec[0]; return true; }},
		{"^maxangle$", [&](const Tok& p) { pp.maxangle = p.rvec[0]; return true; }},
		{"^maxcurvature$", [&](const Tok& p) { pp.maxcurvature = p.rvec[0]; return true; }},
		{"^maxiter$", [&](const Tok& p) { pp.maxiter = p.ivec[0]; return true; }},
		{"^maxsteps$", [&](const Tok& p) { pp.maxsteps = p.rvec[0]; return true; }},
		{"^maxstepsize$", [&](const Tok& p) { pp.maxstepsize = p.rvec[0]; return true; }},
		{"^minstepsize$", [&](const Tok& p) { pp.minstepsize = p.rvec[0]; return true; }},
		{"^reductionfac$", [&](const Tok& p) { pp.reductionfac = p.rvec[0]; return true; }},
	});
	if (!unrec.empty())
		throw runtime_error(vastr("unrecognized Pac spec(s): ",unrec));
	sp.pacspecs = new Pacspecs(pp);
	return true;
}

static bool
param_f(const Tok& opt) {
// all other keywords have been filtered.
// Check for a parameter declaration in the form
//		name(desc)[ min : max ]<conv PU/CU> = value
// or
//		name(desc)[ min : max ]<conv PU/CU> = (implicit list)
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
	T_(Trace trc(1,"param_f");)
	ostringstream os;
	Specs& sp = flutspecs();

	T_(trc.dprint("opt<",opt,"> ",opt.svec.size()," rhs");)

	Par* par(nullptr);
	try {
		if (opt.svec.size() == 1 && opt.svec[0] == "fixed") {
			par = new Par(opt.lhs, "0.0");
			par->set_fixed();
			par->pref = true;
		} else {
			par = new Par(opt.lhs, opt.svec);
			par->pref = 1;
			// if it only has limit (no rhs) set limits targets
			if (opt.svec.empty()) {
				if (par->has_min())
					sp.targets.push_back(new Target(par->name, par->min(), false, "lower"));
				if (par->has_max())
					sp.targets.push_back(new Target(par->name, par->max(), false, "upper"));
			} else {
				if (opt.is_int() || opt.is_real() || opt.is_complex())
					par->set_fixed();
				if (!par->equation.empty())
					par->set_derived();
			}
		}
		// some parameters are "known" and Aux
		if (par->name == "milepost" || par->name == "coord")
			par->set_aux();
	} catch (runtime_error& s) {
		string exc = vastr("illegal parameter definition: ",s.what());
		T_(trc.dprint("throwing exception: ",exc);)
		throw runtime_error(exc);
	}

	// add this parameter to gpset or upgrade the existing parameter
	// XXX should gpset::fetch set all states to nostate?
	Par* existing = gpset::find(par->name);
	if (existing == nullptr)
		return false;
		//!! throw runtime_error(vastr("parameter \"",par->name,"\" is undefined"));
	if (existing->is_indep() && !opt.svec.empty())
		throw runtime_error(vastr(opt,": ",par->name,
			" was already declared ",existing->state));

	par = gpset::get().add(par);

	T_(trc.dprint("got parameter ",par->longsummary());)
	return true;
}
/*------------------------------------------------------------------
 * End of parser functions: names end in _f
 *------------------------------------------------------------------*/

void
parser (string const& prog) {
// Parse the user's options by calling flaps::lexer with handler functions
	T_(Trace trc(1,"parser");)

	Specs& sp = flutspecs();

	// set some default cu2eu's
	sp.cu2eu.push_back({"alt",10000.0});
	sp.cu2eu.push_back({"sigma",10.0});

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

	// first get the "aid" and "source" prefs...
	// parse these two options taking the input from cin (empty string)
	vector<Tok*> unrec = flaps::lexer("", {
		{"^(a)?id", [&](const Tok& p) { sp.aid = p.srhs; return true; }},
		{"^so(urce)?", [&](const Tok& p) { sp.sid = p.srhs; return true; }},
	});

	// read source data
	if (!sp.sid.empty())
		get_source(sp.sid);

	// parse the remaining options taking the input from "unrec", the
	// options not recognized by previous calls to flaps::lexer
	unrec = flaps::lexer(unrec, {
		{"^bif(urcation)?", [&](const Tok& p) { return sp.bifurcation = true; }},
		{"^bifstep", [&](const Tok& p) { sp.bif_init = p.rvec[0]; return true; }},
		{"^checkjac(obian)?", checkjac_f },
		{"^checklooping", [&](const Tok& p) { return sp.checklooping = true; }},
		{"^compare_det", [&](const Tok& p) { return sp.compare_det = true; }},
		{"^constrained", [&](const Tok& p) { return sp.constrained = true; }},
		{"^cu2eu", [&](const Tok& p) {
			if (p.lopt.empty()) {
				flaps::warning("the cu2eu option has the form \"cu2eu{p1=s1,..}\"");
				return false;
			}
			flaps::lexer(p.lopt, {
				{".*",[&](const Tok& q) {
					sp.update(q.lhs, q.rvec[0]);
					return true; }}});
			return true;
		}},
		//!! {"^d(ebug)?$", debug_f },
		{"^decrease$", optimize_f },
		{"^dmatrix$", dmatrix_f },
		{"^ecc$", [&](const Tok& p) {
			if (p.srhs == "no" || p.srhs == "off" || p.srhs == "false")
				sp.ecc = -1;
			else if (p.svec.empty())
				sp.ecc = 1;
			else
				return false;
			return true;
		}},
		{"^freevib(ration)?", [&](const Tok& p) {return sp.free_vibration = true; }},
		{"^gmethod", [&](const Tok& p) {return sp.gmethod = true; }},
		{"^increase$", optimize_f },
		{"^indep(endent)?", indep_f },
		{"^lcostab(ility)?", [&](const Tok& p) {return sp.lcostab = true; }},
		{"^linearize", [&](const Tok& p) {
			Settings::defaults.set("linearize",true);
			return sp.linearize = true; }},
		{".*", matrix_f },  // matrix name?
		{"^opt(imize)?", optimize_f },
		{"^pac", [&](const Tok& p) {return pac_f(p.lopt); }},
		{"^plotfile", [&](const Tok& p) {sp.plotfile = p.srhs; return true; }},
		{"^plot$", [&](const Tok& p) {
			// only plot options are "stepsize" & "homotopy" ("iterations"?)
			for (auto pi : p.svec) {
				if (pi == "stepsize") {
					e9n_create();
					sp.e9n = true;
				} else if (pi == "homotopy")
					sp.plot_homotopy = true;
				else
					return false;
			}
			return true;
		}},
		{"^reflen$", [&](const Tok& p) {sp.reflen = p.rvec[0]; return true; }},
		{"^rflimits$", [&](const Tok& p) { return sp.rflimits = true; }},
		{"^seq(uential)?", [&](const Tok& p) {
			concurrently(0);
			return sp.sequential = true; }},
		{"^start(region)?", start_f },
		{"^svd", [&](const Tok& p) { QR::usesvd(1); return true; }},
		{"^tar(get)?", target_f },
		{"^title", [&](const Tok& p) {sp.title = p.srhs; return true; }},
		{"^print", toprint_f },
		{"^vacm$|^cmvd$", [&](const Tok& p) { return vacm_f(p.rvec[0]); }},
		{"^vvca$", [&](const Tok& p) { return vvca_f(p.rvec[0]); }},
		{"^vzid|amkit", [&](const Tok& p) {
			T_(Trace trc(1,"vzid");)
			sp.vzid = p.srhs;
			T_(trc.dprint("got vzid \"",sp.vzid,"\"");)
			return true; }},
		{".*", param_f },    // parameter name?
	});
		
	if (!unrec.empty())
		throw runtime_error(vastr(unrec.size()," options were not recognized: ",unrec));

	T_(trc.dprint("specs: {\n",sp,"}");)
	return;
}

bool
get_source(string const& sid) {
// set up to run from a previous flut run with id "sid"
// Fetch the parameters used in the source analysis and
// 1) mark them as coming from the source analysis,
// 2) change the state from "indep" to "derived"
// 3) update the global parameter with it, mainly to set fixed states
// 4) leave fixed parameters fixed
// Most will already be in the global pset but some members may be
// taken from the source parameter. Keeping a separate pset in
// addition to the gpset allows passing Fixed states, and ??.
	T_(Trace trc(1,"get_source");)

	string mid{vastr("parameters,",sid)};
	try {
		pset* source = pset::fetch(mid);
		T_(trc.dprint("gpset before updating with source parameters:\n",gpset::get());)
		for (auto& pairi : source->pmap()) {
			Par* pp = pairi.second;
			Par* gpsp = gpset::find(pp->name);
			if (gpsp == nullptr) continue;  // XXX add it to gpset?
			// if Fixed change gpsp to Fixed with pp value
			if (pp->is_fixed()) {
				gpsp->valuef(pp->value());
				gpsp->set_fixed();
				T_(trc.dprint("set ",pp->name," to Fixed ",pp->value());)
			}
			// if Indep change gpsp to Derived if it has candidates,
			// nostate otherwise
			if (pp->is_indep()) {
				if (!gpsp->candidates.empty())
					gpsp->set_derived();
				else
					gpsp->set_nostate();
			}
			// Clear the "solns" arrays
			gpsp->clear_solns();
		}
		T_(trc.dprint("gpset after updating with source parameters:\n",gpset::get());)
	} catch (runtime_error& s) {
		string exc{vastr("source analysis (",sid,") is not available: ",s.what())};
		throw runtime_error(exc);
	}

	// matrices used in source analysis
	try {
		vector<pair<string,string> > mids = Matrix::fetch_mids(sid);
		if (mids.empty()) {
			string exc{vastr("data from source analysis (",sid,
				") is incomplete: matrix list is missing")};
			T_(trc.dprint("throwing exception: ",exc);)
			throw runtime_error(exc);
		}
		for (auto& s : mids) {
			string des = s.first;
			string mid = s.second;
			Matrix* mp = new Matrix(mid);
			if (mp == nullptr)
				throw runtime_error(vastr("matrix \"",mid,"\" was not saved in ",sid));
			mp->desc(des);
			Matrix::insert(mp);
			T_(trc.dprint("fetched source matrix ",mid);)
		}
	} catch (runtime_error& s) {
		throw runtime_error(vastr("data from source analysis (",sid,
			") is incomplete: ", s.what()));
	}
	return true;
}
