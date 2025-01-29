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
#include "fio.h"
#include "matrix.h"
#include "pset.h"
#include "process.h"
#include "specs.h"
#include "startpts.h"  // for startregions()
#include "target.h"
#include "trace.h"

using namespace std;

static bool
post_preferences_defaults(string const& prog);
static void
check_param();
static void
print_prefs();
static pset*
seek_target(vector<Flutcurve*>& curves, Target* target, vector<string>& leaders);
double
precession(pset& plt, Matrix* gyro);

void
pre_process(string const& prog) {
// do some pre-processing assuming the user's preferences have
// been parsed::
//   1) create parameters for each eigenvector component
//   1b) if the user did not specify a D matrix, create a default
//   2) set up evaluation equations for all Derived Par
//   3) check the dependent parameter names for each Matrix and save
//      any that are eigenvector components to add to the AD derivatives.
//   4) set the automatic differentiation derivative names
//   5) re-allocate each gpset parameter's autodiff members to get the
//      correct autodiff parameters
//   6) set some defaults & print some user info
//   7) create targets for independent variable limits
	Trace trc(1,"pre_process");
	vector<Matrix*>& stdmat = Matrix::collection();
	ostringstream os;
	Specs& sp = flutspecs();

#ifdef NEVER // do this in Flutcurve::plot
	// did the user specify a vzid?
	Matrix* gct = Matrix::find_desc("gct");
	if (gct == nullptr) {
		try {
			Matrix* nodes{nullptr};
			Matrix* coords{nullptr};
			Matrix* conn{nullptr};
			Matrix* gct{nullptr};
			// vzmatrices will throw an exception of any are missing
			Matrix::vzmatrices(sp.vzid, nodes, coords, conn, gct);
			Matrix::insert(nodes);
			Matrix::insert(coords);
			Matrix::insert(conn);
			Matrix::insert(gct);
		} catch (const std::exception& s) {
			flaps::warning("cannot save amvz data in the plot file - "
				"do you need to specify vzid?");
		}
	}
#endif // NEVER // do this in Flutcurve::plot

	// check that the user specified at least a mass & stif
	// matrix or the "source" option
	if (stdmat.empty()) {
		os << "no input matrices were specified";
		trc.dprint("throwing exception: ",os.str());
		throw runtime_error(os.str());
	}
	// Mass, stiffness are mandatory, missing aerodynamics
	// just gets a warning
	if (Matrix::find_desc("mass") == nullptr) {
		throw runtime_error("no mass matrix was specified");
	}
	if (Matrix::find_desc("stif") == nullptr) {
		throw runtime_error("no stiffness matrix was specified");
	}
	if (Matrix::find_desc("gaf") == nullptr) {
		flaps::warning("no GAF matrix was specified");
	}
	// create parameter "precession" if a gyro matrix was included
	// XXX mv to pset.cpp, let dependson() eliminate it?
	if (Matrix::find_desc("gyro") != nullptr) {
		Par precession("precession(+/- fwd/backward)","0.0");
		Par* pp = gpset::get().add(&precession);
		pp->set_aux();
	}

	// 1) create a parameter for each real/imag part of each
	//    eigenvector component; to do this it is necessary
	//    to determine the matrix order (get_order())
	size_t ncmplx = get_order();
	gpset::get().make_eigv(ncmplx);

	// 1b) create a D matrix with a default eval function if necessary
	Matrix* dp = Matrix::find_desc("dmatrix");
	if (dp == nullptr)
		dp = new Matrix("dmatrix","dmatrix", ncmplx, ncmplx, true);

#ifdef NEVER // custom dmatrix needs work
	// delete all parameterizations
	dp->pz.clear();
	// add the default parameterization function
	dp->pz.push_back(new CustomPz(default_dmatrix, ncmplx, ncmplx));
	//!! string froot{getEnv("FROOT")};
	//!! paths.push_back(vastr(froot,"/share/flut/default_dmatrix.cpp"));
	// user-written function also?
	//!! string dmatrix_s;
	//!! Settings::defaults.get("dmatrix", dmatrix_s);
	if (!dmatrix.empty())
		dp->pz.push_back(new CustomPz(dmatrix, ncmplx, ncmplx));
	Matrix::insert(dp);
	dp->store();
#endif // NEVER // custom dmatrix needs work

	// set any point start Bndrys in the first Region so that those
	// parameters get changed to Fixed if they are not already Indep
	vector<Region>& srlist = startregions();
	if (!srlist.empty())
		setPointBndrys(srlist[0]);

	// 2) assign equations to derived parameters
	gpset::get().equations();

	// 3) get a list of parameters the dynamic matrix depends on;
	//    if a parameter has not been declared Fixed or Indep and has
	//    no equations, set it to Fixed.
	vector<string> dep = Flutcurve::dependson(gpset::get());
	trc.dprint("dmatrix dependencies: ",dep);

	// get the nonlinear eigenvector components for automatic differentiation
	// and check for "nostate"
	for (auto di : dep) {
		Par* pp = gpset::find(di);
		if (pp == nullptr) {
			flaps::warning(di," (used by D matrix) has not been defined");
		} else {
			if (pp->is_nostate()) {
				if (pp->equation.empty())
					pp->set_fixed();
				else
					pp->set_derived();
				trc.dprint("set ",pp->name," to ",pp->state);
			}
			// add to the list of nonlinear eigenvector components
			if (pp->is_eigv())
				sp.nlev.push_back(di);
		}
	}

	// 4) make a list of the parameter names to be used as AD derivatives
	vector<string> adnames;
	//    indep parameters ...
	vector<Par*> indep = gpset::get().get_indep();
	for (auto pp : indep)
		adnames.push_back(pp->name);
	//    ... nonlinear ev components and gcnorm, unless gcnorm is 0 (fixed)...
	Par* gcnormp = gpset::find("gcnorm");
	if (!(gcnormp->is_fixed() && gcnormp->value() == 0.0)) {
		for (auto nl : sp.nlev)
			adnames.push_back(nl);
		if (!sp.nlev.empty() &&
			find(adnames.begin(), adnames.end(), "gcnorm") == adnames.end())
			adnames.push_back("gcnorm");
	}

	// ... if this is a VOE analysis or if the lcostab option was included
	// add sigma to AD derivatives, also create parameter lcostab.
	if (sp.lcostab) {
		if (gpset::find("lcostab") == nullptr) {
			// create new parameter
			Par pp("lcostab(neg:stable)", "0");
			pp.set_aux();
			pp.pref = true;
			gpset::get().add(&pp);
			trc.dprint("VOE process: added lcostab");
		}
		// turn on lco stability calc
		adnames.push_back("sigma");
		flaps::info("LCO stability will be included in the plot file as \"lcostab\"");
	}

	// create parameter "coord" if not already there
	Par* coord = gpset::find("coord");
	if (coord == nullptr) {
		Par pp("coord(Steps from origin)","0");
		pp.set_aux();
		pp.pref = true;
		gpset::get().add(&pp);
	}
	// create parameter "detsign" if sp.bifurcation & not already there
	if (sp.bifurcation) {
		Par* detsign = gpset::find("detsign");
		if (detsign == nullptr) {
			Par pp("detsign(determinant sign)","0");
			pp.set_aux();
			pp.pref = true;
			gpset::get().add(&pp);
		}
	}
	// create parameter "rfrange": use w=rfrange in vz to show rf out-of-range
	Matrix* gaf = Matrix::find_desc("gaf");
	Par* rf = gpset::find("rf");
	if (gaf != nullptr && rf != nullptr) {
		Par* rfrange = gpset::find("rfrange");
		if (rfrange == nullptr && (rf->has_max() && rf->has_min())) {
			string rhs{vastr("(",rf->max(),"-rf)*(rf-",rf->min(),")")};
			Par pp("rfrange(neg: rf out of range)",rhs);
			rfrange = gpset::get().add(&pp);
		}
		if (rfrange != nullptr) {
			rfrange->set_derived();
			rfrange->value(0.0);
			rfrange->pref = 1;
			rfrange->fresh = true;
		}
	}

	// 5) set the ad parameter names; this will lock the names to
	//    future additions, reallocate all parameter autodiff members,
	//    and set derivatives to 1 in the proper locations
	Ad::initialize(adnames);
	gpset::get().realloc();

	trc.dprint("AD parameters: ",Ad::toString());

	 // 6) do some checking of the parameters
	check_param();

	// 7) add targets to the independent variable limits...
	for (auto pp : indep) {
		if (pp->has_min())
			sp.targets.push_back(new Target(pp->name, pp->min(),false,"lower"));
		if (pp->has_max())
			sp.targets.push_back(new Target(pp->name, pp->max(),false,"upper"));
		// if rf is independent always check it's limits
		if (pp->name == "rf")
			sp.rflimits = true;
	}
	// ... and rf if requested
	if (sp.rflimits) {
		Par* rf = gpset::find("rf");
		assert(rf != nullptr);
		if (rf->has_min())
			sp.targets.push_back(new Target(rf->name, rf->min(),false,"lower"));
		if (rf->has_max())
			sp.targets.push_back(new Target(rf->name, rf->max(),false,"upper"));
	}

	// 8) compute nx: the number of independent variables
	sp.nx = 2*ncmplx + indep.size();

	if (indep.size() > 3) {
		// print some user info re optimization...
		os.str("");
		os << "optimizing: ";
		string sep;
		for (size_t i=0; i<sp.project.size(); i++) {
			os << sep << sp.project[i] << " = " << sp.projdir[i];
			sep = ", ";
		}
		flaps::info(os.str());
	}

	// ... and targets
	if (sp.targets.size() == 1)
		flaps::info("target: ",sp.targets);
	else if (sp.targets.size() > 1)
		flaps::info(sp.targets.size()," targets:\n",sp.targets);

	// check some defaults
	post_preferences_defaults(prog);

	// Store standard parameters, rset, matrix
	// list for subsequent runs
	pset mygpset{gpset::get()};
	if (sp.aid.empty())
		sp.aid = "flut";
	os.str("");
	os << "parameters," << sp.aid;
	mygpset.clear_solns();
	mygpset.store(os.str());

	// store the mids of all matrices in Matrix::collection()
	// for subsequent "source" runs
	Matrix::store_mids (sp.aid);

	// print the preferences
	print_prefs();
}  // pre_process

static void
print_prefs() {
// Print execution options

	for (auto& sr : startregions())
		flaps::info ("start region: ", sr);

	// print summaries of equation matrices...
	cout << separator() << endl;
	cout << stringIndent(16, "Equation Matrices") << endl;
	cout << stringIndent(16, "-----------------") << endl;
	cout << Matrix::inventory();
	cout << separator() << endl;
	// ... and the parameters
	gpset::get().table(cout);
}

static bool
post_preferences_defaults(string const& prog) {
	Trace trc(1,"post_preferences_defaults");
	bool rval = true;
	size_t i, j;
	vector<Par*> indep = gpset::get().get_indep();
	ostringstream os;
	Specs& sp = flutspecs();

	if (sp.title.empty()) {
		os.str("");
		for (i=0; i<indep.size(); i++) {
			os << indep[i]->name;
			if (i < indep.size()-1)
				os << "-";
		}
		sp.title = os.str();
	}

	// Check that the user has not done something like
	// freq[0:10], startregion{freq[0:200]}, 
	vector<Region>& srlist = startregions();
	for (i=0; i<srlist.size(); i++) {
		for (j=0; j<srlist[i].size(); j++) {
			// ignore if setting parameter to a constant
			if (srlist[i][j].min == srlist[i][j].max)
				continue;
			Par* pp = gpset::find(srlist[i][j].parname);
			if (pp != nullptr) {
				double cf = pp->convFactor();
				if (pp->has_min() && pp->min() > srlist[i][j].max) {
					os.str("");
					os << "start region " << pp->name << "["
						<< srlist[i][j].min*cf << ':'
						<< srlist[i][j].max*cf <<
						"] is disjoint from the specified limits " << pp->stringLimits();
					throw runtime_error(os.str());
				}
				if (pp->has_max() && pp->max() < srlist[i][j].min) {
					os.str("");
					os << "start region " << pp->name << "["
						<< srlist[i][j].min*cf << ':'
						<< srlist[i][j].max*cf <<
						"] is disjoint from the specified limits " << pp->stringLimits();
					throw runtime_error(os.str());
				}
				if (pp->has_min() && pp->min() > srlist[i][j].min) {
					flaps::warning("start region min (", srlist[i][j].min*cf,
						") is less than specified min for ", pp->name,
						": using ", pp->premin());
					srlist[i][j].min = pp->min();
				}
				if (pp->has_max() && pp->max() < srlist[i][j].max) {
					flaps::warning("start region max (", srlist[i][j].max*cf,
						") is greater than specified max for ", pp->name,
						": using ", pp->premax());
					srlist[i][j].max = pp->max();
				}
			}
		}
	}

	// Ignore determinant sign changes if > 3 indep
	if (indep.size() > 3) {
		sp.bifurcation = false;
		flaps::info("ignoring determinant sign changes: ",
				indep.size()," indep");
	}

	// default parameters to print: multi-valued fixed, indep,
	//!! if (sp.toprint.empty()) {
		vector<Par*> fixed = gpset::get().get_fixed();
		// multi-valued fixed
		vector<Par*> mfv;
		for (auto pp : fixed)
			if (!pp->altval.empty())
				mfv.push_back(pp);
		for (auto pp : indep)
			add2vector(pp->name, sp.toprint);
			//!! sp.toprint.push_back(pp->name);
		// ... also include mfv parameters
		for (auto pp : mfv)
			add2vector(pp->name, sp.toprint);
			//!! sp.toprint.push_back(pp->name);
		// ... and eigenvector
		add2vector(string("eigenvector"), sp.toprint);
		//!! sp.toprint.push_back("eigenvector");
	//!! }

	// default plot file
	if (sp.plotfile.empty())
		sp.plotfile = sp.aid + ".pf";

	return rval;
}	// post_preferences

static void
check_param() {
// Do various checks on parameters:
//   - must have at least 3 indep
//   - check that Fixed parameters are in range
//   - if gcnorm is indep it must have limits
	Trace trc(1,"check_param");
	ostringstream os;
	vector<Par*> indep = gpset::get().get_indep();
	size_t nindep = indep.size();
	Specs& sp = flutspecs();

	trc.dprint(nindep," indep");

	if (nindep == 0) {
		trc.dprint("throwing exception: no indep");
		throw runtime_error("no indep parameters were specified");
	} else if (nindep < 3) {
		trc.dprint("throwing exception: 3 indep needed");
		throw runtime_error("at least 3 parameters must be indep");
	} else if (nindep > 3 && sp.project.empty()) {
		throw runtime_error("the optimize option must be included when"
			" there are more than 3 independents");
	}

	// Iterates contain eigenvectors, not gc. They will be
	// normalized to 2-norm 1.0 and gc's returned from parval() will
	// be scaled by gcnorm.
	// gcnorm must be given limits if it is indep
	Par* gcnorm = gpset::find("gcnorm");
	if (gcnorm->is_indep()) {
		// gcnorm is indep - set to it's minimum if no source run
		if (!gcnorm->has_min() || !gcnorm->has_max()) {
			trc.dprint("throwing exception: gcnorm indep without limits");
			throw runtime_error("gcnorm is indep with no limits");
		}
	}

	// check that fixed parameters are in-range and
	// if any parameters remain in "nostate"
	for (auto& par : gpset::get().pmap()) {
		Par* pp = par.second;
		if (pp->is_fixed()) {
			if (pp->has_min() && pp->value() < pp->min()) {
				flaps::warning(pp->summary()," is below its minimum (",pp->premin(),")");
			}
			if (pp->has_max() && pp->value() > pp->max()) {
				flaps::warning(pp->summary()," is above its maximum (",pp->premax(),")");
			}
		}
	}
	trc.dprint("returning:\n",gpset::get());
}	//check_param

void
print_start(vector<Flutcurve*>& startpts) {
	Trace trc(1,"print_start");
	ostringstream os;
	Specs& sp = flutspecs();

	os << startpts.size() << " Start Points: ";
	if (!sp.sid.empty()) {
		os << "from analysis \"" << sp.sid << '\"';
	} else {
		os << "low-speed modes";
	}
	string title(os.str());

	int pgwidth = page_width();
	cout << stringCenter(pgwidth, title) << endl;
	string ul(title.size(), '-');
	cout << stringCenter(pgwidth, ul) << endl;

	// get the width of the largest leader: labeloffset
	size_t lwidth{0};
	for (const auto cp : startpts) {
#ifdef NEVER // new cid format
		string lead{vastr(cp->aid(), ": ", cp->cid(), " ")};
#else // NEVER // new cid format
		string lead{vastr(cp->cid(), ": ")};
#endif // NEVER // new cid format
		lwidth = std::max(lwidth, lead.size());
	}
	string labeloffset(lwidth, ' ');
	cout << labeloffset << startpts[0]->make_header() << endl;

	// print a summary of each start point
	for (auto& cp : startpts)
#ifdef NEVER // new cid format
		cout << cp->aid() << ": " << cp->cid() << " " << cp->current_values() << endl;
#else // NEVER // new cid format
		cout << cp->cid() << ":  " << cp->current_values() << endl;
#endif // NEVER // new cid format
}

void
post_process(vector<Flutcurve*>& curves) {
// Do some post-processing on the solution curves:
//   1) search each curve for bifurcation & compute the bifurcation curves
//   2) store the results
//   3) print the results
//   4) plot the results
//   5) search each curve for the requested targets and print them
	Trace trc(1,"postProcess");
	string runid;
	int number = 0;
	Specs& sp = flutspecs();

	if (curves.empty())
		return;

	// 1) check for bifurcation
	if (sp.bifurcation) {
		//!! Trace trc(2,2,"postProcess");
		Fstack bifstack;
		for (auto curve : curves) {
			int i{1};
			for (auto& bp : curve->bif_origins) {
				string cid = curve->cid() + "_bif";
				if (curve->bif_origins.size() > 1) cid += to_string(i++);
				Flutcurve* bif = new Flutcurve(curve->aid(), cid, curve->vzid(), bp);
				bif->specs.bifurcation = false;	// no need to look for bif
				// save the bif origin as coord 0
				bif->origin.coord = 0;
				bif->update(bif->origin);
				bif->back_solns();
				bifstack.push_back(bif);
			}
		}
		flaps::info(bifstack.size()," bifurcations found");
		if (!bifstack.empty()) {
			concurrent(&bifstack, trace_curves);
			for (auto fp : bifstack)
				curves.push_back(fp);
		}
	}

	// 2) store all curves
	for (auto cp : curves)
		cp->store();

	// 3) print a summary of each curve, with footnotes
	for (auto& curve : curves) {
		string curveid = curve->cid();
		string error = curve->error;
		string finished = curve->finished;
		size_t npts = curve->nsolns();
		trc.dprint(curve->aid(),": ",curveid," has ",npts," data points");
		if (npts == 0) {
			string msg{vastr("no solutions computed for ", curveid)};
			if (!error.empty()) {
				msg += ": " + error;
			} else if (!finished.empty()) {
				msg += ": " + finished;
			}
			flaps::warning(msg);
			continue;
		}
		// ... and print with title = curveid
		string title{vastr("(", ++number, ") ", curve->aid(), ": ", curveid,
			" (thread ", curve->thread, ")")};

		string ul(title.size(), '-');
		cout << endl;
		cout << stringCenter(page_width(), title) << endl;
		cout << stringCenter(page_width(), ul) << endl;
		// combine curve (and targets?)
		vector<string> leaders;

		// this does the actual printing:
		cout << curve->print_solns(leaders, true) << endl;

		// footnote to print: how the tracking finished and error message
		if (!finished.empty())
			cout << "  Finished tracking: " << finished << endl;

		if (!error.empty())
			flaps::warning(error);

		cout << separator() << endl;
	}

	// 4) Create plotfile using Flutcurve::plot...
	string aid = curves[0]->aid();
	string plotfile{curves[0]->specs.plotfile};
	if (!plotfile.empty())
		Flutcurve::plot(aid, plotfile, curves);

	// 5) print the targets found
	int (*target_fcn)(pset& plt){nullptr};
	if (!sp.process_targets.empty())
		Custom::get(sp.process_targets, target_fcn);
	for (auto target : sp.targets) {
		if (target->is_lowerlimit || target->is_upperlimit)
			continue;
		vector<string> leaders;
		pset* found = seek_target(curves, target, leaders);
		if (found == nullptr) {
			flaps::info("no targets found for ",*target);
		} else {
			string title = vastr(found->nsolns(),
				" target(s) found for ", *target);
			string ul(title.size(), '-');
			cout << stringCenter(page_width(), title) << endl;
			cout << stringCenter(page_width(), ul) << endl;
			string header{curves[0]->make_header()};
			cout << found->summary_solns(sp.toprint, header, leaders, false) << endl;
			// pass to CustomPz target-processing function
			if (target_fcn != nullptr)
				target_fcn(*found);
		}
	}
} // postProcess()


// class targetcmp is for sorting Targets by increasing values of
// a specified parameter. A vector of targetcmp is created with
// the value of the sort parameter and it's location (index in solns).
// Then this vector is sorted by the sort parameter's values and
// each parameter in the pset of targets has it's solns re-arranged
// according to "loc"
class targetcmp {
public:
	size_t loc;
	double val;
	targetcmp() : loc(0), val(0.0) {}
	targetcmp(size_t i, double v) : loc(i), val(v) {}
	bool operator() (const targetcmp& a, const targetcmp& b) {
		return (a.val < b.val);
	}
};

static
pset*
seek_target(vector<Flutcurve*>& curves, Target* target, vector<string>& leaders) {
// search each curve for "target", put each found target into the solns
// of rval and the curve id into "leaders" for printing.
// Returns: a pset* with targets found in solns, or nullptr if none
	Trace trc(1,"seek_target ",*target);
	pset* rval{nullptr};
	Specs& sp = flutspecs();

	string pname{target->parname};
	double val{target->value};

	// save data for sorting
	string sortname = target->sortpar;
	vector<targetcmp> sortstuff;

	for (auto curve : curves) {
		Par* pp = curve->params.findp(pname);
		if (pp == nullptr)
			continue;
		Par* sortpar{nullptr};
		if (!sortname.empty())
			sortpar = curve->params.findp(sortname);
		vector<pair<size_t,double> > locs = pp->find_solns(val);
		trc.dprint(curve->cid()," has ",locs.size()," ",*target);
		for (auto loc : locs) {
			curve->params.interp(loc);
			// check target windows
			bool oor{false};
			for (auto& window : target->windows) {
				if (!window.inrange(curve->params)) {
					trc.dprint("reject target: out of window");
					oor = true;
					break;
				}
			}
			if (oor)
				continue;
			// evaluate all curve parameters
			curve->params.eval();
			// create pset* rval
			if (rval == nullptr) {
				rval = new pset(curve->params);
				rval->clear_solns();
			}
			if (sortpar != nullptr)
				sortstuff.push_back(targetcmp(leaders.size(), sortpar->value()));
			rval->back_solns(curve->params);
			leaders.push_back(vastr("  ",curve->cid(),": "));
		}
	}

	// sort if requested
	if (!sortname.empty() && sortstuff.size() > 2) {
		trc.dprint("before sorting:\n",rval->summary_solns(sp.toprint,"",leaders,false));
		sort (sortstuff.begin(), sortstuff.end(), targetcmp());
		for (auto& par : rval->pmap()) {
			Par* pp = par.second;
			vector<double> newval;
			for (auto& si : sortstuff) {
				newval.push_back(pp->solns[si.loc]);
			}
			pp->solns.clear();
			for (auto& vi : newval)
				pp->solns.push_back(vi);
		}
		// eigv
		for (auto& pp : rval->eigv()) {
			vector<double> newval;
			for (auto& si : sortstuff) {
				newval.push_back(pp->solns[si.loc]);
			}
			pp->solns.clear();
			for (auto& vi : newval)
				pp->solns.push_back(vi);
		}
		vector<string> newlead;
		for (auto& si : sortstuff) {
			newlead.push_back(leaders[si.loc]);
		}
		leaders = newlead;
		trc.dprint("sorted:\n",rval->summary_solns(sp.toprint,"",leaders,false));
	}

	trc.dprint("returning ",rval==nullptr?0:rval->nsolns()," found");
	return rval;
}

size_t
get_order(size_t* basic) {
// Returns the size of the largest matrix. All matrices must have
// the same basic size but some may have user-subroutine
// parameterizations with "extra" equations
	Trace trc(2,"get_order");
	vector<Matrix*>& stdmat = Matrix::collection();
	static size_t rval{0};
	static size_t smallest{std::numeric_limits<int>::max()};

	if (rval > 0) {
		trc.dprint("already called: matrix order=",rval,", smallest = ",smallest);
		if (basic != nullptr) *basic = smallest;
		return rval;
	}

	if (stdmat.empty())
		throw runtime_error("no input matrices were specified");

	for (auto mp : stdmat) {
		// only consider matrices in the std flutter eqn (mass, stif, etc)
		if (!is_flutmat(mp->desc()))
			continue;
		if (mp->rsize() != mp->csize()) {
			trc.dprint(mp->desc()," is not square: ",mp->summary());
			continue;
		}
		size_t order = mp->rsize();
		// get the largest (rval) and smallest matrix size
		smallest = std::min(smallest, order);
		rval = std::max(rval, order);
	}
	trc.dprint("returning ",rval,", smallest ",smallest);
	if (basic != nullptr) *basic = smallest;
	return rval;
}
