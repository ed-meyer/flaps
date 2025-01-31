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

#include "exim.h"
#include "message.h"
#include "plotcurve.h"
#include "settings.h"
#include "specs.h"
#include "trace.h"

using namespace std;

static void extract_points(deque<double>& solns);
static string load_plotfile (string const& file);
static string load_apffile (string const& file);
static string load_others (string const& file);
static vector<Plotcurve*> getcids(string const& aid);
static void diffcurve (Plotcurve* a, Plotcurve* b);

Plotcurve::
Plotcurve(Curve* curve, const string& yname) : Curve(*curve) {
// Plotcurve constructor. The yname argument is to allow
// for multiple y's from each plotfile; or think of it as
// allowing for multiple views of the same curve
	Specs& sp = specs();
	string exc;
	// get yparam first in case there is no xname
	yparam = params.findp(yname);
	if (yparam == nullptr)
		throw runtime_error(vastr(yname," is not a parameter in ",curve->aid()));

	// if xname was not specified create a parameter of ordinal point numbers
	if (sp.xname.empty()) {
		xparam = new Par("point(ordinal point number)","1");
		for (size_t i=1; i<=yparam->nsolns(); i++)
			xparam->solns.push_back(i);
	} else {
		xparam = params.findp(sp.xname);
		if (xparam == nullptr) {
			// the name may be prefixed by @
			string altxname = "@" + sp.xname;
			xparam = params.findp(altxname);
			if (xparam == nullptr)
				throw runtime_error(vastr(sp.xname," is not a parameter in ",curve->aid()));
			sp.xname = altxname;
		}
	}

	string zname = sp.zname;;
	string wname = sp.wname;

	if (!sp.zname.empty()) {
		zparam = params.findp(sp.zname);
		if (zparam == nullptr) {
			exc = vastr(sp.zname," is not a parameter in ",curve->aid());
			flaps::warning(exc);
		}
	} else {
		zparam = nullptr;
	}
	if (!sp.wname.empty()) {
		wparam = params.findp(sp.wname);
		if (wparam == nullptr) {
			exc = vastr(sp.wname," is not a parameter in ",curve->aid());
			flaps::warning(exc);
		}
	} else {
		wparam = nullptr;
	}

	// set policies xunittext, yunittext
	string xunittext{sp.xname};
	if (xparam != nullptr && !xparam->conv.punits.empty())
		xunittext += " (" + xparam->conv.punits + ")";
	Settings::defaults.set("xunittext", xunittext);
	if (!Settings::defaults.is_defined("yunittext")) {
		string yunittext{yname};
		if (!yparam->conv.punits.empty())
			yunittext += " (" + yparam->conv.punits + ")";
		Settings::defaults.set("yunittext", yunittext);
	}
	// Convert all "solns" arrays to presentation units...
	// and remove their conversion factors
	size_t nv = params.nsolns();
	for (auto& par : params.pmap()) {
		Par* pp = par.second;
		assert(pp->solns.size() == nv);
		extract_points(pp->solns);   // extract [pmin:pmax]
		// set fresh to true, i.e. take all values as-is
		pp->fresh = true;
		if (pp->has_conv()) {
			double fac = pp->conv.factor;
			if (fac != 1.0) {
				for (auto& ap : pp->solns)
					ap *= fac;
				pp->conv = {};
				if (pp->has_min())
					pp->min(fac*pp->min());
				if (pp->has_max())
					pp->max(fac*pp->max());
				// set current value to presentation-unit value
				pp->valuef(pp->solns[0]);
			}
		}
	}
	// set tagged to false
	for (size_t i=0; i<nv; i++)
		tagged.push_back(false);

	// get eigenvectors if available
	evmatrix = params.get_evmatrix();
}

static void
extract_points(deque<double>& solns) {
// replace solns with only points between pmin:pmax if
// they were requested
	Specs& sp = specs();
	if (sp.pmin == 0 && sp.pmax == 0)
		return;
	deque<double> newvals;
	for (int i=sp.pmin; i<=sp.pmax; i++)
		newvals.push_back(solns[i]);
	solns = newvals;
}

void
add_legends(vector<Plotcurve*> curves) {
	T_(Trace trc(1,"add_legends");)

	if (curves.size() == 1) {
		curves[0]->legend = curves[0]->cid();
		return;
	}

	// replace plotfile names with the difference between names
	string ref{curves[0]->plotfile};
	bool use_plotfile{false};
	if (!ref.empty()) {
		curves[0]->plotfile = string_diff(ref, curves[1]->plotfile);
		for (size_t i=1; i<curves.size(); i++) {
			curves[i]->plotfile = string_diff(curves[i]->plotfile, ref);
		}
		use_plotfile = true;
	}

	bool use_aid{false};
	bool use_cid{false};
	bool use_yname{false};
	string aid{curves[0]->aid()};
	string cid{curves[0]->cid()};
	string yname{curves[0]->yparam->name};
	T_(trc.dprint("model legend: ",aid,':',cid,':',yname);)
	for (auto pc : curves) {
		if (pc->aid() != aid)
			use_aid = true;
		if (pc->cid() != cid)
			use_cid = true;
		if (pc->yparam->name != yname)
			use_yname = true;
	}
	T_(trc.dprint("use aid? ",use_aid,", cid? ",use_cid,", yname? ",use_yname);)
	for (auto pc : curves) {
		if (use_plotfile)
			pc->legend = pc->plotfile;
		if (use_aid) {
			if (!pc->legend.empty())
				pc->legend += ":";
			pc->legend = pc->aid();
		}
		if (use_cid) {
			if (!pc->legend.empty())
				pc->legend += ":";
			pc->legend += pc->cid();
		}
		if (use_yname) {
			if (!pc->legend.empty())
				pc->legend += ":";
			pc->legend += pc->yparam->name;
		}
		T_(trc.dprint("legend: ",pc->legend);)
	}
}

bool
has_extension(const string& s, const string& ext) {
// does string "s" have an extension "ext", e.g. ".pf"
	string::size_type idx = s.rfind('.');
	if (idx == string::npos)
		return false;
	if (s.substr(idx) == ext)
		return true;
	return false;
}

static string
accessible_file(const string& base) {
// is there a readable file "base", or "base.pf", or "base.apf",
// if so, return the name that is accessible
	bool haspf = has_extension(base, ".pf");
	bool hasapf = has_extension(base, ".apf");
	if (access(base.c_str(), R_OK) != -1) {
		return base;
	} else if (haspf || hasapf) {
		throw runtime_error(vastr(base," is not available"));
	} else {
		string pffile = vastr(base,".pf");
		if (access(pffile.c_str(), R_OK) != -1) {
			return pffile;
		}
		string apffile = vastr(base,".apf");
		if (access(apffile.c_str(), R_OK) != -1) {
			return apffile;
		}
		throw runtime_error(vastr(base," is not available (also tried ",
					pffile," and ",apffile,')'));
	}
	return "";
}

int Ngc{0};

vector<Plotcurve*>
Plotcurve::
getcurves() {
// read data from one or both possible sources:
// - Curve objects in the flaps data dir from a previous execution
//   of flut
// - plotfiles specified on the command line (interactive usage)
//   or in the vz command in a Flaps control program.
	T_(Trace trc(1,"getcurves");)
	vector<Plotcurve*> rval;
	Specs& sp = specs();

	// process any curves in the data dir first - they may
	// get overwritten by plotfiles below
	for (auto& aid : sp.aids) {
		try {
			vector<Plotcurve*> curves = getcids(aid);
			if (curves.empty())
				flaps::warning("no curves found from analysis \"",aid,"\"");
			for (auto cp : curves)
				rval.push_back(cp);
		} catch (runtime_error& s) {
			// cerr << "caught exception: " << s << endl;
			continue;
		}
	}

	// now load any plotfiles into the Flaps datadir
	for (auto& pf : sp.plotfiles) {
		// check if it is accessible, or if it needs a .pf or .apf extension
		string file = accessible_file(pf);
		string aid;
		// read the plotfile and write the data to the flaps data directory
		if (rsubstr(file, 4) == ".apf")
			aid = load_apffile(file);
		else if (rsubstr(file, 3) == ".pf")
			aid = load_plotfile(file);
		else {
			aid = load_others(file);
		}

		// read requested curves from the data dir, create Plotcurves
		try {
			vector<Plotcurve*> curves = getcids(aid);
			for (auto cp : curves) {
				cp->plotfile = pf;   // save the plotfile name to put in legend
				rval.push_back(cp);
			}
		} catch (runtime_error& s) {
			// cerr << "caught exception: " << s << endl;
			flaps::warning("processing ",pf,": ",s.what());
			continue;
		}
	}

	if (rval.empty())
		return rval;

	// try to shorten the legends by eliminating common components
	add_legends(rval);

	if (rval.empty())
		flaps::warning("no runs match the requested analysis ids");

	// get the order of the generalized coordinates for amvz
	size_t nrev{0};
	for (auto cp : rval) {
		vector<double> ev = cp->params.getev();
		nrev = std::max(nrev, ev.size());
	}
	Ngc = nrev/2;

	T_(trc.dprint("returning ",rval.size()," curves");)
	return rval;
}

int
Plotcurve::
get_order() { return Ngc; }

static string
load_plotfile (string const& file) {
// Restore a plotfile to the datadir, figure out what analysis ids
// (aid) it contains, and add those to default policy "aids"
	T_(Trace trc(1,"load_plotfile ", file);)
	
	// check for existence of plotfile
	int mode = R_OK;
	string actual{file};
	if (access(actual.c_str(), mode) != 0) {
		actual += ".pf";
		if (access(actual.c_str(), mode) != 0) {
			string exc = vastr(file," is not available");
			throw runtime_error(exc);
		}
	}
	// try restoring the plotfile
	vector<string> mats = fio::restore(actual);

	// no matrices loaded? must not be a savefile
	if (mats.empty()) {
		string exc = vastr("could not restore ",actual,": not a plotfile?");
		throw runtime_error(exc);
	}
	// determine the analysis id (aid) assuming curves are stored
	// with names like "curve,aid,cid"
	smatch mch;
	string aid;
	regex re(R"(curve[,.]([^,.]*)[,.].*)");
	for (auto& mp : mats) {
		if (regex_match(mp, mch, re)) {
			aid = mch[1];
			break;
		}
	}
	if (aid.empty()) {
		ostringstream os;
		os << actual << " does not appear to be a plotfile:\n";
		for (auto& mp: mats)
			os << mp << endl;
		throw runtime_error(os.str());
	}
	// check for amvz data, load into Matrix::collection
	vector<string> comp{"nodes", "coords", "conn", "gct"};
	for (auto& s : mats) {
		for (auto& t : comp) {
			if (s.compare(0,t.size(),t) == 0) {
				Matrix* mp = new Matrix(s);
				Matrix::insert(mp);
			}
		}
	}
	return aid;
}

static
vector<Plotcurve*>
getcids(string const& aid) {
// fetch Curves with analysis id (aid) from the Flaps data dir,
// create Plotcurve's from them and return them
	T_(Trace trc(1,"getcids ",aid);)
	vector<pset> ml;
	vector<Plotcurve*> rval;
	vector<Plotcurve*> diff;
	Specs& sp = specs();

	vector<string> curve_mids;
	// no requested curve id's? get all of them, use regex
	string mid = Curve::midrx(aid, ".*");
	curve_mids = fio::catalog(mid);
	if (curve_mids.empty()) {
		string exc = vastr("no curves available for ",aid);
		T_(trc.dprint(exc);)
		return rval;
	}

	// the user may specify pairs of curves to diff
	vector<regex> diffre;
	for (auto& dp : sp.diffids)
		diffre.push_back(regex(dp));
	
	// compile each "cid" regular expression
	// may throw an exception if the regular expression is bad
	vector<regex> cidrx;
	for (auto ci : sp.cids) {
		// anchor the regex
		string cirx = vastr(ci,"$");
		cidrx.push_back(make_regex(cirx));
	}
	
	// compile each "ignore" regular expression
	// may throw an exception if the regular expression is bad
	vector<regex> ignrx;
	for (auto ig : sp.ignore)
		ignrx.push_back(make_regex(ig));

	// fetch each curve, check mode numbers, etc
	vector<string> descriptions;
	for (auto& si : curve_mids) {
		// ignore this curve?
		bool skip{false};
		for (auto& re : ignrx)
			if (regex_search(si, re)) skip = true;
		if (skip) continue;
		// requested?
		if (!cidrx.empty()) {
			skip = true;
			for (auto& re : cidrx)
				if (regex_search(si, re)) skip = false;
			if (skip) continue;
		}

		try {
			// fetch the curve
			Curve* curve = Curve::fetch(si);
			if (curve == nullptr) {
				flaps::warning("curve \"",si,"\" is not available");
				continue;
			}

			string curveid = curve->cid();
			// diffing curves?
			if (!diffre.empty()) {
				for (auto& pat : diffre) {
					if (regex_match(curveid, pat)) {
						for (auto& yname : sp.ynames) {
							diff.push_back(new Plotcurve(curve, yname));
							T_(trc.dprint("added diff curve ",curveid,", y ",yname);)
						}
					}
				}
			} else {
				// for each requested y create a new curve
				// watch out for regular expressions
				for (auto& yname : sp.ynames) {
					T_(trc.dprintn("looking for yname<",yname,">...");)
					vector<Par*> pars;
					Par* par = curve->params.findp(yname);
					if (par == nullptr) {
						pars = curve->params.findRe(yname);
					} else {
						pars.push_back(par);
					}
					T_(trc.dprint("found ",pars.size());)
					for (auto pp : pars) {
						Plotcurve* pc = new Plotcurve(curve, pp->name);
						rval.push_back(pc);
					}
				}
			}
		} catch (runtime_error& s) {
			throw runtime_error(vastr("reading curve ",si,": ",s.what()));
		}
	}

	// were any requested curves identified?
	if (rval.empty() && diff.empty()) {
		if (!sp.cids.empty()) {
			vector<string> all_cids = fio::catalog(Curve::midrx(aid,".*"));
			cerr << "curves: " << all_cids << endl;
			throw runtime_error("none of the requested curve ids found");
		} else {
			flaps::warning("data were not found for analysis ", aid);
		}
	}

	T_(trc.dprint("got ",rval.size()," curves, ",diff.size()," diffs");)

	// if there are diffids create new curves with the diff
	if (!diff.empty()) {
		for (size_t i=0; i<diff.size(); i+= 2) {
			diffcurve(diff[i], diff[i+1]);
			rval.push_back(diff[i]);
			T_(trc.dprint("added diff curve ",diff[i]->params.desc());)
		}
	}

	// set the color & style for each curve
	string coloropt;

	static int colorno{-1};
	colorno++;
	for (size_t i=0; i<rval.size(); i++) {
		rval[i]->color = colorno%MAXATTR;
		if (sp.coloropt.empty())
			colorno++;
		// rval[i]->style = i%MAXATTR;
		rval[i]->style = 0;
	}
	
	T_(trc.dprint("returning ",rval.size()," curves");)
	return rval;
}  // getcids

static
void
diffcurve (Plotcurve* a, Plotcurve* b) {
	T_(Trace trc(1,"diffcurve");)
	size_t j;

	T_(trc.dprint("diffing ",a->params.desc()," and ",b->params.desc());)
	for (auto& ami : a->params.pmap()) {
		string name = ami.first;
		if (name == "coord")
			continue;
		Par* ai = ami.second;
		Par* bi = b->params.findp(name);
		// they may not have the same number of solns - just compare
		// up to the smallest number
		size_t na = ai->nsolns();
		size_t nb = bi->nsolns();
		size_t n = std::min(na, nb);
		T_(trc.dprint(name," has ",na," and ",nb," solns");)
		for (j=0; j<n; j++) {
			ai->solns[j] = abs(ai->solns[j] - bi->solns[j]);
		}
		if (na > n) {
			for (j=n; j<na; j++)
				ai->solns[j] = 0.0;
		}
	}
}

static string
load_apffile (string const& path) {
// Read an ascii plot file (.apf), create a pset for each run
	T_(Trace trc(1,"load_apffile ", path);)
	ostringstream os;

	// the analysis id is the file name minus the .apf extension & directory
	string aid{stringBasename(path)};

	vector<Curve*> curves;
	try {
		curves = Apf::importer(path);
	} catch (runtime_error const& s) {
		flaps::warning("reading ",path,": ",s.what());
		return "";
	}
	// store the curves
	for (auto cp : curves) {
		cp->store();
	}
	return aid;
}

static string
load_others (string const& file) {
// Read a data file with 1 or more items per line, separated by
// blanks, tabs, or a comma. create a pset with a new parameter for
// each column of data (x, y, z, or x1, x2, ... xn). If there is only one
// item per line call it y and x will be the ordinal number (1,2,3..)
	T_(Trace trc(1,"load_others ", file);)
	ostringstream os;
	Specs& sp = specs();

	// open the file
	ifstream input(file);
	if (!input.good()) {
		// try it with a .pf extension
		string pf = vastr(file,".pf");
		if (access(pf.c_str(), R_OK) != -1) {
			return load_plotfile(pf);
		}
		string exc = vastr(file," is not available");
		throw runtime_error(exc);
	}

	// the aid & cid must have been defined
	vector<string> aids;
	vector<string> cids;

	string aid = sp.aids[0];
	// string cid = cids[0];
	string cid = aid;
	// was a skip factor given?
	int skip = sp.skip;

	// get the first line & tokenize it to determine the # of columns (m)
	// if there is only one column add a first column: point number
	string line;
	while (line.empty()) {
		std::getline(input, line);
	}
	int lineno{1};
	vector<string> toks = string2tok(line, " \t,");
	size_t ntok = toks.size();
	// make sure it contains numeric data
	double x;
	if (!str2double(toks[0], x)) {
		string exc = vastr(file," does not contain numeric data.",
				" First line: ",line);
		throw runtime_error(exc);
	}

	// now allocate vectors to hold the data
	size_t ndata = ntok;
	vector<vector<double> > data(ndata);
	// and insert the first line
	size_t nc;
	for (size_t j=0; j<toks.size(); j++)
		data[j].push_back(stod(toks[j], &nc));
	
	// now do the rest of the file
	while(input.good()) {
		std::getline(input, line);
		if (line.empty())
			break;
		lineno++;
		if (skip > 0 && lineno%skip != 0)
			continue;
		toks = string2tok(line, " \t,");
		if (toks.size() != ntok) {
			string exc = vastr("reading ",file," line ",lineno," has ",
				toks.size(), " items, others have ", ntok);
			throw runtime_error(exc);
		}

		for (size_t j=0; j<toks.size(); j++)
			data[j].push_back(stod(toks[j], &nc));
	}

	// the x parameter is the first one, call it "x" unless there is
	// only one column of data - then call it "point" and create
	pset ps;
	// check for xname == "points": create a parameter with pt numbers
	if (sp.xname == "points") {
		int npts = data[0].size();
		Par newpar(sp.xname,"0");
		Par* pp = ps.add(&newpar);
		for (int i=1; i<= npts; i++)
			pp->solns.push_back((double)i);
	}

	// create the new parameter with plotfile data
	Par newpar(file,"0");
	Par* pp = ps.add(&newpar);
	for (size_t i=0; i<data[0].size(); i++)
		pp->solns.push_back(data[0][i]);

	// add this new parameter to the pset of the Curve if
	// it exists
	string mid = Curve::mid(aid, cid);
	Curve* cp = Curve::fetch(mid);
	if (cp == nullptr) {
		// create a Curve with the gpset...
		cp = new Curve(aid, cid, "");
		// ... replace gpset params with ps ...
		cp->params = ps;
	} else {
		for (auto pp : ps.pmap())
			cp->params.add(pp.second);
		for (auto pp : ps.eigv())
			cp->params.add(pp);
	}
	// ... and store it
	cp->store();
	T_(trc.dprint("ps (curve) (",ps.desc(), ") has ",ps.size()," parameters");)

	return aid;
}
