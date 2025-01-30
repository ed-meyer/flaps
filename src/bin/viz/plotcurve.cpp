//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <map>
#include <string>
#include <unistd.h>
#include <vector>

#include "exim.h"
#include "message.h"
#include "plotcurve.h"
#include "settings.h"
#include "specs.h"

using namespace std;

static string load_pffile (string const& file);
static string load_apffile (string const& file);
static string load_others (string const& file);
static vector<Plotcurve*> getcids(string const& aid);
static void diffcurve (Plotcurve* a, Plotcurve* b);

// Legend implementation {

Legend::
Legend(Plotcurve* pc, const string& ynm) {
// Legend ctor
// set *some* members of Legend:
//   lead - either aid or plotfile
//   num  - from pc.cid
//   ext  - ditto
//   yname  ynm
// others cannot be set until we have a collection
// of all curves to plot
	T_(Trace trc(2,"Legend ctor ",ynm);)
	// "lead" is either aid or plotfile (XXX is aid always set to plotfile?)
	lead = pc->aid();
	if (lead.empty())
		lead = fio::shortenpath(pc->plotfile);
	// split cid into num.ext
	string ci = pc->cid();
	auto i = ci.find('.');
	if (i == string::npos)
		num = ci;
	else {
		num = ci.substr(0,i);
		ext = ci.substr(i+1);
	}
	yname = ynm;
}

std::ostream&
operator<<(std::ostream& s, const Legend& t) {
	s << "id \"" << t.id_ << "\", lead \"" << t.lead;
	s << "\", num \"" << t.num << "\", ext \"" << t.ext;
	s << "\", sym \"" << t.sym << "\", key \"" << t.ext << "\"";
	return s;
}

map<string,wxPen> Legend::pens;
map<string,Dottype> Legend::dots;

void
Legend::
drawDot(wxDC& dc, wxCoord ix, wxCoord iy, int size) {
// draw the current "dot" at (ix,iy) assuming "dc" has the desired pen
// XXX how to tell if this point is inside clipping region?
#ifdef NEVER // new size arg
	int dotwidth{6};			// width of any dottype
#else // NEVER // new size arg
	int dotwidth = size;
#endif // NEVER // new size arg
	int dw2 = dotwidth/2;	// radius of any dottype
	Dottype dottype = this->dot;
	if (dottype == Dottype::dot) {
		// draw a dot as a zero-length fat line
		//dotpen.SetWidth(2*dotwidth);
		//dc.SetPen(dotpen);
#ifdef NEVER // try ellipse
		dc.DrawLine(ix, iy, ix, iy);
#else // NEVER // try ellipse
		dc.DrawEllipse(ix, iy, dw2, dw2);
#endif // NEVER // try ellipse
	} else if (dottype == Dottype::x) {
		wxPen xpen = dc.GetPen();		// copy
		xpen.SetWidth(1);
		dc.SetPen(xpen);
#ifdef NEVER // too small
		dc.DrawLine(ix-dw2,iy-dw2,ix+dw2,iy+dw2);
		dc.DrawLine(ix-dw2,iy+dw2,ix+dw2,iy-dw2);
#else // NEVER // too small
		dc.DrawLine(ix-size,iy-size,ix+size,iy+size);
		dc.DrawLine(ix-size,iy+size,ix+size,iy-size);
#endif // NEVER // too small
	} else if (dottype == Dottype::triangle) {
		int r{dw2};
		vector<wxPoint> points{{ix-r,iy+r},{ix,iy-r},{ix+r,iy+r}};
		dc.DrawPolygon(3,points.data());
	} else if (dottype == Dottype::invtriangle) {
		int r{dw2};
		vector<wxPoint> points{{ix-r,iy-r},{ix+r,iy-r},{ix,iy+r}};
		dc.DrawPolygon(3,points.data());
	} else if (dottype == Dottype::square) {
		int r{dw2};
		vector<wxPoint> points{{ix-r,iy-r},{ix+r,iy-r},{ix+r,iy+r},{ix-r,iy+r}};
		dc.DrawPolygon(4,points.data());
	} else if (dottype == Dottype::diamond) {
		int r{dw2};
		vector<wxPoint> points{{ix-r,iy},{ix,iy-r},{ix+r,iy},{ix,iy+r}};
		dc.DrawPolygon(4,points.data());
	} else if (dottype == Dottype::star) {
		int r{dw2*4};
		vector<pair<int,int>> s{{-10,-6}, {-3,-6}, {0,-11},{3,-6},
			{10,-6}, {7,0}, {10,6}, {3,6}, {0,11},
			{-3,6}, {-10,6}, {-7,0}, {-10,-6}};
		vector<wxPoint> pts;
		for (auto& p : s)
			pts.push_back(wxPoint(ix+p.first*r/10,iy+p.second*r/10));
		dc.DrawLines(pts.size(),pts.data());
	}
}

// Legend implementation }

		

Plotcurve::
Plotcurve(Curve* curve) : Curve(*curve) {
// Plotcurve constructor.
	T_(Trace trc(1,"Plotcurve constructor");)
	Specs& sp = specs();

	// Set vzid in Specs from Curve::vzid to pass to amviz
	string vid = curve->vzid();
	if (!vid.empty() && sp.vzid.empty()) {
		sp.vzid = vid;
		T_(trc.dprint("taking vzid \"",sp.vzid,"\" from input curve");)
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
	if (access(base.c_str(), R_OK) != -1) {	// try the raw input name first...
		return base;
	} else if (haspf || hasapf) {
		throw runtime_error(vastr(base," is not available"));
	} else {
		string pffile = vastr(base,".pf");	// ... then with .pf extension ...
		if (access(pffile.c_str(), R_OK) != -1) {
			return pffile;
		}
		string apffile = vastr(base,".apf");	// ... then .apf
		if (access(apffile.c_str(), R_OK) != -1) {
			return apffile;
		}
		throw runtime_error(vastr(base," is not available (also tried ",
					pffile," and ",apffile,')'));
	}
	return "";
}

int Ngc{0};		// XXX use a static in get_order?

size_t
findmin(double p, deque<double> solns) {
// find the first point that is >= p, return it's index
	for (size_t i=0; i<solns.size(); i++) {
		if (solns[i] >= p)
			return i;
	}
	return 0;
}

size_t
findmax(double p, deque<double> solns) {
// find the first point that is <= p, return it's index
	size_t i{solns.size()-1};
	while (i >= 0) {
		if (solns[i] <= p)
			return i;
		if (i == 0)
			break;
		i--;
	}
	return 0;
}

vector<Plotcurve*>
Plotcurve::
getcurves(const vector<string>& aids, const vector<string>& plotfiles) {
// read data from one or both possible sources:
// - Curve objects in the flaps data dir from a previous execution
//   of flut
// - plotfiles specified on the command line (interactive usage)
//   or in the vz command in a Flaps control program.
// Convert the "solns" array to presentation units
	T_(Trace trc(1,"getcurves");)
	vector<Plotcurve*> rval;
	vector<string> paths;

	// process any curves in the data dir first - they may
	// get overwritten by plotfiles below
	for (auto& aid : aids) {
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
	for (auto& pf : plotfiles) {
		// check if it is accessible, or if it needs a .pf or .apf extension
		string file = accessible_file(pf);
		paths.push_back(file);		// full path
		string aid;
		// read the plotfile and write the data to the flaps data directory
		if (rsubstr(file, 4) == ".apf")
			aid = load_apffile(file);
		else if (rsubstr(file, 3) == ".pf")
			aid = load_pffile(file);
		else {
			aid = load_others(file);
		}

		// load_* stores the data in the Flaps data dir
		// read requested curves from the data dir, create Plotcurves
		try {
			vector<Plotcurve*> curves = getcids(aid);
			for (auto cp : curves) {
				cp->plotfile = pf;   // save the plotfile name to put in legend
				cp->path = file;   	// save the path to put in title bar
				cp->aid("");			// Ignore aid for plotfiles
				rval.push_back(cp);
			}
		} catch (runtime_error& s) {
			// cerr << "caught exception: " << s << endl;
			flaps::warning("processing ",pf,": ",s.what());
			continue;
		}
	}

	if (rval.empty()) {
		if (plotfiles.empty())
			flaps::warning("no aids nor plotfiles specified");
		else
			flaps::warning("no runs match the requested analysis ids");
		return rval;
	}

	// Convert all solns arrays to presentation units
	for (auto cp : rval) {
		for (auto pp : cp->params.all()) {
			double eu2pu = pp->conv_factor();
			if (eu2pu != 1.0) {
				for (auto& dp : pp->solns)
					dp *= eu2pu;
			}
		}
	}

	// get the order of the generalized coordinates for amviz
	size_t nrev{0};
	for (auto cp : rval) {
		vector<double> ev = cp->params.getev();
		nrev = std::max(nrev, ev.size());
	}
	Ngc = nrev/2;

	// replace cid's like mode_3 with just 3
	for (auto cp : rval) {
		string nm = cp->cid();
		string::size_type idx = nm.find('_');
		if (idx != string::npos) {
			cp->cid(nm.substr(idx+1));
		}
	}

	T_(trc.dprint("returning ",rval.size()," curves");)
	return rval;
}

int
Plotcurve::
get_order() { return Ngc; }

static string
load_pffile (string const& file) {
// Restore a plotfile to the datadir, figure out what analysis ids
// (aid) it contains, and add those to default policy "aids"
	T_(Trace trc(1,"load_pffile ", file);)
	
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
	regex re(R"(curve[,\.]([^,\.]*)[,\.].*)");
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
#ifdef NEVER // Fma
	// check for amviz data, load into Matrix::collection
	vector<string> comp{"nodes", "coords", "conn", "gct"};
	for (auto& s : mats) {
		for (auto& t : comp) {
			if (s.compare(0,t.size(),t) == 0) {
				Matrix* mp = new Matrix(s);
				Matrix::insert(mp);
			}
		}
	}
#endif // NEVER // Fma
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
	// get all curves with this aid
	curve_mids = fio::catalog(Curve::midrx(aid, ".*"));
	if (curve_mids.empty()) {
		string exc = vastr("no curves available for ",aid);
		T_(trc.dprint(exc);)
		// return rval;
		throw runtime_error(exc);
	}

	// the user may specify pairs of curves to diff
	vector<regex> diffre;
	for (auto& dp : sp.diffids)
		diffre.push_back(regex(dp));
	
	// compile each "cid" regular expression
	// may throw an exception if the regular expression is bad
	vector<regex> cidrx;
	for (auto ci : sp.cids)
		cidrx.push_back(make_regex(ci));
	
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
			// diffing curves? XXX needs work
			if (!diffre.empty()) {
				for (auto& pat : diffre) {
					if (regex_match(curveid, pat)) {
						for (auto& yname : sp.ynames) {
							diff.push_back(new Plotcurve(curve));
							T_(trc.dprint("added diff curve ",curveid,", y ",yname);)
						}
					}
				}
			} else {
				// create a new curve using the first y
				rval.push_back(new Plotcurve(curve));
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
			vector<string> all_cids = fio::catalog(Curve::midrx(aid,".*"));
			//!! flaps::warning("data were not found for analysis ", aid);
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

	// if a vzid was not specified, try to get it from a curve
	if (sp.vzid.empty() && !rval[0]->vzid().empty()) {
		sp.vzid = rval[0]->vzid();
		T_(trc.dprint("taking vzid \"",sp.vzid,"\" from the first curve");)
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
load_others (string const& path) {
// Read an ASCII data file with 1 or more items per line, separated by
// blanks, tabs, or a comma. create a pset with a new parameter for
// each column of data (x, y, z, or x1, x2, ... xn). If there is only one
// item per line call it y and x will be the ordinal number (1,2,3..)
	T_(Trace trc(1,"load_others ", path);)
	ostringstream os;
	Specs& sp = specs();
	string rval;
	string file = stringBasename(path);

	Ascii input;
	try {
		Ascii inp(path);
		input = move(inp);
	} catch (std::exception& s) {
		flaps::error(s.what());
		exit(1);
	}

	// the aid & cid are just the file name
	string aid = file;
	// string cid = cids[0];
	string cid = aid;
	// was a skip factor given?

	// get the first line & tokenize it to determine the # of columns (m)
	// if there is only one column add a first column: point number
	string line;

	input(line, "#");

	vector<string> toks = string2tok(line, " \t,");
	size_t ncol = toks.size();		// number of columns of data
	// make sure it contains numeric data
	double x;
	if (!str2double(toks[0], x)) {
		string exc = vastr(file," does not contain numeric data.",
				" First line: ",line);
		throw runtime_error(exc);
	}

	// now allocate vectors to hold the data
	vector<vector<double> > data(ncol);
	// and insert the first line
	size_t nc;
	for (size_t j=0; j<toks.size(); j++)
		data[j].push_back(stod(toks[j], &nc));
	
	// now do the rest of the file
	while(input.good()) {
		input(line, "#");
		if (line.empty())		// Ascii returns empty on eof
			break;
		toks = string2tok(line, " \t,");
		if (toks.size() != ncol) {
			string exc = vastr("reading ",file," line ",input.line_number()," has ",
				toks.size(), " items, others have ", ncol);
			throw runtime_error(exc);
		}

		for (size_t j=0; j<toks.size(); j++)
			data[j].push_back(stod(toks[j], &nc));
	}

	// the x parameter is the first one, call it "x" unless there is
	// only one column of data - then call it "point" and create
	// a column of point numbers
	pset ps;
#ifdef NEVER // disallow points?
	if (sp.xname == "points") {
		int npts = data[0].size();
		Par newpar(sp.xname,"0");
		Par* pp = ps.add(&newpar);
		for (int i=1; i<= npts; i++)
			pp->solns.push_back((double)i);
	}
#else // NEVER // disallow
	string name{"x"};
	if (ncol == 1) {
		name = "points";
		int npts = data[0].size();
		Par newpar(name,"0");
		Par* pp = ps.add(&newpar);
		for (int i=1; i<= npts; i++)
			pp->solns.push_back((double)i);
	}
	sp.xname = name;
#endif // NEVER // disallow

	// create 1 or 2 new parameters with plotfile data
	if (ncol == 1)
		name = "y";
	Par newpar(name,"0");
	Par* pp = ps.add(&newpar);
	for (size_t i=0; i<data[0].size(); i++)
		pp->solns.push_back(data[0][i]);
	if (ncol > 1) {
		name = "y";
		Par newpar(name,"0");
		Par* pp = ps.add(&newpar);
		for (size_t i=0; i<data[0].size(); i++)
			pp->solns.push_back(data[1][i]);
	}
	sp.ynames.push_back(name);

	// add to sp.views XXX shouldn't be necessary
	sp.views.push_back(vastr(sp.xname,":",sp.ynames[0]));

	// create a Curve with the new pset...
	Curve* cp = new Curve(aid, cid, sp.vzid, &ps);
	// ... and store it
	cp->store();
	T_(trc.dprint("ps (curve) (",ps.desc(), ") has ",ps.size()," parameters");)

	return aid;
}
