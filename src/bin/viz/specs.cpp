//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <sys/types.h>
#include <sys/stat.h>  // stat()

#include "config.h"
#include "lexer.h"
#include "plotcurve.h"
#include "settings.h"
#include "specs.h"
#include "text.h"
#include "trace.h"
#undef Complex

using namespace std;

// return a reference to the Specs: an automatic Singleton
Specs& specs() {
	static Specs thespecs;
	return thespecs;
}

static void
usage();

/*------------------------------------------------------------------
 * Globals: set the defaults for these guys here
 *------------------------------------------------------------------*/
bool Amviz = false;

bool LogX = false;
bool LogY = false;

bool Uncertainty = false;
bool NoUncertainty = false;

string
cmdline (int argc, char** argv) {
/*------------------------------------------------------------------
 * Convert command-line options like "-x veas -y growth" to
 * a string like "x=veas, y=growth"
 *------------------------------------------------------------------*/
	T_(Trace trc(1,"cmdline");)
	ostringstream os;
	int i = 1;   // first arg of interest

	while (i < argc) {
		if (argv[i][0] == '-') {
			string argvi(&argv[i][1]);
			if (argvi == "color") {
				os << "color\n";
				i++;
			} else if (argvi == "uf") {
				os << "uf=" << argv[i+1] << endl;
				i += 2;
			} else if (argvi == "title") {
				os << "title=" << argv[i+1] << endl;
				i += 2;
			} else if (argvi == "ign") {
				os << "ignore=" << argv[i+1] << endl;
				i += 2;
			} else if (argvi == "threads") {
				os << "threads\n";
				i++;
			} else if (argvi == "k") {
				os << "keep\n";
				i++;
			} else if (argvi == "skip") {
				if (i+1 < argc && argv[i+1][0] == '-') {
					os << "skip" << endl;
					i++;
				} else {
					os << "skip=" << argv[i+1] << endl;
					i += 2;
				}
			} else if (argvi == "vzid") {
				os << "vzid=" << argv[i+1] << endl;
				i += 2;
			} else if (compare(argvi, "xlog")) {
				os << "xlog\n";
				i++;
			} else if (compare(argvi, "ylog")) {
				os << "ylog\n";
				i++;
			// Check for an option following the current option; if
			// it starts with a '-' and the second char is not number
			// or period (e.g. -xmin=-.2) it is an option, otherwise it is a rhs
			} else if (i+1 < argc) {
				if (argv[i+1][0] != '-' || isdigit(argv[i+1][1]) ||
						argv[i+1][1] == '.') {
					os << &argv[i][1] << '=' << argv[i+1] << endl;
					i += 2;
				} else {
					i++;
				}
			// handle options like -color
			} else {
				os << &argv[i][1] << endl;
				i += 1;
			}
		// non-options: plotfile names
		} else {
			os << argv[i] << endl;
			i += 1;
		}
	}
	T_(trc.dprint("returning ",os.str());)
	return os.str();
}

bool
no_plotfiles(vector<string> plotfiles) {
// Do any of the strings in "plotfiles" have either
// .apf or .pf extensions? if so return false, otherwise
// return true
	for (auto pi : plotfiles) {
		if ((rsubstr(pi,4) == ".apf") ||
			(rsubstr(pi,3) == ".pf"))
			return false;
	}
	return true;
}

bool
parse_specs (string const& options, Specs& sp) {
	T_(Trace trc(1,"parse_specs");)
	bool rval = true;
	ostringstream os;
	vector<string> ignore;

	vector<Tok*> unrec = flaps::lexer(options, {
			{"^(a)?id$", [&](const Tok& p) {
				sp.aids.insert(sp.aids.end(),p.svec.begin(),p.svec.end());
				return true;
			}},
			{"^amviz", [&](const Tok& p) { return sp.amviz = true; }},
			{"^c(id)?$", [&](const Tok& p) {
				sp.cids.insert(sp.cids.end(),p.svec.begin(),p.svec.end());
				return true;
			}},
			{"^color$", [&](const Tok& p) { sp.colorby = PerCurve; return true;}},
			{"^d(iff)?$", [&](const Tok& p) {
				sp.diffids.insert(sp.diffids.end(),p.svec.begin(),p.svec.end());
				return true;
			}},
			{"^h(elp)?$", [&](const Tok& p) { usage(); exit(0); return true; }},
			{"^k(eep)?$", [&](const Tok& p) { putEnv("KEEPFTMP=1"); return true; }},
			{"ign(ore)?", [&](const Tok& p) {
				sp.ignore.insert(sp.ignore.end(),p.svec.begin(),p.svec.end());
				return true; }},
			{"plot(file)?", [&](const Tok& p) {
				if (p.svec.empty())
					sp.plotfiles.push_back(p.lhs);
				else
					sp.plotfiles.insert(sp.plotfiles.end(),p.svec.begin(),p.svec.end());
				return true; }},
			{"^pmin$", [&](const Tok& p) { sp.pmin = p.ivec[0]; return true; }},
			{"^pmax$", [&](const Tok& p) { sp.pmax = p.ivec[0]; return true; }},
			{"save(clicks)?", [&](const Tok& p) {
				sp.saveclicksfile = p.srhs;
				if (sp.saveclicksfile.empty())
					sp.saveclicksfile = ".vz";
				return true; }},
			{"skip", [&](const Tok& p) {sp.skip = p.ivec[0]; return true; }},
			{"^threads$", [&](const Tok& p) { sp.threads = true; return true; }},
			{"^title$", [&](const Tok& p) { sp.title = p.srhs; return true; }},
			{"uf", [&](const Tok& p) { sp.ufname = p.srhs; return true; }},
			{"^v$", [&](const Tok& p) { sp.views = p.svec; return true; }},
			{"^vzid$", [&](const Tok& p) {
				// append vzid to the default names of the 4 arrays:
				sp.vzid = p.srhs;
#ifdef NEVER // replaced with Fma?
				string ext = vastr("_",sp.vzid);
				sp.nodes += ext; sp.coords += ext; sp.conn += ext; sp.gct += ext;
#endif // NEVER // replaced with Fma?
				return true; }},
			{"^w$", [&](const Tok& p) { sp.wname = p.srhs; return true; }},
			{"^x$", [&](const Tok& p) { sp.xname = p.srhs; return true; }},
			{R"(^x\[.*\])", [&](const Tok& p) {
				sp.xname = p.srhs;
				string min, max;
				Par::parseLimits (p.lhs.substr(1), min, max);
				if (!min.empty()) {
					sp.has_xmin = true;
					str2double(min, sp.xmin);
				}
				if (!max.empty()) {
					sp.has_xmax = true;
					str2double(max, sp.xmax);
				}
				return true; }},
			{"^y$", [&](const Tok& p) {
				sp.ynames.insert(sp.ynames.end(),p.svec.begin(),p.svec.end());
				return true; }},
			{R"(^y\[.*\])", [&](const Tok& p) {
				sp.ynames.insert(sp.ynames.end(),p.svec.begin(),p.svec.end());
				string min, max;
				Par::parseLimits (p.lhs.substr(1), min, max);
				if (!min.empty()) {
					sp.has_ymin = true;
					str2double(min, sp.ymin);
				}
				if (!max.empty()) {
					sp.has_ymax = true;
					str2double(max, sp.ymax);
				}
				return true; }},
			{"^xmin$", [&](const Tok& p) {
				str2double(p.srhs, sp.xmin);
				return sp.has_xmin = true; }},
			{"^xlog$", [&](const Tok& p) { return sp.xlog = true; }},
			{"^xmax$", [&](const Tok& p) {
				str2double(p.srhs, sp.xmax);
				return sp.has_xmax = true; }},
			{"^ymin$", [&](const Tok& p) {
				str2double(p.srhs, sp.ymin);
				return sp.has_ymin = true; }},
			{"^ymax$", [&](const Tok& p) {
				str2double(p.srhs, sp.ymax);
				return sp.has_ymax = true; }},
			{"^ylog$", [&](const Tok& p) { return sp.ylog = true; }},
			{"^z$", [&](const Tok& p) { flaps::warning("the z option is ignored");
				return true; }},
			{"^.*", [&](const Tok& p) { sp.plotfiles.push_back(p.lhs); return true; }}
		}
	);

	if (!unrec.empty()) {
		if (unrec.size() == 1)
			flaps::warning("the following option was not recognized: ",*unrec[0]);
		else
			flaps::warning("the following options were not recognized: ",unrec);
	}

#ifdef NEVER // use only views
	if (!sp.xname.empty() && !sp.ynames.empty()) {
		for (auto yname : sp.ynames)
			sp.views.push_back(vastr(sp.xname,":",yname));
		T_(trc.dprint("sp views: ",sp.views);)
	}
#else // NEVER // use only views
	// if x, y names not given but v=(list) was, take the 1st view
	if (sp.xname.empty() && sp.ynames.empty() && !sp.views.empty()) {
		vector<string> toks = string2tok(sp.views[0], ":");
		sp.xname = toks[0];
		sp.ynames.push_back(toks[1]);
	}
#endif // NEVER // use only views

/*----------------  end of parse_specs ----------------------------------------*/

#ifdef NEVER // move to vis.cpp
	if (sp.aids.empty() && sp.plotfiles.empty()) {
		flaps::error("no analysis id or plot file specified");
		exit(1);
	}

	// check plotfile extensions: if it is .mm fire up matview
	vector<string> notmm;
	for (auto& pf : sp.plotfiles) {
		if (rsubstr(pf, 3) == ".mm") {
			string cmd = vastr("matview -F ", pf,"&");
			T_(trc.dprint("starting matview with ",cmd);)
			system(cmd.c_str());
		} else {
			notmm.push_back(pf);
		}
	}
	// re-set the list of plotfiles without the .mm files
	sp.plotfiles = notmm;
#endif // NEVER // move to vis.cpp

	// if there are 2 plotfiles and they are not .apf or .pf create
	// an aid and cid
#ifdef NEVER // needs work: bad message: no curves found for "text file "
	string xname;
	Settings::defaults.get("xname", xname);
	if (no_plotfiles(notmm) && (notmm.size() == 2 || notmm.size() == 1)) {
		string aid = "text file";
		vector<string> aids;
		aids.push_back(aid);
		Settings::defaults.set("aids", aids);
	}
#endif // NEVER : needs work: bad message: no curves found for "text file "

	return rval;
}

bool
updateBoundingBox (double x, double xmin, double xmax,
		double y, double ymin, double ymax,
		double& llx, double& urx, double& lly, double& ury) {
	Specs& sp = specs();

	if (sp.has_xmin)
		llx = sp.xmin;
	if (sp.has_xmax)
		urx = sp.xmax;
	if (sp.has_ymin)
		lly = sp.ymin;
	if (sp.has_ymax)
		ury = sp.ymax;

	if (!sp.has_xmin && xmin < llx) {
		if (sp.has_ymin && sp.has_ymax) {
			if (y >= sp.ymin && y <= sp.ymax)
				llx = xmin;
		} else {
			llx = xmin;
		}
	}

	if (!sp.has_xmax && xmax > urx) {
		if (sp.has_ymin && sp.has_ymax) {
			if (y >= sp.ymin && y <= sp.ymax)
				urx = xmax;
		} else {
			urx = xmax;
		}
	}

	if (!sp.has_ymin && ymin < lly) {
		if (sp.has_xmin && sp.has_xmax) {
			if (x >= sp.xmin && x <= sp.xmax)
				lly = ymin;
		} else {
			lly = ymin;
		}
	}

	if (!sp.has_ymax && ymax > ury) {
		if (sp.has_xmin && sp.has_xmax) {
			if (x >= sp.xmin && x <= sp.xmax)
				ury = ymax;
		} else {
			ury = ymax;
		}
	}
	return true;
}

static void
usage() {
	cout << "Usage:\n";
	cout << " vz [options] file(s)\n";
	cout << " Options:\n";
	cout << "   -x xparam       x parameter name\n";
	cout << "                   default: first parameter in the file\n";
	cout << "   -y yparam       y parameter name. This option may be repeated\n";
	cout << "                   to specify multiple y parameters\n";
	cout << "                   default: second parameter in the file\n";
	cout << "   -z zparam       z parameter name. If a z parameter is specified,\n";
	cout << "                   a 3D plot is created using gnuplot\n";
	cout << "                   default: a 2D plot is created\n";
	cout << "   -xmin min       Limits are specified with these options.\n";
	cout << "   -xmax max       If a min is specified, a max must also be specified\n";
	cout << "   -ymin min\n";
	cout << "   -ymax max\n";
	cout << "   -color=oneper   If multiple files are given on the command\n";
	cout << "                   line this option causes all curves in a file\n";
	cout << "                   to be the same color. The default is to assign\n";
	cout << "                   colors according to runid, so that if a runid\n";
	cout << "                   appears in more than one file, those curves\n";
	cout << "                   will have the same color.\n";
	cout << "  file(s)          The path names of one or more plotfiles.\n";
}

std::ostream&
operator<<(std::ostream& s, const Specs& t) {
	s << "aids: " << t.aids << endl;
	s << "title: " << t.title << endl;
	s << "cids: " << t.cids << endl;
	//!! std::string coloropt;
	//!! std::vector<std::string> diffids;
	//!! bool gnuplot{false};
	//!! std::vector<std::string> ignore;
	s << "plotfiles: " << t.plotfiles << endl;
	//!! int pmin{0};		// min point #
	//!! int pmax{0};		// max point #
	// x stuff
	s << "x name: " << t.xname << endl;
	s << "has xmin: " << t.has_xmin << endl;
	//!! double xmin;
	//!! bool has_xmax{false};
	//!! double xmax;
	//!! bool xlog;
	// multiple y stuff
	s << "y names: " << t.ynames << endl;
	//!! bool has_ymin{false};
	//!! double ymin;
	//!! bool has_ymax{false};
	//!! double ymax;
	//!! bool ylog;
	//!! std::string saveclicksfile;
	//!! int skip{0};	// skip data points
	//!! std::string ufname;		// universal file path
	s << "views: " << t.views << endl;
	return s;
}
