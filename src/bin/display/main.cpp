//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <regex>
#include <stdlib.h> // system()

#include "config.h"
#include "Ad.h"
#include "exim.h"
#include "fio.h"
#include "main.h"
#include "matrix.h"
#include "specs.h"
#include "pset.h"
#include "settings.h"
#include "version.c"

using namespace std;

int display(const string& mid, pset& ps, Specs& sp);

// if parameters are specified as matrix options, create Altpset
// from gpset and set it's values and states instead of gpset, since
// gpset gets saved and we don't want it changed

int
main (int argc, char** argv) {
// Print one or more matrices to cout if small enough, or export to a Matrix
// Market file for visualization in matview if too large for cout.
	Trace trc(1,"main");
	string progname{argv[0]};

	try {
		version(progname);
		vector<pair<string,Specs>> matrices = parser();
		//!! Specs& sp = specs();

		for (auto& mi : matrices) {
			string mid{mi.first};
			Specs sp{mi.second};
			// if derivs or parameters requested create a copy of gpset to modify
			if (!sp.derivs.empty()) {
				vector<string> adpar;
				for (auto& di : sp.derivs)
					if (gpset::find(di) != nullptr)
						adpar.push_back(di);
				Ad::initialize(adpar);
				// and set the corresponding AD deriv to 1 in gpset
				gpset::get().realloc();
			}
			// create a copy of gpset
			pset altpset{gpset::get()};
			// set all parameters to Fixed
			for (auto& p : altpset.pmap())
				p.second->set_fixed();
					
			// if derivatives requested set the AD derivative parameter names
			if (!sp.derivs.empty()) {
				// set the state of derivs to Indep
				for (auto& par : sp.derivs) {
					Par* pp = altpset.findp(par);
					if (pp != nullptr)
						pp->set_indep();
				}
			}
			// if parameter values requested set them and create equations
			if (!sp.params.empty()) {
				for (auto& ti : sp.params) {
					Par* pp = altpset.findp(get<0>(ti));
					if (pp != nullptr)
						pp->valuef(get<1>(ti));
				}
				// now create equations for altpset
				altpset.equations();
				// set the AD parameter derivatives in altpset
				altpset.set_adpar_derivs();
				// mark all derived parameters !fresh
				altpset.touch();
			}
			// Do the actual printing
			display(mid, altpset, sp);
		}
	// exceptions just get a warning printed - no need to stop the program
	} catch (runtime_error& s) {
		flaps::warning(s.what());
	} catch (const std::exception& s) {
		flaps::warning(s.what());
	}
}

int
display(const string& mid, pset& altpset, Specs& sp) {
// fetch matrix 'mid', evaluate it, and display
	Trace trc(2,"display");
	ostringstream os;
	int rval{0};

	// fetch constructor
	Matrix mat(mid);

	// default path for output file: the matrix name
	string path{sp.outfile};
	if (path.empty())
		path = mat.mid();

	// extract all requested derivatives
	if (sp.derivs.empty()) sp.derivs.push_back("0"); // default: only the values
	for (auto& dername : sp.derivs) {
		int nr = mat.rsize();
		int nc = mat.csize();
		// extract the values or derivative
		vector<double> val = mat.values(dername, altpset);
		// create a header
		string header;
		if (dername == "0")
			header = vastr(mat);
		else
			header = vastr(mat.mid(), " derivative wrt ",dername);
		// print on cout if narrow enough...
		if (!sp.matview && mat.is_complex() && nc < 6) {
			//!! Dprint::dprintm(nr, nc, nr, val.data(), mat.mid());
			cout << endl << header << endl << strArray(nr, nc, nr,
				(complex<double>*)val.data());
		} else if (!sp.matview && nc < 7) {
			//!! Dprint::dprintm(nr,nc,nr,(complex<double>*)val.data(), mat.mid());
			cout << endl << header << endl << strArray(nr, nc, nr, val.data());
		// ... otherwise create a matrix market file
		} else {
			// give the .mm file a name including the deriv parameter
			string filename{path};
			if (dername == "0")
				filename = path + ".mm";
				//!! filename = vastr(stringBasename(path),".mm");
			else
				filename = path + "." + dername + ".mm";
				//!! filename = vastr(stringBasename(path),".",dername,".mm");
			// create the .mm file
			if (mat.is_complex())
				MM::exporter(filename, "display", (complex<double>*)&val[0],
						mat.rsize(), mat.csize());
			else
				MM::exporter(filename, "display", &val[0], mat.rsize(), mat.csize());
			// fire up matview in the background, unless NOVZ is in the environment
			if (sp.matview && getenv("NOVZ") == nullptr)
					rval = system(vastr("matview ",filename,"&").c_str());
		}
	}
	return rval;
}
