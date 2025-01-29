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

#include "config.h"
#include "exim.h"
#include "fma.h"
#include "lexer.h"
#include "matrix.h"
#include "main.h"
#include "message.h"
#include "run_program.h"
#include "specs.h"

using namespace std;

static bool do_matlab (const string& cmd, const string& script,
		const vector<string> matrices, const vector<string> results);

//--------------------------------------------------------------------
// octlab {i=(file, name1, name2, ... ), o=(out1, out2,...)}
//-------------------------------------------------------------------
bool
octlab () {
	Trace trc(2,"octlab");
	vector<string> matrices;
	string script;
	Specs& sp = specs();

	trc.dprint("specs: ",sp);
	
	for (auto& s : sp.input) {
		if (rsubstr(s, 2) == ".m")
			script = s;
		else
			matrices.push_back(s);
	}
	trc.dprint("input matrices: ",matrices);

	// check the user's PATH for either matlab or octave:
	// take the first one found unless use_* set
	int mplace;
	int oplace;
	string path;
	string matlab_path = which("matlab", &mplace);
	string octave_path = which("octave", &oplace);
	if (mplace != -1 && oplace != -1) {
		if (mplace < oplace)
			path = matlab_path;
		else
			path = octave_path;
	} else if (mplace != -1) {
		path = matlab_path;
	} else if (oplace != -1) {
		path = octave_path;
	} else {
		throw runtime_error("neither matlab nor octave are in the PATH");
	}
	if (sp.use_octave) {
		if (oplace == -1)
			throw runtime_error("octave is not available");
		path = octave_path;
	}
	if (sp.use_matlab) {
		if (mplace == -1)
			throw runtime_error("matlab is not available");
		path = matlab_path;
	}

	string the_cmd = rsubstr(path,6);

	// default: all matrices
	if (matrices.empty())
		matrices = fio::catalog(".*");

	// run Matlab
	return do_matlab(the_cmd, script, matrices, sp.results);
}

void
evaluate(vector<Matrix*> ml) {
// evaluate all matrices in "ml", putting the evaluated matrix into
// the Matrix.data_ array. Use gpset as the pset
	Trace trc(2,"evaluate");
	for (auto mp : ml) {
		trc.dprint("working on ",mp->mid());
		if (mp->pz.empty())
			continue;
		int nr = mp->rsize();
		int nc = mp->csize();
		vector<complex<Ad>> result(nr*nc);
		mp->eval(gpset::get(), result);
		if (mp->is_complex()) {
			complex<double>* cp = mp->celem();
			for (int i=0; i<nr*nc; i++) {
				double re = real(result[i]).value();
				double im = imag(result[i]).value();
				cp[i] = complex(re,im);
			}
			trc.dprintm(nr,nc,nr,cp,mp->mid());
		} else {
			double* rp = mp->elem();
			for (int i=0; i<nr*nc; i++) {
				double re = real(result[i]).value();
				rp[i] = re;
			}
			trc.dprintm(nr,nc,nr,rp,mp->mid());
		}
	}
}

static bool
do_matlab (const string& matlab_cmd, const string& script,
		const vector<string> matrices, const vector<string> results) {
// run matlab:
// - export the specified matrices to a .mat file
// - create a Matlab .m script file with the user's commands
//   and commands to import the .mat file and export all results
//   XXX to the same .mat file??
// - run Matlab to execute the script
// - import the .mat file
	Trace trc(2,"do_matlab");

	// change to temp directory using RAII class Chdir
	string ftmp = getftmp();
	Chdir tmp(ftmp);

	// export the requested matrices to temp.mat
	// matlab/octave uses . to mean element-wise mult so change to _
	ostringstream os;
	os << "o=temp.mat";
	Fma* fma{nullptr};
	vector<Matrix*> ml;
	for (auto mi : matrices) {
		// watch out for a gct.vzid: fetch the Fma, create a Matrix
		// and store it with the mid "mi"
		if (mi.substr(0,3) == "gct") {
			string vzid;
			string::size_type idx = mi.find('.');
			if (idx != string::npos)
				vzid = mi.substr(idx+1);
			fma = Fma::fetch(vzid);
			ml.push_back(fma->gct_matrix());
		}
		// only Matrix's are relevant to matlab/octave
		for (auto& re : matrices) {
			vector<Matrix*> mlre = Matrix::fetch(re);
			for (auto mp : mlre)
				ml.push_back(mp);
		}
		// evaluate the matrices
		evaluate(ml);
	}
	Matlab::exporter("temp.mat", ml);

	if (script.empty())
		throw std::runtime_error(vastr("a matlab/octave script was not specified"));

	// open the user's script
	ifstream cmds(script);
	if (!cmds)
		throw std::system_error(errno, std::generic_category(),
			vastr("octlab script \"",script,"\" is not available"));

	// create the script (temp.m): load the .mat file+user's
	//  commands+save results.mat
	ofstream ofs("temp.m", std::ios::trunc);
	ofs << "load temp.mat;\n";
	string line;
	while(!cmds.eof()) {
		std::getline(cmds, line);
		ofs << line << endl;
	}
	// save everything on results.mat
	// the doc doesn't say, but matlab v7 (-7) is always compressed
	// according to octave's load-save.cc; so use v6 if we are using octave
	// XXX how to allow regular-expressions in return matrices?
	if (matlab_cmd == "octave")
		ofs << "save(\"-6\",\"results.mat\"";
	else
		ofs << "save(\"results.mat\",\"-nocompression\"";
	if (!results.empty())
		for (auto& ri : results)
			ofs << ",\"" << ri << '\"';
	ofs << ")\n";
	ofs << "quit\n";
	ofs.close();

	// run matlab - unset LD_PRELOAD while running matlab (libstdc++ conflict)
	string ld_preload = getEnv("LD_PRELOAD");
	if (!ld_preload.empty())
		unsetenv("LD_PRELOAD");
	string cmdline{vastr(matlab_cmd," -batch \"temp\"")};
	if (matlab_cmd == "octave")
		cmdline = vastr(matlab_cmd," temp.m");

	flaps::info("running \"",cmdline,"\"");

	// run the script in matlab/octave
	int status = system(cmdline.c_str());

	if (!ld_preload.empty())
		putEnv(vastr("LD_PRELOAD=",ld_preload));
	if (status != 0) {
		flaps::error(vastr(cmdline," returned status ",status));
		return false;
	}

	// import the results
	string unused;
	ml = Matlab::importer("results.mat",unused);

	for (auto mi : ml) {
		// If a gct was returned, create a new Fma using everything
		// from the old one except the gct: replace it with the imported
		// one, and get the new vzid from the imported gct
		string midi = mi->mid();
		if (fma != nullptr && midi.substr(0,3) == "gct") {
			trc.dprint("read gct: ",mi);
			vector<string> toks = string2tok(midi, ".");
			string newvzid;
			if (toks.size() > 1)
				newvzid = toks[toks.size()-1];
			fma->vzid = newvzid;
			size_t n = mi->rsize()*mi->csize();
			trc.dprint("resize gct to ",n);
			fma->gct.resize(n, 0.0);
			// gct may have gotten converted to real in matlab - convert
			// back to complex
			if (mi->is_complex())
				blas_copy(n, mi->celem(), 1, fma->gct.data(), 1);
			else
				blas_copy(n, mi->elem(), 1, (double*)fma->gct.data(), 2);
			fma->store();
			trc.dprint("saving new Fma:", *fma);
		} else {
			mi->store();
		}
	}
	return true;
}

int
main (int argc, char** argv) {
	try {
		parse_specs("");
		octlab();
	} catch (std::exception& s) {
		flaps::error(s.what());
	}
}
