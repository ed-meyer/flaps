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
#include "functions.h"
#include "lexer.h"
#include "matrix.h"
#include "message.h"
#include "settings.h"
#include "trace.h"

using namespace std;

static void doexp (const string& filename, vector<string>& matrices,
	const string& vzid, bool append);

//--------------------------------------------------------------------
// exporter {o=file, append, name1, name2, ... 
//-------------------------------------------------------------------
bool
exporter (string const& options) {
	T_(Trace trc(2,"exporter");)
	string filename;
	bool append{false};
	string vzid;
	vector<string> matrices;
	vector<Par*> params;

	T_(trc.dprint("options: \"",options,"\"");)

	// lambdas for Parsers
	vector<Tok*> unrec = flaps::lexer(options, {
		{"o", [&](const Tok& p) { filename = p.srhs; return true; }},
		{"append", [&](const Tok& p) { return append = true; }},
		{"vzid", [&](const Tok& p) { vzid = p.srhs; return true; }},
		{".*", [&](const Tok& p) {
			bool rval{true};
			if (p.srhs.empty()) {
				matrices.push_back(p.lhs);
				// matrix options: just parameter names,values
				if (!p.lopt.empty()) {
					vector<Tok*> unrec = flaps::lexer(p.lopt, {
						{".*", [&params](const Tok& p) {
							params.push_back(new Par(p.lhs,p.srhs));
							return true;
						}}
					});
				}
			} else {
				rval = false;
			}
			return rval;
		}}
	});

	if (!unrec.empty())
		throw runtime_error(vastr(unrec.size()," unrecognized options: ",unrec));

	// sanity check
	if (filename.empty())
		throw runtime_error("no output file given");

	// any parameters spec'd?
	T_(trc.dprint(params.size()," parameters specified");)
	if (!params.empty()) {
		pset& gp = gpset::get();
		for (auto pp : params) {
			Par* exist = gp.findp(pp->name);
			if (exist != nullptr) {
				exist->value(pp->value());
				T_(trc.dprint("set ",exist->name," to ",exist->value());)
			} else {
				gp.add(pp);
				flaps::info("added new parameter: ",*pp);
			}
		}
	}

	// do the export
	try {
		doexp(filename, matrices, vzid, append);
	} catch (runtime_error& s) {
		flaps::error(s.what());
		return false;
	}

	return true;
}

void
evaluate(vector<Matrix*> ml) {
// evaluate all matrices in "ml", putting the evaluated matrix into
// the Matrix.data_ array. Use gpset as the pset
	T_(Trace trc(2,"evaluate");)
	for (auto mp : ml) {
		T_(trc.dprint("working on ",mp->mid());)
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
			T_(trc.dprintm(nr,nc,nr,cp,mp->mid());)
		} else {
			double* rp = mp->elem();
			for (int i=0; i<nr*nc; i++) {
				double re = real(result[i]).value();
				rp[i] = re;
			}
			T_(trc.dprintm(nr,nc,nr,rp,mp->mid());)
		}
	}
}

static void
doexp (const string& filename, vector<string>& matrices,
		const string& vzid, bool append) {
// export a set of matrices to a file with the format determined
// by the file extension

	string::size_type idx = filename.find_last_of(".");
	if (idx == string::npos) {
		throw runtime_error("no extension on file name - "
				"cannot determine desired format");
	}
	string ext{filename.substr(idx+1)};
	// export according to the file extension:
	if (ext == "fig") {
		// fetch the Fma
		Fma* fma{nullptr};
		if (vzid.empty()) {
			flaps::warning("a vzid must be given for the \"fig\" option");
			return;
		}
		fma = Fma::fetch(vzid);
		Fig::exporter(filename, fma->coords);
	} else if (ext == "op4") {
		if (matrices.empty())
			throw runtime_error("no matrix names to export were given");
		// fetch the matrices
		vector<Matrix*> ml;
		for (auto& mid : matrices) {
			Matrix* mp = new Matrix(mid);
			ml.push_back(mp);
		}
		Output4::exporter(filename, ml, append);
	} else if (ext == "mat") {
		if (matrices.empty())
			throw runtime_error("no matrix names to export were given");
		// only Matrix's are relevant to matlab/octave
		vector<Matrix*> toexp;
		for (auto& re : matrices) {
			vector<Matrix*> ml = Matrix::fetch(re);

			for (auto mp : ml)
				toexp.push_back(mp);
		}
		// evaluate the matrices
		evaluate(toexp);
		Matlab::exporter(filename, toexp, append);
	} else if (ext == "uf") {
#ifdef NEVER // use vzmatrices
		// export "coordinates", "segments" along with matrices[0] (only 1 matrix)
		if (matrices.size() > 1)
			throw runtime_error("only one matrix may be exported to Universal files");
		Matrix* coords_m = new Matrix("coordinates");
		Matrix* segments_m = new Matrix("segments");
		// convert segments to int...
		vector<int> segments;
		for (auto i : segments_m->data())
			segments.push_back(i);
		Matrix* T_m = new Matrix(matrices[0]);
		// ...and T to complex
		vector<complex<double>> T;
		size_t n = T_m->rsize()*T_m->csize();
		if (T_m->is_complex()) {
			complex<double>* cp = T_m->celem();
			for (int i=0; i<n; i++)
				T.push_back(cp[i]);
		} else {
			for (auto i : T_m->data())
				T.push_back(i);
		}
		UF::exporter(filename, coords_m->data(), segments, T);
#else // NEVER // use vzmatrices
	cout << "new vzmatrices scheme not implemented\n";
#endif // NEVER // use vzmatrices
	} else {
		throw runtime_error(vastr("unrecognized file extension \"",ext,"\""));
	}
}

