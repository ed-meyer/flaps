//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include "exim.h"
#include "lexer.h"
#include "lti.h"	// lib/lti.h

using namespace std;

string make_T(string& entrypt,
   Matrix& A, Matrix& B, Matrix& C, Matrix& D,
   vector<Itd>& Atd, vector<double>& itd, vector<double>& otd, vector<int>& Sidx,
   const string& psiname, const string& stifname, vector<pair<int,double>>& ecolk,
   vector<string>& igains, vector<string>& iphases,
   vector<string>& ogains, vector<string>& ophases, string& outputname);

bool
lti (const string& specs) {
// Create a "T" matrix representing a Linear Time-Invariant (LTI) control law

// Parse options for the "lti" option and create the output matrix with an
// LTIPz; options:
//   lti{file=path, psi=name, e=(list{scale})|name, igain=(list)|parname, ogain=(list)|parname}
	string outputname;
	string path;
	// stuff from the file
	Matrix A;
	Matrix B;
	Matrix C;
	Matrix D;
	vector<Itd> Atd;
	vector<double> itd;
	vector<double> otd;
	vector<int> Sidx;
	// stuff from the user
	string psiname;
	string stifname;
	vector<pair<int,double>> ecolk;
	vector<string> igains;
	vector<string> iphases;
	vector<string> ogains;
	vector<string> ophases;
	double smooth{0.1};		// smoothing parameter for the IntPz in A

	// parse the specs - get "smooth" first for LTI::importer...
	vector<Tok*> unrec = flaps::lexer(specs, {
		{"smooth", [&](const Tok& p) {smooth = p.rvec[0]; return true; }},
	});
	// ... then the rest
	unrec = flaps::lexer(specs, {
		{"file|i", [&](const Tok& p) {
			string path = vastr(getftmp(),"/",p.srhs);
			return LTI::importer (path, A, B, C, D, Atd, itd, otd, Sidx, smooth); }},
		{"psi", [&](const Tok& p) { psiname = p.srhs; return true; }},
		{"stif", [&](const Tok& p) { stifname = p.srhs; return true; }},
		{"e", [&](const Tok& p) {
			if (p.srhs.empty()) return false;
			if (isdigit(p.srhs[0])) {	// pairs
				for (size_t i=0; i<p.ivec.size(); i++) {
					double sf{1.0};
					if (!p.roptvec.empty())
						sf = stod(p.roptvec[i]);
					ecolk.push_back({p.ivec[i]-1,sf});		// zero-based
				}
			} else {
				return false;	// e=name disabled for now
			}
			return true;
		}},
		{"igain", [&](const Tok& p) {		// parameter name or value
			igains = p.svec; return true;
		}},
		{"iphase", [&](const Tok& p) {		// parameter name or value
			iphases = p.svec; return true;
		}},
		{"ogain", [&](const Tok& p) {		// parameter name or value
			ogains = p.svec; return true;
		}},
		{"ophase", [&](const Tok& p) {		// parameter name or value
			ophases = p.svec; return true;
		}},
		{"o", [&](const Tok& p) {outputname = p.srhs; return true; }},
	});

	// do some checking
	if (stifname.empty())
		throw runtime_error("a stiffness matrix was not specified");

	if (outputname.empty())
		outputname = "T";

	// save A, B, C, D
	A.mid(vastr(outputname,".A"));
	A.store();
	B.mid(vastr(outputname,".B"));
	B.store();
	C.mid(vastr(outputname,".C"));
	C.store();
	D.mid(vastr(outputname,".D"));
	D.store();

	// read the stiffness matrix (fetch constructor) and get it's size
	Matrix K(stifname);
	int ne = K.rsize();
	
	// dependent parameters: sigma, freq + K->dependson + A->dependson
	vector<string> dep_par;
	dep_par.push_back("sigma");
	dep_par.push_back("freq");
	auto pbegin = dep_par.begin();
	auto pend = dep_par.end();
	vector<string> kdep = K.dependson(gpset::get());
	for (auto& d : kdep) {
		if (find(pbegin,pend,d) == pend)
			dep_par.push_back(d);
	}
	vector<string> adep = A.dependson(gpset::get());
	for (auto& d : adep) {
		if (find(pbegin,pend,d) == pend)
			dep_par.push_back(d);
	}

	// create the output matrix: (ne+ns,ne+ns) complex...
	string entrypt;
	string code_path = make_T(entrypt, A, B, C, D,
		Atd, itd, otd, Sidx, psiname, stifname, ecolk,
		igains, iphases, ogains, ophases, outputname);
	// ... add it to custom.so...
	Custom::add(code_path);

	int ns = A.rsize();
	size_t m = ne + ns;
	Matrix output_matrix(outputname, "LTI T Matrix", m, m, true);
	// ... and add the CustomPz XXX is nr==m or ne?
	//!! ltifunctor fo(specs);
	//!! output_matrix.pz.push_back(new CustomPz(fo, m, m, ns));
	output_matrix.pz.push_back(new CustomPz(entrypt, m, m, ns));
	output_matrix.store();

	return true;
}

string
make_T(string& entrypt,
   Matrix& A,
   Matrix& B,
   Matrix& C,
   Matrix& D,
   vector<Itd>& Atd,
   vector<double>& itd,
   vector<double>& otd,
   vector<int>& Sidx,
   // user-specified:
   const string& psiname,
   const string& stifname,
   vector<pair<int,double>>& ecolk,
   vector<string>& igains,
   vector<string>& iphases,
   vector<string>& ogains,
   vector<string>& ophases, string& outputname) {

   // create the file on the cwd, return the full path
   entrypt = outputname;
   string rval{vastr(get_cwd(),"/",entrypt,".cpp")};
   ofstream ofs(rval, std::ios::out | std::ios::trunc);
   if (!ofs)
      throw runtime_error(vastr("cannot open ",rval));

	ofs << "int " << entrypt <<
		" (pset& plt, int nr, int nc, vector<complex<Ad>>& result) {\n";
	ofs << "\tMatrix A(\"" << A.mid() << "\");\n";
	ofs << "\tMatrix B(\"" << B.mid() << "\");\n";
	ofs << "\tMatrix C(\"" << C.mid() << "\");\n";
	ofs << "\tMatrix D(\"" << D.mid() << "\");\n";
	string sep;
	// internal time delays
	ofs << "\tvector<Itd> Atd{";
	for (auto& ai : Atd) {
		ofs << sep << "{{";
		string se;
		for (auto& ei : ai.elem) {
			ofs << se << "{" << ei.row << "," << ei.col << "}";
			se = ",";
		}
		ofs << "}," << ai.deltat << "}";
		sep = ",";
	}
	ofs << "};\n";
	// input time delays
	sep = "";
	ofs << "\tvector<double> itd{";
	for (auto& i : itd) {
		ofs << sep << i;
		sep = ",";
	}
	ofs << "};\n";
	// output time delays
	sep = "";
	ofs << "\tvector<double> otd{";
	for (auto& i : otd) {
		ofs << sep << i;
		sep = ",";
	}
	ofs << "};\n";
	// Sidx
	sep = "";
	ofs << "\tvector<int> Sidx{";
	for (auto& i : Sidx) {
		ofs << sep << i;
		sep = ",";
	}
	ofs << "};\n";
	
	ofs << "\tstring psiname{\"" << psiname << "\"};\n";
	ofs << "\tstring stifname{\"" << stifname << "\"};\n";
	// ecolk
	sep = "";
	ofs << "\tvector<pair<int,double>> ecolk{";
	for (auto& i : ecolk) {
		ofs << sep << "{" << i.first << "," << i.second << "}";
		sep = ",";
	}
	ofs << "};\n";
	// igains
	sep = "";
	ofs << "\tvector<string> igains{";
	for (auto& i : igains) {
		ofs << sep << "\"" << i << "\"";
		sep = ",";
	}
	ofs << "};\n";
	// iphases
	sep = "";
	ofs << "\tvector<string> iphases{";
	for (auto& i : iphases) {
		ofs << sep << "\"" << i << "\"";
		sep = ",";
	}
	ofs << "};\n";
	// ogains
	sep = "";
	ofs << "\tvector<string> ogains{";
	for (auto& i : ogains) {
		ofs << sep << "\"" << i << "\"";
		sep = ",";
	}
	ofs << "};\n";
	// ophases
	sep = "";
	ofs << "\tvector<string> ophases{";
	for (auto& i : ophases) {
		ofs << sep << "\"" << i << "\"";
		sep = ",";
	}
	ofs << "};\n";

	// call the evaluator
	ofs << "\treturn LTI::eval(plt, nr, nc, result,\n";
	ofs << "\t\tA, B, C, D, Atd, itd, otd, Sidx,\n";
	ofs << "\t\tpsiname, stifname, ecolk, igains, iphases, ogains, ophases);\n";
		
	ofs << "}\n";

	return rval;
}	// make_T

