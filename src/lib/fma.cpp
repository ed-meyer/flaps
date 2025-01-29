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
#include "fma.h"
#include "matrix.h"
#include "trace.h"

using namespace std;

// constructors
Fma::
Fma (const string& ufname, const string& id) {
// import a Universal file, construct an Fma. If id is empty (default)
// attempt to get it from the ufname: "name.vzid.uf", otherwise the
// vzid is "name", all but the .uf
	Trace trc(2,"Fma ",ufname);
	// connectivity from a uf is in "penlift" format: series of
	// node numbers defining lines, each series separated by a zero
	// (penlift).
	vector<int> penlift;
	UF::importer(ufname, nodes, coords, penlift, gct);
	segidx = UF::penlift2segidx(penlift, nodes);

	// if input "id" is empty take the vzid from the file name
	if (id.empty()) {
		string file{ufname};
		string::size_type idx = ufname.rfind('/');
		// strip off any directories
		if (idx != string::npos)
			file = ufname.substr(idx+1);
		vector<string> toks = string2tok(file, ".");
		size_t n = toks.size();
		trc.dprint("ufname has ",n," toks: ",toks);
		if (n < 2) {
			vzid = toks[0];
		} else if (n == 2) {
			if (toks[1] == "uf")
				vzid = toks[0];
			else
				vzid = toks[1];
		} else if (n > 2) {
			if (toks[n-1] == "uf")
				vzid = toks[n-2];
			else
				vzid = toks[n-1];
		}
		trc.dprint("defaulting vzid to ",vzid);
	} else {
		vzid = id;
	}

	// Save it in the flaps database XXX really?
	store();
}

#ifdef NEVER // only use fetch
Fma::
Fma(const string& vzid) {
// Fma vzid constructor XXX replace fetch?
// Read Fma named "vzid" from the Flaps data directory,
// throws a runtime_error exception if the file is corrupt
	Trace trc(2,"Fma ctor \"",vzid,"\"");

	string filename{"fma"};
	if (!vzid.empty())
		filename += "." + vzid;

	Receiver r(filename);
	// Fio::get will return an object or throw an exception
	Fio* op = Fio::get(r);
	Fma* fop = dynamic_cast<Fma*>(op);
	if (fop == nullptr) {
		string err = vastr("reading \"",filename,
			"\", expecting an Fma, got a ",op->vid());
		trc.dprint("throwing exception: ",err);
		throw runtime_error(err);
	}
	*this = *fop;
	trc.dprint("read ",filename);
}
#endif // NEVER // only use fetch

Matrix*
Fma::
gct_matrix() {
// create a Matrix with the generalized-coordinate transformation (gct)
	Trace trc(2,"Fma::gct_matrix");
	int m = this->coords.size();
   int n = this->gct.size()/m;
	string gctmid{"gct"};
	if (!vzid.empty())
		gctmid += "." + vzid;
	Matrix* rval = new Matrix(gctmid,"gct",m,n,true);
	blas_copy(m*n, this->gct.data(), 1, rval->celem(), 1);
	trc.dprint("returning ",*rval);
	return rval;
}

string
Fma::
mid(string const& vzid) {
// Return the name an Fma is stored under
	if (vzid.empty())
		return "fma";
	return vastr("fma.",vzid);
}

// fio stuff
bool Fma::regd = Fio::register_class(Fma::id(), Fma::get);
Fma::
Fma (Receiver& r) {
	Trace trc(2,"Fma Receiver constructor");
	r.serialize(nodes);
	r.serialize(coords);
	r.serialize(segidx);
	r.serialize(gct);
	r.serialize(vzid);
}

void
Fma::
put(Sender& s) const {
	Trace trc(2,"Fma::put");
	s.serialize(nodes);
	s.serialize(coords);
	s.serialize(segidx);
	s.serialize(gct);
	s.serialize(vzid);
}

void
Fma::
store() {
	Trace trc(2,"Fma::store");
	string filename{"fma"};
	if (!vzid.empty())
		filename += vastr(".",vzid);
	Sender file(filename);
	Fio::put(file, *this);
	trc.dprint("stored \"",filename,"\"");
}

Fma*
Fma::
fetch(const string& vzid) {
// Read Fma named "vzid" from the Flaps data directory,
// throws a runtime_error exception if the file is corrupt
	Trace trc(2,"Fma::fetch \"",vzid,"\"");
	Fma* rval{nullptr};

	string filename{"fma"};
	if (!vzid.empty())
		filename += "." + vzid;

	Receiver r(filename);
	// Fio::get will return an object or throw and exception
	Fio* op = Fio::get(r);
	rval = dynamic_cast<Fma*>(op);
	if (rval == nullptr) {
		string err = vastr("reading \"",filename,
			"\", expecting an Fma, got a ",op->vid());
		trc.dprint("throwing exception: ",err);
		throw runtime_error(err);
	}
	trc.dprint("read ",filename);
	return rval;
}

ostream&
operator<<(ostream& s, const Fma& f) {
	s << Fma::mid(f.vzid) << ":\n";
	s << "   " << f.nodes.size() << " nodes\n";
	s << "   " << f.coords.size() << " coords\n";
	s << "   " << f.segidx.size() << " segidx\n";
	s << "   " << f.gct.size() << " gct";
	return s;
}

#ifdef MAIN
#include "main.h"

int
main (int argc, char** argv) {
	Ftmpdir ftmp;
	if (argc < 2) {
		cerr << "Usage: fma ufname\n";
		exit(1);
	}
	string ufname{argv[1]};
	Fma fma(ufname);
	cout << "created:\n" << fma << endl;
	fma.store();
	Fma fetched;
	fetched.fetch(fma.vzid);
	cout << "fetched:\n" << fma << endl;
}
#endif // MAIN
