//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Implementation of class Curve, a pset with each Par containing
// the same number of elements in the "values" array, describing the Curve

#include "config.h"
#include "Curve.h"
#include "trace.h"

using namespace std;

Curve::
Curve(const string& analysis_id, const string& curveid,
	const string& vzid, const pset* ps) {
// main Curve constructor: set aid, cid,
// initialize the pset to the current values in gpset
	T_(Trace trc(1,"Curve constructor ",curveid);)
	_aid = analysis_id;
	_cid = curveid;
	_vzid = vzid;
	T_(trc.dprint("aid \"",_aid,"\", cid \"",_cid,"\", vzid \"",_vzid,"\"");)
	// set member "params" to a copy of the gpset if "ps" not included
	if (ps != nullptr) {
		params = *ps;
		T_(trc.dprint("using non-global parameter set");)
	} else {
		params = gpset::get();
		// sanity check on pset assignment operator
		assert(params.findp("alt") != gpset::find("alt"));
		// delete all solns
		params.clear_solns();
	}
	// ...and set it's description to curveid
	params.desc(curveid);
	// evaluate params and check for constant Derived parameters
	T_(trc.dprint("created curve aid<",_aid,"> cid<",_cid,">");)
}

// Destructor: stuff required for the Fio vtable
Curve::
~Curve() {
// virtual destructor body required in the .c file for Fio vtable
	T_(Trace trc(1,"Curve destructor ",_cid);)
}

// register class Curve for fio serialization
bool Curve::regd = Fio::register_class(Curve::id(), Curve::get);


Fio*
Curve::
get(Receiver& s) {
// serialize a Curve
	T_(Trace trc(2,"Curve::get");)

	Curve* rval = new Curve();

	s.serialize(rval->_aid);
	s.serialize(rval->_cid);
	s.serialize(rval->_vzid);
	T_(trc.dprint("got aid \"",rval->_aid,"\", cid \"",rval->_cid,"\", vzid \"",rval->_vzid,"\"");)
	Fio* op = pset::get(s);
	pset* pp = dynamic_cast<pset*>(op);
	if (pp == nullptr) {
		throw runtime_error(vastr("reading a Curve, expecting a pset, got a ",
					op->vid()));
	}
	rval->params = *pp;   // need move assignment op
	delete pp;

	s.serialize(rval->finished);
	s.serialize(rval->error);
	return rval;
}

void
Curve::
put(Sender& s) const {
// serialize a Curve
	T_(Trace trc(2,"Curve::put");)
	T_(trc.dprint("putting aid \"",_aid,"\", cid \"",_cid,"\", vzid \"",_vzid,"\"");)
	s.serialize(_aid);
	s.serialize(_cid);
	s.serialize(_vzid);
	params.put(s);
	s.serialize(finished);
	s.serialize(error);
}

string
Curve::
mid(const string& analysis_id, const string& curveid) {
// returns the mid of a Curve with analysis id and curve id
	ostringstream os;
	os << "curve." << analysis_id << "." << curveid;
	return os.str();
}

string
Curve::
midrx(const string& analysis_id, const string& curveid) {
// returns a regular expression of Curve mids with analysis id,
// and curve id, possibly regular expressions
// XXX allow for either commas or period separators
	ostringstream os;
	os << "curve[,\\.]" << analysis_id << "[,\\.]" << curveid;
	return os.str();
}

void
Curve::
mid2comp(const string& mid, string& aid, string& cid) {
// given a Curve mid in the form "curve.aid.cid (usually created
// by Curve::mid()), split it into 2 components: aid and cid
	// XXX allow for either comma or period separators
	vector<string> toks = string2tok(mid, ".,");
	if (toks.size() < 3)
		throw runtime_error(vastr("invalid mid passed to Curve::mid2comp: ",mid));
	aid = toks[1];
	// cid is all the remaining toks
	cid = toks[2];
	for (size_t i=3; i<toks.size(); i++)
		cid += '.' + toks[i];
}

Curve*
Curve::
fetch (const string& mid) {
// read a Curve from the flaps database, return nullptr if it
// does not exist.
// Note: Curve::mid(aid,cid) can be used to create argument "mid"
	T_(Trace trc(1,"Curve::fetch ",mid);)
	Curve* rval{nullptr};
	try {
		Receiver file(mid);
		if (!file.good() || file.eof()) {
			T_(trc.dprint("returning empty: ",mid," does not exist");)
			return rval;
		}
		Fio* op = Fio::get(file);
		rval = dynamic_cast<Curve*>(op);
		if (rval == nullptr) {
			throw runtime_error(vastr("reading ", mid, ", expecting Curve, got ",
					op->vid()));
		}
	} catch (std::exception& s) {
		T_(trc.dprint("returning empty: ",mid," does not exist");)
		return rval;
	}
	T_(trc.dprint("returning ",rval->cid(),", ",rval->params.nsolns()," points");)
	return rval;
}

void
Curve::
store () {
//------------------------------------------------------------------
// Store a Curve in the flaps database with a mid given by Curve::mid()
//------------------------------------------------------------------
	T_(Trace trc(1,"Curve::store ",this->aid());)

	const string mid = Curve::mid(this->aid(), this->cid());
	// Solutions are stored using Curve::store() with a
	// a filename "mid"
	Sender file(mid);
	Fio::put(file, *this);
}


void
Curve::
front_solns() {
// push the current values in "params" onto the front of the "solns" array
	params.front_solns(params);
}

void
Curve::
back_solns() {
// push the current values in "params" onto the back of the "solns" array
	params.back_solns(params);
}

void
Curve::
plot() {
// simplified plotting of a Curve: file name, title, and
// runid are all the Curve::id
	params.plot();
}

void
Curve::
plot (std::string const& path, std::string const& title,
		string const& runid,
		vector<string> const& toPlot, bool append) {
	T_(Trace trc(1,"Curve::plot");)
	T_(trc.dprint(params.nsolns()," points to plot");)

	params.plot (path, title, runid, toPlot, append);
}

#ifdef MAIN

int
main (int argc, char** argv) {
	try {
		string aid("vso");
		string cid("5.3.4");
		string mid(Curve::mid(aid, cid));
		Curve::mid2comp(mid, aid, cid);
		cout << "mid \"" << mid << "\", aid: \"" << aid <<
			"\", cid: \"" << cid << "\"\n";
	} catch (std::exception& s) {
		cerr << "caught exception: " << s.what() << endl;
	}
}

#endif // MAIN
