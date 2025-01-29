//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef CURVE_H
#define CURVE_H

#include <iostream>
#include <string>

#include "pset.h"
//------------------------------------------------------------------
// a Curve consists of a pset with each Par containing a set of
// values in the "solns" member, one per point on the curve
//------------------------------------------------------------------
class Curve : public Fio {
	std::string _aid;
	std::string _cid;
	std::string _vzid;
public:
	pset params;         // a copy of the global parameter list
	std::string finished;	// reason(s) for quitting
	std::string error;	   // if tracking stopped early

	// constructors:
	Curve() {}
	Curve(const std::string& aid, const std::string& cid,
				const std::string& vzid, const pset* ps=nullptr);
	~Curve();
	// default copy constructor & assignment op ok

	// get/set the curve id...
	std::string cid() const { return _cid; }
	void cid(std::string const& id) { _cid = id; }
	// get/set the analysis id
	std::string aid() const { return _aid; }
	void aid(std::string const& id) { _aid = id; }
	// get/set the visualization id
	std::string vzid() const { return _vzid; }
	void vzid(std::string const& id) { _vzid = id; }

	// how many data points in this curve?
	size_t nsolns() { return params.nsolns(); }

	// push the current param values onto the front or back of param.solns
	void front_solns();
	void back_solns();

#ifdef NEVER // rflimits: let the caller decide
	// check that all params (except rf) are in range
	std::string inrange() { return params.inRange(); }
#else // NEVER // rflimits: let the caller decide
	std::vector<std::string> inrange() { return params.inrange(); }
#endif // NEVER // rflimits: let the caller decide
	// find parameter named "name"
	Par* findp(std::string const& name) { return params.findp(name); }
	// get a vector of pointers to all parameters in this Curve
	std::vector<Par*> allpar() { return params.allpar(); }

	void print (std::ostream& stream);

	void plot ();
	void plot (std::string const& path, std::string const& title,
		std::string const& runid,
		std::vector<std::string> const& toPlot, bool append=false);

	// Fio/serialization
	static std::string id() { return std::string("Curve"); }
	std::string vid() const { return id(); }
	void put (Sender&) const;
	static Fio* get (Receiver&);
	static bool regd;  // has been registered?

	// create the mid of a Curve for store/fetch; if the
	// mid must be a regular expression use midrx
	static std::string mid(const std::string& aid, const std::string& cid);
	static std::string midrx(const std::string& aid, const std::string& cid);

	// split a Curve mid into components: aid and cid
	static void mid2comp(const std::string& mid,
		std::string& aid, std::string& curveid);

	// fetch/store a curve
	static Curve* fetch(const std::string& mid);
	void store();

};  // class Curve


std::ostream&
operator<<(std::ostream& s, const Curve& pl);

#endif  // CURVE_H
