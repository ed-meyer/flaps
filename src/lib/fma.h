//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Flutter mode animation data

#include  <complex>
#include  <string>
#include  <vector>
#include "fio.h"
#include "matrix.h"

class Fma : public Fio {
public:
	std::vector<int> nodes;
	std::vector<double> coords;
	std::vector<int> segidx;
	std::vector<std::complex<double>> gct;
	std::string vzid;

	// constructors: the only way to construct is to import a universal file
	// a universal file may have the form "name.vzid.uf", i.e. the vzid
	// is part of the name. To fetch an Fma from the Flaps data directory
	// use fetch()
	Fma() = default;
	Fma(const std::string& ufname,			// read a UF ctor
			const std::string& vzid="");
#ifdef NEVER // only use fetch
	Fma(const std::string& vzid);				// "fetch" constructor
#endif // NEVER // only use fetch
	Fma(const Fma&) = default;					// copy constructor
	Fma& operator=(const Fma&) = default;	// assignment op
	static Fma* fetch(const std::string& vzid);	//  read from database

	// fio stuff
	Fma(Receiver& r);
	static bool regd;			// registered?
	static std::string id() { return "Fma"; }
	std::string vid() const { return id(); }
	static std::string mid(std::string const& vzid="");
	// { return vastr("fma.",vzid); }
	static Fio* get(Receiver& r) { return new Fma(r); }
	void put (Sender& s) const;
	void store();
	// create a Matrix called "gct.vzid" from the gct
	Matrix* gct_matrix();
};

std::ostream& operator<<(std::ostream& s, const Fma& f);
