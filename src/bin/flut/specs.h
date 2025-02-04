//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef specs_h
#define specs_h

#include <string>
#include <vector>

class Specs;

#include "matrix.h"
#include "pac.h"
#include "Par.h"
#include "region.h"
#include "target.h"

class Cu2eu {
public:
	std::string name;
	double value;
	Cu2eu(std::string n, double v) : name(n), value(v) {}
};
std::ostream& operator<<(std::ostream& s, const Cu2eu& t);

class Specs {
public:
	std::string aid;						// analysis id
	std::string sid;						// source id

	double alt_scale{2000.0};			// see Flutcurve::cu2eu()
	std::string apodo;					// nickname: vacm = variable alt, const Mach
	bool bifurcation{false};			// check for bifurcation?
	double bif_init{0.00001};			// initial step for bifurcation
	int checkjacobian{0};
	std::vector<std::string> checkjacobian_par;
	bool checklooping{false};
	bool compare_det{false};
	bool constrained{false};
	std::vector<Cu2eu> cu2eu;			// pre-defined or user-spec
	int direction{1};						// direction the curve is being traced
	std::string dmatrix;
	int ecc{0};								// use ECC? -1/0/1: never/default/always
	bool e9n{false};						// include exploration data in plotfile?
	bool free_vibration{false};
	bool lcostab{false};
	bool linearize{false};
	std::vector<std::string> nlev;
	int nx;									// # of independent variables
	Pacspecs* pacspecs{nullptr};		// specs for Pac
	std::string plotfile;
	bool plot_homotopy{false};
	bool printmatrices{false};
	bool printsummary{true};
	std::string process_targets;		// name of a Custom function
	std::vector<std::string> project;// parameter names
	std::vector<double> projdir;		// -1:1 for each project parameter
	std::vector<double> projectee;	// *the* vector that guides optimization
	double reflen{1.0};
	bool rflimits{false};				// honor rf limits?
	bool sequential{false};
	double sigma_scale{10.0};			// see Flutcurve::cu2eu()
	bool stopatdsc{false};
	bool svd{false};
	std::vector<Target*> targets;
	std::string title;
	std::vector<std::string> toplot;
	std::vector<std::string> toprint;
	std::string vzid;

	void update(const std::string& name, double value);
};

Specs& flutspecs();	// ref to the global specs

std::ostream& operator<<(std::ostream& s, const Specs& t);

void parser (std::string const& prog);

void print_prefs(std::string const& prog);

// is "desc" a matrix in the flutter equation (mass, stif, gaf, etc)?
bool is_flutmat (std::string const& desc);

// increment the multi-valued fixed parameters and return
// a tag. To initialize, call "incr_mvf(0)", and increment
// subsequently, "incr_mvf();"
std::string incr_mvf(int init=-1);

int
is_voe(int state=-1);


#endif // specs_h
