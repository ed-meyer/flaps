//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Specs for viz and amviz
#ifndef Specs_h
#define Specs_h

#include <string>
#include <map>
#include <vector>

constexpr int PerCurve{0};
constexpr int PerCurveNum{1};
constexpr int PerFile{2};
constexpr int PerExt{3};
constexpr int PerY{4};

class Specs {
public:
	std::vector<std::string> aids;
	std::string title;
	bool amviz{false};
	bool noviz{true};
	std::vector<std::string> cids;
	bool threads{false};		// use threads or timer?
	//!! std::string coloropt;
	int colorby{-1};
	std::string connectivityName;
	std::vector<std::string> diffids;
	bool gnuplot{false};
	std::vector<std::string> ignore;
	std::vector<std::string> plotfiles;
	std::vector<std::string> paths;	// plotfiles with .pf or .apf
	int pmin{0};		// min point #
	int pmax{0};		// max point #
	std::vector<std::string> views;	// (vtas:sigma,vtas:freq,...)
	std::string wname;	// change curve color if w changes sign
	// x stuff from the options
	std::string xname;
	bool has_xmin{false};	// x,y,z limits from the options
	double xmin;
	bool has_xmax{false};
	double xmax;
	bool xlog;
	// multiple y stuff
	std::vector<std::string> ynames;
	bool has_ymin{false};
	double ymin;
	bool has_ymax{false};
	double ymax;
	bool ylog;
	std::string saveclicksfile;
	int skip{0};	// skip data points
	// amviz stuff:
	// visualization data comes from either a Universal file:
	std::string ufname;		// universal file path
	// or the Flaps database
	std::string vzid;
};

std::ostream& operator<<(std::ostream& s, const Specs& t);

// convert options from command-line (argc,argv) to flaps options
std::string cmdline (int argc, char** argv);
// get specs, set various policies
bool parse_specs(std::string const& optionList, Specs& sp);

// return a reference to the Specs: an automatic Singleton
Specs& specs();

bool updateBoundingBox (double x, double xmin, double xmax,
		double y, double ymin, double ymax,
		double& llx, double& urx, double& lly, double& ury);

bool getLogx();
bool getLogy();

bool getUncertainty();
bool getNoUncert();

std::string const& getConnectivityName();
std::map<std::string,std::string> const& getAttributes();

#endif // Specs_h
