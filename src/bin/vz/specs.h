//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef Specs_h
#define Specs_h

#include <string>
#include <vector>

class Specs {
public:
	std::vector<std::string> aids;
	std::string title;
	bool amvz{false};
	std::vector<std::string> cids;
	std::string coloropt;
	std::string connectivityName;
	std::vector<std::string> diffids;
	bool gnuplot{false};
	std::vector<std::string> ignore;
	std::vector<std::string> plotfiles;
	int pmin{0};		// min point #
	int pmax{0};		// max point #
	std::string wname;
	std::string xname;
	bool has_xmin{false};
	double xmin;
	bool has_xmax{false};
	double xmax;
	bool xlog;
	std::vector<std::string> ynames;
	bool has_ymin{false};
	double ymin;
	bool has_ymax{false};
	double ymax;
	bool ylog;
	std::string zname;
	std::string saveclicksfile;
	int skip{0};	// skip data points
};

// convert options from command-line (argc,argv) to flaps options
std::string convertOptions (int argc, char** argv);
// get specs, set various policies
bool parse_specs(std::string const& optionList, Specs& sp);
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
