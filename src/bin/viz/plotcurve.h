//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#ifndef PLOTCURVE_H
#define PLOTCURVE_H

#include <cstdio>
#include <string>
#include <vector>
#include <wx/dc.h>
#include <wx/pen.h>

#include "Curve.h"
#include "matrix.h"
#include "Par.h"
#include "specs.h"

class Plotcurve;

enum class Dottype {dot, x, triangle, invtriangle, square, diamond, star};

// class Legend: identify a curve uniquely, primarily in a plot legend panel
// A legend panel has 3 sections:
// 1) an info box showing the coordinates of the cursor
// 2) 1 line for each curve:
//       sym num.ext yname
// 3) footnotes identify the different aid/plotfiles: sym lead
class Legend {
public:
	std::string id_;		// lead.num: unique id for a curve - access with id() XXX used where?
	std::string sid_;		// short version of id_ for CurvesSelect access with sid()
	std::vector<std::string> toks;	// id_ tokenized by .
	// a legend appears as: "sym num" with footnote "sym lead" iff
	// there are multiple "lead"; otherwise as "num" and no footnote
	// A curve cid is num.ext, e.g. 3.4.3 -> num=3, ext=4.3
	std::string lead;		// aid or filename
	std::string num;		// cid without ext: always a number?
	std::string ext;		// cid = num.ext, can be numbers or letters, e.g. 3.c
	std::string sym;		// plot symbol
	std::string key;		// key into the colormap: one of lead, num, or ext
	std::string yname;	// y parameter name if multiple y's: set in myFXYVector ctor
	wxPen pen;
	Dottype dot; 			// a dot, square, diamond, etc
	static std::map<std::string, wxPen> pens;		// all possible pens
	static std::map<std::string, Dottype> dots;	// all possible dot types

	// Legend ctor: all data is added where?
	Legend() {}
	Legend(Plotcurve* pc, const std::string& ynm);

	// unique idendifier for a curve: lead.num.ext
	const std::string& id() const { return id_; }
	// shortened version of id()
	const std::string& sid() const { return sid_; }

	// draw the current "dot" with an even-numbered size ranging from 1-10?
	void drawDot(wxDC& dc, wxCoord ix, wxCoord iy, int size);
	// for comparisons
	bool operator==(const Legend& b) { return lead == b.lead; }
};

std::ostream&
operator<<(std::ostream& s, const Legend& t);

class Plotcurve : public Curve {
public:
	std::string plotfile;	// XXX rm? use aid?
	std::string path;			// plotfile with .pf or .apf
	//!! Legend legend;	XXX unused?
	double xmin{0.0};
	double ymin{0.0};
	double xmax{0.0};
	double ymax{0.0};
	// inherited from Curve:
	// string aid()
	// string cid()
	// string vzid()
	// pset params

	// the one-and-only constructor XXX make private
	Plotcurve(Curve*);

	// Public interface:
	// load plotfiles and read from flaps database
	static std::vector<Plotcurve*> getcurves(
			const std::vector<std::string>& aids,
			const std::vector<std::string>& plotfiles);
	static int get_order();  // set in getcurves
};


#endif // PLOTCURVE_H
