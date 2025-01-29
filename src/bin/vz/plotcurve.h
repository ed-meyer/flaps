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
#include <X11/Xos.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#undef Complex

#undef index
#include "Curve.h"
#include "matrix.h"
#include "Par.h"
#include "xgout.h"
#include "xtb.h"


#define MAXKEYS		50
#define MAXATTR 	8
#define MAXBUFSIZE 	120
#define MAXLS		50

typedef unsigned long Pixel;

/* Globally accessible values */
extern Display *disp;			/* Open display            */
extern Visual *visual;			/* Standard visual         */
extern Colormap cmap;			/* Standard colormap       */
extern int screen;			/* Screen number           */
extern int depth;			/* Depth of screen         */

void do_hardcopy( char *prog, xtb_data info,
	int (*init_fun)(FILE*,int,int,char*,double,char*,double,int,xgOut*,char*),
	const char *dev_spec, char *file_or_dev, double maxdim,
	char *ti_fam,
	double ti_size,
	char* ax_fam,
	double ax_size,
	int doc_p);
void ho_dialog( Window parent, char *prog, xtb_data cookie);
void set_X( Window new_win, xgOut *out_info);
// void do_error(char const* err_text);		/* dialog.c */
void do_error(std::string const& err_text);		/* dialog.c */
void msg_box(char const* title, char const* text); /* dialog.c */

class Plotcurve : public Curve {
public:
	Par* xparam{nullptr};
	Par* yparam{nullptr};
	Par* zparam{nullptr};
	Par* wparam{nullptr};
	std::vector<bool> tagged;  // size altval.size()
	std::vector<std::vector<double> > evmatrix;  // real rep of complex eigenvectors
	std::string plotfile;
	int color{0};
	int style{0};
	std::string legend;

	// the one-and-only constructor
	Plotcurve(Curve*, const std::string& yname);
	// load plotfiles and read from flaps database
	static std::vector<Plotcurve*> getcurves();
	static int get_order();  // set in getcurves
};


#endif // PLOTCURVE_H
