//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
//////////////////////////////////////////////////////////////////////////////
// Animated Modes VIsualiZation
/////////////////////////////////////////////////////////////////////////////

#include <cerrno>
#include <chrono>
#include <ctime>
#include <exception>
#include <future>
#include <mqueue.h>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <iosfwd>

#include "wx/string.h"
#include "wx/prntbase.h"

#include "wx/wx.h"

#if !wxUSE_GLCANVAS
    #error "OpenGL required: set wxUSE_GLCANVAS to 1 and rebuild the library"
#endif
#include <GL/glu.h>
// #include <GL/glut.h>

#undef Complex    // X11 macro
#include "config.h"
#include "amviz.h"
#include "exim.h"
//!! #include "gl2ps.h"
#include "matrix.h"
#include "specs.h"
#include "trace.h"
#include "util.h"

using namespace std;

static string
helpmsg() {
	string rval{"This is a work in progress but here are some features:\n"
			 " - left mouse button will rotate the figure\n"
			 " - middle mouse button translates the figure\n"
			 " - right mouse button zooms:\n"
			 "     - moving it towards the center zooms out,\n"
			 "     - moving it away from the center zooms in\n"
			 " - File|Video creates an mp4 file\n"
			 " - File|Print creates a .ps file\n"
			 " - Edit|Preferences allows you to change some behavior\n"
			 " - The View menu has options for changing orientation and color\n"
			 " - There are 4 sliders on the left:\n"
			 "     - the top one controls the amplitude\n"
			 "     - the next one will stop the motion and adjust the phase\n"
			 "     - the third one will adjust the frequency\n"
			 "     - the fourth one will adjust weight the color scale\n"
			 " - there is a button on the left labeled either \"start\"\n"
			 "   or \"stop\" which will start or stop the motion\n"
			 " - another button on the left toggles the undeformed grid\n"
			 " - at the bottom left there is a box for choosing the mode\n"
			 "   to display from a list of previously-chosen modes\n"
			 " - at the bottom is a status bar\n"};
	rval += buildinfo() + "\n" + OGLVersion();
	Specs& sp = specs();
	if (sp.threads) {
		const int hwthreads = std::thread::hardware_concurrency()-1;
		rval += vastr(" (",hwthreads," threads)");
	}
	return rval;
}

static string
aboutmsg() {
	string rval{"Flaps Animated Mode Visualizer\n"};
	rval += buildinfo() + "\n" + OGLVersion();
	if (isWSL())
		rval += " (WSL)";
	rval += "\nRunning under ";
	if (isWayland())
		rval += "Wayland";
	else
		rval += "X11";

	Specs& sp = specs();
	if (sp.threads) {
		const int hwthreads = std::thread::hardware_concurrency()-1;
		rval += vastr(" (",hwthreads," threads)");
	}
	return rval;
}

void ani();

#ifdef NEVER // mv to util
static void
errormsg(const string& msg, bool quit) {
	if (Amviz::instance().frame != nullptr)
		wxMessageBox(msg, wxT("Error"), wxICON_INFORMATION);
	else
		flaps::error(msg);
	if (quit)
		exit(1);
}
#endif // NEVER // mv to util

// class Amviz implementation

Amviz&
Amviz::
instance() {
// the first time this is called the static instance of Amviz
// will be constructed
	static Amviz instance_;
	return instance_;
}

Amviz::
Amviz() {
// Amviz constructor (private): just check that either a vzid
// or a Universal file name are in Specs - data is loaded in kit()
	T_(Trace trc(1,"Amviz constructor");)
	Specs& sp = specs();

	if (sp.ufname.empty() && sp.vzid.empty())
		error = "neither a vzid nor a uf file name were specified";
	return;
}

Amviz::
~Amviz() {
	T_(Trace trc(2,"Amviz dtor");)
	// tell animate() to quit XXX necessary?
	quit = true;
}

void
amviz_kit() {
// Read an amviz kit: data necessary to create animated modes aka Fma
// 1) Either read a universal file (sp.ufname) or read an Fma from the
//    Flaps data directory with id "sp.vzid".
// 2) normalize the coordinates and shift to center in the bounding box
// 3) If this is a standalone run load the modes from the Fma::gct to
//    visualize
	T_(Trace trc(1,"amviz_kit");)
	Specs& sp = specs();
	Amviz& amviz{Amviz::instance()};

	// Get the 4 arrays necessary for visualization...
	try {
		// ...import a universal file...
		if (!sp.ufname.empty()) {
			amviz.fma = new Fma(sp.ufname);
		// ... or fetch from the flaps data directory
		} else {
			if (sp.vzid.empty()) {
				amviz.error = "neither a vzid nor a uf file name were specified";
			} else {
				amviz.fma = Fma::fetch(sp.vzid);
			}
		}
	} catch (std::exception& e) {
		amviz.error = vastr("amviz kit is not available: ",e.what());
	}

	if (!amviz.error.empty()) {
		T_(trc.dprint("quick return: ",amviz.error);)
		return;
	}
	T_(trc.dprint("got fma:\n",*amviz.fma);)

	// normalize the coordinates and shift the origin to the bounding-box center
	amviz.coord_norm_ = amviz.normalize_coord();
	T_(trc.dprint("normalized & shifted coords: ",amviz.fma->coords);)
	T_(trc.dprint("coord norm ",amviz.coord_norm_);)

	amviz.colormode_ = "abs value";
	amviz.colorskew_ = 5;  // no skewing
	// nodal displacements will be scaled by amplitude_, which is
	// changed in AmvizFrame::OnAmpSlider (amviz.cpp)
	amviz.nominal_scale_ = 0.005;
	amviz.amplitude_ = amviz.nominal_scale_*pow(2.0, double(Init_amp));
	amviz.drawundeformed_ = 0;

	// if this is a uf-file only run put the modes from the
	// uf file on the list of nodal_disp_
	//!! if (sp.respqname.empty())
	if (amviz.standalone) {
		size_t nr = 3*amviz.fma->nodes.size();
		size_t nc = amviz.fma->gct.size()/nr;
		vector<complex<double>> disp(nr);
		complex<double>* cele = amviz.fma->gct.data();
		for (size_t j=0; j<nc; j++) {
			//!! blas::copy(nr, &cele[IJ(0,nc-j-1,nr)], 1, &disp[0], 1);
			blas::copy(nr, &cele[IJ(0,j,nr)], 1, &disp[0], 1);
			double disp_norm = blas::nrm2(nr, &disp[0], 1);
			if (disp_norm == 0.0)
				flaps::warning("mode ",nc-j," is zero");
			amviz.nodal_disp_.push_back(disp);
			amviz.picked_names_.Add(wxString(vastr("mode ",j+1)));
			T_(trc.dprint("added mode ",j+1," to picked_names");)
		}
		// start visualizing on the first one
		amviz.picked_point_ = 0;
	}
}

void
amviz_run() {
// Create the AmvizFrame then call animate()
	T_(Trace trc(2,"amviz_run");)
	Amviz& amviz{Amviz::instance()};
	Specs& sp = specs();

	// did the Amviz ctor and kit() complete?
	if (!amviz.error.empty()) {
		// for now, but how to pop up a message if running from viz?
		cerr << amviz.error << endl;
		return;
	}

	try {
		// startup creates the frame...
		amviz.startup();
		// ... then run animate() either in a thread or with a timer
		if (sp.threads) {
			amviz.ani_future = std::async(ani);
		} else {
			amviz.frame->timer = new myTimer();
			// run the timer a little longer that the default to
			// allow Notify() to initialize OGL
			amviz.frame->timer->Start(5*Default_period);
		}
	} catch (std::exception& s) {
		amviz.error = vastr("could not start amviz: ",s.what());
		return;
	}
}	// amviz_run

int
Amviz::
add(const vector<complex<double>>& ev, const string& name) {
// Multiply the gct by the eigenvector (ev) and add the
// result to the list of nodal displacements; return
// the 0b index of the result in the nodal_disp_ array.
	T_(Trace trc(2,"Amviz::add");)
	complex<double> alpha{1.0};
	complex<double> beta{0.0};
	int m{3*nnodes()};
	int n{ngc()};
	// ev may be larger than the number of columns but not smaller - e.g.
	// if there are extra controls equations
	if ((int)ev.size() < n) {
		string exc{vastr("mismatch: ",ev.size()," gc but ",n," transform columns")};
		T_(trc.dprint("throwing exception: ",exc);)
		//!! throw runtime_error(exc);
		errormsg(exc, false);
	}
	complex<double>* gcxf = fma->gct.data();
	vector<complex<double>> disp(m);
	T_(trc.dprint("eigenvector:",ev);)
	blas::gemv ("n", m, n, alpha, gcxf, m, &ev[0], 1, beta, &disp[0], 1);

	//!! MM::exporter("nodal_disp","nodal_disp",disp.data(),m,1);
	//!! MM::exporter("gct","gct",gcxf,m,n);
	int rval = nodal_disp_.size();
	picked_point_ = rval;     // default picked point
	nodal_disp_.push_back(disp);
	//!! add_picked_name(name);
	picked_names_.Add(wxString(name));
	return rval;
}

vector<double>
Amviz::
coord_displacements(int omegat, int phase, vector<double>& colorvalue) {
// compute the (real) displacement at each node at time omegat
//   x[k] = Real(cos(omegat)+i*sin(omegat)*nodal_disp[k] + coord[k]
//       k = 0:n3-1
	T_(Trace trc(1,"Amviz::coord_displacements");)
	constexpr double deg2rad = atan(1.0)/45.0;
	constexpr double pi = 4.0*atan(1.0);
	static int visit{0};
	int n3 = 3*nnodes();

	vector<double> rval(n3, 0.0);
	Amviz& amviz{Amviz::instance()};

	if ((int)colorvalue.size() != nnodes())
		throw runtime_error(vastr("colorvalue wrong size (",colorvalue.size(),
				") should be ",nnodes()));

	// check for no modes yet
	vector<complex<double>> disp;
	if (nodal_disp_.empty())
		errormsg("no displacements have been computed", true);

	// copy the current picked nodal_disp_
	if (picked_point_ >= 0 && picked_point_ < (int)nodal_disp_.size())
		disp = nodal_disp_[picked_point_];
	else
		throw runtime_error(vastr("picked point ",picked_point_," is out of the range 0-",
			nodal_disp_.size()-1));

	if (visit == 0)
		T_(trc.dprint("omegat ",omegat,", disp = ",disp);)
	visit++;

	// normalize the complex nodal displacements to amplitude_*coord_norm_
	// use the infinity norm and normalize colorvalues with it below
	double dispnorm = abs(disp[blas::icamax(n3, &disp[0], 1)-1]);
	double eps{1.0e-6};
	// return zero if tiny displacements
	if (dispnorm*amplitude_ <= eps) {
		return rval;
	}
	complex<double> ct = amplitude_*coord_norm_/dispnorm;
	blas::scal(n3, ct, &disp[0], 1);

	double omtp{(omegat + phase)*deg2rad};
	double scr{cos(omtp)};
	double sci{sin(omtp)};
	complex<double> csc(scr, sci);

	T_(trc.dprint("omegat+phase ",omtp," rad, (cos,sin) ",csc);)

	vector<complex<double>> cd(3);
	for (int i=0; i<nnodes(); i++) {
		int k = 3*i;
		cd[0] = disp[k]*csc;
		cd[1] = disp[k+1]*csc;
		cd[2] = disp[k+2]*csc;
		double x = cd[0].real();
		double y = cd[1].real();
		double z = cd[2].real();
		rval[k] = x + fma->coords[k];
		rval[k+1] = y + fma->coords[k+1];
		rval[k+2] = z + fma->coords[k+2];
		double normd = sqrt(x*x+y*y+z*z);
		int idx = blas::icamax(3,&cd[0],1)-1;
		// set a color value for this node
		if (amviz.colormode_ == "abs value") {
			colorvalue[i] = normd;
		} else if (amviz.colormode_ == "phase") {
			colorvalue[i] = atan2(cd[idx].imag(), cd[idx].real())/(2.0*pi) + 0.5;
		} else if (Amviz::colormode() == "abs phase") {
			colorvalue[i] = abs(atan2(cd[idx].imag(), cd[idx].real())/pi);
		}
	}

	// scale colorvalues so that the largest colorvalue over an
	// oscillation cycle is 1.0. The largest (normalized) displacement 
	// over a cycle is amplitude_*coord_norm_
	blas::scal(nnodes(), 1.0/(amplitude_*coord_norm_), &colorvalue[0], 1);

	// skewed color variation
	if (amviz.colorskew_ != 5) {
		double cs = ((double)amviz.colorskew_)/5.0;
		for (int i=0; i<nnodes(); i++) {
			double ci = colorvalue[i];
			if (ci > 0.0)
				ci = pow(abs(10.0*ci), cs)/10.0;
			else if (ci < 0.0)
				ci = -pow(abs(10.0*ci), cs)/10.0;
			colorvalue[i] = ci;
		}
	}

	T_(trc.dprint("nodal disp:",rval);)
	return rval;
}

static void
color_values (double t, float& red, float& green, float& blue) {
// given a number (t) between 0 and 1, set the rgb values
// so that they vary linearly with t=0 blue, t=0.5 green,
// and t=1 red
	T_(Trace trc(3,"color_values");)

	if (t > 1.0)
		t = 1.0;
	if (t < -1.0)
		t = -1.0;

	if (t < 0.5) {
		blue = 1.0 - 2.0*t;
		green = 2.0*t;
		red = 0.0;
	} else {
		blue = 0.0;
		green = -2.0*t + 2.0;
		red = 2.0*t - 1.0;
	}
}

void
Amviz::
render(int omegat, int phase) {
	T_(Trace trc(2,"render ",omegat);)

	// if colormode == "none" use the reverse of the bg color
	float red = 1.0 - BGRed;
	float green = 1.0 - BGGreen;
	float blue = 1.0 - BGBlue;

    /* clear color and depth buffers */
	glClearColor(BGRed, BGGreen, BGBlue, BGAlpha);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
#ifdef NEVER // necessary?
	// antialiasing
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// glBlendFunc(1.0, 0.0);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);
#endif // NEVER // necessary?

	// colorvalue is a vector of color values ranging from 0-1 to
	// be used by color_values to give rgb values for GL
	vector<double> colorvalue(fma->nodes.size(), 0.0);

	// compute the coordinates of the nodes at time omegat+phase
	vector<double> dispnodes = this->coord_displacements(omegat,
			phase, colorvalue);

	//! static double* vertexEnergy = 0;

	// we always have segidx which are pairs of 0b indices into the
	// nodenumbers_ array, i.e. the ordinal numbers of nodes.
	// The start of (x,y,z) data for ordinal node j in coordinates,
	// gct, and dispnodes arrays is thus 3*j
	// This is a more useful format for connectivity than either Universal
	// File's penlift format or GL's segment format.
	glLineWidth(2.0);
	glBegin(GL_LINES);
	for (size_t i=0; i<fma->segidx.size(); i++) {
		int k = 3*fma->segidx[i]; // start of this node
		// the Nodes vector is assumed to be in the same order as dispnodes
		float x = dispnodes[k];
		float y = dispnodes[k+1];
		float z = dispnodes[k+2];

		// set the color for this node and the next
		color_values(colorvalue[k/3], red, green, blue);

		glColor3f(red, green, blue);
		glVertex3f(x, y, z);
		i++;
		k = 3*fma->segidx[i];
		x = dispnodes[k];
		y = dispnodes[k+1];
		z = dispnodes[k+2];
		glVertex3f(x, y, z);
	}
	// draw the coord system if requested
	if (this->draw_coord_sys()) {
		vector<float> x;
		for (size_t i=0; i<3; i++) {
			x = origin_;
			color_values(1.0-i*0.5, red, green, blue);
			glColor3f(red, green, blue);
			glVertex3f(x[0], x[1], x[2]);
			x[i] += 0.1;
			glVertex3f(x[0], x[1], x[2]);
		}
	}
	glEnd();
	// draw the undeformed grid if requested
	if (this->drawundeformed_) {
		double* coord = fma->coords.data();
		// red = green = blue = 0.5;  // gray
		red = 1.0 - BGRed;
		green = 1.0 - BGGreen;
		blue = 1.0 - BGBlue;
		glColor3f(red, green, blue);
		glLineWidth(0.2);
		glBegin(GL_LINES);
		for (size_t i=0; i<fma->segidx.size(); i++) {
			int k = 3*fma->segidx[i];
			float x = coord[k];
			float y = coord[k+1];
			float z = coord[k+2];
			glVertex3f(x, y, z);
			i++;
			k = 3*fma->segidx[i];
			x = coord[k];
			y = coord[k+1];
			z = coord[k+2];
			glVertex3f(x, y, z);
		}
		glEnd();
		// glDisable(GL_LINE_STIPPLE);
	}

	// glFlush();
	// wxGLCanvas::SwapBuffers();
} // Amviz::render()

string
Amviz::
colormode(const string& newmode) {
	string rval = instance().colormode_;
	if (!newmode.empty())
		instance().colormode_ = newmode;
	return rval;
}

int
Amviz::
colorskew(int newskew) {
	int rval = instance().colorskew_;
	if (newskew >= 0)
		instance().colorskew_ = newskew;
	return rval;
}

int
Amviz::
npicked() {
	return instance().nodal_disp_.size();
}

int
Amviz::
picked_point(int newpoint) {
	int rval = instance().picked_point_;
	if (newpoint >= 0)
		instance().picked_point_ = newpoint;
	return rval;
}

double
Amviz::
amplitude(int newamp) {
	Amviz& amviz{Amviz::instance()};
	double rval = amviz.amplitude_;
	if (newamp >= 0)
		amviz.amplitude_ = amviz.nominal_scale_*pow(2.0, double(newamp));
	return rval;
}

int
Amviz::
draw_undeformed(int newund) {
// change drawundeformed_ setting, return current
	int rval = instance().drawundeformed_;
	if (newund >= 0)
		instance().drawundeformed_ = newund;
	return rval;
}

int
Amviz::
draw_coord_sys(int newcs) {
// change drawcs_ setting, return current
	int rval = instance().drawcs_;
	if (newcs >= 0)
		instance().drawcs_ = newcs;
	return rval;
}

double
Amviz::
normalize_coord() {
// Normalize the nodal coordinates by scaling so that the largest
// dimension of the enclosing box is 1.0 and the enclosing box
// is centered on the origin
	//!! vector<double>& coord = coord_->data();
	double rval{1.0};

	double big = std::numeric_limits<double>::max();
	double xmax=-big, xmin=big;
	double ymax=-big, ymin=big;
	double zmax=-big, zmin=big;
	for (int i=0; i<nnodes(); i++) {
		int k = 3*i;
		xmin = std::min(xmin, fma->coords[k]);
		ymin = std::min(ymin, fma->coords[k+1]);
		zmin = std::min(zmin, fma->coords[k+2]);
		xmax = std::max(xmax, fma->coords[k]);
		ymax = std::max(ymax, fma->coords[k+1]);
		zmax = std::max(zmax, fma->coords[k+2]);
	}
	double scale = abs(xmax - xmin);
	scale = std::max(scale, abs(ymax - ymin));
	scale = std::max(scale, abs(zmax - zmin));
	scale *= 2.0;
	double xc = (xmax + xmin)/2.0;
	double yc = (ymax + ymin)/2.0;
	double zc = (zmax + zmin)/2.0;
	for (int i=0; i<nnodes(); i++) {
		int k = 3*i;
		fma->coords[k] -= xc;
		fma->coords[k+1] -= yc;
		fma->coords[k+2] -= zc;
		fma->coords[k] /= scale;
		fma->coords[k+1] /= scale;
		fma->coords[k+2] /= scale;
	}
	// save the origin for drawing the coordinate system
	origin_.clear();
	origin_.push_back(-xc/scale);
	origin_.push_back(-yc/scale);
	origin_.push_back(-zc/scale);
	return rval;
}

// end class Amviz implementation
class Quaternion {
public:
	string id;
	vector<float> quat;
	Quaternion(const string& ident, vector<float> q) : id(ident), quat(q) {}
};

static vector<float> Default_quat;
static vector<float> airplane_quat{0.396197, -0.133608, -0.491376, 0.76402};
static vector<float> Isoquat{0.180731, 0.356605, 0.827727, 0.393747};
static vector<float> Xquat{0.5, 0.5, 0.5, 0.5};
static vector<float> Yquat{0.0, 0.707106781, 0.707106781, 0.0};
static vector<float> Zquat{0.0, 0.0, 1.0, 0.0};

// the sizes of the frame, canvas, toolbar and menubar
AmvizSizes amsizes;

AmvizSizes::
AmvizSizes() {
// AmvizSizes ctor
	display = getDisplaySize();
	int side = 0.8*std::min(display.x, display.y);
#ifdef SQUARE_CANVAS
	// square canvas
	canvas = wxSize(side,side);
	//!! canvas = wxSize(0.7*display.x, 0.7*display.y);
	menu = wxSize(canvas.x, 27);
	toolbar = wxSize(canvas.x, 70);
	frame = wxSize(canvas.x, canvas.y+toolbar.y+menu.y);
#else // SQUARE_CANVAS or Frame?
	// square frame
	frame = wxSize(side, side);
	menu = wxSize(frame.x, 27);
	toolbar = wxSize(frame.x, 70);
	canvas = wxSize(frame.x, frame.y - menu.y - toolbar.y);
#endif // SQUARE_CANVAS
}

float OriginX = 0.0;
float OriginY = 0.0;

// Default viewing parameters for glOrtho, changed by mouse
// motion with right mouse button down
// The nodes have been scaled so the bounding box is 1
// and centered, so make the viewing area a little larger
float ViewLeft{-0.6};
float ViewRight{0.6};
float ViewBottom{-0.6};
float ViewTop{0.6};
float ViewNear{-3.0};
float ViewFar{3.0};

bool ViewDefault{true};
bool ViewQuad = false;
bool ViewX = false;
bool ViewY = false;
bool ViewZ = false;

string Printer;

int Omegat = 0;
int Phase = 0;

wxArrayString ModeStrings;

bool
Amviz::
startup() {
// Create an AmvizFrame with a GLCanvas
// Takes the place of an OnInit() when run in viz
	T_(Trace trc(2,"Amviz::startup");)
	Specs& sp = specs();

	wxSize dpy = getDisplaySize();
	wxPoint pos(dpy.x-amsizes.frame.x, dpy.y-amsizes.frame.y);

	// title for the frame
	string title{"Animated Modes Visualization"};
	if (!sp.ufname.empty())
		title = sp.ufname;
	if (sp.threads)
		title += " (threads)";
	else
		title += " (timer)";
   // Create the main frame window...
	this->frame = new AmvizFrame(nullptr, title, pos, amsizes.frame);
	// ... and display it
	this->frame->Show(true);

    return true;
}	// Amviz::startup

// AmvizFrame

BEGIN_EVENT_TABLE(AmvizFrame, wxFrame)
    EVT_MENU(wxID_OPEN, AmvizFrame::OnMenuFileOpen)
    EVT_MENU(wxID_PRINT, AmvizFrame::OnMenuFilePrint)
    EVT_MENU(wxID_EXIT, AmvizFrame::OnMenuFileExit)
    // EVT_MENU(ID_EDIT_PREF, AmvizFrame::OnMenuEditPref)
    // EVT_MENU(ID_HELP_ABOUT, AmvizFrame::OnMenuHelpAbout)
    // EVT_MENU(ID_HELP_HELP, AmvizFrame::OnMenuHelpHelp)
    EVT_MENU(ID_VIEW_DEFAULT, AmvizFrame::OnMenuViewDefault)
    EVT_MENU(ID_VIEW_X, AmvizFrame::OnMenuViewX)
    EVT_MENU(ID_VIEW_Y, AmvizFrame::OnMenuViewY)
    EVT_MENU(ID_VIEW_Z, AmvizFrame::OnMenuViewZ)
    EVT_MENU(ID_VIEW_I, AmvizFrame::OnMenuViewI)
    EVT_MENU(ID_VIEW_QUAD, AmvizFrame::OnMenuViewQuad)
    EVT_MENU(ID_VIEW_COLOR_ABS, AmvizFrame::OnMenuViewColorAbs)
    EVT_MENU(ID_VIEW_COLOR_VAL, AmvizFrame::OnMenuViewColorVal)
    EVT_MENU(ID_VIEW_COLOR_PHASE, AmvizFrame::OnMenuViewColorPhase)
    EVT_MENU(ID_VIEW_COLOR_PHASE_ABS, AmvizFrame::OnMenuViewColorPhaseAbs)
    EVT_MENU(ID_VIEW_COLOR_ENERGY, AmvizFrame::OnMenuViewColorEnergy)
    EVT_MENU(ID_VIEW_COLOR_NONE, AmvizFrame::OnMenuViewColorNone)

	EVT_COMMAND_SCROLL(ID_AMP_SLIDER, AmvizFrame::OnAmpSlider)
	EVT_COMMAND_SCROLL(ID_FREQ_SLIDER, AmvizFrame::OnFreqSlider)
	EVT_COMMAND_SCROLL(ID_COLOR_SLIDER, AmvizFrame::OnColorSlider)

	EVT_BUTTON(ID_START_BUTTON, AmvizFrame::OnStartButton)
	EVT_BUTTON(ID_UNDEF_BUTTON, AmvizFrame::OnUndefButton)
	EVT_BUTTON(ID_CS_BUTTON, AmvizFrame::OnCSButton)
	EVT_BUTTON(ID_MODE_NUMBER, AmvizFrame::OnModeNumber)
END_EVENT_TABLE()

// AmvizFrame constructor
#include "amviz.xpm"
AmvizFrame::AmvizFrame(wxWindow *frame, const wxString& title, const wxPoint& pos,
    const wxSize& size, long style)
    : wxFrame(frame, wxID_ANY, title, pos, size, style) {

	T_(Trace trc(2,"AmvizFrame ctor");)

	// an icon for the app bar
	SetIcon(wxIcon(amviz_xpm));

	// add menubar
	createControls();


	wxBoxSizer* topsizer = new wxBoxSizer(wxVERTICAL);
	SetSizer(topsizer);

	wxPoint tbpos{wxDefaultPosition};
	wxPanel* toolbar = new wxPanel(this, wxID_ANY, tbpos, amsizes.toolbar);
	wxBoxSizer* toolbarsizer = new wxBoxSizer(wxHORIZONTAL);
	toolbar->SetSizer(toolbarsizer);

	int prop{0};	// proportion

	// amplitude & freq sliders need to be in boxes in order to
	// put titles on them
	int sstyle{wxEXPAND|wxTOP|wxBOTTOM};
	int smarg{5};
	// amplitude slider: values from 0 - 10
	wxPanel* ampbox = new wxPanel(toolbar);
	wxBoxSizer* ampsizer = new wxBoxSizer(wxVERTICAL);
	ampbox->SetSizer(ampsizer);
	ampSlider = new wxSlider(ampbox, ID_AMP_SLIDER, Init_amp, 0, 10);
	wxStaticText* label = new wxStaticText(ampbox, wxID_ANY, "Amplitude");
	ampSlider->SetToolTip(wxT("amplitude"));
	// ampsizer->Add(ampSlider, prop, wxEXPAND|wxALL, 5);
	ampsizer->Add(ampSlider, prop, wxEXPAND);
	ampsizer->Add(label, prop, wxALIGN_CENTER_HORIZONTAL);
	toolbarsizer->Add(ampbox, prop, sstyle, smarg);

	// frequency slider
	wxPanel* freqbox = new wxPanel(toolbar);
	wxBoxSizer* freqsizer = new wxBoxSizer(wxVERTICAL);
	freqbox->SetSizer(freqsizer);
	// freq slider values: 5:200
	int freq{50};	// initial freq slider
	freqSlider = new wxSlider(freqbox, ID_FREQ_SLIDER, freq, 5, 200);
	label = new wxStaticText(freqbox, wxID_ANY, "Frequency");
	freqSlider->SetToolTip(wxT("frequency"));
	freqsizer->Add(freqSlider, prop, wxEXPAND);
	freqsizer->Add(label, prop, wxALIGN_CENTER_HORIZONTAL);
	toolbarsizer->Add(freqbox, prop, sstyle, smarg);

	// Buttons: stop/start, undeformed, coord sys, modes, quit
	int buttonstyle{wxEXPAND|wxTOP|wxBOTTOM};
	int bmarg{20};

	// start/stop button
	startButton = new wxButton(toolbar,
		myBind(wxEVT_BUTTON, &AmvizFrame::OnStartButton, this), "Stop");
	startButton->SetToolTip(wxT("start/stop oscillations"));
	toolbarsizer->Add(startButton, prop, buttonstyle, bmarg);
	toolbarsizer->AddStretchSpacer();

	// toggle undeformed grid
	undeformedButton = new wxButton(toolbar,
		myBind(wxEVT_BUTTON, &AmvizFrame::OnUndefButton, this), "undef");
	undeformedButton->SetToolTip(wxT("toggle undeformed grid"));
	toolbarsizer->Add(undeformedButton, prop, buttonstyle, bmarg);
	toolbarsizer->AddStretchSpacer();

	// toggle the coordinate system
	CSButton = new wxButton(toolbar,
		myBind(wxEVT_BUTTON, &AmvizFrame::OnCSButton, this), "coord sys");
	CSButton->SetToolTip(wxT("toggle coordinate system"));
	toolbarsizer->Add(CSButton, prop, buttonstyle, bmarg);
	toolbarsizer->AddStretchSpacer();

	// mode list
	wxButton* modeNumber = new wxButton(toolbar,
		myBind(wxEVT_BUTTON, &AmvizFrame::OnModeNumber, this), "modes");
	modeNumber->SetToolTip(wxT("which mode to display"));
	toolbarsizer->Add(modeNumber, prop, buttonstyle, bmarg);
	toolbarsizer->AddStretchSpacer();

	// quit button
	wxButton* quitbutton = new wxButton(toolbar,
		myBind(wxEVT_BUTTON, &AmvizFrame::OnMenuFileExit, this), "quit");
	quitbutton->SetToolTip(wxT("exit amviz"));
	toolbarsizer->Add(quitbutton, prop, buttonstyle, bmarg);
	toolbarsizer->AddStretchSpacer();

	toolbar->Layout();
	toolbar->Show(true);
	topsizer->Add(toolbar);

	// the OpenGL canvas
	long glstyle{wxBORDER_SIMPLE};
	if (isWSL()) {
		// XXX this may be the wsl problem
		int attlist[3];
		int i = 0;
		attlist[i++] = WX_GL_DOUBLEBUFFER;
		//!! attlist[i++] = WX_GL_RGBA;
		attlist[i++] = 0;
		m_canvas = new FGLCanvas(this, wxID_ANY, attlist, wxDefaultPosition,
			amsizes.canvas, glstyle, wxT("Amviz"));
		flaps::info("running under WSL");
	} else {
		wxGLAttributes dispAttrs;
		dispAttrs.PlatformDefaults().DoubleBuffer().EndList();
		m_canvas = new FGLCanvas(this, dispAttrs, wxID_ANY, wxDefaultPosition,
				amsizes.canvas, glstyle, wxT("Amviz"));
	}

	topsizer->Add(m_canvas);

	topsizer->Layout();

	m_canvas->Show(true);
	Show(true);

	T_(trc.dprint("topsizer after layout: ",sizer2string(*topsizer));)
}	// AmvizFrame ctor

bool
AmvizFrame::
isRunning() {
	Specs& sp = specs();
	if (sp.threads)
		return true;	// XXX needs work
	else
		return (timer->IsRunning());
}

void
AmvizFrame::
createControls() {
	T_(Trace trc(2,"AmvizFrame::createControls");)

    // Make the "File" menu
    wxMenu *fileMenu = new wxMenu;
    //fileMenu->Append(wxID_OPEN, wxT("&Open..."));
    //fileMenu->AppendSeparator();
    fileMenu->Append(wxID_PRINT, wxT("&Print\tALT-P"));
	// using lambdas
	fileMenu->Append(funBind(wxEVT_MENU, [&](wxCommandEvent& ev) {
		// just call the FGLCanvas xfig creator
		m_canvas->xfig(); }, this), "&Xfig\tALT-X");

	fileMenu->Append(funBind(wxEVT_MENU, [&](wxCommandEvent& ev) {
		// just call the FGLCanvas video creator
		m_canvas->video(); }, this), "&Video\tALT-v");
    fileMenu->AppendSeparator();
	 // using std wx ID
    fileMenu->Append(wxID_EXIT, wxT("E&xit\tq"));
#ifdef NEVER // empty
	 // make the "Edit" menu
    wxMenu *editMenu = new wxMenu;
#endif // NEVER // empty
	 // make the "View" menu
    wxMenu *viewMenu = new wxMenu;
    viewMenu->Append(ID_VIEW_DEFAULT, wxT("&Default\td"));
    viewMenu->Append(ID_VIEW_X, wxT("&X\tx"));
    viewMenu->Append(ID_VIEW_Y, wxT("&Y\ty"));
    viewMenu->Append(ID_VIEW_Z, wxT("&Z\tz"));
    viewMenu->Append(ID_VIEW_I, wxT("&I\ti"));
    viewMenu->Append(-1, wxT("-----"));   // supposed to be a separator
	 wxMenu* colorStyleMenu = new wxMenu;
    viewMenu->Append(ID_VIEW_COLOR, wxT("&Color\tALT-C"), colorStyleMenu);
	 colorStyleMenu->Append(ID_VIEW_COLOR_VAL, wxT("value"));
	 colorStyleMenu->Append(ID_VIEW_COLOR_ABS, wxT("abs value"));
	 colorStyleMenu->Append(ID_VIEW_COLOR_PHASE, wxT("phase"));
	 colorStyleMenu->Append(ID_VIEW_COLOR_PHASE_ABS, wxT("abs phase"));
	 colorStyleMenu->Append(ID_VIEW_COLOR_ENERGY, wxT("energy"));
	 colorStyleMenu->Append(ID_VIEW_COLOR_NONE, wxT("none"));

    // Make the "Help" menu
    wxMenu *helpMenu = new wxMenu;
#ifdef NEVER // use lambdas & funBind
    helpMenu->Append(ID_HELP_HELP, wxT("&Help..."));
    helpMenu->Append(ID_HELP_ABOUT, wxT("&About..."));
#else // NEVER // use lambdas & funBind
    helpMenu->Append(funBind(wxEVT_MENU, [](wxCommandEvent& ev) {
	 	wxMessageBox(helpmsg(), "Help"); }, this), "Help");
    helpMenu->Append(funBind(wxEVT_MENU, [](wxCommandEvent& ev) {
	 	wxMessageBox(aboutmsg(), "About"); }, this), "About");
#endif // NEVER // use lambdas & funBind

	 // Add these menus to the menu bar
    wxMenuBar *menuBar = new wxMenuBar;
    menuBar->Append(fileMenu, wxT("&File"));
#ifdef NEVER // empty
    menuBar->Append(editMenu, wxT("&Edit"));
#endif // NEVER // empty
    menuBar->Append(viewMenu, wxT("&View"));
    menuBar->Append(helpMenu, wxT("&Help"));
    SetMenuBar(menuBar);

	 statusBar = new wxStatusBar(this, wxID_ANY, wxST_SIZEGRIP);
	 this->SetStatusBar(statusBar);
	 statusBar->SetStatusText(wxT("Ready"), 0);

} // AmvizFrame::createControls()

bool
AmvizFrame::
ShowToolTips() { return true; }

void
AmvizFrame::
OnAmpSlider(wxScrollEvent& event) {
	if (this->ampSlider) {
		//! DispScale = NomDispScale*this->ampSlider->GetValue();
		int newscale = this->ampSlider->GetValue();
		// double newscale = NomDispScale*this->ampSlider->GetValue();
		Amviz& amviz{Amviz::instance()};
		amviz.amplitude(newscale);
	}
}

void
AmvizFrame::
OnFreqSlider(wxScrollEvent& event) {
// Respond to movement of the freq slider; if the new frequency
// is zero just stop the timer. Values go from 5-200
	T_(Trace trc(2,"OnFreqSlider");)
	Specs& sp = specs();
	int freq = this->freqSlider->GetValue();
	if (sp.threads)
		omegat_incr = freq/5;
	else {
		int period = Default_period;
		timer->Stop();
		if (freq > 0) {
			period = (10*period)/freq;
			T_(trc.dprint("new period: ",period);)
		}
		timer->Start(period);
	}
}

void
AmvizFrame::
OnColorSlider(wxScrollEvent& event) {
	if (this->colorSlider) {
		//! ColorSkew = (float)this->colorSlider->GetValue()/5.0;
		Amviz& amviz{Amviz::instance()};
		amviz.colorskew(this->colorSlider->GetValue());
	}
}

void
AmvizFrame::
OnStartButton(wxCommandEvent& event) {
	Specs& sp = specs();
	// toggle the oscillations
	if (sp.threads) {
		Amviz& amviz{Amviz::instance()};
		amviz.stop = !amviz.stop;
		if (amviz.stop)
			startButton->SetLabel(wxT("start"));
		else
			startButton->SetLabel(wxT("stop"));
	} else {
		static int interval{Default_period};
		if (isRunning()) {
			interval = timer->GetInterval();
			timer->Stop();
			startButton->SetLabel(wxT("start"));
		} else {
			timer->Start(interval);
			startButton->SetLabel(wxT("stop"));
		}
	}
}

void
AmvizFrame::
OnUndefButton(wxCommandEvent& event) {
	//! DrawUndeformed = !DrawUndeformed;
	Amviz& amviz{Amviz::instance()};
	int current = amviz.draw_undeformed();
	if (current == 0)
		amviz.draw_undeformed(1);
	else
		amviz.draw_undeformed(0);
}


void
AmvizFrame::
OnCSButton(wxCommandEvent& event) {
	//! Drawcs = !Drawcs;
	Amviz& amviz{Amviz::instance()};
	int current = amviz.draw_coord_sys();
	if (current == 0)
		amviz.draw_coord_sys(1);
	else
		amviz.draw_coord_sys(0);
}


// File|Open... command
void
AmvizFrame::
OnMenuFileOpen( wxCommandEvent& WXUNUSED(event) ) {
    wxString filename = wxFileSelector(wxT("Choose DXF Model"), wxT(""), wxT(""), wxT(""),
#if wxUSE_ZLIB
        wxT("DXF Drawing (*.dxf;*.dxf.gz)|*.dxf;*.dxf.gz|All files (*.*)|*.*"),
#else
        wxT("DXF Drawing (*.dxf)|*.dxf)|All files (*.*)|*.*"),
#endif
        wxFD_OPEN);
    if (!filename.IsEmpty())
    {
        m_canvas->Refresh(false);
    }
}

// File|Print... command
void
AmvizFrame::
OnMenuFilePrint( wxCommandEvent& WXUNUSED(event) ) {
	errormsg("printing has not been implemented yet",false);
	// print() uses wxPrinter stuff, printPS uses gl2ps
	// print();
	// printPS();
}

wxPrintDialogData g_printDialogData;

void
AmvizFrame::
print() {
#ifdef NEVER // no timer
	int period = m_timer.GetInterval();
	// stop the oscillations
	m_timer.Stop();
#endif // NEVER // no timer
	
	wxPrinter printer (&g_printDialogData);
	MyPrintout printout(wxT("My Printout"));

	if (!printer.Print(this, &printout, true)) {
		if (wxPrinter::GetLastError() == wxPRINTER_ERROR)
			wxMessageBox(wxT("printing error"), wxT("Printing"), wxOK);
		else
			wxMessageBox(wxT("you cancelled printing"),
				wxT("Printing"), wxOK);
	} else {
		g_printDialogData = printer.GetPrintDialogData();
	}

#ifdef NEVER // no timer
	// and re-start
	m_timer.Start(period);
#endif // NEVER // no timer
}

void
AmvizFrame::
printPS() {
#ifdef NEVER // no timer
	int period = m_timer.GetInterval();
	// stop the oscillations
	m_timer.Stop();
	// call the FGLCanvas printer
	m_canvas->print();
	// and re-start
	m_timer.Start(period);
#else // NEVER // no timer
	// call the FGLCanvas printer
	m_canvas->print();
#endif // NEVER // no timer
}

// File|Exit command or simply type "q"
void
AmvizFrame::
OnMenuFileExit( wxCommandEvent& WXUNUSED(event) ) {
	T_(Trace trc(2,"AmvizFrame::OnMenuFileExit");)
	Specs& sp = specs();
	// write out the current quaternion to use as the initial
	// position in FGlCanvas constructor
	// cerr << "current quaternion: " << m_canvas->gldata()->quat << endl;

	if (sp.threads) {
		// tell animate() to return
		Amviz& amviz{Amviz::instance()};
		amviz.quit = true;
		snooze(5.0e+8);
	} else
		timer->Stop();

#ifdef NEVER // try Destroy()
	// true is to force the frame to close
	if (!closed)
		Close(true);
	closed = true;
#else // NEVER // try Destroy()
	Destroy();
#endif // NEVER // try Destroy()
	// kill viz also?
	//exit(0);
}

#ifdef NEVER // needs work
void
AmvizFrame::
OnMenuEditPref( wxCommandEvent& WXUNUSED(event) ) {
	PrefDialog dialog(this);
	if (dialog.ShowModal() == wxID_OK) {
		ostringstream os;
		os << dialog.m_printer;
		Printer = os.str();
		// cout << "got printer = " << Printer << endl;
		if (dialog.m_print == wxT("To File")) {
			Printer.clear();
			// cout << "not sending to printer\n";
		}
	}
}
#endif // NEVER // needs work

void
AmvizFrame::
OnMenuViewDefault(wxCommandEvent& WXUNUSED(event)) {
	ViewX = false;
	ViewY = ViewZ = false;
	ViewQuad = false;
#ifdef NEVER // unnecessary to stop/start
	int period = m_timer.GetInterval();
	m_timer.Stop();
	m_canvas->gldata()->quat = Default_quat;
	m_timer.Start(period);
#else // NEVER // unnecessary to stop/start
	m_canvas->gldata()->quat = Default_quat;
#endif // NEVER // unnecessary to stop/start
	m_canvas->Refresh(false);
}

void
AmvizFrame::
OnMenuViewX(wxCommandEvent& WXUNUSED(event)) {
	m_canvas->gldata()->quat = Xquat;
	m_canvas->Refresh(false);
}

void
AmvizFrame::
OnMenuViewY(wxCommandEvent& WXUNUSED(event)) {
	m_canvas->gldata()->quat = Yquat;
	m_canvas->Refresh(false);
}

void
AmvizFrame::
OnMenuViewZ(wxCommandEvent& WXUNUSED(event)) {
	m_canvas->gldata()->quat = Zquat;
	m_canvas->Refresh(false);
}

void
AmvizFrame::
OnMenuViewI(wxCommandEvent& WXUNUSED(event)) {
	m_canvas->gldata()->quat = Isoquat;
	m_canvas->Refresh(false);
}

void
AmvizFrame::
OnMenuViewQuad(wxCommandEvent& WXUNUSED(event)) {
	if (!ViewQuad) {
		ViewQuad = true;
		m_canvas->viewQuad(Omegat, Phase);
	}
}

void
AmvizFrame::
OnMenuViewColorAbs(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_ABS;
	Amviz& amviz{Amviz::instance()};
	amviz.colormode("abs value");
}

void
AmvizFrame::
OnMenuViewColorVal(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_VAL;
	Amviz& amviz{Amviz::instance()};
	amviz.colormode("value");
}

void
AmvizFrame::
OnMenuViewColorPhase(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_PHASE;
	Amviz& amviz{Amviz::instance()};
	amviz.colormode("phase");
}

void
AmvizFrame::
OnMenuViewColorPhaseAbs(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_PHASE_ABS;
	Amviz& amviz{Amviz::instance()};
	amviz.colormode("abs phase");
}

void
AmvizFrame::
OnMenuViewColorEnergy(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_ENERGY;
	Amviz& amviz{Amviz::instance()};
	amviz.colormode("energy");
}

void
AmvizFrame::
OnMenuViewColorNone(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_NONE;
	Amviz& amviz{Amviz::instance()};
	amviz.colormode("none");
}

void
AmvizFrame::
OnModeNumber( wxCommandEvent& event ) {
// popup a dialog for selecting which mode to display
	T_(Trace trc(2,"OnModeNumber");)
	Amviz& amviz{Amviz::instance()};
	wxArrayString& picks = amviz.picked_names();
	wxSingleChoiceDialog modedlg(this,"Mode","choose",picks);
	// mark the current pick
	int current = amviz.picked_point();
	if (current >= 0)
		modedlg.SetSelection(current);
	// show the dialog, get selection
	int stat = modedlg.ShowModal();	// ??
	T_(trc.dprint("ShowModal returned ",stat);)
	if (stat == wxID_CANCEL)
		return;
	int sel = modedlg.GetSelection();
	amviz.picked_point(sel);
}

// Help|Help... command
void
AmvizFrame::
OnMenuHelpHelp( wxCommandEvent& WXUNUSED(event) ) {
    wxMessageBox(wxT("This is a work in progress but here are some features:\n"
			 " - left mouse button will rotate the figure\n"
			 " - middle mouse button translates the figure\n"
			 " - right mouse button zooms:\n"
			 "     - moving it towards the center zooms out,\n"
			 "     - moving it away from the center zooms in\n"
			 " - File|Video creates an mp4 file\n"
			 " - File|Print creates a .ps file\n"
			 " - Edit|Preferences allows you to change some behavior\n"
			 " - The View menu has options for changing orientation and color\n"
			 " - There are 4 sliders on the left:\n"
			 "     - the top one controls the amplitude\n"
			 "     - the next one will stop the motion and adjust the phase\n"
			 "     - the third one will adjust the frequency\n"
			 "     - the fourth one will adjust weight the color scale\n"
			 " - there is a button on the left labeled either \"start\"\n"
			 "   or \"stop\" which will start or stop the motion\n"
			 " - another button on the left toggles the undeformed grid\n"
			 " - at the bottom left there is a box for choosing the mode\n"
			 "   to display from a list of previously-chosen modes\n"
			 " - at the bottom is a status bar\n"
		), wxT("Help"), wxICON_INFORMATION);

}

// Help|About... command
void
AmvizFrame::
OnMenuHelpAbout( wxCommandEvent& WXUNUSED(event) ) {
    wxMessageBox(wxT("   Animated Mode Visualizer\n"
			"intended to be run either by right-clicking on\n"
			"curves in viz or as a standalone program (amviz)\n"
			"taking input from Universal-formatted files.\n"
			"In either case a set of nodes, connectivities,\n"
			"and modes must be available\n"
		), wxT("About"), wxICON_INFORMATION);

}

// FGLCanvas

BEGIN_EVENT_TABLE(FGLCanvas, wxGLCanvas)
    EVT_SIZE(FGLCanvas::OnSize)
    EVT_PAINT(FGLCanvas::OnPaint)
    EVT_ERASE_BACKGROUND(FGLCanvas::OnEraseBackground)
    EVT_MOUSE_EVENTS(FGLCanvas::OnMouse)
END_EVENT_TABLE()

void
ani() {
	T_(Trace trc(2,"ani");)
	Amviz::instance().frame->GetCanvas()->animate();
}
void
FGLCanvas::
animate() {
	T_(Trace trc(2,"FGLCanvas::animate");)
	static string last_colormode;
	static int omegat = 0;
	int phase = 0;
	Amviz& amviz = Amviz::instance();

#ifdef NEVER // OnPaint not necessary?
	// pause until InitGL() finishes in OnPaint
	T_(trc.dprint("m_init = ",m_init);)
	while (!m_init)
		std::this_thread::sleep_for(std::chrono::seconds{1});
		// snooze(1.0e+8);
#endif // NEVER // OnPaint not necessary?

	if (!glcontext->IsOK())
		errormsg(vastr("could not create an OpenGL context: ",OGLVersion()),true);
// #ifdef NEVER // where should this go
	snooze(5.0e+8);
	Show(true);
	glcontext->SetCurrent(*this);
	InitGL();
// #endif // NEVER // where should this go

	while(!amviz.quit) {
		omegat += amviz.frame->omegat_incr;
		if (omegat >= 360)
			omegat = 0;

		// pause if we are stopped
		while(amviz.stop)
			std::this_thread::sleep_for(std::chrono::milliseconds{10});

		// check to see if an eigenvector has appeared
		if (!amviz.ev.empty()) {
			static int visit{0};
			T_(trc.dprint("new eigenvector, visit ",visit);)
			amviz.add(amviz.ev, to_string(++visit));
			amviz.ev.clear();
		}

		amviz.render(omegat, phase);
		glFlush();
		wxGLCanvas::SwapBuffers();
	}
	T_(trc.dprint("returning: got quit");)
	amviz.quit = false;	// for the next time
}
void
myTimer::
Notify() {
	T_(Trace trc(2,"myTimer::Notify");)
	Amviz& amviz{Amviz::instance()};
	static bool visit{false};
	int interval = GetInterval();
	if (!visit) {
		FGLCanvas* cp = amviz.frame->GetCanvas();
		cp->glcontext->SetCurrent(*cp);
		cp->InitGL();
		visit = true;
		interval = Default_period;
	}
	Stop();
	// render one frame...
	amviz.frame->GetCanvas()->oneframe();
	// ... then re-start the timer
	Start(interval);
}


void
FGLCanvas::
oneframe() {
	T_(Trace trc(2,"FGLCanvas::oneframe");)
	int incr{4};
	Amviz& amviz{Amviz::instance()};
	int phase{0};

	static int omegat{0};
	omegat += incr;
	if (omegat >= 360)
		omegat = 0;

	// check to see if an eigenvector has appeared
	if (!amviz.ev.empty()) {
		static int visit{0};
		T_(trc.dprint("new eigenvector, visit ",visit);)
		amviz.add(amviz.ev, to_string(++visit));
		amviz.ev.clear();
	}

	amviz.render(omegat, phase);
	glFlush();
	wxGLCanvas::SwapBuffers();
}


FGLCanvas::FGLCanvas(wxWindow *parent, wxWindowID id, const int* attlist,
		const wxPoint& pos, const wxSize& size, long style, const wxString& name) :
			wxGLCanvas(parent, id, attlist, pos, size, style, name, wxNullPalette) {
// FGlCanvas constructor for WSL using the old GL attributes and GL context
// XXX combine the 2 ctors somehow
	T_(Trace trc(2,"WSL FGLCanvas ctor");)

   m_gldata.initialized = false;
	m_init = false;

   m_gldata.beginx = 0.0f;
   m_gldata.beginy = 0.0f;
   m_gldata.zoom   = 1.0f;
   // Set the default quat, initialize view matrix with it
	Default_quat = airplane_quat;
	m_gldata.quat = Default_quat;
	float qn = 1.0/blas::nrm2(4, &m_gldata.quat[0], 1);
	if (abs(qn-1.0) > 8.0*std::numeric_limits<float>::epsilon()) {
		blas::scal(4, qn, &m_gldata.quat[0], 1);
		cerr << "inital quat scaled by " << qn << ": " << m_gldata.quat << endl;
	}
	// create a glContext but don't SetCurrent yet
	if (isWSL()) {
		glcontext = new wxGLContext(this,nullptr);
	} else {
		vector<wxGLContextAttrs> catts;
		wxGLContextAttrs c0;
		c0.CoreProfile().PlatformDefaults().EndList();
		catts.push_back(c0);
		wxGLContextAttrs c1;
		c1.CoreProfile().EndList();
		catts.push_back(c1);
		wxGLContextAttrs c2;
		c2.CoreProfile().MajorVersion(1).EndList();
		catts.push_back(c2);
		wxGLContextAttrs c3;
		c3.CompatibilityProfile().EndList();
		catts.push_back(c3);
		// wxGLContextAttrs c4;
		// c4.CoreProfile().OGLVersion(4,6).EndList();
		// catts.push_back(c4);
		// for (auto& ci : catts)
		size_t i{0};
		for (i=0; i<catts.size(); i++) {
			glcontext = new wxGLContext(this,nullptr,&catts[i]);
			if (glcontext->IsOK())
				break;
		}
		T_(trc.dprint("using glcontext attr ",i);)
	}
	if (!glcontext->IsOK())
		throw runtime_error("could not create an OpenGL context");

	// don't show the canvas until animate()
	Show(false);
}	// FGLCanvas WSL ctor

FGLCanvas::FGLCanvas(wxWindow *parent, wxGLAttributes& att, wxWindowID id,
		const wxPoint& pos, const wxSize& size, long style, const wxString& name) :
			wxGLCanvas(parent, att, id, pos, size, style, name, wxNullPalette) {
// FGlCanvas constructor using the newer GL attributes and GL context
	T_(Trace trc(2,"FGLCanvas ctor");)

   m_gldata.initialized = false;
	m_init = false;

   m_gldata.beginx = 0.0f;
   m_gldata.beginy = 0.0f;
   m_gldata.zoom   = 1.0f;
   // Set the default quat, initialize view matrix with it
	Default_quat = airplane_quat;
	m_gldata.quat = Default_quat;
	float qn = 1.0/blas::nrm2(4, &m_gldata.quat[0], 1);
	if (abs(qn-1.0) > 8.0*std::numeric_limits<float>::epsilon()) {
		blas::scal(4, qn, &m_gldata.quat[0], 1);
		cerr << "inital quat scaled by " << qn << ": " << m_gldata.quat << endl;
	}
	// create a glContext but don't SetCurrent yet
	if (isWSL()) {
		glcontext = new wxGLContext(this,nullptr);
	} else {
		vector<wxGLContextAttrs> catts;
		wxGLContextAttrs c0;
		c0.CoreProfile().PlatformDefaults().EndList();
		catts.push_back(c0);
		wxGLContextAttrs c1;
		c1.CoreProfile().EndList();
		catts.push_back(c1);
		wxGLContextAttrs c2;
		c2.CoreProfile().MajorVersion(1).EndList();
		catts.push_back(c2);
		wxGLContextAttrs c3;
		c3.CompatibilityProfile().EndList();
		catts.push_back(c3);
		// wxGLContextAttrs c4;
		// c4.CoreProfile().OGLVersion(4,6).EndList();
		// catts.push_back(c4);
		// for (auto& ci : catts)
		size_t i{0};
		for (i=0; i<catts.size(); i++) {
			glcontext = new wxGLContext(this,nullptr,&catts[i]);
			if (glcontext->IsOK())
				break;
		}
		T_(trc.dprint("using glcontext attr ",i);)
	}
	if (!glcontext->IsOK())
		throw runtime_error("could not create an OpenGL context");

	// don't show the canvas until animate()
	Show(false);
#ifdef NEVER // where should this go
	Show(true);
	glcontext->SetCurrent(*this);
	InitGL();
#endif // NEVER // where should this go
}	// FGLCanvas ctor

FGLCanvas::~FGLCanvas() { }

void
FGLCanvas::
OnPaint( wxPaintEvent& event) {
// This doesn't do any painting, just adjusts the view to
// account for trackball motion  XXX maybe there's a better place for it
	T_(Trace trc(2,"FGLCanvas::OnPaint");)

	glcontext->SetCurrent(*this);
#ifdef NEVER // where should this go
	wxPaintDC dc(this);
    // Init OpenGL once, but after SetCurrent
    if (!m_init) {
        InitGL();
        m_init = true;
    }
#endif // NEVER // where should this go

    // Transformations
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glTranslatef( OriginX, OriginY, -2.0f );

	 // viewing area (for zooming):
	glOrtho(ViewLeft, ViewRight, ViewBottom, ViewTop, ViewNear, ViewFar);

    GLfloat m[4][4];
	 // trackball:
    build_rotmatrix( m, &m_gldata.quat[0] );
    glMultMatrixf( &m[0][0] );

} // FGLCanvas::OnPaint()

void
FGLCanvas::
OnSize(wxSizeEvent& event) {
// this is also necessary to update the context on some platforms
// replacing wxGLCanvas::OnSize() with event.Skip()
// per:
// http://groups.google.com/group/wx-users/browse_thread/thread/cfa90147ee800f0c?pli=1
    //wxGLCanvas::OnSize(event);
	 event.Skip();
    // Reset the OpenGL view aspect
    ResetProjectionMode();
}

void
FGLCanvas::
OnEraseBackground(wxEraseEvent& WXUNUSED(event)) {
    // Do nothing, to avoid flashing on MSW
}

void
screen2world(const wxPoint& pt,
		GLdouble& worldx,
		GLdouble& worldy,
		GLdouble& worldz) {
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];
	GLdouble winx = pt.x;
	GLdouble winy = pt.y;
	GLdouble winz{0.0};

	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);
	gluUnProject(winx, winy, winz, modelview, projection, viewport,
			&worldx, &worldy, &worldz);
}

void
world2screen(GLdouble worldx, GLdouble worldy,
		GLdouble worldz, GLdouble& winx, GLdouble& winy, GLdouble& winz) {
// convert from world to screen coord
	GLint viewport[4];
	GLdouble modelview[16];
	GLdouble projection[16];

	glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
	glGetDoublev(GL_PROJECTION_MATRIX, projection);
	glGetIntegerv(GL_VIEWPORT, viewport);
	// convert from world to screen coord
	int stat = gluProject(worldx, worldy, worldz,
			modelview, projection, viewport, &winx, &winy, &winz);
	if (stat == GL_FALSE) {
		cerr << "gluProject failed\n";
	}
}

void
FGLCanvas::
display_node(double x, double y) {
	T_(Trace trc(1,"display_node");)
	Amviz& amviz{Amviz::instance()};

	if (amviz.frame->isRunning()) {
	// if (DispScale != 0.0) 
		wxMessageBox("node id only with motion stopped (press stop button)");
		return;
	}

	wxSize siz = GetClientSize();
	T_(trc.dprint("window width ", siz.GetWidth(),", height ",siz.GetHeight());)
	double h = siz.GetHeight();
	y = h - y;  // transform from wx to gl coord
	GLdouble winx;
	GLdouble winy;
	GLdouble winz;
	int node{0};
	double mindist{std::numeric_limits<double>::max()};
	vector<double>& coord = amviz.coordinates();
	// find the closest node
	for (size_t i=0; i<coord.size(); i += 3) {
		GLdouble worldx{coord[i]};
		GLdouble worldy{coord[i+1]};
		GLdouble worldz{coord[i+2]};
		// transform this node's coords to screen
		world2screen(worldx, worldy, worldz, winx, winy, winz);
		T_(trc.dprint("node ",i/3," world (",worldx,',',worldy, ',',worldz,") screen (", (int)winx, ',', (int)winy,',',(int)winz,")\n");)
		// compute distance to mouse click
		double dist = sqrt((winx-x)*(winx-x) + (winy-y)*(winy-y));
		if (dist < mindist) {
			node = i/3 + 1;
			mindist = dist;
		}
	}
	T_(trc.dprint("got alt-left at (", x, ", ", y,") placing ",node,"\n");)
	string nodestr = vastr(node);
	wxMessageBox(nodestr);
}

void
FGLCanvas::
OnMouse(wxMouseEvent& event) {
	double x = event.GetX();
	double y = event.GetY();
	// alt-left: display closest node number
	if (event.LeftIsDown() && event.AltDown()) {
		display_node(x,y);
		return;
	}

	if (event.Dragging()) {
		float x = event.GetX();
		float y = event.GetY();
		// dragging the mouse turns off Quad, X, Y, and Z view modes:
		ViewX = ViewY = ViewZ = false;
		if (ViewQuad) {
			ViewQuad = false;
			viewX(Omegat, Phase);
		}

		wxSize sz(GetClientSize());
		if (event.LeftIsDown()) {
			// drag in progress, simulate trackball to rotate the figure
			vector<float> spin_quat(4,0);
			trackball(&spin_quat[0],
            (2.0*m_gldata.beginx - sz.x) / sz.x,
            (sz.y - 2.0*m_gldata.beginy) / sz.y,
            (2.0*x - sz.x)    / sz.x,
            (sz.y - 2.0*y)    / sz.y);
			// add spin_quat to the existing quat, overwriting it
			add_quats(&spin_quat[0], &m_gldata.quat[0], &m_gldata.quat[0]);
			Refresh(false);
		} else if (event.RightIsDown()) {
			// scale the figure
			float xc = x - sz.x/2.0;
			float yc = y - sz.y/2.0;
			float radius = sqrt(xc*xc + yc*yc);
			xc = m_gldata.beginx - sz.x/2.0;
			yc = m_gldata.beginy - sz.y/2.0;
			float lastr = sqrt(xc*xc + yc*yc);
			float scale = 1.0 + (radius - lastr)/std::max(radius, lastr);
			ViewLeft /= scale;
			ViewRight /= scale;
			ViewBottom /= scale;
			ViewTop /= scale;
			Refresh(false);
		} else if (event.MiddleIsDown()) {
			// translation
				// XXX change to beginx, remove lastX from amviz.h:
			OriginX += (event.GetX() - m_gldata.beginx)/1200.0;
			OriginY -= (event.GetY() - m_gldata.beginy)/1200.0;
			Refresh(false);
		}
	}
		m_gldata.beginx = event.GetX();
		m_gldata.beginy = event.GetY();
}

void
FGLCanvas::
InitGL() {
	T_(Trace trc(2,"InitGL");)

    // set viewing projection XXX necessary?
	glMatrixMode(GL_PROJECTION);
	// size of the dots:
	glPointSize(2.0);
	// the nodes have been scaled so the largest coord is 1
	// and centered, so make the viewing area a little larger
	// glOrtho(left, right, bottom, top, near, far)
	glOrtho(ViewLeft, ViewRight, ViewBottom, ViewTop, ViewNear, ViewFar);

#ifdef NEVER // necessary?
    /* position viewer */
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -2.0f);
#endif // NEVER // necessary?

#ifdef NEVER // overrides quat?
    /* position object */
	 glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
	 glRotatef(135.0f, 0.0f, 0.0f, 1.0f);
#endif // NEVER // overrides quat?
}

void
FGLCanvas::
ResetProjectionMode() {
    int w, h;
    GetClientSize(&w, &h);
}

void
FGLCanvas::
print() {

	errormsg("printing has not been implemented yet",false);
#ifdef NEVER // need local/wxpdfdoc library
	wxPdfDocument pdfDoc;
	if (!pdfDoc.Create(wxT("output.pdf"))) {
		wxLogError("Failed to create PDF document.");
		return;
	}

	wxGraphicsContext* gc = wxGraphicsContext::CreateFromNativeWindow(this);
	if (!gc) {
		wxLogError("Failed to create graphics context from wxGLCanvas.");
		return;
	}

	wxPdfPage* page = pdfDoc.AddPage();

	gc->SetPDFContext(page->GetContext());
	// Trigger a redraw of the wxGLCanvas to ensure the latest content is captured
	myGLCanvas->Refresh(false);
	myGLCanvas->Update();
	gc->SetPDFContext(nullptr);

	if (!pdfDoc.Save()) {
		wxLogError("Failed to save PDF document.");
	}
	delete gc; // Clean up the graphics context
#endif // NEVER // need local/wxpdfdoc library

#ifdef NEVER // use wxPDF

// #define XFIG 1
// #ifdef XFIG
// 	FGLCanvas::xfig();
// #else
	GLint viewport[4];
	GLint options = GL2PS_USE_CURRENT_VIEWPORT;
	GLint format = GL2PS_PS;
	GLint sort = GL2PS_NO_SORT;
	GLint colormode = GL_RGB;
	GLint colorsize = 0;   // using GL_RGB instead of color map
	GL2PSrgb* colortable = 0;   // using GL_RGB instead of color map
	GLint nr = 0;
	GLint ng = 0;
	GLint nb = 0;
	GLint buffersize = 1200*900*8;
	char const* filename = "amviz.ps";
	FILE* stream = fopen(filename, "w");
	gl2psBeginPage("amviz", "amviz", viewport, format,
			sort, options, colormode, colorsize, colortable,
			nr, ng, nb, buffersize, stream, filename);
	Render(Omegat, Phase);
	gl2psEndPage();

	ostringstream os;
	os << "Created postscript file " << filename;
	wxChar wst[256];
	for (size_t i=0; i<os.str().size(); i++)
		wst[i] = os.str()[i];
	wxString txt(wst, os.str().size());
	wxMessageBox(txt, _("printPS"), wxOK);
// #endif // XFIG
#endif // NEVER // use wxPDF
}

void
FGLCanvas::
print90270() {
}

struct point {
	int x;
	int y;
	point() : x(0), y(0) {}
	point(GLfloat* buf) {
		x = (int)buf[0];
		y = (int)buf[1];
	}
	point(int a, int b) : x(a), y(b) {}

	bool operator==(const point& b) const {
		return (x == b.x && y == b.y);
	}
};

static int
xfig_ymax(vector<GLfloat>& buffer) {
	int rval{std::numeric_limits<int>::lowest()};
	for( size_t i=0; i<buffer.size(); i++ ) {
		if( buffer[i] == GL_POINT_TOKEN ) {
			i += 7;
		} else if ( buffer[i] == GL_LINE_TOKEN ||
				buffer[i] == GL_LINE_RESET_TOKEN) {
			point p1 = point(&buffer[i+1]);
			point p2 = point(&buffer[i+3]);
			rval = std::max(rval, p1.y);
			rval = std::max(rval, p2.y);
			i += 4;
		}
	}
	return rval;
}

static void
xfig_write(ofstream& stream, vector<point>& polyline, int ymax,
		int style, int width, int color) {
	if (!polyline.empty()) {
		stream << "2 1 " << style << " " << width << " " << color
			<< " 7 50 -1 -1 0.000 0 0 -1 0 0 " << polyline.size() << endl;
		// transform to xfig coords: y <- ymax - y
		for (size_t j=0; j<polyline.size(); j++) {
			int x = polyline[j].x;
			int y = ymax - polyline[j].y;
			stream << x << " " << y << endl;
		}
		polyline.clear();
	}
}

static void
xfig_create(const string& path, vector<GLfloat>& buffer, bool append,
		int style, int width, int color) {
	vector<point> polyline;
	static int ymax{0};
	static int frameno{0};

	std::ios_base::openmode mode = std::ios::trunc;
	if (append) {
		mode = std::ios::app;
		frameno++;
	}
	ofstream stream(path, mode);

	if (!append) {
		// determine the largest y value for transforming coord systems
		// in xfig_write (only if not appending)
		ymax = xfig_ymax(buffer);

		stream << "#FIG 3.2  Produced by amviz\n";
		stream << "Landscape\n";
		stream << "Center\n";
		stream << "Inches\n";
		stream << "Letter\n";
		stream << "100.00\n";
		stream << "Single\n";
		stream << "-2\n";
		stream << "100 2\n";
		frameno = 0;
	}

	stream << "# Frame " << frameno << endl;

	point p1;
	point p2;
	int big{numeric_limits<int>::max()};
	point lastp2(big,big);
	for( size_t i=0; i<buffer.size(); i++ ) {
		if( buffer[i] == GL_POINT_TOKEN ) {
			// points.push_back( FlashPoint( &buffer[i+1]) );
			i += 7;

		} else if ( buffer[i] == GL_LINE_TOKEN ) {
			p1 = point(&buffer[i+1]);
			// p2 = point(&buffer[i+8]);
			p2 = point(&buffer[i+3]);
			if (p1 == lastp2) {
				polyline.push_back(p2);
			} else {
				xfig_write(stream, polyline, ymax, style, width, color);
				polyline.push_back(p1);
				polyline.push_back(p2);
			}
			lastp2 = p2;
			// i += 14;
			i += 4;

		} else if ( buffer[i] == GL_LINE_RESET_TOKEN ) {
			// penlift?
			// lines.push_back( FlashLine( &buffer[i+1]) );
			polyline.push_back(point( &buffer[i+1]) );
			// polyline.push_back(point( &buffer[i+8]) );
			polyline.push_back(point( &buffer[i+3]) );
			xfig_write(stream, polyline, ymax, style, width, color);
			// i += 14;
			i += 4;

		} else if (buffer[i] == GL_PASS_THROUGH_TOKEN) {
			break;
		}
	}
	xfig_write(stream, polyline, ymax, style, width, color);
}

void
FGLCanvas::
xfig() {
	Amviz& amviz{Amviz::instance()};
	int omegat = 0;
	int phase = 0;
	// render into a gl feedback buffer
	int buffersize = 1000000;
	wxString wxpath = wxGetTextFromUser(wxT("Please enter a filename for xfig output "
							  "(include a *.fig extenstion)"), 
							  wxT("xfig filename"), wxT("amviz.fig"));
	string pathname = wxpath.ToStdString();
	wxString wxsteps = wxGetTextFromUser(wxT("How many steps between 0-360 deg?"),
										  wxT("# steps"), wxT("5"));
	string nstepstr = wxsteps.ToStdString();
	int nstep{5};
	bool ok = str2int(nstepstr, nstep);
	if (!ok)
		cerr << "not an int: " << nstepstr << endl;

	vector<GLfloat> feedbuf(buffersize);

	bool append{false};
	int delta{0};
	int start{0};
	int style{0};
	int color{0};
	// nstep points between -90:90
	if (nstep > 1) {
		delta = 180/(nstep-1);
		start = -90;
	}
	for (int i=0; i<nstep; i++) {
		int width{1};
		if (i == 0 || i == nstep-1) {
			width = 2;
			color = 0;
		} else {
			width = 1;
			// color = i;
			color = 0;
		}
		omegat = start + i*delta;
		if (omegat == 0)
			color = 4;		// red for neutral
		glFeedbackBuffer(feedbuf.size(), GL_2D, &feedbuf[0]);
	
		(void)glRenderMode(GL_FEEDBACK);
		// Render(omegat, phase);
		amviz.render(omegat, phase);
		glFlush();
		wxGLCanvas::SwapBuffers();
		glRenderMode(GL_RENDER);
		xfig_create(pathname, feedbuf, append, style, width, color);
		append = true;
	}
}

void
FGLCanvas::
video() {
// create an mp4 video file of a few cycles (ncycles) of animation
	T_(Trace trc(2,"FGLCanvas::video");)
	Amviz& amviz{Amviz::instance()};
	Specs& sp = specs();

	if (sp.threads)
		amviz.stop = true;
	else
		amviz.frame->timer->Stop();

	// pause to let it stop
	std::this_thread::sleep_for(std::chrono::milliseconds{10});

	int phase = 0;
	wxString pathname = wxGetTextFromUser("Please enter a filename for video output "
		"(include an *.mp4 extenstion)", "mp4 filename", "amviz.mp4");

	string path = pathname.ToStdString();
	string basefile = stringBasename(path);
	// the extension determines the output format: mp4 or png
	string ext = rsubstr(path, 4);

	// get window size
	wxSize siz = GetClientSize();
	GLsizei width = siz.GetWidth();
	GLsizei height = siz.GetHeight();
	T_(trc.dprint("window width ", width,", height ",height,", canvas: ",amsizes.canvas);)

	// start ffmpeg
	FILE* ffmpeg{nullptr};
	if (ext == ".mp4") {
		int framerate{60};
		string resolution{vastr(width,"x",height)};
		//!! string file{"amviz.mp4"};
		string cmd = vastr("ffmpeg -r ",framerate," -f rawvideo -pix_fmt yuv420p -s ",
			resolution, " -i - -threads 0 -preset fast -y -crf 21 -vf vflip ",path);
		cout << cmd << endl;
		// open a write pipe to ffmpeg using RAII class Fpipe
		// FILE* ffmpeg = popen(cmd.c_str(), "w");
		try {
			Fpipe ffmpegpipe(cmd, "w");
			ffmpeg = ffmpegpipe();
		} catch (std::exception& s) {
			errormsg(vastr("cannot start ffmpeg: ",s.what()), false);
		}
		// if (ffmpeg == nullptr)
		// 	throw std::system_error(errno, std::generic_category(),
		// 		"cannot start ffmpeg");
	}

	glRenderMode(GL_RENDER);

	size_t bufsize = 3*width*height*sizeof(unsigned char);
	//!! vector<unsigned char> buffer(width*height*3);
	vector<float> buffer(width*height*3);
	int ncycles{4};
	if (ext == ".png")
		ncycles = 1;
	int frameno{1};
	// 360 degrees for ncycles cycles and ncycles*40 frames total.
	for (int j=0; j<ncycles; j++) {
		for (size_t i=0; i < 360; i+=9) {
			amviz.render(i, phase);
			glFlush();
			wxGLCanvas::SwapBuffers();
			//!! glReadPixels(0,0,width,height,GL_RGB, GL_UNSIGNED_BYTE, buffer.data());
			glReadPixels(0,0,width,height,GL_RGB, GL_FLOAT, buffer.data());
			if (ext == ".mp4")
				fwrite(&buffer[0], bufsize, 1, ffmpeg);
			else {
				// append the base filename with, e.g. 0001, 0002, etc
#ifdef NEVER // frameno
				ostringstream os;
				os << basefile << setw(4) << setfill('0') << i;
				os << ".png";
				pngwrite(os.str(), buffer.data(), width, height);
#else // NEVER // frameno
				pngwrite(vastr(basefile,frameno++,".png"),buffer.data(),width,height);
#endif // NEVER // frameno
			}
			std::this_thread::sleep_for(std::chrono::milliseconds{100});
		}
	}

	// close file XXX unecessary: RAII
	//!! pclose(ffmpeg);

	// finished filling feedbuf, return to GL_RENDER mode...
	glRenderMode(GL_RENDER);
	// ... and restart the animation
	if (sp.threads)
		amviz.stop = false;
	else
		amviz.frame->timer->Start(Default_period);
}

void
FGLCanvas::
viewX(int omegat, int phase) {
	int width, height;
	GetClientSize(&width, &height);
	ViewQuad = false;
}

void
FGLCanvas::
viewY(int omegat, int phase) {
	int width, height;
	GetClientSize(&width, &height);
	ViewQuad = false;
}

void
FGLCanvas::
viewZ(int omegat, int phase) {
	int width, height;
	GetClientSize(&width, &height);
	ViewQuad = false;
}


void
FGLCanvas::
viewQuad(int omegat, int phase) {
//    int width = glutGet(GLUT_WINDOW_WIDTH);
//    int height = glutGet(GLUT_WINDOW_HEIGHT);
	// GLint params[4];
	// glGetIntegerv(GL_VIEWPORT, params);
	// int width = params[2];
	// int height = params[3];
	int width, height;
	GetClientSize(&width, &height);
    width = (width+1)/2;
    height = (height+1)/2;
	glEnable(GL_SCISSOR_TEST);

#ifdef NEVER // needs work
	bottom_left;
	bottom;
	Render(omegat, phase);

	bottom_right;
	left;
	Render(omegat, phase);

	top_left;
	front;
	Render(omegat, phase);

	top_right;
	perspective;
	Render(omegat, phase);
#endif // NEVER // needs work
}

