//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
/////////////////////////////////////////////////////////////////////////////
// Flaps 3D plots
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx.h".

#include <cerrno>
#include <ctime>
#include <exception>
#include <stdexcept>
#include <iosfwd>

#include "specs.h"

#include "wx/string.h"
#include "wx/tglbtn.h"
#include "wx/toolbar.h"
#include "wx/prntbase.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#if !defined(__WXGTK__)
	#error "rebuild with WXGTK"
#endif

#if !wxUSE_GLCANVAS
    #error "OpenGL required: set wxUSE_GLCANVAS to 1 and rebuild the library"
#endif
#include <GL/glu.h>
// #include <GL/glut.h>

#undef Complex    // X11 macro
#include "config.h"
#include "vz3d.h"
#include "exim.h"
#include "gl2ps.h"
#include "matrix.h"
#include "prefdialog.h"
#include "swf.h"
#include "trace.h"

using namespace std;

// for a discussion of Quaternions as rotations, see
//   en.wikipedia.org/wiki/Quaternions_and_spatial_rotation
class Quaternion {
public:
	string id;
	vector<float> quat;
	Quaternion(const string& ident, vector<float> q) : id(ident), quat(q) {}
};

vector<float> airplane_quat{0.396197, -0.133608, -0.491376, 0.76402};
vector<float> Isoquat{0.180731, 0.356605, 0.827727, 0.393747};
vector<float> Xquat{-0.5, -0.5, -0.5, -0.5};
vector<float> Yquat{-0.707106781, 0.0, 0.0, -0.707106781};
vector<float> Zquat{0.0, 0.0, 0.0, 1.0};
//!! vector<float> Default_quat{0.274,0.52399,0.707106781,0.364944};
vector<float> Default_quat{Isoquat};

// overall window size in pixels
int WindowSizeX = 1100;
int WindowSizeY = 900;

float OriginX = 0.0;
float OriginY = 0.0;

float OrigViewLeft = -0.7;
float OrigViewRight = 0.7;
float OrigViewBottom = -0.7;
float OrigViewTop = 0.7;
float OrigViewNear = -3.0;
float OrigViewFar = 3.0;

float ViewLeft = OrigViewLeft;
float ViewRight = OrigViewRight;
float ViewBottom = OrigViewBottom;
float ViewTop = OrigViewTop;
float ViewNear = OrigViewNear;
float ViewFar = OrigViewFar;

bool ViewDefault{true};
bool ViewQuad = false;
bool ViewX = false;
bool ViewY = false;
bool ViewZ = false;

string Printer;

static void projection (int width, int height, int perspective);

// for front, ...
float spin_x = 0.0;
float spin_y = 0.0;
#define bottom_left                   \
    glViewport(0, 0, width, height);  \
    glScissor(0, 0, width, height)

#define bottom_right                     \
    glViewport(width, 0, width, height); \
    glScissor(width, 0, width, height)

#define top_left                          \
    glViewport(0, height, width, height); \
    glScissor(0, height, width, height)

#define top_right                             \
    glViewport(width, height, width, height); \
    glScissor(width, height, width, height)

#define front                         \
    projection(width, height, 0);     \
    glRotatef(spin_y, 1.0, 0.0, 0.0); \
    glRotatef(spin_x, 0.0, 1.0, 0.0)

#define back                          \
    projection(width, height, 0);     \
    glRotatef(180.0, 0.0, 1.0, 0.0);  \
    glRotatef(spin_y, 1.0, 0.0, 0.0); \
    glRotatef(spin_x, 0.0, 1.0, 0.0)

#define right                            \
    projection(width, height, 0);        \
    glRotatef(90.0, 0.0, 1.0, 0.0);      \
    glRotatef(spin_y, 1.0, 0.0, 0.0);    \
    glRotatef(spin_x, 0.0, 1.0, 0.0)

#define left                             \
    projection(width, height, 0);        \
    glRotatef(-90.0, 0.0, 1.0, 0.0);      \
    glRotatef(spin_y, 1.0, 0.0, 0.0);    \
    glRotatef(spin_x, 0.0, 1.0, 0.0)

#define top                               \
    projection(width, height, 0);         \
    glRotatef(90.0, 1.0, 0.0, 0.0);       \
    glRotatef(spin_y, 1.0, 0.0, 0.0);     \
    glRotatef(spin_x, 0.0, 1.0, 0.0)

#define bottom                            \
    projection(width, height, 0);         \
    glRotatef(-90.0, 1.0, 0.0, 0.0);      \
    glRotatef(spin_y, 1.0, 0.0, 0.0);     \
    glRotatef(spin_x, 0.0, 1.0, 0.0)

#define perspective                           \
    projection(width, height, 1);             \
    glRotatef(40.0, 0.0, 1.0, 0.0);           \
    glRotatef(-40.0, 1.0, 0.0, 0.0);           \
    glRotatef(20.0, 0.0, 0.0, 1.0);           \
    glRotatef(-10.0, 0.0, 1.0, 0.0);           \
    glRotatef(spin_y, 1.0, 0.0, 0.0);         \
    glRotatef(spin_x, 0.0, 1.0, 0.0)


void
color_values (double t, float& red, float& green, float& blue) {
// given a number (t) between 0 and 1, set the rgb values
// so that they vary linearly with t=0 blue, t=0.5 green,
// and t=1 red

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

Vz3d* Vz3d::instance_{nullptr};

Vz3d* 
Vz3d::
instance (const vector<Plotcurve*>& curves, Specs& specs) {
	instance_ = new Vz3d(curves, specs);
	return instance_;
}

Vz3d::
Vz3d (const vector<Plotcurve*>& cv, Specs& specs) : curves(cv) {
	// copy x, y, z values from curves, converting to PU
	T_(Trace trc(1,"Vz3d constructor");)

	// x, y, z parameter names may not have been set in options
	xname_ = specs.xname;
	yname_ = specs.ynames[0];
	zname_ = specs.zname;
	// get a list of all the parameter names from the first curve (parnames)
	// XXX what if the curves have diff parameters?
	vector<Par*> theparams = cv[0]->allpar();
	for (auto& p : theparams)
		parnames_.push_back(p->name);
	
#ifdef NEVER // use load_data
	double big{std::numeric_limits<double>::max()};
	double xmax{-big};
	double ymax{-big};
	double zmax{-big};
	// XXX require all curves to have the same x,y,z parameters?
	for (auto& cj : curves) {
		vector<double> ci;
		Par* xp = cj->xparam;
		xname_ = cj->xparam->name;
		Par* yp = cj->yparam;
		yname_ = cj->yparam->name;
		Par* zp = cj->zparam;
		zname_ = cj->zparam->name;
		T_(trc.dprint("working on ",cj->cid(),", nsolns ",cj->nsolns());)
		for (size_t i=0; i<xp->nsolns(); i++) {
			double x = xp->solns[i]*xp->conv.factor;
			double y = yp->solns[i]*yp->conv.factor;
			double z = zp->solns[i]*zp->conv.factor;
			//!! xmin = std::min(xmin, x);
			xmax = std::max(xmax, std::abs(x));
			//!! ymin = std::min(ymin, y);
			ymax = std::max(ymax, std::abs(y));
			//!! zmin = std::min(zmin, z);
			zmax = std::max(zmax, std::abs(z));
			ci.push_back(x);
			ci.push_back(y);
			ci.push_back(z);
		}
		T_(trc.dprintm(3,ci.size()/3,3,ci.data(), cj->cid());)
		data.push_back(ci);
	}
	// scale each curve to fit in a t-cube
	double t{0.3};
	//!! double xsf = abs(xmax-xmin)/2.0;
	//!! double ysf = abs(ymax-ymin)/2.0;
	//!! double zsf = abs(zmax-zmin)/2.0;
	T_(trc.dprint("x scale: ",xmax,", y ",ymax,", z ",zmax);)
	for (size_t i=0; i<data.size(); i++) {
		int n = data[i].size()/3;
		blas_scal(n, t/xmax, &data[i][0], 3);
		blas_scal(n, t/ymax, &data[i][1], 3);
		blas_scal(n, t/zmax, &data[i][2], 3);
		T_(trc.dprintm(3,n,3,data[i].data(), vastr("scaled",i));)
	}
#else // NEVER // use load_data
	// if x, y, z parameter names were spec'd load data
	if (!(xname_.empty() || yname_.empty() || zname_.empty()))
		load_data();
#endif // NEVER // use load_data
	T_(trc.dprint("got ",data.size()," curves");)
}

void
Vz3d::
load_data() {
// load the data array from current x, y, & z parameters
	T_(Trace trc(2, "load_data");)
	double big{std::numeric_limits<double>::max()};
	double xmax{-big};
	double ymax{-big};
	double zmax{-big};

	T_(trc.dprint("x: ",xname_,", y: ",yname_,", z: ",zname_);)

	// clear the current data array
	data.clear();

	// XXX require all curves to have the same x,y,z parameters?
	for (auto& cj : curves) {
		vector<double> ci;
		Par* xp = cj->findp(xname_);
		if (xp == nullptr)
			throw runtime_error(vastr(cj->cid()," does not have parameter \"",xname_,"\""));
		Par* yp = cj->findp(yname_);
		if (yp == nullptr)
			throw runtime_error(vastr(cj->cid()," does not have parameter \"",yname_,"\""));
		Par* zp = cj->findp(zname_);
		if (zp == nullptr)
			throw runtime_error(vastr(cj->cid()," does not have parameter \"",zname_,"\""));
		int n = cj->nsolns();
		T_(trc.dprint("working on ",cj->cid(),", nsolns = ",n);)
		double xcf = xp->conv.factor;
		double ycf = yp->conv.factor;
		double zcf = zp->conv.factor;
		T_(trc.dprint("xcf: ",xcf,", ycf: ",ycf,", zcf: ",zcf);)
		for (int i=0; i<n; i++) {
			double x = xp->solns[i]*xcf;
			double y = yp->solns[i]*ycf;
			double z = zp->solns[i]*zcf;
			xmax = std::max(xmax, std::abs(x));
			ymax = std::max(ymax, std::abs(y));
			zmax = std::max(zmax, std::abs(z));
			ci.push_back(x);
			ci.push_back(y);
			ci.push_back(z);
		}
		T_(trc.dprintm(3,ci.size()/3,3,ci.data(), cj->cid());)
		data.push_back(ci);
	}
	// scale each curve to fit in a t-cube
	double t{0.3};
	//!! double xsf = abs(xmax-xmin)/2.0;
	//!! double ysf = abs(ymax-ymin)/2.0;
	//!! double zsf = abs(zmax-zmin)/2.0;
	T_(trc.dprint("x scale: ",xmax,", y ",ymax,", z ",zmax);)
	for (size_t i=0; i<data.size(); i++) {
		int n = data[i].size()/3;
		blas_scal(n, t/xmax, &data[i][0], 3);
		blas_scal(n, t/ymax, &data[i][1], 3);
		blas_scal(n, t/zmax, &data[i][2], 3);
		T_(trc.dprintm(3,n,3,data[i].data(), vastr("scaled",i));)
	}
	T_(trc.dprint("got ",data.size()," curves");)
}

void
Vz3d::
render() {
	T_(Trace trc(1,"render");)

	// if colormode == "none" use the reverse of the bg color
	float red = 1.0 - BGRed;
	float green = 1.0 - BGGreen;
	float blue = 1.0 - BGBlue;

	if (data.empty())
		throw runtime_error("no data to plot");

    /* clear color and depth buffers */
	glClearColor(BGRed, BGGreen, BGBlue, BGAlpha);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	// antialiasing
	glEnable(GL_LINE_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	// glBlendFunc(1.0, 0.0);
	glHint(GL_LINE_SMOOTH_HINT, GL_DONT_CARE);


	// colorvalue is a vector of color values ranging from 0-1 to
	// be used by color_values to give rgb values for GL
	int ncurve = data.size();
	vector<double> colorvalue(ncurve, 0.0);
	double intvl{1.0/ncurve};
	for (int i=0; i<ncurve; i++)
		colorvalue[i] = (i+0.5)*intvl;

	glLineWidth(2.0);
	glBegin(GL_LINES);
	for (int j=0; j<ncurve; j++) {
		if (!toplot.empty() &&
			std::find(toplot.begin(), toplot.end(), j) == toplot.end())
				continue;
		// set the color for this curve
		color_values(colorvalue[j], red, green, blue);
		glColor3f(red, green, blue);
		int nm = data[j].size()/3 - 1;
		for (int i=0; i<nm; i++) {
			double* p = &data[j][3*i];
			float x = *p++;
			float y = *p++;
			float z = *p++;
			//!! glVertex3f(float(*p++), float(*p++), float(*p++));
			glVertex3f(x, y, z);
			//!! glVertex3f(float(*p++), float(*p++), float(*p++));
			x = *p++; y = *p++; z = *p++;
			glVertex3f(x, y, z);
		}
	}
	// draw the coord system if requested
	vector<float> t;
	vector<float> origin{0.0, 0.0, 0.0};
	for (int i=0; i<3; i++) {
		t = origin;
		color_values(1.0-i*0.5, red, green, blue);
		glColor3f(red, green, blue);
		glVertex3f(t[0], t[1], t[2]);
		t[i] += 0.3;
		glVertex3f(t[0], t[1], t[2]);
	}
	glEnd();

	// glFlush();
	// wxGLCanvas::SwapBuffers();
} // Vz3d::render()

void
Vz3d::
run() {
	T_(Trace trc(1,"Vz3d::run");)
	int argc = 1;
	char* argv[3];
	char prog[8];
	strcpy(prog, "vz3d");
	argv[0] = prog;
	argv[1] = nullptr;
	try {
		wxEntry(argc, argv);
	} catch (runtime_error const& s) {
		cerr << s.what() << std::endl;
	} catch (std::exception const& t) {
		cerr << t.what() << std::endl;
	}
}

static void
projection (int width, int height, int persp) {
    float ratio = (float)width/height;
	 // float vr = 0.6*ratio;
	 float vr = ratio;

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
	if (persp) {
		gluPerspective(25, ratio, 1, 256);
	} else {
      glOrtho(-vr, vr, -vr, vr, 1, 256);
	}
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
}


// MyApp

// `Main program' equivalent, creating windows and returning main app frame

MyFrame* Theframe{nullptr};

bool
MyApp::
OnInit() {
	T_(Trace trc(1,"OnInit");)
    // Create the main frame window...
	Theframe = new MyFrame(nullptr, wxT("Flaps 3D plot"),
			 wxDefaultPosition, wxSize(WindowSizeX,WindowSizeY));
	Theframe->Show(true);

	//!! Theframe->glcanvas()->setcurrent();
	//!! Theframe->glcanvas()->Render();

	return true;
}

IMPLEMENT_APP_NO_MAIN(MyApp)

// MyFrame

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(wxID_OPEN, MyFrame::OnMenuFileOpen)
    EVT_MENU(wxID_PRINT, MyFrame::OnMenuFilePrint)
    EVT_MENU(ID_XFIG, MyFrame::OnMenuFileXfig)
    //!! EVT_MENU(wxID_VIDEO, MyFrame::OnMenuFileVideo)
    EVT_MENU(wxID_EXIT, MyFrame::OnMenuFileExit)
    EVT_MENU(ID_EDIT_PREF, MyFrame::OnMenuEditPref)
    //!! EVT_CHOICE(ID_MODE_NUMBER, MyFrame::OnMenuMode)
    EVT_MENU(ID_HELP_ABOUT, MyFrame::OnMenuHelpAbout)
    EVT_MENU(ID_HELP_HELP, MyFrame::OnMenuHelpHelp)
    EVT_MENU(ID_VIEW_DEFAULT, MyFrame::OnMenuViewDefault)
    EVT_MENU(ID_VIEW_X, MyFrame::OnMenuViewX)
    EVT_MENU(ID_VIEW_Y, MyFrame::OnMenuViewY)
    EVT_MENU(ID_VIEW_Z, MyFrame::OnMenuViewZ)
    EVT_MENU(ID_VIEW_I, MyFrame::OnMenuViewI)
    EVT_MENU(ID_VIEW_QUAD, MyFrame::OnMenuViewQuad)
    //!! EVT_MENU(ID_VIEW_COLOR_ABS, MyFrame::OnMenuViewColorAbs)
    //!! EVT_MENU(ID_VIEW_COLOR_VAL, MyFrame::OnMenuViewColorVal)
    //!! EVT_MENU(ID_VIEW_COLOR_PHASE, MyFrame::OnMenuViewColorPhase)
    //!! EVT_MENU(ID_VIEW_COLOR_PHASE_ABS, MyFrame::OnMenuViewColorPhaseAbs)
    //!! EVT_MENU(ID_VIEW_COLOR_ENERGY, MyFrame::OnMenuViewColorEnergy)
    //!! EVT_MENU(ID_VIEW_COLOR_NONE, MyFrame::OnMenuViewColorNone)

	//!! EVT_COMMAND_SCROLL(ID_AMP_SLIDER, MyFrame::OnAmpSlider)
	//!! EVT_COMMAND_SCROLL(ID_PHASE_SLIDER, MyFrame::OnPhaseSlider)
	//!! EVT_COMMAND_SCROLL(ID_FREQ_SLIDER, MyFrame::OnFreqSlider)
	//!! EVT_COMMAND_SCROLL(ID_COLOR_SLIDER, MyFrame::OnColorSlider)

	//!! EVT_TIMER(ID_TIMER, MyFrame::OnTimer)

	//!! EVT_BUTTON(ID_START_BUTTON, MyFrame::OnStartButton)
	//!! EVT_BUTTON(ID_UNDEF_BUTTON, MyFrame::OnUndefButton)
	//!! EVT_BUTTON(ID_CS_BUTTON, MyFrame::OnCSButton)
	// top tools
	//!! EVT_CHOICE(ID_XPAR, MyFrame::OnXChoice)
	EVT_MENU_RANGE(ID_XPAR, ID_YPAR-1, MyFrame::OnXChoice)
	EVT_MENU_RANGE(ID_YPAR, ID_ZPAR-1, MyFrame::OnYChoice)
	EVT_MENU_RANGE(ID_ZPAR, ID_DONE-1, MyFrame::OnZChoice)
	EVT_MENU(ID_DONE, MyFrame::OnDone)
	// Curves menu
	//!! EVT_MENU_RANGE(ID_INCLUDE, ID_EXCLUDE-1, MyFrame::OnInclude)
	//!! EVT_MENU_RANGE(ID_EXCLUDE, ID_CURVES_DONE-1, MyFrame::OnExclude)
	EVT_MENU(ID_CURVES_SELECT, MyFrame::OnCurvesSelect)
	// EVT_MENU(ID_CURVES_DONE, MyFrame::OnCurvesDone)
	// left side toolbar
	EVT_CHOICE(ID_CURVE_NUMBER, MyFrame::OnCurveNumber)
	//!! EVT_RADIOBOX(ID_RADIO_VIEW, MyFrame::OnRadioView)
	EVT_RADIOBUTTON(ID_Z_VIEW, MyFrame::OnZView)
	EVT_RADIOBUTTON(ID_Y_VIEW, MyFrame::OnYView)
	EVT_RADIOBUTTON(ID_X_VIEW, MyFrame::OnXView)
	EVT_RADIOBUTTON(ID_DEFAULT_VIEW, MyFrame::OnDefaultView)
	// XXX can I use ID_CURVES_SELECT in 2 places?
	EVT_BUTTON(ID_CURVES_SELECT, MyFrame::OnCurvesSelect)
	// EVT_BUTTON(ID_PARAM_SELECT, MyFrame::OnParamSelect)
	EVT_BUTTON(ID_XPARAM_SELECT, MyFrame::OnXParamSelect)
	EVT_BUTTON(ID_YPARAM_SELECT, MyFrame::OnYParamSelect)
	EVT_BUTTON(ID_ZPARAM_SELECT, MyFrame::OnZParamSelect)
	EVT_BUTTON(ID_PLOT, MyFrame::OnPlot)
	EVT_BUTTON(ID_QUIT_BUTTON, MyFrame::OnMenuFileExit)
END_EVENT_TABLE()

// MyFrame constructor
#include "vz3d.xpm"
MyFrame::MyFrame(wxFrame *frame, const wxString& title, const wxPoint& pos,
    const wxSize& size, long style) : wxFrame(frame, wxID_ANY, title, pos, size, style) {
	 T_(Trace trc(2,"MyFrame constructor");)

	// an icon for the app bar
	SetIcon(wxIcon(vz3d_xpm));

	 // add menubar, toolbar + controls
	 createControls();

	wxGLAttributes dispAttrs;
	// dispAttrs.PlatformDefaults().MinRGBA(8,8,8,8).DoubleBuffer().Depth(32).EndList();
	// dispAttrs.PlatformDefaults().RGBA().DoubleBuffer().Depth(32).EndList();
	dispAttrs.Defaults().EndList();
	long sty{style|wxSUNKEN_BORDER|wxFULL_REPAINT_ON_RESIZE};	// XXX deprecate
	// XXX only need (this,dispAttrs)
#ifdef NEVER // try size
	m_canvas = new FGLCanvas(this, dispAttrs);
#else // NEVER // try size
	m_canvas = new FGLCanvas(this, dispAttrs,wxID_ANY, wxDefaultPosition, wxSize(900,600));
#endif // NEVER // try size
	//!! wxID_ANY, wxDefaultPosition, GetClientSize(), sty, wxT("vz3d"));
}

#ifdef NEVER // no timer
bool
MyFrame::
isRunning() {
	return (m_timer.IsRunning());
}
#endif // NEVER // no timer

void
MyFrame::
createControls() {
// create all controls for this frame:
	wxPoint pos{wxDefaultPosition};
	wxSize sz{wxDefaultSize};
	Vz3d* vz3d{Vz3d::instance()};

    // Make the "File" menu
    wxMenu *fileMenu = new wxMenu;
    //fileMenu->Append(wxID_OPEN, wxT("&Open..."));
    //fileMenu->AppendSeparator();
    fileMenu->Append(wxID_PRINT, wxT("Print\tALT-P"));
    fileMenu->Append(ID_XFIG, wxT("Xfig\tALT-X"));
    fileMenu->AppendSeparator();
    fileMenu->Append(wxID_EXIT, wxT("&Exit\tq"));
	 // make the "Edit" menu
    wxMenu *editMenu = new wxMenu;
    editMenu->Append(ID_EDIT_PREF, wxT("&Preferences\tALT-E"));

	// make the "View" menu
   wxMenu *viewMenu = new wxMenu;
	// choose X submenu
	wxMenu* xparmenu = new wxMenu;
	int npar = vz3d->parnames().size();
	for (int i=0; i<npar; i++)
		xparmenu->AppendRadioItem(ID_XPAR+i, wxString(vz3d->parnames()[i]));
	viewMenu->AppendSubMenu(xparmenu, wxT("X param"));
	// choose Y submenu
	wxMenu* yparmenu = new wxMenu;
	for (int i=0; i<npar; i++)
		yparmenu->AppendRadioItem(ID_YPAR+i, wxString(vz3d->parnames()[i]));
	viewMenu->AppendSubMenu(yparmenu, wxT("Y param"));
	// choose Z submenu
	wxMenu* zparmenu = new wxMenu;
	for (int i=0; i<npar; i++)
		zparmenu->AppendRadioItem(ID_ZPAR+i, wxString(vz3d->parnames()[i]));
	viewMenu->AppendSubMenu(zparmenu, wxT("Z param"));
	// "done" button
	//!! wxButton* doneButton = new wxButton(viewMenu, ID_DONE, wxT("Done"));
	//!! viewMenu->AddControl(doneButton, wxT("Done"));
	viewMenu->Append(ID_DONE, wxT("Done"));

	// curve pick menu
   curvesMenu = new wxMenu;
	// select menu item
	curvesMenu->Append(ID_CURVES_SELECT, wxT("Select"));

    // Make the "Help" menu
    wxMenu *helpMenu = new wxMenu;
    helpMenu->Append(ID_HELP_HELP, wxT("&Help..."));
    helpMenu->Append(ID_HELP_ABOUT, wxT("&About..."));

	 // Add these menus to the menu bar
    wxMenuBar *menuBar = new wxMenuBar;
    menuBar->Append(fileMenu, wxT("&File"));
    menuBar->Append(editMenu, wxT("&Edit"));
    menuBar->Append(viewMenu, wxT("Parameters"));
    menuBar->Append(curvesMenu, wxT("Curves"));
    menuBar->Append(helpMenu, wxT("&Help"));
    SetMenuBar(menuBar);


	//------------------------------------------------------------------
	// toolbar on the top
	//------------------------------------------------------------------
#ifdef NEVER // use includemenu
	wxToolBar* toptools = new wxToolBar(this, ID_TOOLBAR);
	wxArrayString par_names;
	for (auto& pc : vz3d->parnames)
		par_names.Add(pc);
#ifdef NEVER // try listbox
	wxChoice* xchoice = new wxChoice(toptools, ID_XPAR, pos, sz, par_names);
#else // NEVER // try listbox
	wxListBox* xchoice = new wxListBox(toptools, ID_XPAR, pos, sz, par_names, wxLB_SINGLE);
#endif // NEVER // try listbox
	xchoice->SetToolTip(wxT("choose X parameter"));
	toptools->AddControl(xchoice, wxT("X param"));
#endif // NEVER // use includemenu

	//------------------------------------------------------------------
	// toolbar on the left side
	//------------------------------------------------------------------
	//!! long style{wxTB_VERTICAL|wxTB_TEXT};
	long style{wxTB_VERTICAL};
	 //!! wxToolBar* toolbar = CreateToolBar(style, ID_TOOLBAR);
	wxToolBar* toolbar = new wxToolBar(this, ID_TOOLBAR, pos, sz, style);
	// toolbar->SetToolPacking(5);
	// toolbar->SetToolSeparation(10); // XXX doesn't work?
	// toolbar->SetMargins(2,4);
	// toolbar->SetMargins(2,50);

	// curves select button
	toolbar->AddSeparator(); toolbar->AddSeparator(); toolbar->AddSeparator();
	wxButton* curvesButton = new wxButton(toolbar, ID_CURVES_SELECT, wxT("Curves"));
	curvesButton->SetToolTip(wxT("select curves"));
	toolbar->AddControl(curvesButton, wxT("Select Curves"));

	toolbar->AddSeparator(); toolbar->AddSeparator(); toolbar->AddSeparator();
	// x parameter select button
	wxButton* xparButton = new wxButton(toolbar, ID_XPARAM_SELECT, wxT("X param"));
	xparButton->SetToolTip(wxT("choose X param"));
	toolbar->AddControl(xparButton, wxT("X par"));
	// y parameter select button
	wxButton* yparButton = new wxButton(toolbar, ID_YPARAM_SELECT, wxT("Y param"));
	yparButton->SetToolTip(wxT("choose Y param"));
	toolbar->AddControl(yparButton, wxT("Y par"));
	// z parameter select button
	wxButton* zparButton = new wxButton(toolbar, ID_ZPARAM_SELECT, wxT("Z param"));
	zparButton->SetToolTip(wxT("choose Z param"));
	toolbar->AddControl(zparButton, wxT("Z par"));
	// create the plot button
	wxButton* plotButton = new wxButton(toolbar, ID_PLOT, wxT("plot"));
	plotButton->SetToolTip(wxT("draw the plot"));
	toolbar->AddControl(plotButton, wxT("Plot"));

	// make a radio box for choosing view
	toolbar->AddSeparator(); toolbar->AddSeparator(); toolbar->AddSeparator();
#ifdef NEVER // try radioButtons
	wxArrayString axes;
	axes.Add("Iso");
	axes.Add(vastr(vz3d->xname(),"-",vz3d->yname())); // Zquat
	axes.Add(vastr(vz3d->xname(),"-",vz3d->zname())); // Yquat
	axes.Add(vastr(vz3d->yname(),"-",vz3d->zname())); // Xquat
	radioView = new wxRadioBox(toolbar, ID_RADIO_VIEW, wxString("Views"),
		pos, wxDefaultSize, axes, 0, wxRA_SPECIFY_ROWS );
	toolbar->AddControl(radioView, wxString("Views"));
#else // NEVER // try radioButtons
	string defaultlabel = "Default"; // Default_quat
	long radiostyle{wxRB_GROUP};
	defaultbutton = new wxRadioButton(toolbar, ID_DEFAULT_VIEW, defaultlabel, pos, sz, radiostyle);

	toolbar->AddControl(defaultbutton, wxString("Default view"));
	radiostyle = 0;
	// string zview = vastr(vz3d->xname(),"-",vz3d->yname()); // Zquat
	string zview;
	zbutton = new wxRadioButton(toolbar, ID_Z_VIEW, zview, pos, sz, radiostyle);
	//!! zview = vastr("<span foreground='red'>",vz3d->xname(),
		//!! "</span>-<span foreground='green'>",vz3d->yname(),"</span>");
	zview = vastr("<b>",vz3d->xname(),
		"</b>-<span foreground='green'>",vz3d->yname(),"</span>");
	zbutton->SetLabelMarkup(zview);
	toolbar->AddControl(zbutton, wxString("Z view"));
	string yview = vastr(vz3d->xname(),"-",vz3d->zname()); // Yquat
	ybutton = new wxRadioButton(toolbar, ID_Y_VIEW, yview);
	toolbar->AddControl(ybutton, wxString("Y view"));
	string xview = vastr(vz3d->yname(),"-",vz3d->zname()); // Xquat
	xbutton = new wxRadioButton(toolbar, ID_X_VIEW, xview);
	toolbar->AddControl(xbutton, wxString("X view"));
#endif // NEVER // try radioButtons

	// quit button
	toolbar->AddSeparator(); toolbar->AddSeparator(); toolbar->AddSeparator();
	wxButton* quitButton = new wxButton(toolbar, ID_QUIT_BUTTON, wxT("&Quit"));
	quitButton->SetToolTip(wxT("quit"));
	// toolbar->AddStretchableSpace();
	toolbar->AddControl(quitButton, wxT("Quit"));

	toolbar->Realize();
	this->SetToolBar(toolbar);

	 statusBar = new wxStatusBar(this, wxID_ANY, wxST_SIZEGRIP);
	 this->SetStatusBar(statusBar);
	 statusBar->SetStatusText(wxT("Ready"), 0);

	 // SetSizer(topsizer);
} // MyFrame::createControls()


bool
MyFrame::
ShowToolTips() { return true; }

#ifdef NEVER // no csbutton
void
MyFrame::
OnCSButton(wxCommandEvent& event) {
	//! Drawcs = !Drawcs;
	int current = Amdata::instance()->draw_coord_sys();
	if (current == 0)
		Amdata::instance()->draw_coord_sys(1);
	else
		Amdata::instance()->draw_coord_sys(0);
}
#endif // NEVER // no csbutton


// File|Open... command
void
MyFrame::
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
MyFrame::
OnMenuFilePrint( wxCommandEvent& WXUNUSED(event) ) {
	// print() uses wxPrinter stuff, printPS uses gl2ps
	// print();
	printPS();
}

// File|Xfig... command
void
MyFrame::
OnMenuFileXfig( wxCommandEvent& WXUNUSED(event) ) {
	//!! int period = m_timer.GetInterval();
	// stop the oscillations
	//!! m_timer.Stop();
	// call the FGLCanvas xfig creator
	m_canvas->xfig();
	// and re-start
	//!! m_timer.Start(period);
}

wxPrintDialogData g_printDialogData;

void
MyFrame::
print() {
	//!! int period = m_timer.GetInterval();
	// stop the oscillations
	//!! m_timer.Stop();
	
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

	// and re-start
	//!! m_timer.Start(period);
}

void
MyFrame::
printPS() {
	//!! int period = m_timer.GetInterval();
	// stop the oscillations
	//!! m_timer.Stop();
	// call the FGLCanvas printer
	m_canvas->print();
	// and re-start
	//!! m_timer.Start(period);
}

// File|Exit command or simply type "q"
void
MyFrame::
OnMenuFileExit( wxCommandEvent& WXUNUSED(event) ) {
	// write out the current quaternion to use as the initial
	// position in FGlCanvas constructor
	cerr << "current quaternion: " << m_canvas->gldata()->quat << endl;
	float theta = acos(m_canvas->gldata()->quat[0])*2.0;
	float st = sin(theta/2.0);
	float* v = &m_canvas->gldata()->quat[1];
	cerr << "angle: " << theta*57.295779 << " deg\n";
	cerr << "vector: " << v[0]*st << ", " << v[1]*st << ", " << v[2]*st << endl;
    // true is to force the frame to close
    Close(true);
}

void
MyFrame::
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

void
MyFrame::
OnMenuViewDefault(wxCommandEvent& WXUNUSED(event)) {
	ViewX = false;
	ViewY = ViewZ = false;
	ViewQuad = false;
	// int period = m_timer.GetInterval();
	// m_timer.Stop();
	m_canvas->gldata()->quat = Default_quat;
	// m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewX(wxCommandEvent& event) {
	m_canvas->gldata()->quat = Xquat;
	m_canvas->Refresh(false);
}

// top tools:
void
MyFrame::
OnXChoice(wxCommandEvent& event) {
// just get the new xname
	T_(Trace trc(1,"OnXChoice");)
	Vz3d* vz3d{Vz3d::instance()};
	int k = event.GetId()-ID_XPAR;		// 0b parameter number
	vz3d->xname(vz3d->parnames()[k]);
	T_(trc.dprint("selection: ",k," x: ",vz3d->xname());)
#ifdef NEVER // just get xname
	int ncurves = vz3d->data.size();
	for (int j=0; j<ncurves; j++) {
		Par* xp = vz3d->curves[j]->findp(vz3d->xname());
		double xcf{xp->conv.factor};
		int npt = xp->solns.size();
		for (int i=0; i<npt; i++) {
			double x = xp->solns[i]*xcf;
			vz3d->data[j][3*i] = x;
		}
	}
	m_canvas->Refresh(false);
#endif // NEVER // just get xname
}

void
MyFrame::
OnYChoice(wxCommandEvent& event) {
// just get the new yname
	T_(Trace trc(1,"OnYChoice");)
	Vz3d* vz3d{Vz3d::instance()};
	int k = event.GetId()-ID_YPAR;		// 0b parameter number
	vz3d->yname(vz3d->parnames()[k]);
	T_(trc.dprint("selection: ",k," y: ",vz3d->yname());)
}

void
MyFrame::
OnZChoice(wxCommandEvent& event) {
// just get the new zname
	Vz3d* vz3d{Vz3d::instance()};
	int k = event.GetId() - ID_ZPAR;		// 0b parameter number
	vz3d->zname(vz3d->parnames()[k]);
}

void
MyFrame::
OnDone(wxCommandEvent& event) {
// plot after choosing new x, y, z parameters
	Vz3d* vz3d{Vz3d::instance()};
	vz3d->load_data();
	m_canvas->Refresh(false);
	m_canvas->Render();
}

#ifdef NEVER // unused
void
MyFrame::
OnInclude(wxCommandEvent& event) {
// get a list of curves to include
	int k = event.GetId() - ID_INCLUDE;		// 0b parameter number
	Vz3d* vz3d{Vz3d::instance()};
	if (std::find(vz3d->toplot.begin(), vz3d->toplot.end(), k) == vz3d->toplot.end())
		vz3d->toplot.push_back(k);
}

void
MyFrame::
OnExclude(wxCommandEvent& event) {
// get a list of curves to exclude
}
#endif // NEVER // unused

void
MyFrame::
OnCurvesSelect(wxCommandEvent& event) {
// pop up a dialog for selecting which curves to plot
	T_(Trace trc(1,"OnCurvesSelect");)
	Vz3d* vz3d{Vz3d::instance()};
	wxArrayString curve_names;
	curve_names.Add("All");
	for (auto& pc : vz3d->curves)
		curve_names.Add(pc->cid());
	// create the dialog
	wxMultiChoiceDialog curvestool(this,"Curve Selection", "Curves", curve_names);
	// checkmark selections which are in "toplot" - all if toplot is empty
	T_(trc.dprint("input toplot: ",vz3d->toplot);)
	// put checkmarks on those already in "toplot"
	wxArrayInt checkmarks;
	if (vz3d->toplot.empty()) {
		for (size_t i=1; i<curve_names.size(); i++)
			checkmarks.push_back(i);
	} else {
		for (size_t i=0; i<vz3d->toplot.size(); i++)
			checkmarks.Add(vz3d->toplot[i]+1);
	}
	curvestool.SetSelections(checkmarks);

	// show the dialog, get selections
	int stat = curvestool.ShowModal();
	T_(trc.dprint("dialog returned ",stat);)
	if (stat == wxID_CANCEL)
		return;
	wxArrayInt sel = curvestool.GetSelections();
	if (sel.empty())
		return;
	// empty out toplot 
	vz3d->toplot.clear();
	checkmarks.clear();
	// if sel == 0 => all curves so just return with toplot empty
	if (sel[0] == 0) {
		for (size_t i=1; i<curve_names.size(); i++)
			checkmarks.push_back(i);
		curvestool.SetSelections(checkmarks);
		return;
	}
	// otherwise sel are 1b curve numbers
	for(size_t i=0; i<sel.size(); i++) {
		T_(trc.dprint("got ",sel[i]);)
		vz3d->toplot.push_back(sel[i]-1);
		checkmarks.push_back(sel[i]);
	}
	T_(trc.dprint("output toplot: ",vz3d->toplot);)
	m_canvas->Refresh(false);
	m_canvas->Render();
}

string
MyFrame::
getParChoice(const string& current, const string& axis) {
// popup a dialog with all parameter names to get the choice of one
	T_(Trace trc(1, "getParChoice");)
	string rval;

	Vz3d* vz3d{Vz3d::instance()};
	wxArrayString par_names;
	const vector<string>& parnames = vz3d->parnames();
	for (const auto& pc : parnames)
		par_names.Add(pc);

	// create the dialog with all parameter names
	wxSingleChoiceDialog pardlg(this,vastr(axis," Parameter"), "choose", par_names);

	// checkmark the parameter which is the current name
	if (!current.empty()) {
		T_(trc.dprint("input name: ",current);)
		for (int i=0; i<parnames.size(); i++) {
			if (parnames[i] == current) {
				pardlg.SetSelection(i);
				break;
			}
		}
	}

	// show the dialog, get selections
	int stat = pardlg.ShowModal();
	T_(trc.dprint("dialog returned ",stat);)
	if (stat == wxID_CANCEL)
		return rval;
	int sel = pardlg.GetSelection();

	rval = parnames[sel];
	T_(trc.dprint("output parameter name: ",rval);)
	return rval;
}

void
MyFrame::
OnXParamSelect(wxCommandEvent& event) {
// pop up a dialog for selecting which parameter is x
	T_(Trace trc(1,"OnXParamSelect");)
	Vz3d* vz3d{Vz3d::instance()};
	string newname = getParChoice(vz3d->xname(), "Z");
	if (!newname.empty())
		vz3d->xname(newname);
}

void
MyFrame::
OnYParamSelect(wxCommandEvent& event) {
// pop up a dialog for selecting which parameter is y
	T_(Trace trc(1,"OnYParamSelect");)
	Vz3d* vz3d{Vz3d::instance()};
	string newname = getParChoice(vz3d->yname(), "Y");
	if (!newname.empty())
		vz3d->yname(newname);
}

void
MyFrame::
OnZParamSelect(wxCommandEvent& event) {
// pop up a dialog for selecting which parameter is z
	T_(Trace trc(1,"OnZParamSelect");)
	Vz3d* vz3d{Vz3d::instance()};
	string newname = getParChoice(vz3d->zname(), "Z");
	if (!newname.empty())
		vz3d->zname(newname);
}

void
MyFrame::
OnPlot(wxCommandEvent& event) {
// replot after choosing new parameter or curves
	T_(Trace trc(1,"OnPlot");)
	Vz3d* vz3d{Vz3d::instance()};
	std::sort(vz3d->toplot.begin(), vz3d->toplot.end());
	vz3d->load_data();
	// update the perspectives radiobox
	zbutton->SetLabel(vastr(vz3d->xname(),"-",vz3d->yname())); // Zquat
	ybutton->SetLabel(vastr(vz3d->xname(),"-",vz3d->zname())); // Yquat
	xbutton->SetLabel(vastr(vz3d->yname(),"-",vz3d->zname())); // Xquat
	// radioView->SetLabel(2,vastr(vz3d->xname(),"-",vz3d->zname())); // Yquat
	// radioView->SetLabel(3,vastr(vz3d->yname(),"-",vz3d->zname())); // Xquat
	m_canvas->Refresh(false);
	m_canvas->Render();
}

#ifdef NEVER // unused
void
MyFrame::
OnCurvesDone(wxCommandEvent& event) {
// plot after choosing new curves
	Vz3d* vz3d{Vz3d::instance()};
	std::sort(vz3d->toplot.begin(), vz3d->toplot.end());
	vz3d->load_data();
	m_canvas->Refresh(false);
	m_canvas->Render();
}
#endif // NEVER // unused

// Toolbar events
void
MyFrame::
OnCurveNumber(wxCommandEvent& event) {
	int i = event.GetSelection();
	Vz3d* vz3d{Vz3d::instance()};
	vz3d->toplot.clear();	// all curves
	if (i > 0)			// i == 0 => all
		vz3d->toplot.push_back(i);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnRadioView(wxCommandEvent& event) {
	int axis = event.GetSelection();
cerr << "axis: " << axis << endl;
	if (axis == 0)
		m_canvas->gldata()->quat = Default_quat;
	else if (axis == 1)
		m_canvas->gldata()->quat = Zquat;
	else if (axis == 2)
		m_canvas->gldata()->quat = Yquat;
	else
		m_canvas->gldata()->quat = Xquat;
	m_canvas->Refresh(false);
}

void
MyFrame::
OnZView(wxCommandEvent& event) {
	m_canvas->gldata()->quat = Zquat;
	m_canvas->Refresh(false);
}

void
MyFrame::
OnYView(wxCommandEvent& event) {
	m_canvas->gldata()->quat = Yquat;
	m_canvas->Refresh(false);
}

void
MyFrame::
OnXView(wxCommandEvent& event) {
	m_canvas->gldata()->quat = Xquat;
	m_canvas->Refresh(false);
}

void
MyFrame::
OnDefaultView(wxCommandEvent& event) {
	m_canvas->gldata()->quat = Default_quat;
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewY(wxCommandEvent& WXUNUSED(event)) {
	//!! int period = m_timer.GetInterval();
	//!! m_timer.Stop();
	m_canvas->gldata()->quat = Yquat;
	//!! m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewZ(wxCommandEvent& WXUNUSED(event)) {
	//!! int period = m_timer.GetInterval();
	//!! m_timer.Stop();
	m_canvas->gldata()->quat = Zquat;
	//!! m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewI(wxCommandEvent& WXUNUSED(event)) {
	//!! int period = m_timer.GetInterval();
	//!! m_timer.Stop();
	m_canvas->gldata()->quat = Isoquat;
	//!! m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewQuad(wxCommandEvent& WXUNUSED(event)) {
	if (!ViewQuad) {
		ViewQuad = true;
		m_canvas->viewQuad();
	}
}

#ifdef NEVER // no color
void
MyFrame::
OnMenuViewColorAbs(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_ABS;
	Amdata::colormode("abs value");
}

void
MyFrame::
OnMenuViewColorVal(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_VAL;
	Amdata::colormode("value");
}

void
MyFrame::
OnMenuViewColorPhase(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_PHASE;
	Amdata::colormode("phase");
}

void
MyFrame::
OnMenuViewColorPhaseAbs(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_PHASE_ABS;
	Amdata::colormode("abs phase");
}

void
MyFrame::
OnMenuViewColorEnergy(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_ENERGY;
	Amdata::colormode("energy");
}

void
MyFrame::
OnMenuViewColorNone(wxCommandEvent& WXUNUSED(event)) {
	// ViewColor = ID_VIEW_COLOR_NONE;
	Amdata::colormode("none");
}
#endif // NEVER // no color

#ifdef NEVER // no mode selection
void
MyFrame::
OnMenuMode( wxCommandEvent& event ) {
	// the items in the choice box are in reverse order
	// (most recent first) so GetSelection 0 is the last
	// column of Disp
	Amdata::picked_point(event.GetSelection());
}
#endif // NEVER // no mode selection

// Help|Help... command
void
MyFrame::
OnMenuHelpHelp( wxCommandEvent& WXUNUSED(event) ) {
    wxMessageBox(wxT("This is a work in progress but here are some features:\n"
			 " - left mouse button will rotate the figure\n"
			 " - middle mouse button translates the figure\n"
			 " - right mouse button zooms:\n"
			 "     - moving it towards the center zooms out,\n"
			 "     - moving it away from the center zooms in\n"
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
MyFrame::
OnMenuHelpAbout( wxCommandEvent& WXUNUSED(event) ) {
    wxMessageBox(wxT("   Animated Mode Visualizer\n"
			"intended to be run either by right-clicking on\n"
			"curves in vz or as a standalone program (amvz)\n"
			"taking input from Universal-formatted files.\n"
			"In either case a set of nodes, connectivities,\n"
			"and modes must be available\n"
		), wxT("About"), wxICON_INFORMATION);

}

// FGLCanvas implementation

BEGIN_EVENT_TABLE(FGLCanvas, wxGLCanvas)
    EVT_SIZE(FGLCanvas::OnSize)
    EVT_PAINT(FGLCanvas::OnPaint)
    EVT_ERASE_BACKGROUND(FGLCanvas::OnEraseBackground)
    EVT_MOUSE_EVENTS(FGLCanvas::OnMouse)
    // EVT_CHAR(FGLCanvas::OnChar)
END_EVENT_TABLE()

FGLCanvas::FGLCanvas(wxWindow *parent, const wxGLAttributes& dispAttrs,
		wxWindowID id, const wxPoint& pos, const wxSize& size, long style, const wxString& name) :
			wxGLCanvas(parent, dispAttrs) {
// FGlCanvas constructor

   m_gldata.initialized = false;
	m_init = false;

   m_gldata.beginx = 0.0f;
   m_gldata.beginy = 0.0f;
   m_gldata.zoom   = 1.0f;
   // Set the default quat, initialize view matrix with it
	//!! Default_quat = airplane_quat;
	m_gldata.quat = Default_quat;
	float qn = 1.0/blas_snrm2(4, &m_gldata.quat[0], 1);
	if (abs(qn-1.0) > 8.0*std::numeric_limits<float>::epsilon()) {
		blas_scal(4, qn, &m_gldata.quat[0], 1);
		cerr << "inital quat scaled by " << qn << ": " << m_gldata.quat << endl;
	}
}

void
FGLCanvas::
setcurrent() {
	static wxGLContext* context{nullptr};
	// create a wxGLContext the first time this is called
	if (context == nullptr)
		context = new wxGLContext(this, 0);

	// ... and call it's SetCurrent
	bool stat = context->SetCurrent(*this);
	if (!stat)
		throw runtime_error("wxGLContext::SetCurrent failed");
}

void
FGLCanvas::
Render() {
// XXX merge with Vz3d::render()
	T_(Trace trc(1,"Render");)
	this->setcurrent();
	//!! vz3d->render();
	Vz3d::instance()->render();
	glFlush();
	wxGLCanvas::SwapBuffers();

	// Theframe->Refresh(false);
} // Render()

void
FGLCanvas::
OnPaint( wxPaintEvent& WXUNUSED(event) ) {

	wxPaintDC dc(this);

	//! if (!GetContext()) return;

	// direct gl commands to this FGLCanvas
	//! wxGLCanvas::SetCurrent();
	 // wxGLContext* context = new wxGLContext(this, 0);
	 // context->SetCurrent(*this);
	 // SetCurrent(this);
	 setcurrent();

    // Init OpenGL once, but after SetCurrent
    if (!m_init) {
        InitGL();
        m_init = true;
    }

    // Transformations
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
	 //glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
	 //glRotatef(135.0f, 0.0f, 0.0f, 1.0f);
    glTranslatef( OriginX, OriginY, -2.0f );

	 // viewing area (for zooming):
	glOrtho(ViewLeft, ViewRight, ViewBottom, ViewTop, ViewNear, ViewFar);

    GLfloat m[4][4];
	 // trackball:
    build_rotmatrix( m, &m_gldata.quat[0] );
    glMultMatrixf( &m[0][0] );
	// XXX try this:
	Render();
} // OnPaint()

void
FGLCanvas::
OnSize(wxSizeEvent& event) {
    // this is also necessary to update the context on some platforms
	 // replacing wxGLCanvas::OnSize() with event.Skip()
	 // per http://groups.google.com/group/wx-users/browse_thread/thread/cfa90147ee800f0c?pli=1
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
OnMouse(wxMouseEvent& event) {

	if (event.Dragging()) {
		float x = event.GetX();
		float y = event.GetY();
		// dragging the mouse turns off Quad, X, Y, and Z view modes:
		ViewX = ViewY = ViewZ = false;
		ViewQuad = false;

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
				// XXX change to beginx, remove lastX from amvz.h:
			OriginX += (event.GetX() - m_gldata.lastX)/1200.0;
			OriginY -= (event.GetY() - m_gldata.lastY)/1200.0;
			Refresh(false);
		}
	}
	m_gldata.beginx = event.GetX();
	m_gldata.beginy = event.GetY();
	m_gldata.lastX = event.GetX();
	m_gldata.lastY = event.GetY();
	Render();
	
	event.Skip();
}

#ifdef NEVER // unnecessary
void
FGLCanvas::
OnChar(wxKeyEvent& event) {
	setcurrent();
	if (event.GetKeyCode() == wxKeyCode('q'))
		Theframe->Close();
	if (event.GetKeyCode() == wxKeyCode('Q'))
		Theframe->Close();
	if (event.GetKeyCode() == wxKeyCode('x')) {
		gldata()->quat = Xquat;
		Render();
		// Refresh(false);
		event.Skip();
	} else if (event.GetKeyCode() == wxKeyCode('y'))
		gldata()->quat = Yquat;
	else if (event.GetKeyCode() == wxKeyCode('z'))
		gldata()->quat = Zquat;
	else
		event.Skip();
	
	Refresh(false);
	Render();
	//event.Skip();
}
#endif // NEVER // unnecessary

void
FGLCanvas::
InitGL() {

    /* set viewing projection */
	glMatrixMode(GL_PROJECTION);
	// size of the dots:
	glPointSize(2.0);
	// the nodes have been scaled so the largest coord is 1
	// and centered, so make the viewing area a little larger
	// glOrtho(left, right, bottom, top, near, far)
	glOrtho(ViewLeft, ViewRight, ViewBottom, ViewTop, ViewNear, ViewFar);

    /* position viewer */
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(0.0f, 0.0f, -2.0f);

    /* position object */
    // glRotatef(30.0f, 1.0f, 0.0f, 0.0f);
    // glRotatef(30.0f, 0.0f, 1.0f, 0.0f);
	 glRotatef(-90.0f, 1.0f, 0.0f, 0.0f);
	 glRotatef(135.0f, 0.0f, 0.0f, 1.0f);

    //glEnable(GL_DEPTH_TEST);
    //glEnable(GL_LIGHTING);
    //glEnable(GL_LIGHT0);
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
#define XFIG 1
#ifdef XFIG
	FGLCanvas::xfig();
#else
	GLint viewport[4];
	GLint options = GL2PS_USE_CURRENT_VIEWPORT;
	GLint format = GL2PS_PS;
	GLint sort = GL2PS_NO_SORT;
	GLint colormode = GL_RGBA;
	GLint colorsize = 0;   // using GL_RGBA instead of color map
	GL2PSrgba* colortable = 0;   // using GL_RGBA instead of color map
	GLint nr = 0;
	GLint ng = 0;
	GLint nb = 0;
	GLint buffersize = 1200*900*8;
	char const* filename = "amvz.ps";
	FILE* stream = fopen(filename, "w");
	gl2psBeginPage("amvz", "amvz", viewport, format,
			sort, options, colormode, colorsize, colortable,
			nr, ng, nb, buffersize, stream, filename);
	Render();
	gl2psEndPage();

	ostringstream os;
	os << "Created postscript file " << filename;
	wxChar wst[256];
	for (size_t i=0; i<os.str().size(); i++)
		wst[i] = os.str()[i];
	wxString txt(wst, os.str().size());
	wxMessageBox(txt, _("printPS"), wxOK);
#endif // XFIG
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

		stream << "#FIG 3.2  Produced by amvz\n";
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
	// render into a gl feedback buffer
	int buffersize = 1000000;
	wxString wxpath = wxGetTextFromUser(wxT("Please enter a filename for xfig output "
							  "(include a *.fig extenstion)"), 
							  wxT("xfig filename"), wxT("amvz.fig"));
	string pathname = wxpath.ToStdString();
#ifdef NEVER // no steps
	wxString wxsteps = wxGetTextFromUser(wxT("How many steps between 0-360 deg?"),
										  wxT("# steps"), wxT("5"));
	string nstepstr = wxsteps.ToStdString();
	int nstep{5};
	bool ok = str2int(nstepstr, nstep);
	if (!ok)
		cerr << "not an int: " << nstepstr << endl;
#endif // NEVER // no steps

	vector<GLfloat> feedbuf(buffersize);

	bool append{false};
	int style{0};
	int color{0};
	int width{2};

	glFeedbackBuffer(feedbuf.size(), GL_2D, &feedbuf[0]);

	(void)glRenderMode(GL_FEEDBACK);
	Render();
	glRenderMode(GL_RENDER);
	xfig_create(pathname, feedbuf, append, style, width, color);
}

void
FGLCanvas::
viewQuad() {
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

	bottom_left;
	bottom;
	Render();

	bottom_right;
	left;
	Render();

	top_left;
	front;
	Render();

	top_right;
	perspective;
	Render();
}

#ifdef STANDALONE

#include <cassert>
// #include "main.h"

double XParValue;
double YParValue;
string ModeName;

int
main(int argc, char** argv) {
	char* dbg = getenv("DEBUG");
	if (dbg != nullptr) Trace::debug(std::stoi(dbg));
	T_(Trace trc(1,"vz3d");)
	Ftmpdir ftmp;

	if (argc < 1) {
		cerr << "Usage: vz3d\n";
		exit(1);
	}

	// parse the specs
	Specs& sp = specs();
	string options = cmdline(argc, argv);
	if (!parse_specs(options, sp))
		return 1;
	try {
		// read the curves
		vector<Plotcurve*> curves = Plotcurve::getcurves(sp);
		// create a Vz3d...
		Vz3d* vz3d = Vz3d::instance(curves, sp);
		// ... and run it
		vz3d->run();

	} catch (runtime_error& s) {
		flaps::error("caught exception: ",s.what());
	} catch (const std::exception& t) {
		flaps::error("caught exception: ",t.what());
	}
}
#endif // STANDALONE
