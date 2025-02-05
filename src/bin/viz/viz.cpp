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
// Flaps plots
/////////////////////////////////////////////////////////////////////////////

#include <cassert>
#include <cerrno>
#include <ctime>
#include <exception>
#include <fnmatch.h>
#include <stdexcept>
#include <sys/wait.h>
#include <thread>
#include <utility>
#include <iosfwd>
#include <wx/accel.h>
#include <wx/wx.h>
#include <wx/msgdlg.h>
#include <wx/scrolwin.h>
#include <wx/sizer.h>
#include <wx/colour.h>
#include <wx/filename.h>
#include <wx/string.h>
#include <wx/tglbtn.h>
#include <wx/toolbar.h>
#include <wx/prntbase.h>

#include "config.h"
#include "amviz.h"
#include "exim.h"
#include "fio.h"
#include "lexer.h"
#include "mathplot.h"
#include "matrix.h"
#include "plotcurve.h"
#include "specs.h"
#include "run_program.h"
#include "trace.h"
#include "viz.h"
#include "util.h"

using namespace std;

VizSizes vizsizes;
VizSizes::
VizSizes() {
// VizSizes ctor
	T_(Trace trc(2,"VizSizes ctor");)
	// display = wxGetDisplaySize();
	//!! wxDisplay dpy;
	//!! display = dpy.GetGeometry().GetSize();
	display = getDisplaySize();
#ifdef SQUARE
	// square plot
	int side = 0.7*std::min(display.x, display.y);
	plot = wxSize(side,side);
	panel = wxSize(140, plot.y);
	infobox = wxSize(panel.x, 50);
	legendbox = wxSize(panel.x, 550);
	footnotes = wxSize(panel.x, 200);
	int framex = plot.x + panel.x;
	menu = wxSize(framex, 27);
	toolbar = wxSize(framex, 42);
	frame = wxSize(framex, plot.y+toolbar.y+menu.y);
#else // SQUARE
	// 60/70% display but watch out for double-monitors
	if (display.x > 2*display.y)
		display.x = 1.2*display.y;
	frame = wxSize(0.6*display.x, 0.7*display.y);
	menu = wxSize(frame.x, 27);
	toolbar = wxSize(frame.x, 42);
	int panelx = 140;
	plot = wxSize(frame.x - panelx, frame.y - menu.y - toolbar.y);
	panel = wxSize(panelx, plot.y);
	infobox = wxSize(panel.x, 50);
	footnotes = wxSize(panel.x, 200);
	legendbox = wxSize(panel.x, plot.y - infobox.y - footnotes.y);
#endif // SQUARE
	T_(trc.dprint("display ", display);)
	T_(trc.dprint("plot ", plot);)
	T_(trc.dprint("panel ", panel);)
	T_(trc.dprint("legendbox ", legendbox);)
	T_(trc.dprint("infobox ", infobox);)
	T_(trc.dprint("footnotes ", footnotes);)
	T_(trc.dprint("menu ", menu);)
	T_(trc.dprint("toolbar ", toolbar);)
	T_(trc.dprint("frame ", frame);)
}

ostream&
operator<<(ostream& s, const VizSizes& t) {
	s << "display " << t.display << endl;
	s << "plot " << t.plot << endl;
	s << "panel " << t.panel << endl;
	s << "legendbox " << t.legendbox << endl;
	s << "infobox " << t.infobox << endl;
	s << "footnotes " << t.footnotes << endl;
	s << "menu " << t.menu << endl;
	s << "toolbar " << t.toolbar << endl;
	s << "frame " << t.frame << endl;
	return s;
}
	
#if !defined(__WXGTK__)
	#error "rebuild with WXGTK"
#endif

ostream&
operator<<(ostream& s, const vector<Vizplot*>& t) {
	for (auto& ti : t)
		s << ti << endl;
	return s;
}
 
void
VizFrame::
set_title(const string& t) {
// set the frame title and status bar to a list of plotfiles
// or aids or a user supplied title
	Specs& sp = specs();
	string title{t};
	if (frameno != 0)
		title = vastr("(",frameno,") ",t," ");

	// user-supplied title?
	if (!sp.title.empty())
		title += " " + sp.title;

	// file names or aids
	Viz& viz{Viz::instance()};
	vector<string> leaders;
	for (auto pc : viz.curves) {
		string li{pc->aid()};		// li: either aid or plotfile
		if (li.empty()) li = pc->path;
		add2vector(li, leaders);
	}
	string sep;
	for (auto& li : leaders) {
		title += sep + li;
		sep = ", ";
	}

	statusBar->SetStatusText(title);
	SetTitle(title);
}


static string
aboutmsg() {
	string rval{"Flaps 2D plotter\n"};
	rval +=  buildinfo();
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

static string
helpmsg() {
	string rval{"This is a work in progress but here are some features:\n"
		" - clicking the left mouse button pops up a list of parameter values\n"
		"   at the closest point\n"
		" - dragging with the left mouse button defines an area to zoom in\n"
		" - dragging with the right mouse button pans\n"
		" - middle mouse pops up a menu of plot manipulations\n"
		" - right mouse button starts amviz or adds a mode if already running\n"
		" - the wheel moves the plot up and down, shift-wheel sideways\n"
		" - The top row of buttons labeled:\n"
		"     * Curves: choose which curves to show\n"
		"     * Parameters: choose X and/or Y parameter(s)\n"
		"     * X-Y Views: choose previously shown combinations\n"
		"     * Curve colors: 3 options for Curve colors:\n"
		"       - by curve number\n"
		"       - by plotfile or aid\n"
		"       - by curve id extension, e.g. 1.a, 1.b, ...\n"
		"     * Legend/Values: toggles legend and parameter values \n"};
	rval +=  buildinfo();
	if (isWSL())
		rval += " (WSL)";
	Specs& sp = specs();
	if (sp.threads) {
		const int hwthreads = std::thread::hardware_concurrency()-1;
		rval += vastr(" (",hwthreads," threads)");
	}
	return rval;
}


void
message(const string& msg) {
    wxMessageBox(msg, wxT("Error"), wxICON_INFORMATION);
}

// LeftPanel implementation {
LeftPanel::
LeftPanel(LowerPanel* parent, wxWindowID id, const wxPoint& pos, const wxSize& size) :
		wxPanel(parent, id, pos, size), lowerpanel(parent) {
// A LeftPanel consist of an Infobox, a scrollable legend, and a
// scrollable footnote panel, all contained in a vertical sizer (lpsizer)
	T_(Trace trc(2,"LeftPanel ctor");)

	wxPoint defpos{wxDefaultPosition};

	// create a sizer
	lpsizer = new wxBoxSizer(wxVERTICAL);
	SetSizer(lpsizer);

	// info box
	infobox = new Infobox(this);
	lpsizer->Add(infobox);

	// create a scrollable subpanel: legendbox
	legendbox = new myScroll(this, wxID_ANY, defpos, vizsizes.legendbox);
	lpsizer->Add(legendbox);

	// create a scrollable footnotes panel
	footnotes = new Footnotes(this, wxID_ANY, defpos, vizsizes.footnotes);
	lpsizer->Add(footnotes);

	Bind(wxEVT_PAINT, &LeftPanel::OnPaint, this);

	lpsizer->Layout();

	T_(trc.dprint("LeftPanel sizer: ",sizer2string(*lpsizer));)
}	// LeftPanel ctor

void
LeftPanel::
OnPaint(wxPaintEvent& event) {
// just paint the footnotes if any
	wxPaintDC dc(this);

	Refresh();
}	// LeftPanel::OnPaint

// LeftPanel implementation }

// Infobox implementation {
Infobox::
Infobox(LeftPanel* parent) : wxPanel(parent, wxID_ANY,
	wxDefaultPosition, vizsizes.infobox, wxBORDER_SIMPLE) {
	T_(Trace trc(2,"Infobox ctor");)
	leftpanel = parent;
	Bind(wxEVT_PAINT, &Infobox::OnPaint, this);
}

void
Infobox::
OnPaint(wxPaintEvent& event) {
// paint the current info box coordinates (ibcoords)
// An Infobox's parent is a LeftPanel, who's parent is a LowerPanel,
// who's parent is a Vizframe,
// which contains a Vizplot* (currentview)
	T_(Trace trc(2,"Infobox::OnPaint");)
	wxPaintDC dc(this);
	dc.DrawText(leftpanel->lowerpanel->vizframe->currentview->ibcoords, 10, 5);
	Refresh();
}	// Infobox::OnPaint
// Infobox implementation }

// myScroll implementation {
myScroll::
myScroll(LeftPanel* p, wxWindowID id, const wxPoint& pos, const wxSize& size) :
	wxScrolled<wxPanel>(p, id, pos, size, wxBORDER_SIMPLE) {
	T_(Trace trc(2,"myScroll ctor");)

	leftpanel = p;

	Bind(wxEVT_PAINT, &myScroll::OnPaint, this);
	SetScrollRate(5,5);

	int width, height;
	this->GetSize(&width, &height);
	// everything must be written to "canvas" XXX even though I don't actually
	// draw on "canvas" the size determine scrolling size
	canvas = new wxPanel(this, wxID_ANY, wxPoint(0,0), wxSize(250,1200));
	//!! canvas = new wxPanel(this, wxID_ANY, wxPoint(0,0), wxSize(150,height+100));
	canvas->SetVirtualSize(250,900);

	// the drawing canvas needs to be in a sizer
	csizer = new wxBoxSizer(wxVERTICAL);
	csizer->Add(canvas, 1, wxEXPAND);
	SetSizer(csizer);
	T_(trc.dprint("sizer: ",sizer2string(*csizer));)
}

void
myScroll::
OnPaint(wxPaintEvent& event) {
// myScroll parent is a LeftPanel, who's parent is a LowerPanel,
// who's parent is a Vizframe, which has Vizplot* currentview
	//!! Viz& viz{Viz::instance()};
	//!! Vizplot* mp = viz.frame->plot();
	VizFrame* frame = leftpanel->lowerpanel->vizframe;
	Vizplot* mp = frame->currentview;
	wxPaintDC dc(this);		// must use a PaintDC
	PrepareDC(dc);

	if (mp == nullptr)
		return;

	if (mp->showLegend) {
		frame->legendButton->SetLabel("Values");
		paintLegend(dc, mp);
	} else {
		frame->legendButton->SetLabel("Legend");
		paintPset(dc, frame);
	}
	//!! LeftPanel->Refresh();
	Refresh();
	mp->Refresh();		// XXX I don't understand why this is needed
}
// myScroll implementation }

// Footnotes implementation {
Footnotes::
Footnotes(LeftPanel* p, wxWindowID id, const wxPoint& pos, const wxSize& size) :
	wxScrolled<wxPanel>(p, id, pos, size, wxBORDER_SIMPLE) {
	T_(Trace trc(2,"Footnotes ctor");)

	leftpanel = p;

	Bind(wxEVT_PAINT, &Footnotes::OnPaint, this);
	SetScrollRate(5,5);

	int width, height;
	this->GetSize(&width, &height);
	// everything must be written to "canvas" XXX even though I don't actually
	// draw on "canvas" the size determine scrolling size
	canvas = new wxPanel(this, wxID_ANY, wxPoint(0,0), wxSize(250,200));
	//!! canvas = new wxPanel(this, wxID_ANY, wxPoint(0,0), wxSize(150,height+100));
	canvas->SetVirtualSize(250,900);

	// the drawing canvas needs to be in a sizer
	csizer = new wxBoxSizer(wxVERTICAL);
	csizer->Add(canvas, 1, wxEXPAND);
	SetSizer(csizer);
	T_(trc.dprint("sizer: ",sizer2string(*csizer));)
}

void
Footnotes::
OnPaint(wxPaintEvent& event) {
// Footnotes parent is a LeftPanel, who's parent is a LowerPanel,
// who's parent is a Vizframe, which has Vizplot* currentview
	VizFrame* frame = leftpanel->lowerpanel->vizframe;
	Vizplot* mp = frame->currentview;
	// the footnotes are kept in the singlton Viz
	Viz& viz{Viz::instance()};
	wxPaintDC dc(this);		// must use a PaintDC
	PrepareDC(dc);

	if (mp == nullptr)
		return;

	int ix{20};
	int iy{20};
	for (auto& fi : viz.footnotes) {
		string label = fi.second.lead;
		int width, height;
		dc.GetTextExtent(label, &width, &height);
		wxPen lpen = fi.second.pen;
		// the dot...
		wxPen dotpen = lpen;
		dotpen.SetWidth(5);
		dotpen.SetColour(lpen.GetColour());
		dc.SetPen(dotpen);
		fi.second.drawDot(dc,ix,iy, 4);
		// ... and it's label
		dc.DrawText(label, ix+10, iy-height/2);
		iy += height + 4;
	}

	Refresh();
	mp->Refresh();		// XXX I don't understand why this is needed
}
// Footnotes implementation }

// Viz implementation {

vector<Plotcurve*>
Viz::
loadcurves(const vector<string>& aids, const vector<string>& plotfiles) {
// read Plotcurves from given plotfile/aid names, add to existing
	T_(Trace trc(2,"Viz::loadcurves");)

	// read the aid/plotfiles...
	vector<Plotcurve*> newcurves = Plotcurve::getcurves(aids, plotfiles);
	// ... add to existing curves
	for (auto nc : newcurves)
		curves.push_back(nc);

	// if Specs does not have a vzid, get one from a curve
	Specs& sp = specs();
	if (sp.vzid.empty() && !curves.empty()) {
		sp.vzid = curves[0]->vzid();
	}

	// make a list of all the parameter names from all the curves
	vector<Par*> theparams;
	for (auto ci : curves) {
		vector<Par*> theparams = ci->allpar();
		for (auto pp : theparams)
			add2vector(pp->name, parnames_);
	}
	
	return newcurves;
}	// Viz::loadcurves

void
Vizplot::
set_legends() {
// create a legend for each curve with name, color, linestyle, and dotstyle
// chosen to make identification easier
	T_(Trace trc(2,"Viz::set_legends");)
	string leadsep{"."};		// separator between leader and cid
	Specs& sp = specs();
	Viz& viz{Viz::instance()};

	// get all the plot layers
	vector<myFXYVector*> layers = plotlayers();

	// set up some data in each curve's Legend
	// check to see if there are multiple leaders and/or extensions
	vector<string> nums;
	vector<string> leaders;
	vector<string> extensions;
	for (auto ly : layers) {
		T_(trc.dprint("working on lead ",ly->legend.lead,", num ",ly->legend.num,", yname ",ly->legend.yname);)
		add2vector(ly->legend.num, nums);
		add2vector(ly->legend.lead, leaders);
		add2vector(ly->legend.ext, extensions);
	}
	T_(trc.dprint("leaders: ",leaders,", extensions: ",extensions,", colorby ",viz.colorby);)

	// if colorby has not been set, set to a default:
	// - from command line (sp.colorby)
	// - to PerFile if multiple leaders
	// - to PerExt if multiple extensions
	// - to PerY if multiple Y's
	// - to PerCurve otherwise
	if (viz.colorby == -1)
		viz.colorby = sp.colorby;
	if (viz.colorby == -1 && leaders.size() > 1)
		viz.colorby = PerFile;
	if (viz.colorby == -1 && extensions.size() > 1)
		viz.colorby = PerExt;
	if (viz.colorby == -1 && ynames.size() > 1)
		viz.colorby = PerY;
	if (viz.colorby == -1 && nums.size() > 1)
		viz.colorby = PerCurveNum;
	if (viz.colorby == -1)
		viz.colorby = PerCurve;

	// id_ = lead.num.ext yname
	// sid_ = short version of id_: lead, ext, yname only if multiple
	for (auto ly : layers) {
		if (ly->legend.lead.empty()) {
			ly->legend.id_ = ly->legend.num;
			ly->legend.sid_ = ly->legend.num;
		} else {
			ly->legend.id_ = ly->legend.lead + '.' + ly->legend.num;
			if (leaders.size() > 1)
				ly->legend.sid_ = ly->legend.lead + '.' + ly->legend.num;
			else
				ly->legend.sid_ = ly->legend.num;
		}

		if (!ly->legend.ext.empty()) {
			ly->legend.id_ += '.' + ly->legend.ext;
			ly->legend.sid_ += '.' + ly->legend.ext;
		}

		// yname?
		if (ynames.size() > 1) {
			ly->legend.id_ += " " + ly->legend.yname;
			ly->legend.sid_ += "." + ly->legend.yname;
		}

		// tokenize id_ XXX is toks used?
		ly->legend.toks = string2tok(ly->legend.id_, '.');
		T_(trc.dprint("id_: ",ly->legend.id_);)
		T_(trc.dprint("sid_: ",ly->legend.sid_);)
	}

	// set the colorkey for each curve depending on colorby

	T_(trc.dprint("colorby ",viz.colorby);)

	int linewidth{1};	// line width
	// colors
	vector<wxColour> colors;
	colors.push_back(getcolor("black"));
	// colors.push_back(getcolor("blue"));
	colors.push_back(getcolor("slate blue"));
	colors.push_back(getcolor("brown"));
	colors.push_back(getcolor("coral")); 		// "orange"
	colors.push_back(getcolor("cyan"));
	// colors.push_back(getcolor("dark orchid"));	// "purple"
	colors.push_back(getcolor("purple"));
	colors.push_back(getcolor("firebrick"));
	colors.push_back(getcolor("green"));
	colors.push_back(getcolor("magenta"));		// XXX use for am_spot?
	colors.push_back(getcolor("red"));
	//colors.push_back(getcolor("yellow"));	// use yellow for Mark

	// all line types
	vector<wxPenStyle> linetypes{wxPENSTYLE_SOLID, wxPENSTYLE_LONG_DASH,
			wxPENSTYLE_SHORT_DASH, wxPENSTYLE_DOT_DASH};

	// all dot types
	vector<Dottype> dottypes{Dottype::dot,
			Dottype::x,
			Dottype::triangle,
			Dottype::diamond,
			Dottype::invtriangle, Dottype::square};

	// make up a set of wxPens: color, linestyle (7x4=28)
	vector<wxPen> allpens;
	for (auto& line : linetypes)
		for (auto& color : colors)
			allpens.push_back(wxPen(color, linewidth, line));

	// set up a map of lead vs dottype
	int idot{0};
	for (auto& li : leaders)
		Legend::dots[li] = dottypes[idot++];
	if (idot == (int)dottypes.size())
		idot = 0;

	// set up a map of key vs pen: Legend::pens
	int ipen{0};
	int curveno{0};	// for PerCurve
	int npens = allpens.size();
	for (auto ly : layers) {
		string key;
		if (viz.colorby == PerCurve)
			key = to_string(++curveno);
		else if (viz.colorby == PerCurveNum)
			key = ly->legend.num;
		else if (viz.colorby == PerFile)
			key = ly->legend.lead;
		else if (viz.colorby == PerExt)
			key = ly->legend.ext;
		else if (viz.colorby == PerY)
			key = ly->legend.yname;
		// add this pen to Legend::pens if not already there
		if (Legend::pens.find(key) == Legend::pens.end()) {
			Legend::pens[key] = allpens[ipen++];
			if (ipen == npens)
				ipen = 0;
		}
		ly->legend.pen = Legend::pens[key];
		ly->SetPen(Legend::pens[key]);

		ly->legend.dot = Legend::dots[ly->legend.lead];
		// add to footnotes
		if (!ly->legend.lead.empty()) {
			if (viz.footnotes.find(ly->legend.lead) == viz.footnotes.end())
				viz.footnotes[ly->legend.lead] = ly->legend;
		}
	}
	T_(trc.dprint(viz.footnotes.size()," footnotes");)
}

string
Viz::
parTitle(const string& name) const {
// Returns a string which can be used as axis title, e.g.
//   "vtas (Velocity TAS) m/s"
	string rval{name};
	if (!name.empty() && !curves.empty()) {
		Par* par = curves[0]->findp(name);
		if (par == nullptr)
			return rval;
		string des = par->desc();
		if (!des.empty())
			rval += vastr(" (",des,")");
		if (par->has_conv())
			rval += " " + par->presUnits();
	}
	return rval;
}

void
Viz::
run() {
	T_(Trace trc(1,"Viz::run");)

#ifdef wxUSE_STD_IOSTREAM
	// redirect messages from wx to FTMP
	ofstream wxerr(vastr(getftmp(),"/wx.err"));
	wxLogStream wxlog(&wxerr);
#endif // wxUSE_STD_IOSTREAM

	char* argv[2];
	char* prog = strdup("viz");
	argv[0] = prog;
	argv[1] = nullptr;
	int argc = 1;
	try {
		T_(trc.dprint("calling wxEntry");)
		wxEntry(argc, argv);
		T_(trc.dprint("returned from wxEntry");)
	} catch (runtime_error const& s) {
		cerr << s.what() << std::endl;
	} catch (std::exception const& t) {
		cerr << t.what() << std::endl;
	}
}
// end Viz implementation }

// Vizplot implementation {

vector<myFXYVector*>
Vizplot::
plotlayers() {
// Returns a vector of all myFXYVector layers
	T_(Trace trc(2,"plotlayers");)
	vector<myFXYVector*> rval;
	for (unsigned int p2 = 0; p2 < CountAllLayers(); p2++) {
		myFXYVector* ly = dynamic_cast<myFXYVector*>(GetLayer(p2));
		if (ly != nullptr)
			rval.push_back(ly);
	}
	return rval;
}

vector<myFXYVector*>
Vizplot::
visplotlayers() {
// Returns a vector of all visible myFXYVector layers
	T_(Trace trc(2,"visplotlayers");)
	vector<myFXYVector*> rval;
	for (unsigned int p2 = 0; p2 < CountAllLayers(); p2++) {
		myFXYVector* ly = dynamic_cast<myFXYVector*>(GetLayer(p2));
		if (ly != nullptr && ly->IsVisible())
			rval.push_back(ly);
	}
	return rval;
}

long Vizplot::flags{0};

Vizplot::
Vizplot(LowerPanel* parent, wxWindowID id, const wxPoint& pos, const wxSize& size) :
			mpWindow(parent, id, pos, size, Vizplot::flags) {
	T_(Trace trc(2,"Vizplot ctor");)

	lowerpanel = parent;
	vizframe = parent->vizframe;

	SetMPScrollbars(false);
	// axes
	ticks = false;			// start with a grid

	// margins around plot
	SetMargins(vizsizes.top, vizsizes.right, vizsizes.bottom, vizsizes.left);

	EnableDoubleBuffer(false);
	// mouse events:
	Bind(wxEVT_LEFT_DOWN, &Vizplot::OnLeftDown, this);
	Bind(wxEVT_LEFT_UP, &Vizplot::OnLeftUp, this);
	Bind(wxEVT_MIDDLE_DOWN, &Vizplot::OnMiddleDown, this);
	Bind(wxEVT_MIDDLE_UP, &Vizplot::OnMiddleUp, this);
	Bind(wxEVT_RIGHT_DOWN, &Vizplot::OnRightDown, this);
	Bind(wxEVT_RIGHT_UP, &Vizplot::OnRightUp, this);
	Bind(wxEVT_MOTION, &Vizplot::OnMotion, this);
	Bind(wxEVT_PAINT, &Vizplot::OnPaint, this);

	// middle mouse-click menu: use mathplot handlers for fit & lock aspect
	middleMenu.Append(myBind(wxEVT_MENU, &Vizplot::OnOrigView, this),
		"Orig View", "Reset to original view");

#ifdef NEVER // I don't understand why this acts the way it does
	middleMenu.Append(mpID_LOCKASPECT, "Lock aspect",
		"Lock vertical/horizontal aspect ratio");
#endif // NEVER // I don't understand why this acts the way it does
	// use my own set of mouse actions, especially amviz
	middleMenu.Append(myBind(wxEVT_MENU, &Vizplot::OnMouseHelp, this),
			"Mouse actions...","Show help about mouse actions");
}	// Vizplot ctor

void
Vizplot::
loadxy(const string& xnm, const vector<string>& ynms, const vector<int>& toplot,
		double xmin, double xmax, double ymin, double ymax) {
// Second step in constructing a Vizplot: loading the x,y data
// from Viz. If xnm is empty create a set of points
// XXX rm xmin,xmax,ymin,ymax rely on clipping? otherwise the layer numbering
//     get screwed up
	T_(Trace trc(2,"Vizplot::loadxy");)
	Viz& viz{Viz::instance()};

#ifdef NEVER // ignore input xmin, etc: rely on clipping
	if (xmin > xmax) std::swap(xmin, xmax);
	if (ymin > ymax) std::swap(ymin, ymax);
#else // NEVER // ignore input xmin, etc: rely on clipping
	xmin = -std::numeric_limits<double>::max();
	xmax = std::numeric_limits<double>::max();
	ymin = -std::numeric_limits<double>::max();
	ymax = std::numeric_limits<double>::max();
	T_(trc.dprint("xmin,xmax: ",xmin,", ",xmax);)
	T_(trc.dprint("ymin,ymax: ",ymin,", ",ymax);)
#endif // NEVER // ignore input xmin, etc: rely on clipping

	// ynms might contain wildcards - expand them
	vector<string> ys;
	for (auto& pattern : ynms) {
		for (auto& str : viz.parnames()) {
			if (flaps::wildcard(pattern, str))
				ys.push_back(str);
		}
	}

	this->xname = xnm;
	this->ynames = ys;

	// x and y axes
	xaxis = new mpScaleX(viz.parTitle(xname),mpALIGN_BOTTOM,ticks);
	yaxis = new mpScaleY(viz.parTitle(ynames[0]),mpALIGN_LEFT,ticks);

   wxFont graphFont(11,wxFONTFAMILY_DEFAULT,wxFONTSTYLE_NORMAL,wxFONTWEIGHT_NORMAL);
	xaxis->SetFont(graphFont);
	yaxis->SetFont(graphFont);
	xaxis->SetDrawOutsideMargins(false);
	yaxis->SetDrawOutsideMargins(false);
	AddLayer(xaxis);
	AddLayer(yaxis);

	// create a myFXYVector layer for each curve in viz.curves;
	// If there are multiple ynames each viz.curve will have
	// a myFXYVector for each yname
	// get the bounding box of the data, which may be altered by specs
	bbxmin = std::numeric_limits<double>::max();
	bbxmax = -std::numeric_limits<double>::max();
	bbymin = std::numeric_limits<double>::max();
	bbymax = -std::numeric_limits<double>::max();
	int lid{0};
	for (size_t i=0; i<viz.curves.size(); i++) {
		for (auto yname : ynames) {
			Plotcurve* cp = viz.curves[i];
			int n = cp->nsolns();
			Par* xp{nullptr};
			if (!xname.empty()) {
				xp = cp->findp(xname);
			} else {
				// no xname given: create a parameter of point numbers
				xp = new Par("points","0");
				for (int k=1; k<= n; k++)
					xp->solns.push_back(k);
				// add this new Par to curve
				cp->params.add(xp);
			}
			if (xp == nullptr) {
				message(vastr(xname," is not a  parameter in ",cp->cid()));
				continue;
			}
			Par* yp = cp->findp(yname);
			if (yp == nullptr) {
				message(vastr(yname," is not a parameter in ",cp->cid()));
				continue;
			}

			// gather the x,y data from the requested curves within xmin, etc
			// XXX what about when a curve appears multiple times in the box?
			vector<double> xs;
			vector<double> ys;
			bool inside{false};
			for (int k=0; k<n; k++) {
				double x = xp->solns[k];
				double y = yp->solns[k];
				// inside?
				if (x >= xmin && x <= xmax && y >= ymin && y <= ymax) {
					if (k > 0 && !inside) {
						xs.push_back(xp->solns[k-1]);
						ys.push_back(yp->solns[k-1]);
					}
					inside = true;
					xs.push_back(x);
					ys.push_back(y);
				} else {
					// last point was inside?
					if (inside) {
						xs.push_back(x);
						ys.push_back(y);
					}
					inside = false;
				}					
			}
			// XXX should each layer get a unique int id?
			T_(trc.dprint("curve ",lid," has ",xs.size()," points and is ",inside?"":"not"," visible");)

			if (xs.size() > 0) {
				// mpFXYVector layer is the main drawing class
				myFXYVector* m_vector = new myFXYVector(this, cp, yname);
				// let mpFXVector::Plot draw lines, we'll add the dots in my Plot
				m_vector->SetContinuity(true);

				auto xmm = minmax_element(xs.begin(),xs.end());
				auto ymm = minmax_element(ys.begin(),ys.end());
				cp->xmin = *(xmm.first);
				cp->xmax = *(xmm.second);
				bbxmin = std::min(bbxmin, cp->xmin);
				bbxmax = std::max(bbxmax, cp->xmax);
				cp->ymin = *(ymm.first);
				cp->ymax = *(ymm.second);
				bbymin = std::min(bbymin, cp->ymin);
				bbymax = std::max(bbymax, cp->ymax);
				T_(trc.dprint(cp->cid()," x min/max: ",cp->xmin,", ",cp->xmax,", y min/max: ",cp->ymin,", ",cp->ymax);)
			 
				m_vector->SetData(xs, ys);
				m_vector->ShowName(false);		// do not show the curve name
				if (!toplot.empty() && !contains(toplot, lid))
					m_vector->SetVisible(false);
				else
					m_vector->SetVisible(true);	// default: all curves
				AddLayer(m_vector);
				lid++;
			}
		}
	}
	T_(trc.dprint("this Vizplot has ",plotlayers().size(),", ",visplotlayers().size()," visible");)

	// if x,y limits were given call fit to add these
	Specs& sp = specs();
	if (sp.has_xmin && sp.xmin > xmin) bbxmin = sp.xmin;
	if (sp.has_xmax && sp.xmax < xmax) bbxmax = sp.xmax;
	if (sp.has_ymin && sp.ymin > ymin) bbymin = sp.ymin;
	if (sp.has_ymax && sp.ymax < ymax) bbymax = sp.ymax;
	Fit(bbxmin, bbxmax, bbymin, bbymax);
	// save the limits for the "orig view" menu item
	origxmin = bbxmin;
	origxmax = bbxmax;
	origymin = bbymin;
	origymax = bbymax;

	// set up legends for each curve
	set_legends();

	// set the name of each layer to the legend.id()
	vector<myFXYVector*> layers = plotlayers();
	for (auto li : layers) {
		li->SetName(li->legend.id());
		T_(trc.dprint("set layer \"",li->GetName(),"\", inside? ",li->inside);)
	}

}	// Vizplot::loadxy

void
Vizplot::
add(vector<Plotcurve*>& newcurves) {
// add new curves to an existing Vizplot
	T_(Trace trc(2,"Vizplot::add");)

	// create a myFXYVector layer for each new curve
	// get the bounding box of the data, which may be altered by specs
	bbxmin = std::numeric_limits<double>::max();
	bbxmax = -std::numeric_limits<double>::max();
	bbymin = std::numeric_limits<double>::max();
	bbymax = -std::numeric_limits<double>::max();
	//!! for (size_t i=0; i<viz.curves.size(); i++)
	for (auto cp : newcurves) {
		for (auto yname : ynames) {
			Par* xp = cp->findp(xname);
			if (xp == nullptr) {
				message(vastr(xname," is not a valid parameter"));
				continue;
			}
			Par* yp = cp->findp(yname);
			if (yp == nullptr) {
				message(vastr(yname," is not a valid parameter"));
				continue;
			}

			// mpFXYVector layer is the main drawing class
			myFXYVector* m_vector = new myFXYVector(this, cp, yname);

			// gather the x,y data from the requested curves
			vector<double> xs;
			vector<double> ys;
			int n = cp->nsolns();
			for (int i=0; i<n; i++) {
				xs.push_back(xp->solns[i]);
				ys.push_back(yp->solns[i]);
			}
			if (n > 0) {
				auto xmm = minmax_element(xs.begin(),xs.end());
				auto ymm = minmax_element(ys.begin(),ys.end());
				cp->xmin = *(xmm.first);
				cp->xmax = *(xmm.second);
				bbxmin = std::min(bbxmin, cp->xmin);
				bbxmax = std::max(bbxmax, cp->xmax);
				cp->ymin = *(ymm.first);
				cp->ymax = *(ymm.second);
				bbymin = std::min(bbymin, cp->ymin);
				bbymax = std::max(bbymax, cp->ymax);
				T_(trc.dprint(cp->cid()," x min/max: ",cp->xmin,", ",cp->xmax,", y min/max: ",cp->ymin,", ",cp->ymax);)
			}
			 
			m_vector->SetData(xs, ys);
			m_vector->ShowName(false);		// do not show the curve name
			AddLayer(m_vector);
		}
	}
}

void
Vizplot::
feedAmviz(int ix, int iy) {
// given screen coords (ix,iy) of a picked point, find the closest curve,
// interpolate the eigenvector at the 2 enclosing points, send it to
// amviz via evq. Called on right click
	T_(Trace trc(2, "feedAmviz");)
	VizFrame* frame = lowerpanel->vizframe;
	Amviz& amviz{Amviz::instance()};

	if (!amviz.error.empty())
		message(amviz.error);

	// create a mark to show where the eigenvector is from
	static int visit{1};
	Mark newmark(this,ix,iy,to_string(visit++));
	if (frame->am_spot != nullptr)
		frame->am_spot->push_back(newmark);

	Specs& sp = specs();
	if (sp.threads) {
		if (amviz.kit_future.valid())
			amviz.kit_future.get();
	}

	// check that this curve has the vzid that amviz is using
	// XXX amviz needs a separate Fma for each vzid
	vector<myFXYVector*> layers = plotlayers();
	Plotcurve* cp = layers[newmark.layerno]->plotcurve;
	string vzid = cp->vzid();
	if (!vzid.empty() && vzid != amviz.fma->vzid) {
		message(vastr("this curve has vzid ",vzid," but amviz is using ",
			amviz.fma->vzid));
		return;
	}
	// copy the eigenvector XXX need a guard
	try {
		amviz.ev = cp->params.get_eigenvector(newmark.idx);
	} catch (std::exception& s) {
		message(s.what());
		return;
	}
	// start amviz if not running
	if (sp.threads) {
		if (!amviz.run_future.valid())
			amviz.run_future = std::async(amviz_run);
	} else {
		if (visit == 2)
			amviz_run();
	}
	if (!amviz.error.empty())
		message(amviz.error);
}

wxPoint
Vizplot::
transform(Vizplot* from, wxPoint& p) {
// transform a point on "from" to the equivalent point on this
	T_(Trace trc(2,"transform");)
	T_(trc.dprint("input p: ",p);)
	double px = from->p2x(p.x);
	double py = from->p2y(p.y);
	T_(trc.dprint("px ",px," scale ",from->GetXscl(),", py ",py," scale ",from->GetYscl());)
	wxPoint rval(this->x2p(px), this->y2p(py));
	T_(trc.dprint("output rval: ",rval,", scalex ",GetXscl(),", scaley ",GetYscl());)
	return rval;
}


Mark::
Mark(Vizplot* mp, int ix, int iy, string txt) {
// Mark constructor: find the closest solution to (ix,iy), 
// and the corresponding Plotcurve
	T_(Trace trc(2,"Mark ctor");)
	// closestCurve returns (myFXYVector*, idx)
	auto cpt = mp->closestCurve(ix,iy);
	layerno = get<0>(cpt);
	idx = get<1>(cpt);
	if (layerno == -1)	// not found?
		return;
	T_(trc.dprint("layer number \"",layerno,"\", idx ",idx);)
	text = txt;
}

void
Mark::
drawrect(Vizplot* vp, wxDC& dc, int size) {
// Draw a rectangle at the screen coordinates of layer[idx]; if member
// text is not empty put it above the rectangle

	if (layerno == -1)
		return;

	// layer may not be correct if the X-Y View has changed so
	// get the current x,y parameters from the plotcurve, the x,y
	// values from the soln members (index idx), then convert x,y to
	// screen coordinates.
	// XXX if this works, layer is not necessary? also vizplot?
	int ix, iy;
	myFXYVector* lp = vp->plotlayers()[layerno];
	lp->screencoords(idx, ix, iy);

	int w = size/2;
	// XXX allow different colors for pset_spot, am_spot?
	//dc.SetPen(wxPen(wxColour("ORANGE"), 2));
	dc.SetPen(wxPen(getcolor("yellow"), 5));
	dc.DrawRectangle(wxPoint(ix-w,iy-w),wxSize(2*w,2*w));
	if (!text.empty()) {
		int width, height;
		dc.GetTextExtent(text, &width, &height);
		dc.DrawText(text, ix-width/2, iy-height);
	}
}

tuple<int,int>
Vizplot::
closestCurve(int ix, int iy) {	// XXX rename to closest_soln?
// find the closest curve to (x,y), returns 2 ints and 2 doubles:
//   get<0>(rval): 0b curve number
//   get<1>(rval): 0b soln number
//   get<2>(rval): x (plot units)
//   get<3>(rval): y (plot units)
// vector<complex<double>> ev = viz.curves[jmin]->params.get_eigenvector(imin);
	T_(Trace trc(2,"closestCurve");)
	// bounding box in plot units
	double xmin = p2x(GetMarginLeft());
	double ymax = p2y(GetMarginTop());
	double xmax = p2x(GetScrX()-GetMarginRight());
	double ymin = p2y(GetScrY()-GetMarginBottom());
	T_(trc.dprint("bounding box: (",xmin,", ",ymin,") (",xmax,", ",ymax,")");)
	// lambda for checking in range
	auto oor = [&](double x, double y) {
		if (x >= xmin && x <= xmax && y >= ymin && y <= ymax)
			return false;
		return true;
	};

	int dmin{std::numeric_limits<int>::max()};
	int imin{0};	// point number
	int jmin{0};	// solution number
	T_(myFXYVector* closest_layer{nullptr};)
	vector<myFXYVector*> layers = visplotlayers();	// visible plot layers
	for (size_t j=0; j<layers.size(); j++) {
		string name = layers[j]->GetName().ToStdString();
		T_(trc.dprint("working on layer ",name);)
		const vector<double>& xs = layers[j]->getxs();	// XXX return reference?
		const vector<double>& ys = layers[j]->getys();
		int n = xs.size();
		// Must do this in screen coord so that x and y are scaled the same
		for (int i=0; i<n; i++) {
			double xi = xs[i];
			int jxi = x2p(xi);
			double yi = ys[i];
			int jyi = y2p(yi);
			if (oor(xi,yi)) continue;
			int dx = ix - jxi;
			int dy = iy - jyi;
			int d = dx*dx + dy*dy;
			T_(trc.dprint("curve ",j," soln ",i," d ",d);)
			if (d < dmin) {
				imin = i;
				jmin = j;
				dmin = d;
				T_(closest_layer = layers[j];)
			}
		}
	}
	if (jmin == -1) {
		T_(trc.dprint("no layers found in range");)
		return {jmin,0};
	}
	T_(trc.dprint("returning layer \"",closest_layer->GetName(),"\", imin ",imin);)
	return std::tuple<int,int>(jmin,imin);
}

void
myFXYVector::
screencoords(int idx, int& ix, int& iy) {
// Returns the sceen coordinates (ix,iy) of the idx point
// in my (xs, ys)
	double x = getxs()[idx];
	double y = getys()[idx];
	ix = vizplot->x2p(x);
	iy = vizplot->y2p(y);
}

int
myFXYVector::
curveidx(int idx) {
// Returns the index in the soln member of my plotcurve
// corresponding to my index "idx" XXX should always be the same?
	T_(Trace trc(2,"curveidx");)
	Plotcurve* pc = this->plotcurve;
	Par* xp = pc->findp(this->vizplot->xname);
	Par* yp = pc->findp(this->yname);
	int nx = xp->solns.size();
	int ny = yp->solns.size();
	assert(nx == ny);
	const vector<double>& xs = this->getxs();
	const vector<double>& ys = this->getys();
	double x = xs[idx];
	double y = ys[idx];
	double dmin{std::numeric_limits<double>::max()};
	int rval{0};
	for (int i=0; i<nx; i++) {
		double dx = x - xp->solns[i];
		double dy = y - yp->solns[i];
		double d = dx*dx + dy*dy;
		if (d < dmin) {
			rval = i;
			dmin = d;
		}
	}
	T_(trc.dprint("index ",idx," corresponds to ",rval);)
	return rval;
}

void
Vizplot::
OnLeftDown(wxMouseEvent& event) {
// Handle a left mouse press: just grab the current position
	T_(Trace trc(2,"OnLeftDown");)
	this->startPos = event.GetPosition();
	T_(trc.dprint("start pos: ",this->startPos.x,", ",this->startPos.y);)
	this->leftDown = true;
}
void
Vizplot::
OnLeftUp(wxMouseEvent& event) {
// Handle a left-mouse release: might be a zoom (if position changed) or
// just a parameter-list request
	T_(Trace trc(2,"OnLeftUp");)
	VizFrame* frame = lowerpanel->vizframe;
	this->currentPos = event.GetPosition();
	int dist = abs(this->startPos.x-this->currentPos.x) +
		abs(this->startPos.y-this->currentPos.y);
	T_(trc.dprint("current pos: ",this->currentPos.x,", ",this->currentPos.y);)
	if (dist < 5) this->isDragging = false;
	if (this->isDragging) {
		this->isDragging = false;
		this->leftDown = false;
		// open a new zoom frame based on this
		wxPoint pos(10,20);
		// create a new window using the VizFrame(VizFrame*,Vizplot*) ctor
		VizFrame* zoomframe = new VizFrame(this, frame, "", pos, vizsizes.frame);
		zoomframe->Show(true);
		zoomframe->set_title("zoom");
	} else {
		// this was a left click => show parameter list
		this->leftDown = false;
		Viz& viz{Viz::instance()};
		*frame->pset_spot = Mark(this, this->currentPos.x, this->currentPos.y);

		// write the parameter list to member Vizframe::pset
		// first lines identify the curve, e.g.
		//   vso.pf
		//   2.b
		std::ostringstream os;

		vector<myFXYVector*> layers = plotlayers();
		myFXYVector* lp = layers[frame->pset_spot->layerno];
		int lidx = frame->pset_spot->idx;
		int idx = lp->curveidx(lidx);	// XXX should be the same?
		Plotcurve* cp = lp->plotcurve;

		if (!lp->legend.lead.empty())
			os << "   " << lp->legend.lead << endl;
		os << "   " << lp->legend.num;
		if (!lp->legend.ext.empty())
			os << '.' << lp->legend.ext;
		os << endl;
		//!! if (!lp->legend.yname.empty())
		if (ynames.size() > 1)
			os << "   " << lp->legend.yname << endl;
		os << "---------\n";
		//!! os << "Curve " << cp->legend.num << std::endl;
		for (auto nm : viz.parnames()) {
			Par* pp = cp->findp(nm);
			if (pp != nullptr)
				os << "  " << nm << " = " << pp->solns[idx] << std::endl;
		}
		// set my frame's "pset" member
		*frame->pset = os.str();

		this->showLegend = false;	// tell OnPaint to call paintPset

		// update the left-panel Legend/Param button
		frame->legendButton->SetLabel("Legend");
		frame->leftpanel->Refresh();	// seems like this should all that's needed...
		frame->leftpanel->legendbox->Scroll(0,0);
		frame->Refresh();	// seems like this should all that's needed...
		this->UpdateAll();	// XXX ... I don't understand why this is needed
		//Refresh();
	}
}

void
Vizplot::
OnMiddleDown(wxMouseEvent& event) {
	T_(Trace trc(2,"OnMiddleDown");)
	startPos = event.GetPosition();
	T_(trc.dprint("start pos: ",startPos.x,", ",startPos.y);)
	//!! isDragging = true;
	middleDown = true;
}

void
Vizplot::
OnOrigView(wxCommandEvent& WXUNUSED(event)) {
// return to the original view: all curves with original fit
	VizFrame* frame = lowerpanel->vizframe;
	vector<myFXYVector*> layers = plotlayers();
	for (auto ly : layers)
		ly->SetVisible(true);
	Fit(origxmin, origxmax, origymin, origymax);
	UpdateAll();
	frame->Refresh();
}

void
Vizplot::
OnMouseHelp(wxCommandEvent& WXUNUSED(event)) {
    wxMessageBox(_("Supported Mouse Actions:\n \
        - Left button down + Mark area: Rectangular zoom\n \
        - Wheel: Vertical scroll\n \
        - Wheel + SHIFT: Horizontal scroll\n \
        - Wheel + CTRL: Zoom in/out"),_("Mouse Buttons"),wxOK,this);
}

void
Vizplot::
OnMiddleUp(wxMouseEvent& event) {
	T_(Trace trc(2,"OnMiddleUp");)
	currentPos = event.GetPosition();
	int dist = abs(startPos.x-currentPos.x) + abs(startPos.y-currentPos.y);
	T_(trc.dprint("current pos: ",currentPos.x,", ",currentPos.y,", dist = ",dist);)
	if (dist < 5) isDragging = false;
	if (!isDragging) {
		middleDown = false;
		//!! wxMessageBox(vastr("in OnMiddleUp mousehelp: ",mousehelp,", lowest ",wxID_LOWEST,", highest: ",wxID_HIGHEST));
		PopupMenu(&middleMenu, currentPos.x, currentPos.y);
	}
}

void
Vizplot::
OnRightDown(wxMouseEvent& event) {
	T_(Trace trc(2,"OnRightDown");)
	startPos = event.GetPosition();
	T_(trc.dprint("start pos: ",startPos.x,", ",startPos.y);)
	//!! isDragging = true;
	rightDown = true;
}
void
Vizplot::
OnRightUp(wxMouseEvent& event) {
	T_(Trace trc(2,"OnRightUp");)

	currentPos = event.GetPosition();
	int dist = abs(startPos.x-currentPos.x) + abs(startPos.y-currentPos.y);
	T_(trc.dprint("current pos: ",currentPos.x,", ",currentPos.y);)
	if (dist < 5) isDragging = false;
	if (isDragging) {
		event.Skip();
	} else {
		Amviz& amviz{Amviz::instance()};
		if (!amviz.error.empty()) {
			static int visit{0};
			if (visit++ == 0)
				message(amviz.error);
		} else {
			rightDown = false;
			feedAmviz(currentPos.x,currentPos.y);
		}
	}
}

void
Vizplot::
OnMotion(wxMouseEvent& event) {
	T_(Trace trc(2,"OnMotion");)
	currentPos = event.GetPosition();
	if (leftDown) {
		if (currentPos != startPos) isDragging = true;
		if (isDragging) {
			T_(trc.dprint("current pos: ",currentPos.x,", ",currentPos.y);)
			Refresh();
		} else {
			event.Skip();
		}
	} else {
		// write the current coordinate into member "ibcoords"
		ibcoords.Printf(" x = %f\n y = %f", p2x(event.GetX()), p2y(event.GetY()));
		T_(trc.dprint("ibcoords: ",ibcoords);)
		//!! event.Skip();
	}
}

void
Vizplot::
OnPaint(wxPaintEvent& event) {
	T_(Trace trc(2,"Vizplot::OnPaint ",this);)

	// call the base class OnPaint
	mpWindow::OnPaint(event);

	wxClientDC dc(this);
	dc.SetPen(*wxRED_PEN);
	dc.SetBrush(*wxTRANSPARENT_BRUSH);
	if (isDragging) {
		// draw the temporary rectangle here using wxDC
		// or wxGraphicContext
		T_(trc.dprint("drawing rectangle to current pos: ",currentPos.x,", ",currentPos.y);)
		dc.DrawRectangle(startPos.x,startPos.y,currentPos.x-startPos.x,
				currentPos.y-startPos.y);
	}
	T_(trc.dprint("plot x dim: ",GetScrX());)

	// paint the frame::pset_spot...
	VizFrame* frame = lowerpanel->vizframe;
	frame->pset_spot->drawrect(this, dc, 15);
	// ... then frame::am_spot: amviz ellipses
	if (frame->am_spot != nullptr) {
		for (auto& mp : *frame->am_spot)
			mp.drawrect(this, dc, 15);
	}

	// draw target lines
	if (ynames[0] == "sigma" || ynames[0] == "growth") {
		wxPen pen(*wxRED_PEN);
		pen.SetStyle(wxPENSTYLE_DOT);
		pen.SetWidth(3);
		dc.SetPen(pen);
		wxCoord x0{GetMarginLeft()};
		wxCoord x1{GetScrX()-GetMarginRight()};
		wxCoord y0 = y2p(0.0);
		dc.DrawLine(x0, y0, x1, y0);
	}
	if (xname == "sigma" || xname == "growth") {
		wxPen pen(*wxRED_PEN);
		pen.SetStyle(wxPENSTYLE_DOT);
		pen.SetWidth(3);
		dc.SetPen(pen);
		wxCoord y0{GetMarginTop()};
		wxCoord y1{GetScrY()-GetMarginBottom()};
		wxCoord x0 = x2p(0.0);
		dc.DrawLine(x0, y0, x0, y1);
	}

}	// Vizplot::OnPaint
// end Vizplot implementation }

// CurvesDialog implementation {
CurvesDialog::
CurvesDialog(wxWindow* parent, wxArrayString& curvelist, bool clr) :
		wxDialog(parent,wxID_ANY,"Wildcard doc in manual") {
	wxBoxSizer* vsizer = new wxBoxSizer(wxVERTICAL);
	wxPoint defpos{wxDefaultPosition};
	wxSize defsize{wxDefaultSize};

	clear = clr;

	// text input lines
	vsizer->Add(new wxStaticText(this,wxID_ANY,"Include"),0,wxALL,5);
	includes = new wxTextCtrl(this,wxID_ANY);
	vsizer->Add(includes, 0, wxEXPAND|wxALL, 5);

	vsizer->Add(new wxStaticText(this,wxID_ANY,"Exclude"),0,wxALL,5);
	excludes = new wxTextCtrl(this,wxID_ANY);
	vsizer->Add(excludes, 0, wxEXPAND|wxALL, 5);

	// buttons
	wxStdDialogButtonSizer* buttonSizer = new wxStdDialogButtonSizer();
	wxButton* okButton = new wxButton(this, wxID_OK, "Ok");
	okButton->SetDefault();
	buttonSizer->AddButton(okButton);
	buttonSizer->AddButton(new wxButton(this, wxID_CANCEL, "Cancel"));
	buttonSizer->Realize();
	vsizer->Add(buttonSizer, 0, wxALIGN_CENTER|wxTOP|wxBOTTOM, 10);
	
	// checklist of all curves
	// clear the list?
#ifdef NEVER // doesn't work: lambda never called
	clearbox = new wxCheckBox(this, wxID_ANY, "clear?");
	funBind(wxEVT_CHECKBOX, [this](wxCommandEvent& ev) {
		cerr << "clearbox lambda called\n";
		this->clear = true;
		}, this);
	vsizer->Add(clearbox);
#endif // NEVER // doesn't work: lambda never called

	long int style = wxLB_EXTENDED;
	checklist = new wxListBox(this, wxID_ANY,
		wxDefaultPosition, wxDefaultSize, curvelist, style);
	checklist->SetFirstItem(0);
	vsizer->Add(checklist);

	SetSizerAndFit(vsizer);
}

vector<int>
CurvesDialog::
getChecklist() {
	vector<int> rval;
	wxArrayInt sel;
	checklist->GetSelections(sel);
	for (size_t i=0; i<sel.size(); i++)
		rval.push_back(sel[i]);
	return rval;
}

// end CurvesDialog implementation }

// LimitsDialog implementation {
LimitsDialog::
LimitsDialog(VizFrame* parent) : wxDialog(parent,wxID_ANY,"Set Limits") {
// create and display a dialog to get new (x,y) limits
	vizframe = parent;
	Vizplot* vp = vizframe->currentview;
	wxBoxSizer* vsizer = new wxBoxSizer(wxVERTICAL);
	wxPoint defpos{wxDefaultPosition};
	wxSize defsize{wxDefaultSize};
	// text input lines
	vsizer->Add(new wxStaticText(this,wxID_ANY,"xmin"),0,wxALL,5);
#ifdef NEVER // use GetDesired*min
	xminctrl = new wxTextCtrl(this,wxID_ANY, to_string(vp->bbxmin));
#else // NEVER // use GetDesired*min
	xminctrl = new wxTextCtrl(this,wxID_ANY, to_string(vp->GetDesiredXmin()));
#endif // NEVER // use GetDesired*min
	vsizer->Add(xminctrl, 0, wxEXPAND|wxALL, 5);

	vsizer->Add(new wxStaticText(this,wxID_ANY,"xmax"),0,wxALL,5);
#ifdef NEVER // use GetDesired*min
	xmaxctrl = new wxTextCtrl(this,wxID_ANY, to_string(vp->bbxmax));
#else // NEVER // use GetDesired*min
	xmaxctrl = new wxTextCtrl(this,wxID_ANY, to_string(vp->GetDesiredXmax()));
#endif // NEVER // use GetDesired*min
	vsizer->Add(xmaxctrl, 0, wxEXPAND|wxALL, 5);

	vsizer->Add(new wxStaticText(this,wxID_ANY,"ymin"),0,wxALL,5);
	yminctrl = new wxTextCtrl(this,wxID_ANY, to_string(vp->bbymin));
	vsizer->Add(yminctrl, 0, wxEXPAND|wxALL, 5);

	vsizer->Add(new wxStaticText(this,wxID_ANY,"ymax"),0,wxALL,5);
	ymaxctrl = new wxTextCtrl(this,wxID_ANY, to_string(vp->bbymax));
	vsizer->Add(ymaxctrl, 0, wxEXPAND|wxALL, 5);

	// buttons
	wxStdDialogButtonSizer* buttonSizer = new wxStdDialogButtonSizer();
	wxButton* okButton = new wxButton(this, wxID_OK, "Ok");
	okButton->SetDefault();
	buttonSizer->AddButton(okButton);
	buttonSizer->AddButton(new wxButton(this, wxID_CANCEL, "Cancel"));
	buttonSizer->Realize();
	vsizer->Add(buttonSizer, 0, wxALIGN_CENTER|wxTOP|wxBOTTOM, 10);

	SetSizerAndFit(vsizer);
}
// end LimitsDialog implementation }

// ParamsDialog implementation {
ParamsDialog::
ParamsDialog(VizFrame* parent) :
	wxDialog(parent,wxID_ANY,"Select or type X and Y") {
// create and display a dialog to get new (x,y) limits
	T_(Trace trc(2,"ParamsDialog");)
	vizframe = parent;
	Vizplot* vp = vizframe->currentview;

	// set my xsel and ysel to current values
	xsel = vp->xname;
	ysels = vp->ynames;

	T_(trc.dprint("current x ",xsel,", y ",ysels);)

	wxBoxSizer* vsizer = new wxBoxSizer(wxVERTICAL);

	// X button pops up a list of all parameter with the current x highlighted
	wxButton* xButton = new wxButton(this, uniqueId(), "X");
	//!! xButton->Bind(wxEVT_BUTTON, xselect);
	xButton->Bind(wxEVT_BUTTON, [&](wxCommandEvent& ev) {
		xsel = vizframe->getParChoice(xsel, "x");
		xctrl->SetValue(xsel);
	});
	xButton->SetToolTip(wxT("select X"));
	vsizer->Add(xButton, 0);
	// ditto y
	wxButton* yButton = new wxButton(this, uniqueId(), "Y");
	//!! yButton->Bind(wxEVT_BUTTON, yselect);
	yButton->Bind(wxEVT_BUTTON, [&](wxCommandEvent& ev) {
		ysels = vizframe->getYParChoice(ysels);
		yctrl->SetValue(vastr(ysels));		// XXX set all ys?
	});
	yButton->SetToolTip(wxT("select Y"));
	vsizer->Add(yButton, 0);

	// text input lines allow the user to type the names of x,y
	vsizer->Add(new wxStaticText(this,wxID_ANY,"x"),0,wxALL,5);
	//!! xctrl = new wxTextCtrl(this,wxID_ANY, vp->xname);
	xctrl = new wxTextCtrl(this,wxID_ANY);
	xctrl->SetValue(xsel);
	vsizer->Add(xctrl, 0, wxEXPAND|wxALL, 5);

	vsizer->Add(new wxStaticText(this,wxID_ANY,"y"),0,wxALL,5);
	//!! yctrl = new wxTextCtrl(this,wxID_ANY, vp->yname);
	yctrl = new wxTextCtrl(this,wxID_ANY);
	yctrl->SetValue(ysels[0]);		// XXX set all ys?
	vsizer->Add(yctrl, 0, wxEXPAND|wxALL, 5);

	// ok, cancel buttons XXX necessary?
	wxStdDialogButtonSizer* buttonSizer = new wxStdDialogButtonSizer();
	wxButton* okButton = new wxButton(this, wxID_OK, "Ok");
	okButton->SetDefault();
	buttonSizer->AddButton(okButton);
	buttonSizer->AddButton(new wxButton(this, wxID_CANCEL, "Cancel"));
	buttonSizer->Realize();
	vsizer->Add(buttonSizer, 0, wxALIGN_CENTER|wxTOP|wxBOTTOM, 10);

	SetSizerAndFit(vsizer);
}
// end ParamsDialog implementation }

bool
VizApp::
OnInit() {
	T_(Trace trc(1,"VizApp::OnInit");)

	try {
		// Create the main frame window...
		wxPoint pos(10,20);
		VizFrame* vizframe = new VizFrame("", pos, vizsizes.frame);
		vizframe->Show(true);
		vizframe->set_title();
	} catch (std::exception& e) {
		flaps::error(e.what());
	}
	return true;
}

// IMPLEMENT_APP_NO_MAIN(VizApp):
	VizApp& wxGetApp() {
	 	T_(Trace trc(2,"wxGetApp");)
	 	return *static_cast<VizApp*>(wxApp::GetInstance());
	}
	wxAppConsole *wxCreateApp() {
		T_(Trace trc(2,"wxCreateApp");)
		wxAppConsole::CheckBuildOptions(WX_BUILD_OPTIONS_SIGNATURE, "your program");
		return new VizApp;
	}
	wxAppInitializer wxTheAppInitializer((wxAppInitializerFunction) wxCreateApp);

// VizFrame implementation {
//!! #include "viz.xpm"
int VizFrame::number_of_frames{0};

VizFrame::VizFrame(const wxString& title,
	const wxPoint& pos, const wxSize& size, long style) :
		wxFrame(nullptr, wxID_ANY, title, pos, size, style) {
// VizFrame ctor
	T_(Trace trc(2,"VizFrame initial ctor");)
	frameno = number_of_frames++;

	// create a top-level vertical Sizer
	topsizer = new wxBoxSizer(wxVERTICAL);
	SetSizer(topsizer);

	// add menubar, toolbar + controls
	createControls();

	Specs& sp = specs();
	statusBar = new wxStatusBar(this);
	string stat("no file");
	if (!sp.aids.empty())
		stat = sp.aids[0];
	else if (!sp.plotfiles.empty())
		stat = sp.plotfiles[0];
	statusBar->SetStatusText(stat);
	SetStatusBar(statusBar);

	// legend panel and plot go in a LowerPanel which has a horizonal sizer
	lowerpanel = new LowerPanel(this);
	lowersizer = new wxBoxSizer(wxHORIZONTAL);
	lowerpanel->SetSizer(lowersizer);

	// legend panel on left: add to lowersizer before m_plot
	leftpanel = new LeftPanel(lowerpanel,wxID_ANY,wxPoint(10,20),vizsizes.panel);
	lowersizer->Add(leftpanel, 0);
	leftpanel->Show(true);

	// Create a mathplot window for each view ...
	if (sp.views.empty()) {
		// no views - just create a blank Vizplot if the user did not specify x, y
		Vizplot* m_plot = new Vizplot(lowerpanel, -1, wxDefaultPosition, vizsizes.plot);
		// ... and load the data if x & y parameters have been specified
		if (!sp.xname.empty() && !sp.ynames.empty()) {
			m_plot->loadxy(sp.xname, sp.ynames);
			// add it to the plotz vector of views
			addplot(m_plot);
			setcurrent(0);
			m_plot->Show(true);
		}
	} else {
		for (auto& vp : sp.views) {
			Vizplot* m_plot = new Vizplot(lowerpanel, -1, wxDefaultPosition, vizsizes.plot);
			// ... and load the data if x & y parameter have been specified
			vector<string> toks = string2tok(vp, ":");
			sp.xname = toks[0];
			sp.ynames.clear();
			sp.ynames.insert(sp.ynames.begin(), toks.begin()+1, toks.end());
			m_plot->loadxy(sp.xname, sp.ynames);
			addplot(m_plot);
			m_plot->Show(false);
		}
		setcurrent(0);
	}
	// add lowerpanel to the topsizer
	topsizer->Add(lowerpanel, 1, wxEXPAND);
	topsizer->Layout();

	Show(true);
	SetPosition(pos);		// XXX ??

	// allocate a Mark & vector<Mark> for pset_spot & am_spot
	am_spot = new vector<Mark>;
	pset_spot = new Mark;
	pset = new string;

	T_(trc.dprint("plotz: ",plotz);)
	T_(trc.dprint("frame size: ",GetSize(),", pos ",GetPosition());)
	T_(trc.dprint("topsizer: ",*topsizer);)

}	// VizFrame initial ctor

// the zoom constructor creates a new VizFrame with a new plot size
VizFrame::VizFrame(Vizplot* vizplot, VizFrame* p, const wxString& title,
	const wxPoint& pos, const wxSize& size, long style) :
		wxFrame(nullptr, wxID_ANY, title, pos, size, style) {
// VizFrame zoom ctor: take toplot,pset_spot,am_spot from the input Vizplot
	T_(Trace trc(2,"VizFrame zoom ctor");)
	frameno = number_of_frames++;

	// save the frame we were created from so we can propagate Marks
	parent = p;

	// Mark pointers
	pset_spot = p->pset_spot;
	pset = p->pset;
	am_spot = p->am_spot;

	// take stuff from vizplot
	toplot = vizplot->vizframe->toplot;
	T_(trc.dprint("toplot from vizplot: ",toplot);)

	// create a top-level vertical Sizer
	topsizer = new wxBoxSizer(wxVERTICAL);
	SetSizer(topsizer);

	// add menubar, toolbar + controls
	createControls();

	Specs& sp = specs();
	statusBar = new wxStatusBar(this);
	string stat("no file");
	if (!sp.aids.empty())
		stat = sp.aids[0];
	else if (!sp.plotfiles.empty())
		stat = sp.plotfiles[0];
	statusBar->SetStatusText(stat);
	SetStatusBar(statusBar);

	// legend panel and plot go in a LowerPanel which has a horizonal sizer
	lowerpanel = new LowerPanel(this);
	lowersizer = new wxBoxSizer(wxHORIZONTAL);
	lowerpanel->SetSizer(lowersizer);

	// legend panel on left: add to lowersizer before m_plot
	leftpanel = new LeftPanel(lowerpanel,wxID_ANY,wxPoint(10,20),vizsizes.panel);
	lowersizer->Add(leftpanel, 0);
	leftpanel->Show(true);

	// create a Vizplot: plot will be clipped to it's size
	Vizplot* m_plot = new Vizplot(lowerpanel, -1, wxDefaultPosition, vizsizes.plot);
	// ... and load the data if x & y parameters have been specified
	// XXX xmin, etc ignored in loadxy: rely on clipping
	double xmin = vizplot->p2x(vizplot->startPos.x);
	double xmax = vizplot->p2x(vizplot->currentPos.x);
	double ymin = vizplot->p2y(vizplot->startPos.y);
	double ymax = vizplot->p2y(vizplot->currentPos.y);
	m_plot->loadxy(vizplot->xname, vizplot->ynames, this->toplot, xmin, xmax, ymin, ymax);
	// add it to the plotz vector of views
	addplot(m_plot);
	setcurrent(0);
	// zoom it: changes the size (hence the clipping region) of m_plot
	wxPoint start = m_plot->transform(vizplot, vizplot->startPos);
	wxPoint end = m_plot->transform(vizplot, vizplot->currentPos);
	m_plot->ZoomRect(start, end);
	m_plot->Show(true);

	// add lowerpanel to the topsizer
	topsizer->Add(lowerpanel, 1, wxEXPAND);
	topsizer->Layout();

	Show(true);
	SetPosition(pos);		// XXX ??

	T_(trc.dprint("plotz: ",plotz);)
	T_(trc.dprint("frame size: ",GetSize(),", pos ",GetPosition());)
	T_(trc.dprint("topsizer: ",*topsizer);)

}	// VizFrame zoom ctor

void
VizFrame::
createControls() {
// create all controls for this frame:
//   menubar
//   toolbar
	T_(Trace trc(2,"createControls");)
	wxPoint pos{wxDefaultPosition};
	wxSize sz{wxDefaultSize};

	// accumulate accelerator key in accels
	vector<wxAcceleratorEntry> accels;

	// Respond to an EVT_CLOSE
	Bind(wxEVT_CLOSE_WINDOW, &VizFrame::OnClose, this);

    // Make the "File" menu
    wxMenu *fileMenu = new wxMenu;
    fileMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnFileOpen, this), "&Open...");
    fileMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnFileAdd, this), "&Add...");
    fileMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnNewWindow, this),"New Window");
#ifdef NEVER // disable: use screenshot
    fileMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnFilePrint, this),"&Print\tALT-P");
#endif // NEVER // disable: use screenshot
    fileMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnScreenshot, this), "Screenshot");
    fileMenu->AppendSeparator();
    fileMenu->Append(funBind(wxEVT_MENU, [this](wxCommandEvent& ev) {
	 	this->Close(true); }, this), "Exit");

	 // make the "Edit" menu
#ifdef NEVER // needs work
    wxMenu *editMenu = new wxMenu;
    editMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnMenuEditPref, this), "&Preferences");
#endif // NEVER // needs work

    // Make the "View" menu
    wxMenu *viewMenu = new wxMenu;
    viewMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnToggleGrid, this),
			"Toggle grid/ticks");
    viewMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnAlignXAxis, this),
	 		"Switch X axis alignment");
    viewMenu->Append(myBind(wxEVT_MENU, &VizFrame::OnAlignYAxis, this),
	 			"Switch Y axis alignment");

    // Make the "Help" menu
    wxMenu *helpMenu = new wxMenu;
	 // about
    helpMenu->Append(funBind(wxEVT_MENU, [](wxCommandEvent& ev) {
	 	wxMessageBox(aboutmsg(), "viz"); }, this), "About");
	 // Help
    helpMenu->Append(funBind(wxEVT_MENU, [](wxCommandEvent& ev) {
	 	wxMessageBox(helpmsg(), "Help"); }, this), "Help");

	 // Add these menus to the menu bar
    menubar = new wxMenuBar;
    menubar->Append(fileMenu, wxT("&File"));
    menubar->Append(viewMenu, wxT("&View"));
#ifdef NEVER // needs work
    menubar->Append(editMenu, wxT("&Edit"));
#endif // NEVER // needs work
    menubar->Append(helpMenu, wxT("&Help"));

	SetMenuBar(menubar);

	T_(trc.dprint("menubar size: ",menubar->GetSize());)
	//------------------------------------------------------------------
	// toolbar on the top
	//------------------------------------------------------------------
	toolbar = new wxPanel(this, wxID_ANY, wxDefaultPosition, vizsizes.toolbar);
	wxBoxSizer* toolsizer = new wxBoxSizer(wxHORIZONTAL);
	toolbar->SetSizer(toolsizer);

	int prop{1};
	int buttonstyle{wxEXPAND|wxTOP|wxBOTTOM};
	int bmarg{10};

	// curves select button
	toolsizer->AddStretchSpacer();
	wxButton* curvesButton = new wxButton(toolbar,
		myBind(wxEVT_BUTTON, &VizFrame::OnCurvesSelect, this),"Curves");
	curvesButton->SetToolTip("wild cards & implicit lists");
	toolsizer->Add(curvesButton, prop, buttonstyle, bmarg);
	toolsizer->AddStretchSpacer();

	// parameters select button
	wxButton* paramsButton = new wxButton(toolbar,
		myBind(wxEVT_BUTTON, &VizFrame::OnParamsSelect, this), "Parameters");
	paramsButton->SetToolTip(wxT("choose x,y params"));
	toolsizer->Add(paramsButton, prop, buttonstyle, bmarg);
	toolsizer->AddStretchSpacer();

	// X-Y Views button: press to choose previous x-y combos...
	wxButton* viewButton = new wxButton(toolbar,
		myBind(wxEVT_BUTTON, &VizFrame::OnViewSelect, this), wxT("X-Y Views"));
	viewButton->SetToolTip(wxT("select an x-y combo"));
	toolsizer->Add(viewButton, prop, buttonstyle, bmarg);
	toolsizer->AddStretchSpacer();

	// Color button: press to choose how to color curves
	int cmdid = myBind(wxEVT_BUTTON, &VizFrame::OnColor, this);
	wxButton* colorButton = new wxButton(toolbar, cmdid, wxT("Curve colors"));
	colorButton->SetToolTip(wxT("color curves by (C)"));
	toolsizer->Add(colorButton, prop, buttonstyle, bmarg);
	toolsizer->AddStretchSpacer();
	accels.push_back({0,'C',cmdid});

	// Limits button: press to bring up limits dialog
	cmdid = myBind(wxEVT_BUTTON, &VizFrame::OnLimits, this);
	limitsButton = new wxButton(toolbar, cmdid, "Limits");
	limitsButton->SetToolTip(wxT("set parameter limits (L)"));
	toolsizer->Add(limitsButton, prop, buttonstyle, bmarg);
	toolsizer->AddStretchSpacer();
	accels.push_back({0,'L',cmdid});

	// Legend button: press to toggle legend/pset
	cmdid = myBind(wxEVT_BUTTON, &VizFrame::OnLegend, this);
	legendButton = new wxButton(toolbar, cmdid, "Param list");
	legendButton->SetToolTip(wxT("toggle legend/parameter values (P)"));
	toolsizer->Add(legendButton, prop, buttonstyle, bmarg);
	toolsizer->AddStretchSpacer();
	accels.push_back({0,'P',cmdid});

	// quit button lambda
	cmdid = funBind(wxEVT_BUTTON, [this](wxCommandEvent& ev) {
		T_(Trace trc(2,"quit button lambda");)
// #ifdef THREADS
		Amviz& amviz{Amviz::instance()};
		amviz.quit = true;
// #endif // THREADS
		this->Close(true);}, this);
	wxButton* quitButton = new wxButton(toolbar, cmdid, "Quit");
	quitButton->SetToolTip("Exit (q)");
	toolsizer->Add(quitButton, prop, buttonstyle, bmarg);
	toolsizer->AddStretchSpacer();
	accels.push_back({0,'Q',cmdid});

	// toolbar->Realize();
	// this->SetToolBar(toolbar);

	T_(trc.dprint("toolbar size: ",toolbar->GetSize());)

	// add the toolbar to the top sizer
	//!! topsizer->Add(toolbar, 0, wxLEFT);
	topsizer->Add(toolbar, 0, wxEXPAND, 10);

	// add the accelerators
	SetAcceleratorTable(wxAcceleratorTable(accels.size(), accels.data()));

} // VizFrame::createControls()

void
VizFrame::
OnClose(wxCloseEvent& ev) {
	T_(Trace trc(2,"VizFrame::OnClose");)
	Destroy();
}


bool
VizFrame::
ShowToolTips() { return true; }


// File|Open... command
void
VizFrame::
OnFileOpen( wxCommandEvent& event ) {
// open a dialog, get the user's plotfile choice
	T_(Trace trc(2,"OnFileOpen");)
	Specs& sp = specs();

	FilterDialog filterdlg(this);
	if (filterdlg.ShowModal() == wxID_CANCEL)
		return;
	wxString filters = filterdlg.getfilters();
	wxFileDialog dlg(this, "Open plot file(s)", "", "", filters, wxFD_MULTIPLE);
		//"Flaps plot files (*.pf)|*.pf|Ascii plot files (*.apf)|*.apf",wxFD_MULTIPLE);
	//!! myFileDialogHook dlghook;
	//!! dlg.SetCustomizeHook(dlghook);
	if (dlg.ShowModal() == wxID_CANCEL)
		return;
	wxArrayString fn;
	dlg.GetFilenames(fn);
	vector<string> filenames;
	// convert to std::string, watch for uf files
	for (size_t i=0; i<fn.size(); i++) {
		string fi = fn[i].ToStdString();
		if (rsubstr(fi, 3) == ".uf") {
			sp.ufname = fi;
			T_(trc.dprint("got uf filename: ",fi);)
			Amviz& amviz{Amviz::instance()};
			if (!amviz.error.empty()) {
				amviz.error = "";
			}
			// ... then load amviz data in a separate thread
			// if the user right-clicks feedAmviz will be called which
			// will call amviz_run
			if (sp.threads)
				amviz.kit_future = async(amviz_kit);
			else
				amviz_kit();
		} else
			filenames.push_back(fi);
	}
	T_(trc.dprint("picked files: ",filenames);)

    if (filenames.empty())
	 	return;

	// grab some stuff from the current Vizplot:
	Vizplot* current = plot();
	string xnm = current->xname;
	vector<string> ynms = current->ynames;
	wxSize plotsize = current->GetSize();
	T_(trc.dprint("current size ",plotsize.x,", ",plotsize.y);)

	// load the new file(s)
	// XXX put all this stuff in a Viz member function so it doesn't need to be public?
	Viz& viz{Viz::instance()};
	vector<string> aids;
	
	// XXX would it be better to just re-allocate viz?
	viz.clear();

	viz.loadcurves(aids, filenames);

	// delete all plots in plotz
	clear();

	// create a new one
	Vizplot* vp = new Vizplot(lowerpanel, -1, wxDefaultPosition, plotsize);

	// load x,y data
	vp->loadxy(xnm, ynms);

	// add it to the list of views...
	int cv = addplot(vp);
	// ... and make it the current view
	setcurrent(cv);

	// update the frame title & status bar
	set_title();

	Layout();
	Refresh();
}

// File|Add... command
void
VizFrame::
OnFileAdd( wxCommandEvent& event) {
// add curves from another plotfile; plot them all together
// with the same parameters used for the current plot
	T_(Trace trc(2,"OnFileAdd");)
	Specs& sp = specs();

	FilterDialog filterdlg(this);
	if (filterdlg.ShowModal() == wxID_CANCEL)
		return;
	wxString filters = filterdlg.getfilters();
	T_(trc.dprint("got filters: ",filters);)
	wxFileDialog dlg(this, "Add plot file(s)", "", "", filters, wxFD_MULTIPLE);
	if (dlg.ShowModal() == wxID_CANCEL)
		return;
	wxArrayString fn;
	dlg.GetFilenames(fn);
	vector<string> filenames;
	// convert to std::string, watch for uf files
	for (size_t i=0; i<fn.size(); i++) {
		string fi = fn[i].ToStdString();
		if (rsubstr(fi, 3) == ".uf") {
			sp.ufname = fi;
			T_(trc.dprint("got uf filename: ",fi);)
			Amviz& amviz{Amviz::instance()};
			if (!amviz.error.empty()) {
				amviz.error = "";
			}
			// ... then load amviz data in a separate thread
			// if the user right-clicks feedAmviz will be called which
			// will call amviz_run in a thread
			if (sp.threads)
				amviz.kit_future = async(amviz_kit);
			else
				amviz_kit();
		} else
			filenames.push_back(fi);
	}
	T_(trc.dprint("picked files: ",filenames);)

    if (filenames.empty())
	 	return;

	Viz& viz{Viz::instance()};

	vector<string> aids;
	
	// load the new aids/plotfiles...
	vector<Plotcurve*> newcurves = viz.loadcurves(aids, filenames);
	// ... then add the new curves to all views (plotz)
	for (auto vp : plotz) {
		vp->add(newcurves);
		vp->set_legends();
	}

	// re-set colorby
	viz.colorby = -1;

	// new frame title
	set_title();

	Refresh();
}

// File|New Window
void
newframe() {
	wxPoint pos(15,30);
	VizFrame* frame = new VizFrame("New Window", pos, vizsizes.frame);
	frame->Show(true);
	frame->set_title();
}

void
VizFrame::
OnNewWindow( wxCommandEvent& event) {
// open a new Vizplot window
// XXX should this new window inherit toplot?
	T_(Trace trc(2,"OnNewWindow");)

	// Create a new (blank) VizFrame in a thread
	std::thread new_thread(newframe);
	new_thread.detach();
}

#ifdef NEVER // disable: use screenshot
// File|Print... command
void
VizFrame::
OnFilePrint( wxCommandEvent& WXUNUSED(event) ) {
}

void
VizFrame::
printPS() {
}
#endif // NEVER // disable: use screenshot


#ifdef NEVER // needs work
void
VizFrame::
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

vector<int>
implist(const string& str) {
// given a string like 1:5 or 1:2:5 convert to a vector of ints
	T_(Trace trc(2,"implist");)
	vector<int> rval;
	vector<string> toks = string2tok(str, ":");
	T_(trc.dprint("toks: ",toks);)
	if (toks.size() < 2)
		return rval;
	int incr{1};
	int first = stoi(toks[0]);
	int last;
	if (toks.size() == 2)
		last = stoi(toks[1]);
	else if (toks.size() == 3) {
		last = stoi(toks[2]);
		incr = stoi(toks[1]);
	} else {
		flaps::error("invalid range: ",str);
		return rval;
	}
	T_(trc.dprint("first ",first,", incr ",incr,", last ",last);)

	for (int i=first; i<=last; i+=incr)
		rval.push_back(i);
	return rval;
}

vector<int>
intList(const string& list) {
// parse a string of 1b ints like 1:2:5,7,8 into a vector
// of 0b ints: {0,2,4,6,7}
	T_(Trace trc(2,"intList");)
	vector<int> rval;

	// split the list into comma-separated tokens
	vector<string> toks = string2tok(list, ",");
	T_(trc.dprint("comma-separated tokens: ",toks);)

	// check each one for colons (1 or 2)
	for (auto& tok : toks) {
		vector<int> bits = implist(tok);	
		if (bits.empty()) {
			rval.push_back(stoi(tok));
		} else {
			for (auto& bit : bits)
				rval.push_back(bit);
		}
	}
	// convert from 1b to 0b
	for (size_t i=0; i<rval.size(); i++)
		rval[i] -= 1;
	return rval;
}

void
VizFrame::
OnCurvesSelect(wxCommandEvent& event) {
// pop up a dialog for selecting which curves to plot
	T_(Trace trc(2,"OnCurvesSelect");)
	Viz& viz{Viz::instance()};
	Vizplot* m_plot = plot();

	if (viz.curves.empty()) {
		message("no curves defined yet");
		return;
	}
	// make up a list of *all* plot layers
	vector<myFXYVector*> layers = m_plot->plotlayers();
	int ncurves = layers.size();

#ifdef NEVER // move
	// make up a list of all curve ids
	vector<string> curve_ids;
	for (auto& cp : layers)
		curve_ids.push_back(cp->legend.sid());
#endif // NEVER // move

	// popup a question
	bool clear{false};
	int answer = wxMessageBox("Clear the existing curves list?",
                            "Question", wxYES_NO | wxICON_QUESTION, this);
	if (answer == wxYES)
		clear = true;

	// create a dialog with these curve ids
	wxArrayString curve_names;
#ifdef NEVER // move
	for (auto& ci : curve_ids)
		curve_names.push_back(ci);
#else // NEVER // move
	for (auto& cp : layers)
		curve_names.push_back(cp->legend.sid());
#endif // NEVER // move
	CurvesDialog cdlg(this, curve_names, clear);

	// mark visible curves (clear=no)?
	if (!cdlg.clear) {
		for (int i=0; i<ncurves; i++)
			if (layers[i]->IsVisible())
				cdlg.checklist->SetSelection(i);
	}

	// pop up the dialog...
	if (cdlg.ShowModal() == wxID_CANCEL)
		return;
	
	// get the user's answers
	string inclstr = cdlg.GetInclText().ToStdString();
	vector<string> includes = string2tok(inclstr, " ,");

	string exclstr = cdlg.GetExclText().ToStdString();
	vector<string> excludes = string2tok(exclstr, " ,");

	vector<int> sel = cdlg.getChecklist();
	T_(trc.dprint("got sel: ",sel);)
	T_(trc.dprint("got includes \"",includes,"\", excludes \"",excludes,"\"");)

	// set up this->toplot according to the 3 answers
	// includes and excludes may be wildcard - test with fnmatch()
	vector<int> newtoplot;
	for (int i=0; i<ncurves; i++) {
		string sid = layers[i]->legend.sid();
		bool bi{false};
		if (!includes.empty()) {
			for (auto& incl : includes)
				if (flaps::wildcard(incl, sid)) bi = true;
		}
		bool be{false};
		if (!excludes.empty()) {
			for (auto& excl : excludes)
				if (flaps::wildcard(excl, sid)) be = true;
		}
		T_(trc.dprint("result: include ",bi,", exclude ",be);)
		if ((contains(sel, i) || bi) && !be)
			newtoplot.push_back(i);
	}
	if (newtoplot.empty()) {
		string msg("no curves match these specs");
		if (cdlg.clear)
			msg += "\nmaybe try without clearing the existing curves?";
		message(msg);
		return;
	} else {
		T_(trc.dprint("new toplot: ",newtoplot);)
		this->toplot = newtoplot;
	}

	// set toplot layers visible
	for (int i=0; i<ncurves; i++)
		layers[i]->SetVisible(contains(newtoplot, i));

	leftpanel->legendbox->Scroll(0,0);
	m_plot->UpdateAll();
	Refresh();
	T_(trc.dprint("after choosing curves ",topsizer->GetItemCount()," in topsizer");)
}

string
VizFrame::
getParChoice(const string& current, const string& axis) {
// popup a dialog with all parameter names to get the choice of one
	T_(Trace trc(2, "getParChoice");)
	string rval;

	Viz& viz{Viz::instance()};
	wxArrayString par_names;
	const vector<string>& parnames = viz.parnames();
	for (const auto& pc : parnames)
		par_names.Add(pc);

	// create the dialog with all parameter names
	wxSingleChoiceDialog pardlg(this,vastr(axis," Parameter"),
		vastr("Choose ",axis), par_names);

	// checkmark the parameter which is the current name
	T_(trc.dprint("current ",axis,": ",current);)
	if (!current.empty()) {
		for (size_t i=0; i<parnames.size(); i++) {
			if (parnames[i] == current) {
				pardlg.SetSelection(i);
				break;
			}
		}
	}

	// show the dialog, get selection
	int stat = pardlg.ShowModal();
	T_(trc.dprint("dialog returned ",stat);)
	if (stat == wxID_CANCEL)
		return rval;
	int sel = pardlg.GetSelection();

	rval = parnames[sel];
	T_(trc.dprint("returning selection ",sel,": ",rval);)
	return rval;
}

vector<string>
VizFrame::
getYParChoice(const vector<string>& current) {
// popup a dialog with all parameter names to get the choice of one or more Y
	T_(Trace trc(2, "getYParChoice");)
	vector<string> rval;

	Viz& viz{Viz::instance()};
	wxArrayString par_names;
	const vector<string>& parnames = viz.parnames();
	for (const auto& pc : parnames)
		par_names.Add(pc);

	// create the dialog with all parameter names
	wxMultiChoiceDialog pardlg(this,"Y Parameters",
		"Choose Y", par_names);

	wxArrayInt sels;

#ifdef NEVER // don't check current
	// checkmark the current parameters
	T_(trc.dprint("current Y: ",current);)
	for (auto yname : current) {
		int i = find(parnames.begin(), parnames.end(), yname) - parnames.begin();
		sels.push_back(i);
	}
	pardlg.SetSelections(sels);
#endif // NEVER // don't check current

	// show the dialog, get selections
	int stat = pardlg.ShowModal();
	T_(trc.dprint("dialog returned ",stat);)
	if (stat == wxID_CANCEL)
		return rval;
	sels = pardlg.GetSelections();

	for (auto i : sels)
		rval.push_back(parnames[i]);
	T_(trc.dprint("returning selections: ",rval);)
	return rval;
}

void
VizFrame::
OnParamsSelect(wxCommandEvent& event) {
// pop up a dialog for selecting which parameters are x, y
// Allow for the situation where
// - viz was started with no x,y parameters
// - viz was started with no x,y parameters or plotfiles/aids
// In either case the plot window (Vizplot) is up but the Vizplot
// is incomplete so we just need to complete it.
	T_(Trace trc(2,"OnParamsSelect");)

	// pop up...
	ParamsDialog dlg(this);
	if (dlg.ShowModal() == wxID_CANCEL)
		return;
	// ... and get the new x & y names...
	wxString xstr = dlg.GetXText();
	vector<string>  ynms = dlg.GetYText();
	string xnm = xstr.ToStdString();
	// XXX new multiple choice
	//vector<string> ynms{ystr.ToStdString()};

	// the current Vizplot:
	Vizplot* current = plot();

	// load the current plot if incomplete...
	if (!current->is_complete()) {
		current->loadxy(xnm, ynms);
		return;
	}
	// ... or create a new one
	wxSize plotsize = current->GetSize();
	T_(trc.dprint("current size ",plotsize.x,", ",plotsize.y);)
	Vizplot* vp = new Vizplot(lowerpanel, -1, wxDefaultPosition, plotsize);
	// load data
	vp->loadxy(xnm, ynms, this->toplot);

	// add it to the list of views...
	int cv = addplot(vp);
	// ... and make it the current view
	setcurrent(cv);


	// set the bounding box for parameters that are the same
	string xcurr = current->xname;
	// XXX need multiple choice
	string ycurr = current->ynames[0];
	if (xcurr == xnm) {
#ifdef NEVER // use GetDesired*min
		vp->bbxmin = current->bbxmin;
		vp->bbxmax = current->bbxmax;
#else // NEVER // use GetDesired*min
		vp->bbxmin = current->GetDesiredXmin();
		vp->bbxmax = current->GetDesiredXmax();
#endif // NEVER // use GetDesired*min
	}
	if (xcurr == ynms[0]) {
#ifdef NEVER // use GetDesired*min
		vp->bbymin = current->bbxmin;
		vp->bbymax = current->bbxmax;
#else // NEVER // use GetDesired*min
		vp->bbymin = current->GetDesiredXmin();
		vp->bbymax = current->GetDesiredXmax();
#endif // NEVER // use GetDesired*min
	}
		
	if (ycurr == ynms[0]) {
#ifdef NEVER // use GetDesired*min
		vp->bbymin = current->bbymin;
		vp->bbymax = current->bbymax;
#else // NEVER // use GetDesired*min
		vp->bbymin = current->GetDesiredYmin();
		vp->bbymax = current->GetDesiredYmax();
#endif // NEVER // use GetDesired*min
	}
	if (ycurr == xnm) {
#ifdef NEVER // use GetDesired*min
		vp->bbxmin = current->bbymin;
		vp->bbxmax = current->bbymax;
#else // NEVER // use GetDesired*min
		vp->bbxmin = current->GetDesiredYmin();
		vp->bbxmax = current->GetDesiredYmax();
#endif // NEVER // use GetDesired*min
	}	
	vp->Fit(vp->bbxmin,vp->bbxmax, vp->bbymin, vp->bbymax);
}

int
VizFrame::
addplot(Vizplot* p) {
// Add p to the list of Vizplots (VizFrame::plotz) if it is not already there
	T_(Trace trc(2,"addplot");)
	bool found{false};
	int rval{0};
	// add2vector(p, plotz);
	for (size_t j=0; j<plotz.size(); j++) {
		Vizplot* pp = plotz[j];
		if (pp->xname != p->xname)
			continue;
		if (pp->ynames.size() != p->ynames.size())
			continue;
		size_t i;
		for (i=0; i<pp->ynames.size(); i++) {
			if (pp->ynames[i] != p->ynames[i])
				break;
		}
		if (i == pp->ynames.size()) {
			found = true;
			rval = j;
			break;
		}
	}

	// if not found add it to plotz
	if (!found) {
		plotz.push_back(p);
		rval = plotz.size() - 1;
		T_(trc.dprint("added ",p,", now have ",plotz.size(),", lowerpanel: ",*lowersizer);)
	} else {
		T_(trc.dprint("x ",p->xname," y ",p->ynames[0]," already in plotz[",rval,"]");)
	}
	return rval;
}

void
VizFrame::
OnViewSelect(wxCommandEvent& event) {
// pop up a dialog for selecting which set of parameters to plot
	T_(Trace trc(2,"OnViewSelect");)
	ostringstream os;

	// create a list of current views XXX maybe make plotz a map?
	wxArrayString xy;
	for (auto& pi : plotz) {
		os.str("");
		os << pi->xname;
		for (auto& yi : pi->ynames)
			os << "-" << yi;
		xy.Add(os.str());
	}
	wxSingleChoiceDialog viewdlg(this, "X-Y View", "choose", xy);
	// mark the current pick
	int current = currentidx();
	T_(trc.dprint("current idx ",current);)
	if (current >= 0)
		viewdlg.SetSelection(current);
	// show the dialog, get selection
	int stat = viewdlg.ShowModal();
	if (stat == wxID_CANCEL)
		return;
	int sel = viewdlg.GetSelection();
	T_(trc.dprint("got sel ",sel);)
	Vizplot* vp = setcurrent(sel);
	vp->UpdateAll();
	Refresh();
}

void
VizFrame::
OnColor(wxCommandEvent& event) {
// pop up a dialog for selecting how to color curves
	T_(Trace trc(2,"OnColor");)
	Viz& viz{Viz::instance()};

	// 3 possibilities:
	wxArrayString choices;
	vector<string> colorbys{"ordinal curve", "curve number", "plotfile/aid", "extension", "Y"};
	for (auto& ci : colorbys)
		choices.Add(ci);

	// dialog
	wxSingleChoiceDialog colordlg(this, "Color by:", "choose", choices);
	// mark the current pick
	int current = viz.colorby;
	if (current >= 0)
		colordlg.SetSelection(current);
	// show the dialog, get selection
	int stat = colordlg.ShowModal();
	if (stat == wxID_CANCEL)
		return;
	viz.colorby = colordlg.GetSelection();
	T_(trc.dprint("got colorby ",viz.colorby,", current: ",current);)

	// re-create legends if colorby has changed
	if (viz.colorby != current) {
		currentview->set_legends();
		Refresh();
	}
}

void
VizFrame::
OnLegend(wxCommandEvent& event) {
// toggle between legend and parameter values on right panel
	T_(Trace trc(1,"OnLegend");)
	Vizplot* m_plot = plot();
	m_plot->showLegend = !m_plot->showLegend;
	if (m_plot->showLegend) {
		legendButton->SetLabel("Values");
	} else {
		legendButton->SetLabel("Legend");
	}
	leftpanel->legendbox->Scroll(0,0);
	leftpanel->Refresh();
	//!! m_plot->Refresh();		// XXX I don't understand why this is needed
}

// Help|Help... command
void
VizFrame::
OnMenuHelpHelp( wxCommandEvent& WXUNUSED(event) ) {
	wxString name{"wxWidgets"};
	wxVersionInfo vi(name);
	wxMessageBox(wxVERSION_STRING, "Help", wxICON_INFORMATION);
	//wxMessageBox(Help, wxT("Help"), wxICON_INFORMATION);
}

void
VizFrame::
OnFit(wxCommandEvent& WXUNUSED(event)) {
	T_(Trace trc(2,"OnFit");)
	Vizplot* m_plot = plot();
	m_plot->Fit();
	vector<double> bb(6,0.0);
	m_plot->GetBoundingBox(bb.data());
	T_(trc.dprint("bounding box: ",bb);)
}

void
VizFrame::
OnAlignXAxis( wxCommandEvent &WXUNUSED(event) ) {
	Vizplot* mp = plot();
    mp->axesPos[0] = (int) (mp->axesPos[0]+1)%5;
    wxString temp;
    //!! temp.sprintf(wxT("axesPos = %d\n"), axesPos);
    //!! m_log->AppendText(temp);
	if (mp->xaxis == nullptr) return;
	if (mp->axesPos[0] == 0) {
            mp->xaxis->SetAlign(mpALIGN_BOTTOM);
            mp->SetMarginTop(vizsizes.top);
            mp->SetMarginBottom(vizsizes.bottom);
	}
	if (mp->axesPos[0] == 1) {
        //((mpScaleX*)(mp->GetLayer(0)))->SetAlign(mpALIGN_BOTTOM);
            mp->xaxis->SetAlign(mpALIGN_BORDER_BOTTOM);
            mp->SetMarginTop(0);
            mp->SetMarginBottom(50);
	}
	if (mp->axesPos[0] == 2) {
        //((mpScaleX*)(mp->GetLayer(0)))->SetAlign(mpALIGN_CENTER);
            mp->xaxis->SetAlign(mpALIGN_CENTER);
            mp->SetMarginTop(0);
            mp->SetMarginBottom(0);
	}
	if (mp->axesPos[0] == 3) {
        //((mpScaleX*)(mp->GetLayer(0)))->SetAlign(mpALIGN_TOP);
            mp->xaxis->SetAlign(mpALIGN_TOP);
            mp->SetMarginTop(50);
            mp->SetMarginBottom(0);
	}
	if (mp->axesPos[0] == 4) {
        ((mpScaleX*)(mp->GetLayer(0)))->SetAlign(mpALIGN_BORDER_TOP);
            mp->xaxis->SetAlign(mpALIGN_BORDER_TOP);
            mp->SetMarginTop(0);
            mp->SetMarginBottom(0);
	}
    mp->UpdateAll();
}

void
VizFrame::
OnAlignYAxis( wxCommandEvent &WXUNUSED(event) ) {
	Vizplot* mp = plot();
    mp->axesPos[1] = (int) (mp->axesPos[1]+1)%5;
    wxString temp;
    //!! temp.sprintf(wxT("axesPos = %d\n"), axesPos);
    //!! m_log->AppendText(temp);
	if (mp->yaxis == nullptr) return;
	if (mp->axesPos[1] == 0) {
        //((mpScaleY*)(mp->GetLayer(1)))->SetAlign(mpALIGN_BORDER_LEFT);
		mp->yaxis->SetAlign(mpALIGN_LEFT);
		mp->SetMarginLeft(vizsizes.left);
		mp->SetMarginRight(vizsizes.right);
	}
	if (mp->axesPos[1] == 1) {
        //((mpScaleY*)(mp->GetLayer(1)))->SetAlign(mpALIGN_LEFT);
		mp->yaxis->SetAlign(mpALIGN_BORDER_LEFT);
		mp->SetMarginLeft(70);
		mp->SetMarginRight(0);
	}
	if (mp->axesPos[1] == 2) {
        //((mpScaleY*)(mp->GetLayer(1)))->SetAlign(mpALIGN_CENTER);
		mp->yaxis->SetAlign(mpALIGN_CENTER);
		mp->SetMarginLeft(0);
		mp->SetMarginRight(0);
	}
	if (mp->axesPos[1] == 3) {
        //((mpScaleY*)(mp->GetLayer(1)))->SetAlign(mpALIGN_RIGHT);
		mp->yaxis->SetAlign(mpALIGN_RIGHT);
		mp->SetMarginLeft(0);
		mp->SetMarginRight(70);
	}
	if (mp->axesPos[1] == 4) {
        //((mpScaleY*)(mp->GetLayer(1)))->SetAlign(mpALIGN_BORDER_RIGHT);
		mp->yaxis->SetAlign(mpALIGN_BORDER_RIGHT);
		mp->SetMarginLeft(0);
		mp->SetMarginRight(0);
	}
    mp->UpdateAll();
}

void
VizFrame::
OnLimits( wxCommandEvent &WXUNUSED(event) ) {
// handle limits button
	T_(Trace trc(2,"OnLimits");)
	Vizplot* mp = plot();

	T_(trc.dprint("plot: ",mp->xname,", ",mp->ynames[0]);)
	LimitsDialog dlg(this);
	if (dlg.ShowModal() == wxID_CANCEL)
		return;
	wxString xminstr = dlg.GetXminText();
	wxString xmaxstr = dlg.GetXmaxText();
	wxString yminstr = dlg.GetYminText();
	wxString ymaxstr = dlg.GetYmaxText();
	T_(trc.dprint("got [",xminstr,", ",xmaxstr,"], [",yminstr,", ",ymaxstr,"]");)
	mp->bbxmin = stod(xminstr.ToStdString());
	mp->bbxmax = stod(xmaxstr.ToStdString());
	mp->bbymin = stod(yminstr.ToStdString());
	mp->bbymax = stod(ymaxstr.ToStdString());
	mp->Fit(mp->bbxmin,mp->bbxmax, mp->bbymin, mp->bbymax);
}

void
VizFrame::
OnToggleGrid( wxCommandEvent &WXUNUSED(event) ) {
	// XXX re-write: set name on axes, getlayer by name, GetTicks
	Vizplot* mp = plot();
    mp->ticks = !mp->ticks;
	mp->xaxis->SetTicks(mp->ticks);
	mp->yaxis->SetTicks(mp->ticks);
    mp->UpdateAll();
}

void
VizFrame::
OnToggleLegend( wxCommandEvent &WXUNUSED(event) ) {
// handler for the View->toggle legend/parameters XXX deprecate?
// switch between the legend and parameter list on right border
	Vizplot* mp = plot();
	mp->showLegend = !mp->showLegend;
	leftpanel->Refresh();
	mp->Refresh();		// XXX I don't understand why this is needed
}

void
VizFrame::
OnScreenshot(wxCommandEvent& event) {
// pop up a dialog to allow the user to name the file; extension determines
// the file format
	Vizplot* mp = plot();

	// enable all image handlers
	wxInitAllImageHandlers();

	// create the dialog...
	wxFileDialog fileDialog(this, _("Save a screenshot"), wxT(""), wxT(""),
		"JPEG (*.jpg) | *.jpeg;*.jpg | GIF (*.gif)|*.gif | BMP (*.bmp) | *.bmp|PNG (*.png) | *.png",
		//!! "BMP image (*.bmp) | *.bmp|JPEG image (*.jpg) | *.jpeg;*.jpg|PNG image (*.png) | *.png",
				wxFD_SAVE|wxFD_OVERWRITE_PROMPT);

	// pop up and get user's response
	if(fileDialog.ShowModal() == wxID_OK) {
		wxFileName namePath(fileDialog.GetPath());
#ifdef NEVER // default jpg
		int fileType = wxBITMAP_TYPE_BMP;
#else // NEVER // default jpg
		int fileType = wxBITMAP_TYPE_JPEG;
		if( namePath.GetExt().CmpNoCase(wxT("bmp")) == 0 ) fileType = wxBITMAP_TYPE_BMP;
#endif // NEVER // default jpg
		if( namePath.GetExt().CmpNoCase(wxT("jpeg")) == 0 ) fileType = wxBITMAP_TYPE_JPEG;
		if( namePath.GetExt().CmpNoCase(wxT("jpg")) == 0 )  fileType = wxBITMAP_TYPE_JPEG;
		if( namePath.GetExt().CmpNoCase(wxT("png")) == 0 )  fileType = wxBITMAP_TYPE_PNG;
		if( namePath.GetExt().CmpNoCase(wxT("gif")) == 0 )  fileType = wxBITMAP_TYPE_GIF;
		wxSize imgSize(500,500);
		mp->SaveScreenshot(fileDialog.GetPath(), fileType, imgSize, false);
	}
	event.Skip();
}
// end VizFrame implementation }

void
myFXYVector::
Plot(wxDC & dc, mpWindow & w) {
// My override of mpFXYVector::Plot to
//  - create a clipping region
//  - draw both lines and dots
	T_(Trace trc(2,"myFXYVector::Plot");)

	if (!m_visible) {
		T_(trc.dprint("quick return: not visible");)
		return;
	}

	dc.SetPen(legend.pen);
 
	// set up a clipping region
	int tlx = w.GetMarginLeft();
	int tly = w.GetMarginTop();
	int width = w.GetScrX()-w.GetMarginRight() - tlx;
	int height = w.GetScrY() - w.GetMarginBottom() - tly;
	T_(trc.dprint("clipping rect: ",tlx,", ",tly,", w ",width,", h ",height);)
	
	// try expanding the box a little so pts on the boundary show
	tlx -= 2;
	tly -= 2;
	width += 4;
	height += 4;
	wxRect plotBox(tlx, tly, width, height);
	dc.SetClippingRegion(plotBox);

	// call the base class mpFXYVector::Plot: it will draw
	// lines assuming SetContinuity(true) has been called
	mpFXYVector::Plot(dc,w);

	wxCoord ix{0};
	wxCoord iy{0};
	double x, y;
	// draw dots at each point
	// check that at least 1 point is in the clipping region
	inside = false;
	Rewind();
	wxPen dotpen = legend.pen;
	dotpen.SetWidth(5);
	dc.SetPen(dotpen);
	while (GetNextXY(x,y)) {
		ix = w.x2p(x);
		iy = w.y2p(y);
		if (plotBox.Contains(ix,iy)) {
			inside = true;
			legend.drawDot(dc, ix, iy, 2);
		}
	}

	dc.DestroyClippingRegion();
}

// myScroll implementation {

void
myScroll::
paintLegend(wxDC& dc, Vizplot* vp) {
// write the current legend to dc
	T_(Trace trc(2,"paintLegend");)

	dc.Clear();
	int width, height;
	wxCoord ix{20};
	wxCoord iy{60};
	// only treat visible myFXYVector layers
	// However, even visible layers might not show up due to clipping
	vector<myFXYVector*> layers = vp->visplotlayers();
	for (auto ly : layers) {
		string label = ly->legend.num;
		if (!ly->legend.ext.empty())
			label += '.' + ly->legend.ext;
		if (vp->ynames.size() > 1)
			label += ' ' + ly->legend.yname;
		if (label.empty()) continue;
		dc.GetTextExtent(label, &width, &height);
		wxPen lpen = ly->legend.pen;
		// the dot...
		wxPen dotpen = lpen;
		dotpen.SetWidth(5);
		if (!ly->inside) {
			wxColour dim = lpen.GetColour();
			dim.MakeDisabled(150);
			dotpen.SetColour(dim);
		} else
			dotpen.SetColour(lpen.GetColour());
		dc.SetPen(dotpen);
		ly->legend.drawDot(dc,ix,iy, 4);
		// ... and it's label, shaded if not visible
		if (!ly->inside)
			dc.SetTextForeground(wxColour(*wxLIGHT_GREY));
		else
			dc.SetTextForeground(wxColour(*wxBLACK));
		dc.DrawText(label, ix+10, iy-height/2);
		iy += height + 4;
	}
}

void
myScroll::
paintPset(wxDC& dc, VizFrame* vizframe) {
// write the string "pset" to dc
	T_(Trace trc(2,"paintPset");)

	dc.Clear();

	wxCoord ix{0};
	wxCoord iy{0};
	dc.DrawText(*vizframe->pset, ix, iy);
}

// rm'd - caused problems with GL?
//!! #include "main.h"

double XParValue;
double YParValue;
string ModeName;

int
main(int argc, char** argv) {
	char* dbg = getenv("DEBUG");
	if (dbg != nullptr) Trace::debug(std::stoi(dbg));
	Ftmpdir ftmp;
	T_(Trace trc(1,"main");)

	T_(trc.dprint("component sizes in pixels:\n",vizsizes);)

		// parse the specs...
	Specs& sp = specs();
	try {
		// ...from the command line ...
		string options = cmdline(argc, argv);
		if (options.empty() && isatty(fileno(stdin))) {
			cout << helpmsg() << endl;
		// ... or flaps options
		} else if (!parse_specs(options, sp)) {
			cout << helpmsg() << endl;
			return 1;
		}
	} catch (const std::exception& t) {
		flaps::error("caught exception: ",t.what());
		exit(1);
	}

	try {
		// get the Viz instance - this will create the singleton...
		Viz& viz{Viz::instance()};
		// ... then load requested files/aids
		viz.loadcurves(sp.aids, sp.plotfiles);
		// Create an Amviz singleton...
		Amviz& amviz{Amviz::instance()};
		if (amviz.error.empty()) {
			// ... then load amviz data in a separate thread
			// if the user right-clicks feedAmviz will be called which
			// will call amviz_run in a thread
			if (sp.threads)
				amviz.kit_future = async(amviz_kit);
			else
				amviz_kit();
		}
		// and create the plot - this calls wxEntry
		viz.run();

	} catch (runtime_error& s) {
		flaps::error("caught exception: ",s.what());
	} catch (const std::exception& t) {
		flaps::error("caught exception: ",t.what());
	}
}	// main
