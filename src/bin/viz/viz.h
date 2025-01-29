//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
////////////////////////////////////////////////////////////////////////////
//  Flaps plotting
/////////////////////////////////////////////////////////////////////////////

#ifndef VIZ_H
#define VIZ_H

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "mathplot.h"
#include "plotcurve.h"
#include "specs.h"

#include "wx/defs.h"
#include "wx/display.h"
#include "wx/utils.h"
#include "wx/app.h"
#include "wx/button.h"
#include "wx/cursor.h"
#include "wx/dc.h"
#include "wx/dcclient.h"
#include "wx/event.h"
#include "wx/menu.h"
#include "wx/wfstream.h"
#if wxUSE_ZLIB
#include "wx/zstream.h"
#endif

#include "wx/wx.h"
#include "wx/dir.h"
#include "wx/filedlg.h"
#include "wx/filedlgcustomize.h"
#include "wx/sizer.h"
#include "wx/gtk/print.h"

#include "trackball.h"

void message(const std::string& msg);

class Vizplot;
class VizFrame;
class myFXYVector;
class LeftPanel;
class LowerPanel;

// Sizes of the various components
class VizSizes {
public:
	VizSizes();

	wxSize display;
	wxSize plot;
	wxSize panel;
	wxSize infobox;
	wxSize legendbox;
	wxSize footnotes;
	wxSize menu;
	wxSize toolbar;
	wxSize frame;
	int top{20};
	int right{20};
	int bottom{60};
	int left{100};
};

// A Mark is used to put a small rectangle around a solution point,
// usually one that is closest to a point clicked on
class Mark {
	std::string text;
public:
	int layerno{-1};
	int idx{0};		// 0b index of the point in layer
	// curve is the Plotcurve corresponding to "layer"
	Plotcurve* plotcurve;
	Mark() {}
	Mark(Vizplot* vp, int ix, int iy, std::string txt="");

	// draw a rectangle at the closest point on vp (may be different than vp
	// in the ctor)
	void drawrect(Vizplot* vp, wxDC& dc, int size);
};

class myFXYVector : public mpFXYVector {
public:
	// ctor
	myFXYVector(Vizplot* vp, Plotcurve* pc, const std::string& ynm) :
			mpFXYVector(""), vizplot(vp), plotcurve(pc) {
	// doesn't do much: data is loaded in SetData()
		m_type = mpLAYER_PLOT;
		legend = Legend(pc,ynm);	// set Plotcurve, yname for set_legends
		yname = ynm;
		// SetName(legend.id());
	 }
	// XXX also include Vizplot*? ie which Vizplot is this associated with
	//     this would allow for resizing the Vizplot
	Vizplot* vizplot{nullptr};
	Plotcurve* plotcurve{nullptr};
	Legend legend;
	std::string yname;
	bool inside{false};	// the clipping region?

	// given the index of a point on this layer, return a pointer to
	// the corresponding Plotcurve and the index of the corresponding point

	// GetMinX, etc are protected in mpFXYVector
	double getminx() { return GetMinX(); }
	double getmaxx() { return GetMaxX(); }
	double getminy() { return GetMinY(); }
	double getmaxy() { return GetMaxY(); }
	void rewind() { return Rewind(); }
	bool getnextxy(double& x, double& y) { GetNextXY(x,y); return true; }

	// return const references
	const std::vector<double>& getxs() { return m_xs; }
	const std::vector<double>& getys() { return m_ys; }

	 // my custom Plot to draw the dots
	 void Plot(wxDC& dc, mpWindow& w);

	 // get the index in the soln member of my plotcurve
	 // corresponding to an index into my xs,ys members
	 int curveidx(int idx);

	 // get the screen coordinates of the idx point in
	 // my (xs, ys)
	 void screencoords(int idx, int& ix, int& iy);
};	// myFXYVector

class VizApp : public wxApp {
	bool OnInit();	// for wxApp
};

// vector<string> specialization is in viz.cpp
template<typename T>
bool contains(const std::vector<T>& v, const T& i) {
	return (std::find(v.begin(), v.end(), i) != v.end());
}

class CurvesDialog : public wxDialog {
public:
	// ctor
	CurvesDialog(wxWindow* parent, wxArrayString&, bool clr);
	// getters for input values
	wxString GetInclText() const {return includes->GetValue(); }
	wxString GetExclText() const {return excludes->GetValue(); }
	std::vector<int> getChecklist();
	wxListBox* checklist;	// a list of all curves
	bool clear{false};		// clear the list?
private:
	wxTextCtrl* includes;
	wxTextCtrl* excludes;
	wxCheckBox* clearbox{nullptr};
};	// class CurvesDialog

class LimitsDialog : public wxDialog {
public:
	// ctor
	LimitsDialog(VizFrame* parent);

	VizFrame* vizframe;

	// getters for input values
	wxString GetXminText() const {return xminctrl->GetValue(); }
	wxString GetXmaxText() const {return xmaxctrl->GetValue(); }
	wxString GetYminText() const {return yminctrl->GetValue(); }
	wxString GetYmaxText() const {return ymaxctrl->GetValue(); }
private:
	wxTextCtrl* xminctrl{nullptr};
	wxTextCtrl* xmaxctrl{nullptr};
	wxTextCtrl* yminctrl{nullptr};
	wxTextCtrl* ymaxctrl{nullptr};
};	// class LimitsDialog

class ParamsDialog : public wxDialog {
	std::string xsel;
	std::vector<std::string> ysels;
public:
	// ctor
	ParamsDialog(VizFrame* parent);

	VizFrame* vizframe;

	// getters for input values
	wxString GetXText() const {
		wxString rval = xctrl->GetValue();
		if (rval.empty())
			rval = xsel;
		return rval;
	}
	std::vector<std::string> GetYText() const {
	// y names may be separated by a space, comma, or colon
		Trace trc(2,"GetYText");
		std::string sel = yctrl->GetValue().ToStdString();
		trc.dprint("typed: \"",sel,"\"");
		std::vector<std::string> rval = string2tok(sel," ,:");
		if (rval.empty())
			rval = ysels;
		trc.dprint("returning ",rval);
		return rval;
	}
private:
	wxTextCtrl* xctrl{nullptr};
	wxTextCtrl* yctrl{nullptr};
};	// class ParamsDialog

// Infobox is a panel within Leftpanel
class Infobox : public wxPanel {
	LeftPanel* leftpanel;
public:
	Infobox(LeftPanel* parent);
	void OnPaint(wxPaintEvent& event);
};

// myScroll is a scrolling panel with drawing area "canvas"
class myScroll : public wxScrolled<wxPanel> {
	wxPanel* canvas{nullptr};
	wxBoxSizer* csizer{nullptr};
public:
	myScroll(LeftPanel* p, wxWindowID id, const wxPoint& pos, const wxSize& size);

	LeftPanel* leftpanel{nullptr};

	void paintLegend(wxDC& dc, Vizplot* mp);
	void paintPset(wxDC& dc, VizFrame* vizframe);

	void OnPaint(wxPaintEvent& event);

};

// Footnotes is a scrolling panel with drawing area "canvas"
class Footnotes : public wxScrolled<wxPanel> {
	wxPanel* canvas{nullptr};
	wxBoxSizer* csizer{nullptr};
public:
	Footnotes(LeftPanel* p, wxWindowID id, const wxPoint& pos, const wxSize& size);

	LeftPanel* leftpanel{nullptr};

	void paintFootnotes(wxDC& dc, Vizplot* mp);

	void OnPaint(wxPaintEvent& event);

};


// vertical panel on the left side
class LeftPanel : public wxPanel {
public:
	LeftPanel(LowerPanel* parent, wxWindowID id=-1,
			const wxPoint& pos=wxDefaultPosition, const wxSize& size=wxDefaultSize);

	wxBoxSizer* lpsizer{nullptr};

	LowerPanel* lowerpanel{nullptr};

	void OnPaint(wxPaintEvent& event);

	// components of Leftpanel:
	Infobox* infobox{nullptr};
	myScroll* legendbox{nullptr};
	Footnotes* footnotes{nullptr};
};

// A Viz contains the data that we wish to plot
// It is a singleton with an automatic instance so that it
// will be deleted on exit
class Viz {
	std::vector<std::string> parnames_;

	// ctor: does not do much: load data with load()
	Viz() {}
	// delete the copy ctor & assignment op
	Viz(Viz&) = delete;
	void operator=(const Viz&) = delete;
public:
	//!! VizFrame* frame{nullptr};		// the main frame XXX NO - multiple frames
	std::map<std::string,Legend> footnotes;
	int colorby{-1};
	static Viz& instance() {
		static Viz instance_;
		return instance_;
	}
	std::vector<Plotcurve*> loaddata(const std::vector<std::string>& aids,
		const std::vector<std::string>& plotfiles);
	std::vector<Plotcurve*> curves;
	void run();

	// all available parameter names
	const std::vector<std::string>& parnames() { return parnames_; }
// XXX use add2vector
	void add_parname(std::string& nm) {
		// add a name to the list of parameters if it is not already there
		if (std::find(parnames_.begin(),parnames_.end(),nm) == parnames_.end())
			parnames_.push_back(nm);
	}

	// parameter title for axes
	std::string parTitle(const std::string& name) const;

	// delete all curves and parnames
	void clear() {
		for (auto cp : curves)
			delete cp;
		curves.clear();
		parnames_.clear();
		footnotes.clear();
	}
};	// Viz

class Vizplot : public mpWindow {
public:
	std::string xname;
	std::vector<std::string> ynames;
	int axesPos[2]{0,0};	// default: bottom,left
	bool ticks{false};
	//!! mpInfoCoords* infocoords{nullptr};	// show current mouse coords
	wxString ibcoords;	// for LeftPanel->infobox
	// show either the legend or pset
	bool showLegend{true};
	wxMenu middleMenu;   // menu shown on middle-click
	double origxmin;		// original limit in screen units
	double origxmax;
	double origymin;
	double origymax;
	double bbxmin;			// current limit in screen units
	double bbxmax;
	double bbymin;
	double bbymax;
	mpScaleX* xaxis{nullptr};
	mpScaleY* yaxis{nullptr};
	std::vector<myFXYVector*> plotlayers();
	std::vector<myFXYVector*> visplotlayers();
	bool is_complete() { return !plotlayers().empty(); }
	LowerPanel* lowerpanel{nullptr};		// parent window
	VizFrame* vizframe{nullptr};						// lowerpanel->vizframe

	// ctor
	Vizplot(LowerPanel* parent, wxWindowID id,
			const wxPoint& pos, const wxSize& size);

	static constexpr double a__ = -std::numeric_limits<double>::max();
	static constexpr double b__ = std::numeric_limits<double>::max();
	void loadxy(const std::string& xnm, const std::vector<std::string>& ynms,
			const std::vector<int>& topl={},
			double xmin=a__, double xmax=b__, double ymin=a__, double ymax=b__);
	void add(std::vector<Plotcurve*>& newcurves);

	// for zooming
	wxPoint startPos;
	wxPoint currentPos;
	wxPoint transform(Vizplot* from, wxPoint& p);
private:
	// send the eigenvector corresponding to the closest point to a right mouse click
	void feedAmviz(int x, int y);

	static long flags;

	 // mouse events
	void OnOrigView(wxCommandEvent& WXUNUSED(event));
	void OnMouseHelp(wxCommandEvent& WXUNUSED(event));
	bool isDragging{false};
	bool leftDown{false};
	bool middleDown{false};
	bool rightDown{false};
	void OnLeftDown(wxMouseEvent& event);
	void OnLeftUp(wxMouseEvent& event);
	void OnMiddleDown(wxMouseEvent& event);
	void OnMiddleUp(wxMouseEvent& event);
	void OnRightDown(wxMouseEvent& event);
	void OnRightUp(wxMouseEvent& event);
	void OnMotion(wxMouseEvent& event);
public:
	// find the closest point to a left mouse click
	// Returns: (curve#, soln#, x, y, myFXYVector*)
	// New: (myFXYVector*,jx,jy)
	std::tuple<int,int> closestCurve(int ix, int iy);

	void OnPaint(wxPaintEvent& event);

	// create a Legend for each curve
	void set_legends();

}; // class Vizplot


// class VisFrame: container for a Vizplot, left panel, tools, and menus
class VizFrame: public wxFrame {
private:
	wxBoxSizer* topsizer{nullptr};	// tools + lowersizer
	wxBoxSizer* lowersizer{nullptr};	// LowerPanel sizer: left panel + plot
	wxMenuBar* menubar{nullptr};
	wxPanel* toolbar{nullptr};
	LowerPanel* lowerpanel{nullptr};
	std::vector<int> toplot;
	int frameno{0};		// frame number
public:
	VizFrame* parent{nullptr};		// only for zoom ctor
	static int number_of_frames;
	// show either the legend or pset
	bool showLegend{true};
	LeftPanel* leftpanel{nullptr};

	Vizplot* currentview{nullptr};
	std::vector<Vizplot*> plotz;
	// 2 types of marks
	Mark* pset_spot{nullptr};
	std::vector<Mark>* am_spot{nullptr};
	std::string* pset{nullptr};

	// initial ctor: no parent
	VizFrame(const wxString& title, const wxPoint& pos,
		const wxSize& size, long style = wxDEFAULT_FRAME_STYLE);
	// zoom ctor: take data from parent & vizplot
	VizFrame(Vizplot* vizplot, VizFrame* parent, const wxString& title,
		const wxPoint& pos, const wxSize& size, long style = wxDEFAULT_FRAME_STYLE);

	// set the frame title (upper border)
	void set_title(const std::string& t="");

	// delete all accumulated plots
	void clear() {
		currentview = nullptr;
		for (auto& pi : plotz)
			delete pi;
		plotz.clear();
	}


	void createControls();

	// add a plot, return it's 0b index into plotz: must call setcurrent next
	// XXX unless this is the first addition to an empty plotz?
	int addplot(Vizplot* p);

	int currentidx() {
	// returns the 0b index in plotz of currentview
		for (size_t i=0; i<plotz.size(); i++) {
			if (plotz[i] == currentview)
				return i;
		}
		return -1;
	}

	Vizplot* setcurrent(int cv) {
	// set plotz[cv] as the current plot
		Trace trc(2,"setcurrent");
		// sanity check: cv in range
		if (cv < 0 || cv >= (int)plotz.size())
			throw std::runtime_error(vastr("attempt to make view ",cv," current: only ",
					plotz.size(), " available"));
		if (currentview == nullptr) {
			currentview = plotz[0];
			lowersizer->Add(currentview, 1, wxEXPAND);
			currentview->Show(true);
			trc.dprint("added initial ",currentview,", lowersizer: ",*lowersizer);
			return currentview;
		}

		// first detach the current plot...
		trc.dprint("current view ",currentview,", new view ",cv,", plotz: ",plotz.size());
		if (plotz[cv] == currentview)
			return currentview;
		// .. and the currentview
		bool ok = lowersizer->Detach(currentview);
		trc.dprint("detached ",currentview,", new plot ",plotz[cv]);
		if (ok)
			currentview->Show(false);
		else
			std::cerr << "detaching current view failed\n";
		// ... then add the new one
		currentview = plotz[cv];
		currentview->Show(true);
		lowersizer->Add(currentview, 1, wxEXPAND);
		topsizer->Layout();
		Layout();
		Refresh();
		return currentview;
	}

	// access the current plot
	Vizplot* plot() { return currentview; }

	wxStatusBar* statusBar;
	wxMenu* curvesMenu{nullptr};	// XXX unused?
	wxMenu* includeMenu{nullptr};	// XXX unused?
	wxRadioBox* radioView{nullptr};	// XXX unused?

	// event handlers
	void OnClose(wxCloseEvent& event);
    void OnFileOpen(wxCommandEvent& event);
    void OnFileAdd(wxCommandEvent& event);
    void OnNewWindow(wxCommandEvent& event);
    void OnFileXfig(wxCommandEvent& event);
    //!! void OnFileVideo(wxCommandEvent& event);
    //!! void OnExit(wxCommandEvent& event);
    void OnEditPref(wxCommandEvent& event);
	 // mathplot stuff: View menu
	 void OnFit(wxCommandEvent& event);
	 void OnAlignXAxis(wxCommandEvent& event);
	 void OnAlignYAxis(wxCommandEvent& event);
	 void OnToggleGrid(wxCommandEvent& WXUNUSED(event));
	 void OnToggleLegend(wxCommandEvent& WXUNUSED(event));
	void OnScreenshot(wxCommandEvent& WXUNUSED(event));
	void OnOrigView(wxCommandEvent& WXUNUSED(event));
    // Help menu
    void OnMenuHelpAbout(wxCommandEvent& event);
    void OnMenuHelpHelp(wxCommandEvent& event);

	 // toolbar:
	 // Curve selection:
	 void OnCurvesSelect(wxCommandEvent& event);

	 // parameter selection:
	 std::string getParChoice(const std::string& current, const std::string& axis);
	 std::vector<std::string> getYParChoice(const std::vector<std::string>& current);
	 void OnParamsSelect(wxCommandEvent& event);

	 // view selection
	 void OnViewSelect(wxCommandEvent& event);
	 // coloring curves
	 void OnColor(wxCommandEvent& event);

	 // set parameter limits XXX rm from top bar?
	 void OnLimits(wxCommandEvent& WXUNUSED(event));
	 wxButton* limitsButton{nullptr};

	 // toggle legend/pset
	 void OnLegend(wxCommandEvent& event);
	 wxButton* legendButton{nullptr};

	 bool ShowToolTips();

};  // VizFrame

// the LowerPanel contains the LeftPanel and the Vizplot
class LowerPanel : public wxPanel {
public:
	LowerPanel(VizFrame* vf) : wxPanel(vf), vizframe(vf) {}

	VizFrame* vizframe;				// parent
	Vizplot* vizplot{nullptr};		// currentview?
};

class VizPrintout : public wxPrintout {
public:
	VizPrintout(wxString const& title) : wxPrintout(title) { }
	void GetPageInfo(int* minPage, int* maxPage,
			int* pageFrom, int* pageTo) {
		*minPage = *maxPage = 1;
		*pageFrom = *pageTo = 1;
	}
	bool HasPage(int pagenum) { return true; }
	bool OnPrintPage(int page) { return true; }
};  // VizPrintout

#endif // #ifndef VIZ_H

