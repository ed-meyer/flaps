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
//  Flaps 3D Visualizer
/////////////////////////////////////////////////////////////////////////////

#ifndef VZ3D_H
#define VZ3D_H

#include <vector>

#include "plotcurve.h"
#include "specs.h"

#include "wx/defs.h"
#include "wx/utils.h"
#include "wx/app.h"
#include "wx/button.h"
#include "wx/cursor.h"
#include "wx/event.h"
#include "wx/menu.h"
#include "wx/dcclient.h"
#include "wx/wfstream.h"
#if wxUSE_ZLIB
#include "wx/zstream.h"
#endif

#include "wx/glcanvas.h"
#include "wx/wx.h"
#include "wx/gtk/print.h"

#include "trackball.h"

class MyApp : public wxApp {
	bool OnInit();	// for wxApp
};

class Vz3d {
	static Vz3d* instance_;
	std::string xname_;
	std::string yname_;
	std::string zname_;
	std::vector<std::string> parnames_;
	// constructor
	Vz3d(const std::vector<Plotcurve*>& curves, Specs& specs);
	// delete the copy constructor & assignment op
	Vz3d(Vz3d&) = delete;
	void operator=(const Vz3d&) = delete;
public:
	std::vector<Plotcurve*> curves;
	std::vector<int> toplot;		// 0b curve number
	std::vector<std::vector<double>> data;

	static Vz3d* instance() { return instance_; }
	static Vz3d* instance(const std::vector<Plotcurve*>& curves, Specs& specs);
	void run();
	void render();
	double amplitude(int newamp);
	void load_data();

	// parameter names
	const std::vector<std::string>& parnames() { return parnames_; }
	const std::string& xname() { return xname_; }
	const std::string& xname(const std::string& newname) {
		xname_ = newname;
		return xname_;
	}
	const std::string& yname() { return yname_; }
	const std::string& yname(const std::string& newname) {
		yname_ = newname;
		return yname_;
	}
	const std::string& zname() { return zname_; }
	const std::string& zname(const std::string& newname) {
		zname_ = newname;
		return zname_;
	}
};

// main entry point to vz3d XXX replace with run()
//!! void vz3d();

// default background: white
constexpr float BGRed{1.0};
constexpr float BGGreen{1.0};
constexpr float BGBlue{1.0};
constexpr float BGAlpha{0.0};

// OpenGL view data
struct GLData {
    bool initialized;           // has OpenGL been initialized?
    float beginx, beginy;       // position of mouse
	 float lastRadius;
    float lastX, lastY;       // position of mouse XXX remove
	 std::vector<float> quat;  // orientation of object
    float zoom;                 // field of view in degrees
};

// Define a new frame type
class FGLCanvas;
//!! class MyTimer;

class MyFrame: public wxFrame {
private:
	FGLCanvas *m_canvas{nullptr};
	wxToggleButton* toggleX{nullptr};
public:
    MyFrame(wxFrame *frame, const wxString& title, const wxPoint& pos,
        const wxSize& size, long style = wxDEFAULT_FRAME_STYLE);

	void createControls();

	FGLCanvas* glcanvas() { return m_canvas; }
	 
	wxStatusBar* statusBar;
	wxMenu* curvesMenu{nullptr};
	wxMenu* includeMenu{nullptr};
	wxRadioBox* radioView{nullptr};
	wxRadioButton* zbutton{nullptr};
	wxRadioButton* ybutton{nullptr};
	wxRadioButton* xbutton{nullptr};
	wxRadioButton* defaultbutton{nullptr};

	// event handlers
    void OnMenuFileOpen(wxCommandEvent& event);
    void OnMenuFilePrint(wxCommandEvent& event);
    void OnMenuFileXfig(wxCommandEvent& event);
    //!! void OnMenuFileVideo(wxCommandEvent& event);
    void OnMenuFileExit(wxCommandEvent& event);
    void OnMenuEditPref(wxCommandEvent& event);
    //!! void OnMenuMode(wxCommandEvent& event);
    void OnMenuHelpAbout(wxCommandEvent& event);
    void OnMenuHelpHelp(wxCommandEvent& event);
    void OnMenuViewDefault(wxCommandEvent& event);
    void OnMenuViewX(wxCommandEvent& event);
	 void OnToggleX(wxCommandEvent& event);
    void OnMenuViewY(wxCommandEvent& event);
    void OnMenuViewZ(wxCommandEvent& event);
    void OnMenuViewI(wxCommandEvent& event);
    void OnMenuViewQuad(wxCommandEvent& event);
	 // top tools:
	 void OnXChoice(wxCommandEvent& event);
	 void OnYChoice(wxCommandEvent& event);
	 void OnZChoice(wxCommandEvent& event);
	 void OnDone(wxCommandEvent& event);
	 // Curve selection:
	 // void OnInclude(wxCommandEvent& event);
	 // void OnExclude(wxCommandEvent& event);
	 void OnCurvesSelect(wxCommandEvent& event);
	 //!! void OnCurvesDone(wxCommandEvent& event);
	 // left side toolbar:
	 std::string getParChoice(const std::string& current, const std::string& axis);
	 void OnCurveNumber(wxCommandEvent& event);
	 void OnXParamSelect(wxCommandEvent& event);
	 void OnYParamSelect(wxCommandEvent& event);
	 void OnZParamSelect(wxCommandEvent& event);
	 void OnPlot(wxCommandEvent& event);
	 void OnRadioView(wxCommandEvent& event);
	 void OnZView(wxCommandEvent& event);
	 void OnYView(wxCommandEvent& event);
	 void OnXView(wxCommandEvent& event);
	 void OnDefaultView(wxCommandEvent& event);

	 bool ShowToolTips();

	 void vibrate(void);

    //!! void SetCanvas( FGLCanvas *canvas ) { m_canvas = canvas; }
    //!! FGLCanvas *GetCanvas() { return m_canvas; }

	 //!! bool isRunning();

	 void print();
	 void printPS();

    DECLARE_EVENT_TABLE()
};  // MyFrame

class MyPrintout : public wxPrintout {
public:
	MyPrintout(wxString const& title) : wxPrintout(title) { }
	void GetPageInfo(int* minPage, int* maxPage,
			int* pageFrom, int* pageTo) {
		*minPage = *maxPage = 1;
		*pageFrom = *pageTo = 1;
	}
	bool HasPage(int pagenum) { return true; }
	bool OnPrintPage(int page) { return true; }
};  // MyPrintout


#if wxUSE_GLCANVAS

class FGLCanvas: public wxGLCanvas {
public:
#ifdef NEVER // new wxGLAttributes
    FGLCanvas(wxWindow *parent, wxWindowID id = wxID_ANY,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize, long style = 0,
        const wxString& name = wxT("FGLCanvas"), const int* attlist = nullptr);
#else // NEVER // new wxGLAttributes
    FGLCanvas(wxWindow *parent, const wxGLAttributes& dispAttrs,
			wxWindowID id = wxID_ANY,
			const wxPoint& pos = wxDefaultPosition,
			const wxSize& size = wxDefaultSize, long style = 0,
			const wxString& name = wxT("FGLCanvas"));
#endif // NEVER // new wxGLAttributes

    void Render();
    void InitGL();
    void Rotate(GLfloat deg);
    static GLfloat CalcRotateSpeed(unsigned long acceltime);
    static GLfloat CalcRotateAngle( unsigned long lasttime,
        unsigned long acceltime );
    void Action( long code, unsigned long lasttime,
        unsigned long acceltime );

	 void print();
	 void xfig();

	 void viewQuad();
    void ResetProjectionMode();
	 void setcurrent();
	 GLData* gldata() { return &m_gldata; }

//!! protected:
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnEraseBackground(wxEraseEvent& event);
    void OnMouse(wxMouseEvent& event);
    void OnChar(wxKeyEvent& event);

private:
    GLData m_gldata;
    bool   m_init;
    GLuint m_gllist;
    long   m_rleft;
    long   m_rright;

    static unsigned long  m_secbase;
    static int            m_TimeInitialized;
    static unsigned long  m_xsynct;
    static unsigned long  m_gsynct;

    long           m_Key;
    unsigned long  m_StartTime;
    unsigned long  m_LastTime;
    unsigned long  m_LastRedraw;


    DECLARE_EVENT_TABLE()
}; // FGLCanvas

#define ID_XFIG 20002

#define ID_TOOLBAR 10001
#define ID_AMP_SLIDER 10003
#define ID_HELP_HELP 10010
#define ID_HELP_ABOUT 10011
#define ID_EDIT_PREF 10012
#define ID_VIEW_DEFAULT 10013
#define ID_VIEW_QUAD 10014
#define ID_VIEW_X 10015
#define ID_VIEW_Y 10016
#define ID_VIEW_Z 10017
#define ID_VIEW_I 10018
#define ID_VIEW_COLOR 10020
#define ID_VIEW_COLOR_NONE 10026

// top toolbar
#define ID_XPAR 11000
#define ID_YPAR 11100
#define ID_ZPAR 11200
#define ID_DONE 11300

// Curves selection dialog
// #define ID_INCLUDE 11400
// #define ID_EXCLUDE 11500
#define ID_CURVES_SELECT 11600
// #define ID_CURVES_DONE 11700

// left side toolbar
// curve numbers start with ID_CURVE_NUMBER
#define ID_CURVE_NUMBER 12000
// #define ID_RADIO_VIEW 12001
#define ID_Z_VIEW 12011
#define ID_Y_VIEW 12012
#define ID_X_VIEW 12013
#define ID_DEFAULT_VIEW 12014

#define ID_XPARAM_SELECT 12002
#define ID_YPARAM_SELECT 12003
#define ID_ZPARAM_SELECT 12004
#define ID_PLOT 12005
#define ID_QUIT_BUTTON 12010

#endif // #if wxUSE_GLCANVAS

#endif // #ifndef VZ3D_H

