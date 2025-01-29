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
//  Flaps Animated Mode Visualizer
/////////////////////////////////////////////////////////////////////////////

#ifndef AMVIZ_H
#define AMVIZ_H

#include <complex>
#include <future>
#include <thread>
#include <vector>

#include "fio.h"
#include "fma.h"
#include "trace.h"
#include "trackball.h"
#include "wx/defs.h"
#include "wx/utils.h"
#include "wx/app.h"
#include "wx/button.h"
#include "wx/cursor.h"
#include "wx/menu.h"
#include "wx/dcclient.h"
#include "wx/wfstream.h"
#if wxUSE_ZLIB
#include "wx/zstream.h"
#endif
#include "wx/glcanvas.h"
#include "wx/wx.h"
#include "wx/timer.h"
#include "wx/gtk/print.h"

// Sizes of the various components
class AmvizSizes {
public:
	AmvizSizes();

	wxSize display;
	wxSize canvas;
	wxSize menu;
	wxSize toolbar;
	wxSize frame;
};

class myTimer : public wxTimer {
public:
	myTimer() : wxTimer() {}
	// override Notify to render one frame
	void Notify() override;
};

// Define a new frame type
class FGLCanvas;

class AmvizFrame: public wxFrame {
private:
	FGLCanvas *m_canvas{nullptr};
	wxSlider* ampSlider;
	//!! wxSlider* phaseSlider;
	wxSlider* freqSlider;
	wxSlider* colorSlider;
	wxButton* startButton;
	wxButton* undeformedButton;
	wxButton* CSButton;
public:
	// constructor
	AmvizFrame(wxWindow *frame, const wxString& title, const wxPoint& pos,
        const wxSize& size, long style = wxDEFAULT_FRAME_STYLE|wxCLOSE_BOX);

	 void createControls();
	 
	wxStatusBar* statusBar;
	// has this frame been Close()'d?
	bool closed{false};


    void OnMenuFileOpen(wxCommandEvent& event);
    void OnMenuFilePrint(wxCommandEvent& event);
    void OnMenuFileXfig(wxCommandEvent& event);
#ifdef NEVER // using funBind
    void OnMenuFileVideo(wxCommandEvent& event);
#endif // NEVER // using funBind
    void OnMenuFileExit(wxCommandEvent& event);
    void OnMenuEditPref(wxCommandEvent& event);
    void OnMenuHelpAbout(wxCommandEvent& event);
    void OnMenuHelpHelp(wxCommandEvent& event);
    void OnMenuViewDefault(wxCommandEvent& event);
    void OnMenuViewX(wxCommandEvent& event);
    void OnMenuViewY(wxCommandEvent& event);
    void OnMenuViewZ(wxCommandEvent& event);
    void OnMenuViewI(wxCommandEvent& event);
    void OnMenuViewQuad(wxCommandEvent& event);
    void OnMenuViewColorAbs(wxCommandEvent& event);
    void OnMenuViewColorVal(wxCommandEvent& event);
    void OnMenuViewColorPhase(wxCommandEvent& event);
    void OnMenuViewColorPhaseAbs(wxCommandEvent& event);
    void OnMenuViewColorEnergy(wxCommandEvent& event);
    void OnMenuViewColorNone(wxCommandEvent& event);

	 void OnAmpSlider(wxScrollEvent& event);
#ifdef NEVER // deprecate
	 void OnPhaseSlider(wxScrollEvent& event);
#endif // NEVER // deprecate
	 void OnFreqSlider(wxScrollEvent& event);
	 int omegat_incr{4};
	 void OnColorSlider(wxScrollEvent& event);

	 void OnTimer(wxTimerEvent& event);

	 void OnStartButton(wxCommandEvent& event);
	 void OnUndefButton(wxCommandEvent& event);
	 void OnCSButton(wxCommandEvent& event);
    void OnModeNumber(wxCommandEvent& event);

	 bool ShowToolTips();

	 void vibrate(void);

    void SetCanvas( FGLCanvas *canvas ) { m_canvas = canvas; }
    FGLCanvas* GetCanvas() { return m_canvas; }

	 bool isRunning();
	myTimer* timer{nullptr};

	 void print();
	 void printPS();

    DECLARE_EVENT_TABLE()
};  // AmvizFrame

// Amviz is an automatic* singleton holding the data necessary
// to visualize animated modes
// * the instance is an automatic, not a pointer
class Amviz {
	// private constructor
	Amviz();
	~Amviz();
	// singleton: delete the copy constructor & assignment operator
	Amviz(Amviz&) = delete;
	void operator=(const Amviz&) = delete;
public:

	// nodal_disp_: a collection of picked complex displacements
	// (dx,dy,dz) at each node: 3*nnode-complex vectors
	std::vector<std::vector<std::complex<double>>> nodal_disp_;
	wxArrayString picked_names_;
	int picked_point_;  // 0b index in nodal_disp_
	Fma* fma{nullptr};
	double coord_norm_;
	double amplitude_;   // displacement scale: slider-value*nominal_scale
	double nominal_scale_;
	std::string colormode_;
	int colorskew_;
	int drawundeformed_;
	// for drawing the coord axes:
	std::vector<float> origin_;
	int drawcs_;
	// XXX initialize all these above


	// if we are running standalone to view a uf file:
	bool standalone{false};

	// the main wxWidgets frame
	AmvizFrame* frame{nullptr};

	std::string error;

	double normalize_coord();

	// eigenvector from viz: check for non-empty periodically
	std::vector<std::complex<double>> ev;

	// access the one-and-only Amviz - an automatic so that it
	// is deleted on exit. This will load data according to Specs
	static Amviz& instance();
// #ifdef THREADS
	std::thread kit_thread;
	std::thread run_thread;
	std::thread ani_thread;
	// async futures
	std::future<void> kit_future;
	std::future<void> run_future;
	std::future<void> ani_future;
	// commands for animate
	bool quit{false};
	bool stop{false};
// #endif // THREADS
	bool startup();

	// POSIX message queue
	//!! Msgq* pmq{nullptr};

	// returns the (0b) disp_coord_ index of this ev
	int add(const std::vector<std::complex<double>>& ev, const std::string& name);
	// compute the coordinates of the displaced structure at time omegat:
	// XXX make this private, called from Amviz::render
	std::vector<double> coord_displacements(int omegat, int phase,
			std::vector<double>& colorvalue);
	int nnodes() { return instance().fma->nodes.size(); }
	int ngc() { return instance().fma->gct.size()/instance().fma->coords.size(); }
	std::vector<double>& coordinates() { return instance().fma->coords; }
	// compute the nodal displacements at time omegat+phase and feed
	// the segments to GL
	void render(int omegat, int phase);
	std::string colormode(const std::string& newmode="");
	int colorskew(int newskew=-1);
	int npicked();
	int picked_point(int newpoint=-1);
	wxArrayString& picked_names() { return picked_names_; }
	double amplitude(int newamp=-1);
	int draw_undeformed(int newund=-1);
	int draw_coord_sys(int newcs=-1);
};

// these may run in threads
void amviz_kit();
void amviz_run();

// default background: white
constexpr float BGRed{1.0};
constexpr float BGGreen{1.0};
constexpr float BGBlue{1.0};
constexpr float BGAlpha{0.0};

constexpr int Init_amp{3};
// these may need to be adjusted if the processor is overwhelmed
constexpr int Init_freq{2};			// initial setting for freq slider
constexpr int Default_period{20};	// milliseconds
constexpr int omegat_incr{4};		// increment omegat by this each timer step

// OpenGL view data
struct GLData
{
    bool initialized;           // has OpenGL been initialized?
    float beginx, beginy;       // position of mouse
	 float lastRadius;
	 std::vector<float> quat;  // orientation of object
    float zoom;                 // field of view in degrees
};

// Define a new application type
class AmvizApp: public wxApp
{
public:
    bool OnInit();
};


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
	// 2 constructors: the old interface with int* attlist is for WSL
    FGLCanvas(wxWindow *parent, wxWindowID id = wxID_ANY, const int* attlist=nullptr,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize, long style = 0,
        const wxString& name = wxT("FGLCanvas"));
    FGLCanvas(wxWindow *parent, wxGLAttributes& att, wxWindowID id = wxID_ANY,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize, long style = 0,
        const wxString& name = wxT("FGLCanvas"));

    ~FGLCanvas();

	void animate();
	void oneframe();
    void Render(int omegat, int phase);
    void InitGL();
    void Rotate(GLfloat deg);
    static GLfloat CalcRotateSpeed(unsigned long acceltime);
    static GLfloat CalcRotateAngle( unsigned long lasttime,
        unsigned long acceltime );
    void Action( long code, unsigned long lasttime,
        unsigned long acceltime );

	 void print();
	 void video();
	 void xfig();
	 void print90270();

	 void viewX(int,int);
	 void viewY(int,int);
	 void viewZ(int,int);
	 void viewQuad(int,int);
    void ResetProjectionMode();
	 //!! void setcurrent();
	 GLData* gldata() { return &m_gldata; }

	// one-and-only gl context: created in ctor
	wxGLContext* glcontext{nullptr};

protected:
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnEraseBackground(wxEraseEvent& event);
    void OnMouse(wxMouseEvent& event);
	void display_node(double x, double y);

private:

    GLData m_gldata;
    bool   m_init{false};
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

#define ID_TOOLBAR 10001
#define ID_AMP_SLIDER 10003
#define ID_FREQ_SLIDER 10005
#define ID_COLOR_SLIDER 10006
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
#define ID_VIEW_COLOR_ABS 10021
#define ID_VIEW_COLOR_VAL 10022
#define ID_VIEW_COLOR_PHASE 10023
#define ID_VIEW_COLOR_PHASE_ABS 10024
#define ID_VIEW_COLOR_ENERGY 10025
#define ID_VIEW_COLOR_NONE 10026
#define ID_START_BUTTON 10100
#define ID_UNDEF_BUTTON 10101
#define ID_CS_BUTTON 10102


// mode numbers start with ID_MODE_NUMBER
#define ID_MODE_NUMBER 11000
// for the wxTimer
#define ID_TIMER 12000

#endif // #if wxUSE_GLCANVAS

// bool startup();
void ani();

#endif // #ifndef AMVIZ_H

