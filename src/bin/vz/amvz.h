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

#ifndef AMVZ_H
#define AMVZ_H

#include <vector>

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
#include "wx/gtk/print.h"

#include "trackball.h"

// main entry point to amvz
void amvz();

// default background: white
constexpr float BGRed{1.0};
constexpr float BGGreen{1.0};
constexpr float BGBlue{1.0};
constexpr float BGAlpha{0.0};

constexpr int Init_amp{3};
// these may need to be adjusted if the processor is overwhelmed
constexpr int Init_freq{2};			// initial setting for freq slider
constexpr int Default_period{100};	// milliseconds
constexpr int omegat_incr{4};		// increment omegat by this each timer step

// OpenGL view data
struct GLData
{
    bool initialized;           // has OpenGL been initialized?
    float beginx, beginy;       // position of mouse
	 float lastRadius;
    float lastX, lastY;       // position of mouse XXX remove
	 std::vector<float> quat;  // orientation of object
    float zoom;                 // field of view in degrees
};

// Define a new application type
class MyApp: public wxApp
{
public:
    bool OnInit();
};

// Define a new frame type
class FGLCanvas;
//!! class MyTimer;

class MyFrame: public wxFrame {
private:
	FGLCanvas *m_canvas;
	//!! MyTimer* m_timer;
	wxTimer m_timer;
	wxSlider* ampSlider;
	wxSlider* phaseSlider;
	wxSlider* freqSlider;
	wxSlider* colorSlider;
	wxButton* startButton;
	wxButton* undeformedButton;
	wxButton* CSButton;
public:
    MyFrame(wxFrame *frame, const wxString& title, const wxPoint& pos,
        const wxSize& size, long style = wxDEFAULT_FRAME_STYLE);

	 void createControls();
	 
	wxStatusBar* statusBar;

    void OnMenuFileOpen(wxCommandEvent& event);
    void OnMenuFilePrint(wxCommandEvent& event);
    void OnMenuFileXfig(wxCommandEvent& event);
    void OnMenuFileVideo(wxCommandEvent& event);
    void OnMenuFileExit(wxCommandEvent& event);
    void OnMenuEditPref(wxCommandEvent& event);
    void OnMenuMode(wxCommandEvent& event);
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
	 void OnPhaseSlider(wxScrollEvent& event);
	 void OnFreqSlider(wxScrollEvent& event);
	 void OnColorSlider(wxScrollEvent& event);

	 void OnTimer(wxTimerEvent& event);

	 void OnStartButton(wxCommandEvent& event);
	 void OnUndefButton(wxCommandEvent& event);
	 void OnCSButton(wxCommandEvent& event);

	 bool ShowToolTips();

	 void vibrate(void);

    void SetCanvas( FGLCanvas *canvas ) { m_canvas = canvas; }
    FGLCanvas *GetCanvas() { return m_canvas; }

	 bool isRunning();

	 void print();
	 void printPS();

    DECLARE_EVENT_TABLE()
};  // MyFrame

// timer
#ifdef NEVER // new timer scheme
class MyTimer : public wxTimer {
	bool running;
public:
	MyTimer(MyFrame* f);
	MyFrame* myframe;
	void Notify();
	void start();
	void stop();
	bool isRunning() { return running; }
};  // MyTimer
#endif // NEVER : new timer scheme

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
    FGLCanvas(wxWindow *parent, wxWindowID id = wxID_ANY,
        const wxPoint& pos = wxDefaultPosition,
        const wxSize& size = wxDefaultSize, long style = 0,
        const wxString& name = wxT("FGLCanvas"), int* attlist = 0);

    ~FGLCanvas();

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
	 void setcurrent();
	 GLData* gldata() { return &m_gldata; }

protected:
    void OnPaint(wxPaintEvent& event);
    void OnSize(wxSizeEvent& event);
    void OnEraseBackground(wxEraseEvent& event);
    void OnMouse(wxMouseEvent& event);
	void display_node(double x, double y);

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

#define wxID_VIDEO 20001
#define wxID_XFIG 20002

#define ID_TOOLBAR 10001
#define ID_AMP_SLIDER 10003
#define ID_PHASE_SLIDER 10004
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

#endif // #ifndef AMVZ_H

