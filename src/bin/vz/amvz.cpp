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
// Animated Modes VisualiZation
/////////////////////////////////////////////////////////////////////////////

// For compilers that support precompilation, includes "wx.h".

#include <cerrno>
#include <ctime>
#include <exception>
#include <stdexcept>
#include <iosfwd>

#include "wx/string.h"
#include "wx/prntbase.h"

#ifndef WX_PRECOMP
#include "wx/wx.h"
#endif

#if !wxUSE_GLCANVAS
    #error "OpenGL required: set wxUSE_GLCANVAS to 1 and rebuild the library"
#endif
#include <GL/glu.h>
// #include <GL/glut.h>

#undef Complex    // X11 macro
#include "config.h"
#include "amdata.h"
#include "amvz.h"
#include "exim.h"
#include "gl2ps.h"
#include "matrix.h"
#include "prefdialog.h"
#include "swf.h"
#include "trace.h"

class Quaternion {
public:
	string id;
	vector<float> quat;
	Quaternion(const string& ident, vector<float> q) : id(ident), quat(q) {}
};

vector<float> Default_quat;
vector<float> airplane_quat{0.396197, -0.133608, -0.491376, 0.76402};
vector<float> Isoquat{0.180731, 0.356605, 0.827727, 0.393747};
vector<float> Xquat{0.5, 0.5, 0.5, 0.5};
vector<float> Yquat{0.0, 0.707106781, 0.707106781, 0.0};
vector<float> Zquat{0.0, 0.0, 1.0, 0.0};

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

int Omegat = 0;
int Phase = 0;

wxChoice* ModeMenu{nullptr};
wxArrayString ModeStrings;

extern double XParValue;
extern double YParValue;
extern string ModeName;

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
amvz() {
	T_(Trace trc(1,"amvz");)
	int argc = 1;
	char* argv[3];
	char prog[8];
	strcpy(prog, "amvz");
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

MyFrame* Frame{nullptr};

bool
MyApp::
OnInit() {
	T_(Trace trc(1,"OnInit");)
    // Create the main frame window...
	Frame = new MyFrame(nullptr, wxT("Animated Mode Visualization"),
			 wxDefaultPosition, wxSize(WindowSizeX,WindowSizeY));
	// ... and display it
    Frame->Show(true);

    return true;
}

IMPLEMENT_APP_NO_MAIN(MyApp)

// MyFrame

BEGIN_EVENT_TABLE(MyFrame, wxFrame)
    EVT_MENU(wxID_OPEN, MyFrame::OnMenuFileOpen)
    EVT_MENU(wxID_PRINT, MyFrame::OnMenuFilePrint)
    EVT_MENU(wxID_XFIG, MyFrame::OnMenuFileXfig)
    EVT_MENU(wxID_VIDEO, MyFrame::OnMenuFileVideo)
    EVT_MENU(wxID_EXIT, MyFrame::OnMenuFileExit)
    EVT_MENU(ID_EDIT_PREF, MyFrame::OnMenuEditPref)
    EVT_CHOICE(ID_MODE_NUMBER, MyFrame::OnMenuMode)
    EVT_MENU(ID_HELP_ABOUT, MyFrame::OnMenuHelpAbout)
    EVT_MENU(ID_HELP_HELP, MyFrame::OnMenuHelpHelp)
    EVT_MENU(ID_VIEW_DEFAULT, MyFrame::OnMenuViewDefault)
    EVT_MENU(ID_VIEW_X, MyFrame::OnMenuViewX)
    EVT_MENU(ID_VIEW_Y, MyFrame::OnMenuViewY)
    EVT_MENU(ID_VIEW_Z, MyFrame::OnMenuViewZ)
    EVT_MENU(ID_VIEW_I, MyFrame::OnMenuViewI)
    EVT_MENU(ID_VIEW_QUAD, MyFrame::OnMenuViewQuad)
    EVT_MENU(ID_VIEW_COLOR_ABS, MyFrame::OnMenuViewColorAbs)
    EVT_MENU(ID_VIEW_COLOR_VAL, MyFrame::OnMenuViewColorVal)
    EVT_MENU(ID_VIEW_COLOR_PHASE, MyFrame::OnMenuViewColorPhase)
    EVT_MENU(ID_VIEW_COLOR_PHASE_ABS, MyFrame::OnMenuViewColorPhaseAbs)
    EVT_MENU(ID_VIEW_COLOR_ENERGY, MyFrame::OnMenuViewColorEnergy)
    EVT_MENU(ID_VIEW_COLOR_NONE, MyFrame::OnMenuViewColorNone)

	EVT_COMMAND_SCROLL(ID_AMP_SLIDER, MyFrame::OnAmpSlider)
	EVT_COMMAND_SCROLL(ID_PHASE_SLIDER, MyFrame::OnPhaseSlider)
	EVT_COMMAND_SCROLL(ID_FREQ_SLIDER, MyFrame::OnFreqSlider)
	EVT_COMMAND_SCROLL(ID_COLOR_SLIDER, MyFrame::OnColorSlider)

	EVT_TIMER(ID_TIMER, MyFrame::OnTimer)

	EVT_BUTTON(ID_START_BUTTON, MyFrame::OnStartButton)
	EVT_BUTTON(ID_UNDEF_BUTTON, MyFrame::OnUndefButton)
	EVT_BUTTON(ID_CS_BUTTON, MyFrame::OnCSButton)
END_EVENT_TABLE()

// MyFrame constructor
#include "amvz.xpm"
MyFrame::MyFrame(wxFrame *frame, const wxString& title, const wxPoint& pos,
    const wxSize& size, long style)
    : wxFrame(frame, wxID_ANY, title, pos, size, style), m_timer(this, ID_TIMER) {

	// an icon for the app bar
	SetIcon(wxIcon(amvz_xpm));

	 // add menubar, toolbar + controls
	 createControls();

	 int attlist[3];
	 int i = 0;
	 attlist[i++] = WX_GL_DOUBLEBUFFER;
	 attlist[i++] = WX_GL_RGBA;
	 attlist[i++] = 0;
    m_canvas = new FGLCanvas(this, wxID_ANY, wxDefaultPosition,
        // wxSize(900, 600),
		  GetClientSize(),
		  wxSUNKEN_BORDER, wxT("Amvz"), attlist);

	int period = Default_period/Init_freq;
	m_timer.Start(period);
}

bool
MyFrame::
isRunning() {
	return (m_timer.IsRunning());
}

void
MyFrame::
createControls() {

	// Top-level sizer
	//! wxBoxSizer* topsizer = new wxBoxSizer(wxVERTICAL);

    // Make the "File" menu
    wxMenu *fileMenu = new wxMenu;
    //fileMenu->Append(wxID_OPEN, wxT("&Open..."));
    //fileMenu->AppendSeparator();
    fileMenu->Append(wxID_PRINT, wxT("&Print\tALT-P"));
    fileMenu->Append(wxID_XFIG, wxT("&Xfig\tALT-X"));
    fileMenu->Append(wxID_VIDEO, wxT("&Video\tALT-v"));
    fileMenu->AppendSeparator();
    fileMenu->Append(wxID_EXIT, wxT("E&xit\tq"));
	 // make the "Edit" menu
    wxMenu *editMenu = new wxMenu;
    editMenu->Append(ID_EDIT_PREF, wxT("&Preferences\tALT-E"));
	 // make the "View" menu
    wxMenu *viewMenu = new wxMenu;
    viewMenu->Append(ID_VIEW_DEFAULT, wxT("&Default\td"));
#ifdef NEVER // disable quad, use quats
    viewMenu->Append(ID_VIEW_QUAD, wxT("&Quad-View\t4"));
#endif // NEVER : disable quad, use quats
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
    helpMenu->Append(ID_HELP_HELP, wxT("&Help..."));
    helpMenu->Append(ID_HELP_ABOUT, wxT("&About..."));

	 // Add these menus to the menu bar
    wxMenuBar *menuBar = new wxMenuBar;
    menuBar->Append(fileMenu, wxT("&File"));
    menuBar->Append(editMenu, wxT("&Edit"));
    menuBar->Append(viewMenu, wxT("&View"));
    menuBar->Append(helpMenu, wxT("&Help"));
    SetMenuBar(menuBar);

	 int mindisp = 0;
	 int maxdisp = 10;
	 int value = 3;
	 wxToolBar* toolbar = CreateToolBar(wxTB_VERTICAL|wxTB_TEXT, ID_TOOLBAR);
	// toolbar->SetToolPacking(5);
	// toolbar->SetToolSeparation(20);
	toolbar->SetMargins(2,4);
	//! wxBoxSizer* toolbarSizer = new wxBoxSizer(wxVERTICAL);

	// amplitude slider
	ampSlider = new wxSlider(toolbar, ID_AMP_SLIDER, Init_amp, mindisp, maxdisp,
			 wxDefaultPosition, wxSize(-1,100), wxSL_VERTICAL|wxSL_LABELS);
	ampSlider->SetToolTip(wxT("amplitude"));
	toolbar->AddControl(ampSlider);
	 // toolbarSizer->Add(ampSlider);

	mindisp = 0;
	maxdisp = 360;
	value = 0;
	phaseSlider = new wxSlider(toolbar, ID_PHASE_SLIDER, value, mindisp, maxdisp,
			 wxDefaultPosition, wxSize(-1,100), wxSL_VERTICAL|wxSL_LABELS);
	phaseSlider->SetToolTip(wxT("phase"));
	toolbar->AddSeparator();
	toolbar->AddControl(phaseSlider);
	 // toolbarSizer->Add(phaseSlider);

	mindisp = 0;
	maxdisp = 10;
	freqSlider = new wxSlider(toolbar, ID_FREQ_SLIDER, Init_freq, mindisp, maxdisp,
			 wxDefaultPosition, wxSize(-1,100), wxSL_VERTICAL|wxSL_LABELS);
	freqSlider->SetToolTip(wxT("frequency"));
	toolbar->AddSeparator();
	toolbar->AddControl(freqSlider);
	 // toolbarSizer->Add(freqSlider);

	mindisp = 0;
	maxdisp = 10;
	value = 5;
	colorSlider = new wxSlider(toolbar, ID_COLOR_SLIDER, value, mindisp, maxdisp,
			 wxDefaultPosition, wxSize(-1,100), wxSL_VERTICAL|wxSL_LABELS);
	colorSlider->SetToolTip(wxT("color gradient"));
	toolbar->AddSeparator();
	toolbar->AddControl(colorSlider);
	 // toolbarSizer->Add(colorSlider);

	startButton = new wxButton(toolbar, ID_START_BUTTON, wxT("stop"),
			 wxDefaultPosition, wxSize(45,30));
	startButton->SetToolTip(wxT("start/stop oscillations"));
	toolbar->AddSeparator();
	toolbar->AddControl(startButton);
	 // toolbarSizer->Add(startButton);

	 undeformedButton = new wxButton(toolbar, ID_UNDEF_BUTTON,
			 wxT("&undef"),
			 wxDefaultPosition, wxSize(50,30));
	undeformedButton->SetToolTip(wxT("toggle undeformed grid"));
	toolbar->AddSeparator();
	toolbar->AddControl(undeformedButton);
	 // toolbarSizer->Add(undeformedButton);

	// button for toggling the coordinate system
	 CSButton = new wxButton(toolbar, ID_CS_BUTTON,
			 wxT("&coord sys"),
			 wxDefaultPosition, wxSize(50,30));
	undeformedButton->SetToolTip(wxT("toggle coordinate system"));
	toolbar->AddSeparator();
	toolbar->AddControl(CSButton);

	wxArrayString empty;
	ModeMenu = new wxChoice(toolbar, ID_MODE_NUMBER,
			wxDefaultPosition, wxDefaultSize, Amdata::instance()->picked_names());
	toolbar->AddSeparator();
	toolbar->AddControl(ModeMenu);
	// toolbarSizer->Add(ModeMenu, wxSizerFlags(0).Center());

	// toolbar->SetSizer(toolbarSizer);
	toolbar->Realize();
	this->SetToolBar(toolbar);

	// topsizer->Add(toolbarSizer);

	 statusBar = new wxStatusBar(this, wxID_ANY, wxST_SIZEGRIP);
	 this->SetStatusBar(statusBar);
	 statusBar->SetStatusText(wxT("Ready"), 0);

	 // SetSizer(topsizer);
} // MyFrame::createControls()

bool
MyFrame::
ShowToolTips() { return true; }

void
MyFrame::
OnAmpSlider(wxScrollEvent& event) {
	if (this->ampSlider) {
		//! DispScale = NomDispScale*this->ampSlider->GetValue();
		int newscale = this->ampSlider->GetValue();
		// double newscale = NomDispScale*this->ampSlider->GetValue();
		Amdata::amplitude(newscale);
	}
}

void
MyFrame::
OnPhaseSlider(wxScrollEvent& event) {
	if (this->ampSlider) {
		// stop the oscillations
		if (m_timer.IsRunning()) {
			wxCommandEvent ev;
			OnStartButton(ev);
		}
		Omegat = 0;
		Phase = this->phaseSlider->GetValue();
		//Refresh(false);
		m_canvas->Render(Omegat, Phase);
	}
}

void
MyFrame::
OnFreqSlider(wxScrollEvent& event) {
// Respond to movement of the freq button; if the new
// frequency is zero just stop the timer
	int freq = this->freqSlider->GetValue();
	m_timer.Stop();
	if (freq > 0) {
		int period = Default_period/freq;
		m_timer.Start(period);
	}
}

void
MyFrame::
OnColorSlider(wxScrollEvent& event) {
	if (this->colorSlider) {
		//! ColorSkew = (float)this->colorSlider->GetValue()/5.0;
		Amdata::colorskew(this->colorSlider->GetValue());
	}
}

void
MyFrame::
OnStartButton(wxCommandEvent& event) {
	static int interval{Default_period};
	if (m_timer.IsRunning()) {
		interval = m_timer.GetInterval();
		m_timer.Stop();
		startButton->SetLabel(wxT("start"));
	} else {
		m_timer.Start(interval);
		startButton->SetLabel(wxT("stop"));
	}
}

void
MyFrame::
OnUndefButton(wxCommandEvent& event) {
	//! DrawUndeformed = !DrawUndeformed;
	int current = Amdata::draw_undeformed();
	if (current == 0)
		Amdata::draw_undeformed(1);
	else
		Amdata::draw_undeformed(0);
}


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
	int period = m_timer.GetInterval();
	// stop the oscillations
	m_timer.Stop();
	// call the FGLCanvas xfig creator
	m_canvas->xfig();
	// and re-start
	m_timer.Start(period);
}

// File|Video... command
void
MyFrame::
OnMenuFileVideo( wxCommandEvent& WXUNUSED(event) ) {
/*
 * Create a flash animation file using the GL feedback
 * mechanism
 */
	//wxMessageBox(wxT("video not yet implemented"),
	//			wxT("Video"), wxOK);
	//return;

	int period = m_timer.GetInterval();
	// stop the oscillations
	m_timer.Stop();
	// call the FGLCanvas video creator
	m_canvas->video();
	// and re-start
	m_timer.Start(period);
}

wxPrintDialogData g_printDialogData;

void
MyFrame::
print() {
	int period = m_timer.GetInterval();
	// stop the oscillations
	m_timer.Stop();
	
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
	m_timer.Start(period);
}

void
MyFrame::
printPS() {
	int period = m_timer.GetInterval();
	// stop the oscillations
	m_timer.Stop();
	// call the FGLCanvas printer
	m_canvas->print();
	// and re-start
	m_timer.Start(period);
}

// File|Exit command or simply type "q"
void
MyFrame::
OnMenuFileExit( wxCommandEvent& WXUNUSED(event) ) {
	// write out the current quaternion to use as the initial
	// position in FGlCanvas constructor
	// cerr << "current quaternion: " << m_canvas->gldata()->quat << endl;
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
	int period = m_timer.GetInterval();
	m_timer.Stop();
	m_canvas->gldata()->quat = Default_quat;
	m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewX(wxCommandEvent& WXUNUSED(event)) {
	int period = m_timer.GetInterval();
	m_timer.Stop();
	m_canvas->gldata()->quat = Xquat;
	m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewY(wxCommandEvent& WXUNUSED(event)) {
	int period = m_timer.GetInterval();
	m_timer.Stop();
	m_canvas->gldata()->quat = Yquat;
	m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewZ(wxCommandEvent& WXUNUSED(event)) {
	int period = m_timer.GetInterval();
	m_timer.Stop();
	m_canvas->gldata()->quat = Zquat;
	m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewI(wxCommandEvent& WXUNUSED(event)) {
	int period = m_timer.GetInterval();
	m_timer.Stop();
	m_canvas->gldata()->quat = Isoquat;
	m_timer.Start(period);
	m_canvas->Refresh(false);
}

void
MyFrame::
OnMenuViewQuad(wxCommandEvent& WXUNUSED(event)) {
	if (!ViewQuad) {
		ViewQuad = true;
		m_canvas->viewQuad(Omegat, Phase);
	}
}

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

void
MyFrame::
OnMenuMode( wxCommandEvent& event ) {
	// the items in the choice box are in reverse order
	// (most recent first) so GetSelection 0 is the last
	// column of Disp
	Amdata::picked_point(event.GetSelection());
}

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
			 " - File|Video creates an animated flash (.swf) file\n"
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

// FGLCanvas

BEGIN_EVENT_TABLE(FGLCanvas, wxGLCanvas)
    EVT_SIZE(FGLCanvas::OnSize)
    EVT_PAINT(FGLCanvas::OnPaint)
    EVT_ERASE_BACKGROUND(FGLCanvas::OnEraseBackground)
    EVT_MOUSE_EVENTS(FGLCanvas::OnMouse)
END_EVENT_TABLE()

void
MyFrame::
OnTimer(wxTimerEvent& event) {
	static string last_colormode;
	static int omegat = 0;
	int phase = 0;
	omegat += omegat_incr;
	if (omegat >= 360)
		omegat = 0;

#ifdef NEVER // disable quad, use quats
	if (ViewQuad) {
		//Use a very low frequency so we don't overwhelm the processor
		// XXX move Stop & Start to where the button is pressed
		m_timer.Stop();
		m_timer.Start(Default_period);
		Frame->GetCanvas()->viewQuad(omegat, phase);
	} else if (ViewX) {
		Frame->GetCanvas()->viewX(omegat, phase);
	} else if (ViewY) {
		Frame->GetCanvas()->viewY(omegat, phase);
	} else if (ViewZ) {
		Frame->GetCanvas()->viewZ(omegat, phase);
	} else {
		Frame->GetCanvas()->Render(omegat, phase);
	}
#else // NEVER : disable quad, use quats
		Frame->GetCanvas()->Render(omegat, phase);
	// XXX try this
	m_timer.Start();
#endif // NEVER : disable quad, use quats
}


FGLCanvas::FGLCanvas(wxWindow *parent, wxWindowID id, const wxPoint& pos,
	const wxSize& size, long style, const wxString& name, int* attlist)
    : wxGLCanvas(parent, id, attlist, pos, size, style|wxFULL_REPAINT_ON_RESIZE, name, wxNullPalette) {
// FGlCanvas constructor

   m_gldata.initialized = false;
	m_init = false;

   m_gldata.beginx = 0.0f;
   m_gldata.beginy = 0.0f;
   m_gldata.zoom   = 1.0f;
   // Set the default quat, initialize view matrix with it
	Default_quat = airplane_quat;
	m_gldata.quat = Default_quat;
	float qn = 1.0/blas_snrm2(4, &m_gldata.quat[0], 1);
	if (abs(qn-1.0) > 8.0*std::numeric_limits<float>::epsilon()) {
		blas_scal(4, qn, &m_gldata.quat[0], 1);
		cerr << "inital quat scaled by " << qn << ": " << m_gldata.quat << endl;
	}
	// Create a GL context (new with 3.0)
	// wxGLContext* context = new wxGLContext(this, 0);
	// context->SetCurrent(*this);
}

FGLCanvas::~FGLCanvas() { }


void
FGLCanvas::
setcurrent() {
	static wxGLContext* context = 0;
	if (!context) {
		context = new wxGLContext(this, 0);
	}
	bool rval = context->SetCurrent(*this);
	if (!rval) {
		cerr << "SetCurrent failed\n";
	}
}

void
FGLCanvas::
Render(int omegat, int phase) {

	this->setcurrent();
	Amdata::instance()->render(omegat, phase);
	glFlush();
	wxGLCanvas::SwapBuffers();

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
display_node(double x, double y) {
	T_(Trace trc(1,"display_node");)

	if (Frame->isRunning()) {
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
	vector<double>& coord = Amdata::coordinates();
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
}

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
	Render(Omegat, Phase);
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

void
FGLCanvas::
print90270() {
#ifdef NEVER // disabled: needs work to convert to Amdata
// XXX why doesn't this just use feedback?
	size_t i, j;
	float normd;
	float maxnormd;

	// gdImagePtr im = gdImageCreate(1200,900);
	// int black = gdImageColorAllocate(im, 0, 0, 0);
	// int white = gdImageColorAllocate(im, 255, 255, 255);

    /* clear color and depth buffers */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	int omegat = 90;
	int phase = 0;

	//! vector<double> coord = Coordinates->data();
	vector<double>& coord = Amdata::coordinates();
	vector<float> colorvalue(coord.size(), 0.0);

	vector<double> dispcoord = omegat_coord(Disp[DispColumn], coord,
		omegat, phase, normd, maxnormd, &colorvalue[0]);

	// double d2p = 500.0;
	float red = 0.0;
	float green = 0.0;
	float blue = 0.0;

	float maxd;

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

	red = 1.0;
	glColor3f(red, green, blue);

	if (!Connectivity.empty()) {
		glBegin(GL_LINES);
		for (i=0; i<Connectivity.size(); i++) {
			j = Connectivity[i];
			float x = dispcoord[j];
			float y = dispcoord[j+1];
			float z = dispcoord[j+2];
			maxd = std::max(std::max(abs(x), abs(y)), abs(z));
			red = maxd/normd;
			blue = 1.0 - red;
			glColor3f(red, green, blue);
			glVertex3f(x, y, z);
			i++;
			j = Connectivity[i];
			x = dispcoord[j];
			y = dispcoord[j+1];
			z = dispcoord[j+2];
			glVertex3f(x, y, z);
		}
		glEnd();
	} else {
		glBegin(GL_POINTS);
		for (i=0; i<dispcoord.size(); i++) {
			glVertex3f(dispcoord[i],dispcoord[i+1],dispcoord[i+2]);
		}
		glEnd();
	}
	// At 270 degrees:
	omegat = 270;
	dispcoord = omegat_coord(Disp[DispColumn], coord,
		omegat, phase, normd, maxnormd, &colorvalue[0]);

	red = 0.0;
	blue = 1.0;
	glColor3f(red, green, blue);
	if (!Connectivity.empty()) {
		glBegin(GL_LINES);
		for (i=0; i<Connectivity.size(); i++) {
			j = Connectivity[i];
			float x = dispcoord[j];
			float y = dispcoord[j+1];
			float z = dispcoord[j+2];
			maxd = std::max(std::max(abs(x), abs(y)), abs(z));
			red = maxd/normd;
			blue = 1.0 - red;
			glColor3f(red, green, blue);
			glVertex3f(x, y, z);
			i++;
			j = Connectivity[i];
			x = dispcoord[j];
			y = dispcoord[j+1];
			z = dispcoord[j+2];
			glVertex3f(x, y, z);
			// gdImageLine(im, x1, y1, x2, y2, white);
		}
		glEnd();
	} else {
		glBegin(GL_POINTS);
		for (i=0; i<dispcoord.size(); i += 3) {
			glVertex3f(dispcoord[i],dispcoord[i+1],dispcoord[i+2]);
		}
		glEnd();
	}
	gl2psEndPage();
#ifdef NEVER
    glFlush();
    SwapBuffers();

	 char* data = new char[2*1200*900];
	 // float* data = new float[8*im->sx*im->sy];
	 // glReadPixels(0, 0, im->sx, im->sy, GL_RGBA, GL_FLOAT, data);
	 glReadPixels(0, 0, im->sx, im->sy, GL_RGB, GL_BYTE, data);
	 wxBitmap bm(data, wxBITMAP_TYPE_XBM, 1200,900);
	 //wxImage* image = new wxImage(&data);
	 wxImage image(bm.ConvertToImage());
	 wxPNGHandler pnghandler;
	 wxFileOutputStream stream(wxT("amvz.png"));
	 // pnghandler.SaveFile(&image, stream);
	 //memcpy (*(im->pixels), data, 1200*900);
	 for (i=0; i<im->sx; i++)
		 for (j=0; j<im->sy; j++)
			 im->pixels[j][i] = data[IJ(i,j,1200)];
			 //im->tpixels[j][i] = data[IJ(i,j,1200)];
	 delete[] data;
	 // write png image file:
	wxMessageBox(wxT("is this ok?"));

	FILE* fp = fopen("amvz.png", "w");
	gdImagePng(im, fp);
	gdImageDestroy(im);
#endif // NEVER
#endif // NEVER : disabled: needs work to convert to Amdata
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
	int omegat = 0;
	int phase = 0;
	// render into a gl feedback buffer
	int buffersize = 1000000;
	wxString wxpath = wxGetTextFromUser(wxT("Please enter a filename for xfig output "
							  "(include a *.fig extenstion)"), 
							  wxT("xfig filename"), wxT("amvz.fig"));
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
		Render(omegat, phase);
		glRenderMode(GL_RENDER);
		xfig_create(pathname, feedbuf, append, style, width, color);
		append = true;
	}
}

void
FGLCanvas::
video() {

	int phase = 0;
	// render into a gl feedback buffer
	// int buffersize = 1000000;
	// swf flashvid;
	// wxString pathname = wxGetTextFromUser(wxT("Please enter a filename for video output "
	// 										  "(include an *.swf extenstion)"), 
	// 									  wxT("Shockwave filename"), wxT("amvz.swf"));

	// get window size
	wxSize siz = GetClientSize();
	GLsizei width = siz.GetWidth();
	GLsizei height = siz.GetHeight();
	// T_(trc.dprint("window width ", siz.GetWidth(),", height ",siz.GetHeight());)

	// start ffmpeg
	int framerate{60};
	string resolution{vastr(width,"x",height)};
	string file{"amvz.mp4"};
	// string cmd = vastr("ffmpeg -r ",framerate," -f rawvideo -pix_fmt rgba -s ",resolution,
	// 	" -i - -threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip ",file);
	string cmd = vastr("ffmpeg -r ",framerate," -f rawvideo -pix_fmt rgba -s ",resolution,
	 	" -i - -threads 0 -preset fast -y -crf 21 -vf vflip ",file);

	cout << cmd << endl;
		
	FILE* ffmpeg = popen(cmd.c_str(), "w");
	if (ffmpeg == nullptr)
		throw std::system_error(errno, std::generic_category(),"cannot start ffmpeg");

	glRenderMode(GL_RENDER);

	size_t bufsize = width*height*sizeof(int);
	vector<int> buffer(width*height);
	// 360 degrees for one cycle and 40 frames total.
	int ncycles{4};
	for (int j=0; j<ncycles; j++) {
		for (size_t i=0; i < 360; i+=9) {
			Render(i, phase);
			wxGLCanvas::SwapBuffers();
			glReadPixels(0,0,width,height,GL_RGBA, GL_UNSIGNED_BYTE, &buffer[0]);
			fwrite(&buffer[0], bufsize, 1, ffmpeg);
			struct timespec reg, rem;
			reg.tv_sec = 0;
			reg.tv_nsec = 1000000;
			nanosleep(&reg, &rem);
		}
	}

	// close file
	// flashvid.closeBuffer();
	pclose(ffmpeg);

	// write File
	// flashvid.writeBuffer( std::string(pathname.mb_str()) );

	// finished filling feedbuf, return to GL_RENDER mode
	glRenderMode(GL_RENDER);
}

void
FGLCanvas::
viewX(int omegat, int phase) {
	int width, height;
	GetClientSize(&width, &height);
	ViewQuad = false;
#ifdef NEVER // deprecated
	if ( GetContext() ) {
		// direct gl commands to this FGLCanvas
		wxGLCanvas::SetCurrent();
		glDisable(GL_SCISSOR_TEST);
		bottom_left;
		left;
		Render(omegat, phase);
	}
#endif // NEVER : deprecated
}

void
FGLCanvas::
viewY(int omegat, int phase) {
	int width, height;
	GetClientSize(&width, &height);
	ViewQuad = false;
#ifdef NEVER // deprecated
	if ( GetContext() ) {
		// direct gl commands to this FGLCanvas
		wxGLCanvas::SetCurrent();
		glDisable(GL_SCISSOR_TEST);
		bottom_left;
		bottom;
		Render(omegat, phase);
	}
#endif // NEVER : deprecated
}

void
FGLCanvas::
viewZ(int omegat, int phase) {
	int width, height;
	GetClientSize(&width, &height);
	ViewQuad = false;
#ifdef NEVER // deprecated
	if ( GetContext() ) {
		// direct gl commands to this FGLCanvas
		wxGLCanvas::SetCurrent();
		glDisable(GL_SCISSOR_TEST);
		bottom_left;
		front;
		Render(omegat, phase);
	}
#endif // NEVER : deprecated
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
}

#ifdef STANDALONE

#include <cassert>

double XParValue;
double YParValue;
string ModeName;

int
main(int argc, char** argv) {

	if (argc < 2) {
		cerr << "Usage: amvz uf_file\n";
		exit(1);
	}

	try {
		string path(argv[1]);
		// if a file with a .uf extension exists, prefer it
		if (rsubstr(path, 3) != ".uf") {
			string ufpath = path + ".uf";
			if (access(ufpath.c_str(), R_OK) != -1)
				path = ufpath;
		}
		Amdata::initialize(path);
		amvz();
	} catch (runtime_error& s) {
		flaps::error("caught exception: ",s.what());
	} catch (const std::exception& t) {
		flaps::error("caught exception: ",t.what());
	}
}
#endif // STANDALONE
