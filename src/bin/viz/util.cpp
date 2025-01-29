//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <iostream>
#include <GL/gl.h>
#include <png.h>
#include <sstream>
#include <sys/stat.h>	// is_wsl
#include <wx/gdicmn.h>	// colourdatabase
#include <X11/Xlib.h>

#include "fio.h"		// File
#include "trace.h"
#include "util.h"
#include "wx/wx.h"
#include "wx/textctrl.h"

using namespace std;

void
errormsg(const string& msg, bool quit) {
	wxMessageBox(msg, wxT("Error"), wxICON_INFORMATION);
	if (quit)
		exit(1);
}

ostream&
operator<<(ostream& s, const wxPoint& t) {
	s << "(" << t.x << ", " << t.y << ")";
	return s;
}

ostream&
operator<<(ostream& s, const wxSize& t) {
	s << "(" << t.x << ", " << t.y << ")";
	return s;
}

ostream&
operator<<(ostream& s, const wxBoxSizer& t) {
	s << sizer2string(t);
	return s;
}

ostream&
operator<<(ostream& s, const wxString& t) {
	s << t.ToStdString();
	return s;
}

string
buildinfo() {
	wxString v = wxVERSION_STRING;
	string rval = vastr("Built with ",v,"\n", wxGetOsDescription());
	return rval;
}

wxColour
getcolor(const string& name) {
	wxColourDatabase* colordb{nullptr};
	if (colordb == nullptr)
		colordb = new wxColourDatabase();
	wxColour rval = colordb->Find(stringUpper(name));
	if (!rval.IsOk()) {
		errormsg(vastr("invalid color: \"",name,"\""), false);
		rval = colordb->Find("BLACK");
	}
	return rval;
}

string
OGLVersion() {
	GLint majv;
	glGetIntegerv(GL_MAJOR_VERSION, &majv);
	int majorv = majv;
	GLint minv;
	glGetIntegerv(GL_MINOR_VERSION, &minv);
	int minorv = minv;
	string rval = vastr("OGL version supported: ", majorv, ".", minorv);
	return rval;
}

bool isWayland() {
	Trace trc(2,"isWayland");
	// only works under MSW:
	// wxDisplay display; // Get the primary display
	// trc.dprint("display name: ",display.GetName());
	//return display.GetName().Contains("Wayland");

	// check WAYLAND_DISPLAY
	const char* w = getenv("WAYLAND_DISPLAY");
	if (w != nullptr)
		return true;
	return false;
}

bool isWSL() {
	Trace trc(2,"isWSL");
	struct stat buffer;
	bool rval{false};
	rval = (stat("/proc/sys/fs/binfmt_misc/WSLInterop", &buffer) == 0);
	trc.dprint("running under WSL? ",rval);
	return rval;
}

wxSize
getDisplaySize() {
// returns the size of the display in pixels by opening
// the default X11 screen and getting its dimensions
	Trace trc(2,"getDisplaySize");

    Display* display = XOpenDisplay(nullptr);
    if (!display) {
        wxLogError("Failed to open X display.");
        return wxSize(0, 0);
    }

    int screenNum = DefaultScreen(display);
    int width = DisplayWidth(display, screenNum);
    int height = DisplayHeight(display, screenNum);

    XCloseDisplay(display);
    return wxSize(width, height);
}


// FilterDialog implementation {

FilterDialog::
FilterDialog(wxWindow* parent) : wxDialog(parent,wxID_ANY,"Filter the file names?") {
	Trace trc(2,"FilterDialog");

	wxBoxSizer* vbox = new wxBoxSizer(wxVERTICAL);
	//SetSizer(vbox);

	// input wildcard filter
	vbox->Add(new wxStaticText(this, wxID_ANY,
		"separate multiple filters with semicolon :"), 0, wxALL, 5);
	long int style = wxTE_PROCESS_ENTER;
	filterctrl = new wxTextCtrl(this,wxID_ANY,"*.pf;*.apf;*.uf",
		wxDefaultPosition,wxDefaultSize,style);
#ifdef NEVER // get value in getfilters()
	filterctrl->Bind(wxEVT_TEXT_ENTER, [&](wxCommandEvent& ev) {
		Trace trc(2,"filter lambda");
		wxString f = filterctrl->GetValue();
	});
#endif // NEVER // get value in getfilters()
	vbox->Add(filterctrl, 0, wxEXPAND|wxALL, 5);

	// Ok and Cancel buttons
	wxBoxSizer* buttonSizer = new wxBoxSizer(wxHORIZONTAL);
	buttonSizer->Add(new wxButton(this, wxID_CANCEL, "Cancel"));
	wxButton* okButton = new wxButton(this, wxID_OK, "Ok");
	okButton->SetDefault();
#ifdef NEVER // get value in ?
	okButton->Bind(wxEVT_BUTTON, [&](wxCommandEvent& ev) {
		Trace trc(2,"ok lambda");
		wxString f = filterctrl->GetValue();
		filterstr = replace_char(f.ToStdString(), ',', '|');
		trc.dprint("got filters: ",filterstr);
		ev.Skip();
	});
#endif // NEVER // get value in ?
	buttonSizer->Add(okButton);
	buttonSizer->Layout();
	vbox->Add(buttonSizer, 0,wxBOTTOM, 10);

	SetSizerAndFit(vbox);
	vbox->Layout();
}
// FilterDialog implementation }

regex
wild2rx(const string& wild) {
// convert a string possibly containing wild-cards (*) to
// a regular-expression; allow non-escaped periods
	Trace trc(2,"wild2rx");
	int n = wild.size();
	string rval;

	for (int i=0; i<n; i++) {
		char c = wild[i];
		if (c == '.') {
			rval += "\\.";
		} else if (c == '*') {
			rval += ".*";
		} else {
			rval += c;
		}
	}
	trc.dprint("creating regex \"",rval,"\"");
	return regex(rval);
}
// vector<string> specialization
bool
containsrx(const vector<regex>& v, const string& s) {
// return true if "s" matches any regex in "v", false otherwise
	Trace trc(2,"containsrx ");
	trc.dprint("checking ",s," against ",v.size()," rx");
	for (auto& vi : v) {
		trc.dprint("testing ",s);
		// use search instead of match to catch, e.g. "1.b absgc1"
		//!! if (regex_search(s, vi))
		if (regex_match(s, vi)) {
			trc.dprint("matches");
			return true;
		}
	}
	trc.dprint("no match");
	return false;
}

int
uniqueId() {
	static int rval{wxID_HIGHEST-1};
	return rval--;
}

void
snooze(double nanosec) {
// sleep "nanosec" nano seconds (1e-9 sec) (1 < microsec < 1e9)
// if sec < 1 use nanosleep, otherwise use sleep()?
	if (nanosec < 0)
		throw runtime_error(vastr("illegal arg to snooze: ",nanosec));
	if (nanosec >= 1.0e+9) {
		unsigned int ns = nanosec/1.0e+9;
		sleep(ns);
	} else {
		struct timespec reg;
		reg.tv_sec = 0;
		reg.tv_nsec = nanosec;
		nanosleep(&reg, nullptr);
	}
}

string
sizer2string(const wxBoxSizer& sizer) {
// print some info about "sizer"
	Trace trc(2,"sizer2string");
	ostringstream os;
	 
	const wxSizerItemList& children = sizer.GetChildren();
	string orient = sizer.GetOrientation() == wxHORIZONTAL?"horizontal":"vertical";
	os << children.size() << " " << orient <<  " items\n";
	int i{0};
	for (const auto& ci : children) {
		wxSize sz = ci->GetSize();
		wxSize minsz = ci->GetMinSize();
		trc.dprint(++i,") item size ",sz.x,", ",sz.y,", min ",minsz.x,", ",minsz.y);
		os << "current size (" << sz.x << "," << sz.y << ") min size ("
			<< minsz.x << "," << minsz.y << ")\n";
		wxSizer* nested = ci->GetSizer();
		wxBoxSizer* nestedbox = dynamic_cast<wxBoxSizer*>(nested);
		if (nestedbox != nullptr) {
			os << "   nested sizer:" << stringIndent(3,sizer2string(*nestedbox));
		}
	}
	return os.str();
}

void
printsizes(const string& title, wxWindow* t) {
	int width, height;
	t->GetSize(&width, &height);
	wxSize min = t->GetMinSize();
	wxSize max = t->GetMaxSize();
	cerr << title << " size: (" << width << ", " << height << ")";
	cerr << " min: (" << min.x << ", " << min.y << ")";
	cerr << " max: (" << max.x << ", " << max.y << ")\n";
}

void
//!! pngwrite(const string& file, unsigned char* pixels, int width, int height)
pngwrite(const string& file, float* pixels, int width, int height) {
	Trace trc(2,"pngwrite ",file);

	trc.dprint("width ",width,", height ",height);

	File thefile(file, "wb");
	FILE* fp = thefile();

	png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
		nullptr, nullptr, nullptr);
	if (png_ptr == nullptr) {
		errormsg("cannot start png: create_write failed", false);
		return;
	}
	png_infop info_ptr = png_create_info_struct(png_ptr);
	if (info_ptr == nullptr) {
		errormsg("cannot start png: create_info failed", false);
		return;
	}
	// error handling
	if (setjmp(png_jmpbuf(png_ptr)) != 0) {
		cerr << "Error writing png file " << file << endl;
		return;
	}

	// initialize io
	png_init_io(png_ptr, fp);

	// set image info
	png_set_IHDR(png_ptr, info_ptr, width, height,
		8,			// bit depth
		PNG_COLOR_TYPE_RGB,	// color type
		PNG_INTERLACE_NONE,	// interlace type
		PNG_COMPRESSION_TYPE_DEFAULT,
		PNG_FILTER_TYPE_DEFAULT);

	// write the image information
	png_write_info(png_ptr, info_ptr);

	// write the image data, converting float to unsigned char
	// glReadPixels reads from top to bottom so we'll reverse it
	vector<unsigned char> row(width*3);
	for (int y=height-1; y>=0; y--) {
		for (int x=0; x < width*3; x++)
			row[x] = static_cast<unsigned char>(pixels[y*width*3+x]*255.0f);
		png_write_row(png_ptr, row.data());
	}

	png_write_end(png_ptr, nullptr);

	png_destroy_write_struct(&png_ptr, &info_ptr);
}

#ifdef MAIN
#undef MAIN

class MyFrame : public wxFrame
{
public:
    MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size);

private:
    void OnPaint(wxPaintEvent& event);
};

class MyApp : public wxApp
{
public:
    virtual bool OnInit();
};

// Implementations

MyFrame::MyFrame(const wxString& title, const wxPoint& pos, const wxSize& size)
    : wxFrame(nullptr, wxID_ANY, title, pos, size)
{
    // Connect the paint event to the OnPaint function
    Connect(wxEVT_PAINT, wxPaintEventHandler(MyFrame::OnPaint));
	// quit when any key is pressed
	// Bind(wxEVT_KEY_DOWN, [this](wxKeyEvent& s){ Close(true); }, wxID_ANY);
	// quit when ESC is pressed
	Bind(wxEVT_CHAR, [this](wxKeyEvent& s){
			if (s.GetKeyCode() == WXK_ESCAPE)
				Close(true);
		}, wxID_ANY);
}

MyFrame *frame{nullptr};

class Color {
public:
	string name;
	wxColour color;
	Color(const string& nm) : name(nm), color(getcolor(nm)) {}
};

vector<Color> thecolors{
{"AQUAMARINE"},
{"BLACK"},
{"BLUE"},
{"BLUE VIOLET"},
{"BROWN"},
{"CADET BLUE"},
{"CORAL"},
{"CORNFLOWER BLUE"},
{"CYAN"},
{"DARK GREY"},
{"DARK GREEN"},
{"DARK OLIVE GREEN"},
{"DARK ORCHID"},
{"DARK SLATE BLUE"},
{"DARK SLATE GREY"},
{"DARK TURQUOISE"},
{"DIM GREY"},
{"FIREBRICK"},
{"FOREST GREEN"},
{"GOLD"},
{"GOLDENROD"},
{"GREY"},
{"GREEN"},
{"GREEN YELLOW"},
{"INDIAN RED"},
{"KHAKI"},
{"LIGHT BLUE"},
{"LIGHT GREY"},
{"LIGHT STEEL BLUE"},
{"LIME GREEN"},
{"MAGENTA"},
{"MAROON"},
{"MEDIUM AQUAMARINE"},
{"MEDIUM BLUE"},
{"MEDIUM FOREST GREEN"},
{"MEDIUM GOLDENROD"},
{"MEDIUM ORCHID"},
{"MEDIUM SEA GREEN"},
{"MEDIUM SLATE BLUE"},
{"MEDIUM SPRING GREEN"},
{"MEDIUM TURQUOISE"},
{"MEDIUM VIOLET RED"},
{"MIDNIGHT BLUE"},
{"NAVY"},
{"ORANGE"},
{"ORANGE RED"},
{"ORCHID"},
{"PALE GREEN"},
{"PINK"},
{"PLUM"},
{"PURPLE"},
{"RED"},
{"SALMON"},
{"SEA GREEN"},
{"SIENNA"},
{"SKY BLUE"},
{"SLATE BLUE"},
{"SPRING GREEN"},
{"STEEL BLUE"},
{"TAN"},
{"THISTLE"},
{"TURQUOISE"},
{"VIOLET"},
{"VIOLET RED"},
{"WHEAT"},
{"WHITE"},
{"YELLOW"},
{"YELLOW GREEN"}};

void MyFrame::OnPaint(wxPaintEvent& event) {
    wxPaintDC dc(this); // Create a device context for drawing

	// Bind(wxEVT_KEY_UP, [&](wxKeyEvent& s){ frame->Close(true); }, wxID_EXIT);

	dc.SetFont(wxFont(13,wxFONTFAMILY_DEFAULT,wxFONTSTYLE_NORMAL,
		wxFONTWEIGHT_BOLD));
	int ix{20};
	int iy{10};
	int w, h;
	dc.GetTextExtent(thecolors[0].name, &w, &h);
	int skip = h + 5;
	int k{0};
	for (int j=0; j<4; j++) {
		iy = 0;
		for (int i=0; i<17; i++) {
			auto c = thecolors[k++];
			dc.SetTextForeground(c.color);
			dc.DrawText(c.name, ix, iy); // Draw the text
			iy += skip;
		}
		ix += 2*w;
	}
}

bool MyApp::OnInit()
{
	wxSize display = getDisplaySize();
	MyFrame *frame = new MyFrame("Available Colors", wxPoint(50, 50),
			wxSize(display.x-600, display.y-300));
	frame->Show(true);
	// Bind(wxEVT_KEY_UP, [&](wxCommandEvent& s){ frame->Close(true); }, wxID_EXIT);
	// frame->Bind(wxEVT_KEY_DOWN, [&](wxKeyEvent& s){ frame->Close(true); }, wxID_EXIT);
	// frame->Bind(wxEVT_KEY_UP, [&](wxKeyEvent& s){ frame->Close(true); }, wxID_EXIT);
	return true;
}

wxIMPLEMENT_APP(MyApp); // Provide দোকানের the entry point

#endif // MAIN
