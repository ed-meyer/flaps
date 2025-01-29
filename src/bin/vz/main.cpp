//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
/*
 * xgraph - A Simple Plotter for X
 *
 * David Harrison
 * University of California,  Berkeley
 * 1986, 1987, 1988, 1989
 *
 * Please see copyright.h concerning the formal reproduction rights
 * of this software.
 */

#define __USE_BSD 			/* caddr_t on LINUX */
#include <fcntl.h>       // mkdir
#include <sys/stat.h>    // mkdir
#include <pwd.h>
#include <cfloat>		/* added (EEM) */
#include <sys/utsname.h>
#include <X11/Xresource.h>

#undef Complex
#include "config.h"
#include "amdata.h"
#include "amvz.h"
#include "hard_devices.h"
// enabling fp exceptions causes amvz to crash in wx (divide by zero?)
// #include "main.h"
#include "params.h"
#include "specs.h"
#include "settings.h"
#include "text.h"

using namespace std;

pthread_mutex_t Disp_mutex = PTHREAD_MUTEX_INITIALIZER;

#define MAXMSG 1024

static bool Hardcopy = false;

// currently picked displacements and the x-y values
// on the plot
vector<vector<complex<double>> > Disp;
int DispColumn = 0;
vector<string> DispTitles;
static vector<double> Freqs;
static vector<complex<double>> Ses;
static vector<double> Veases;

double XParValue = 0.0;
double YParValue = 0.0;
string ModeName;
int PickedPoint = 0;


static void gnuplot();

#define ZOOM
#define TOOLBOX
#define NEW_READER

#define GRIDPOWER 	10
#define INITSIZE 	128

#define CONTROL_D	'\004'
#define CONTROL_C	'\003'
#define TILDE		'~'

#define BTNPAD		1
#define BTNINTER	3

#ifndef MAX					/* added (EEM) */
#define MAX(a,b)	((a) > (b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b)	((a) < (b) ? (a) : (b))
#endif
#ifndef ABS
#define ABS(x)		((x) < 0 ? -(x) : (x))
#endif
#define ZERO_THRES	1.0E-07

/* To get around an inaccurate log */
#define nlog10(x)	(x == 0.0 ? 0.0 : log10(x) + 1e-15)

#define ISCOLOR		(wi->dev_info.dev_flags & D_COLOR)

#define PIXVALUE(set) 	((set) % MAXATTR)

#define LINESTYLE(set) \
(ISCOLOR ?  ((set)/MAXATTR) : ((set) % MAXATTR))

#define MARKSTYLE(set) \
(colorMark ? COLMARK(set) : BWMARK(set))

#define COLMARK(set) \
((set) / MAXATTR)

#define BWMARK(set) \
((set) % MAXATTR)

#define LOG_X	0x01
#define LOG_Y	0x02


/*
 * Default settings for xgraph parameters
 */

#define DEF_BORDER_WIDTH	"2"
#define DEF_BORDER_COLOR	"Black"
#define DEF_TITLE_TEXT		"X Graph"
#define DEF_XUNIT_TEXT		"X"
#define DEF_YUNIT_TEXT		"Y"
#define DEF_TICK_FLAG		"off"
#define DEF_MARK_FLAG		"on"
#define DEF_PIXMARK_FLAG	"off"
#define DEF_LARGEPIX_FLAG	"off"
#define DEF_DIFFMARK_FLAG	"off"
#define DEF_BB_FLAG		"on"
#define DEF_NOLINE_FLAG		"off"
#define DEF_LOGX_FLAG		"off"
#define DEF_LOGY_FLAG		"off"
#define DEF_BAR_FLAG		"off"
#define DEF_BAR_BASE		"0.0"
#define DEF_BAR_WIDTH		"-1.0"
#define DEF_LINE_WIDTH		"1"
#define DEF_GRID_SIZE		"0"
#define DEF_GRID_STYLE		"10"
// #define DEF_LABEL_FONT		"helvetica-12"
// #define DEF_TITLE_FONT		"helvetica-16"
#define DEF_TITLE_FONT		"fixed-16"
#define DEF_LABEL_FONT		"fixed-12"
#define DEF_GEOMETRY		""
#define DEF_REVERSE		"off"
#define DEF_DEVICE		"Postscript"
#define DEF_DISPOSITION		"To Device"
#define DEF_FILEORDEV		"lpr"

#define DEF_MARKER_FLAG		"on"
#define DEF_DIFFMARK_FLAG	"off"
#define DEF_PIXMARK_FLAG	"off"
#define DEF_LARGEPIX_FLAG	"off"

/* Low > High means set it based on the data */
#define DEF_LOW_LIMIT		"1.0"
#define DEF_HIGH_LIMIT		"0.0"

/* Black and white defaults */
#define DEF_BW_BACKGROUND	"white"
#define DEF_BW_BORDER		"black"
#define DEF_BW_ZEROCOLOR	"black"
#define DEF_BW_ZEROWIDTH	"3"
#define DEF_BW_ZEROSTYLE	"1"
#define DEF_BW_FOREGROUND	"black"

/* Color defaults */
/* #define DEF_COL_BACKGROUND	"#ccc" */
#define DEF_COL_BACKGROUND	"white"
#define DEF_COL_BORDER		"black"
#define DEF_COL_ZEROCOLOR	"black"
#define DEF_COL_ZEROWIDTH	"2"
#define DEF_COL_ZEROSTYLE	"1"
#define DEF_COL_FOREGROUND	"black"
#define DEF_COL_FIRSTSTYLE	"1"

/* Default line styles */
static const char *defStyle[MAXATTR] = {
	"1", "10", "11110000", "010111", "1110",
	"1111111100000000", "11001111", "0011000111"
};

/* Default color names */
static const char *defColors[MAXATTR] = {
	"red", "green", "blue",
	"brown", "magenta", "orange", "cyan1", "black"
};



static char *tildeExpand(char*,char*);
static void ReverseIt();
static void Traverse(int);

vector<Plotcurve*> Curves;

// XSegment is a struct defined in /usr/include/X11/Xlib.h:
// typedef struct {
//    short x1, y1, x2, y2;
// } XSegment;


/* Basic transformation stuff */
static double llx, lly, urx, ury; /* Bounding box of all data */

#define HARDCOPY_IN_PROGRESS	0x01

typedef struct local_win {
	double loX, loY, hiX, hiY;	/* Local bounding box of window         */
	int XOrgX, XOrgY;		/* Origin of bounding box on screen     */
	int XOppX, XOppY;		/* Other point defining bounding box    */
	double UsrOrgX, UsrOrgY;	/* Origin of bounding box in user space */
	double UsrOppX, UsrOppY;	/* Other point of bounding box          */
	double XUnitsPerPixel;	/* X Axis scale factor                  */
	double YUnitsPerPixel;	/* Y Axis scale factor                  */
	xgOut dev_info;		/* Device information                   */
	Window close, hardcopy;	/* Buttons for closing and hardcopy     */
	Window about;		/* Version information                  */
	int flags;			/* Window flags                         */
} LocalWin;

#define VISIBLE(ws,userX,userY) \
	((userX) < ws->UsrOppX && (userX) > ws->UsrOrgX && \
	(userY) < ws->UsrOppY && (userY) > ws->UsrOrgY)

#define SCREENX(ws, userX) \
	(((int) (((userX) - ws->UsrOrgX)/ws->XUnitsPerPixel + 0.5)) + ws->XOrgX)
#define SCREENY(ws, userY) \
	(ws->XOppY - ((int) (((userY) - ws->UsrOrgY)/ws->YUnitsPerPixel + 0.5)))

#define USERX(ws, screenX) \
	((double)(screenX - ws->XOrgX)*ws->XUnitsPerPixel + (double)ws->UsrOrgX)
#define USERY(ws, screenY) \
	(ws->UsrOrgY - (double)(screenY - ws->XOppY)*ws->YUnitsPerPixel)

static XContext win_context = (XContext) 0;

/* Other globally set defaults */

Display *disp{nullptr};			/* Open display            */
Visual *visual{nullptr};			/* Standard visual         */
Colormap cmap;			/* Standard colormap       */
int screen;			/* Screen number           */
int depth;			/* Depth of screen         */


/* Total number of active windows */
static int Num_Windows = 0;
static char *Prog_Name;


static int XErrHandler( Display *disp_ptr, XErrorEvent *evt);

static int InitSets();
static int ReadDefaults();
static void set_mark_flags( int *markFlag, int *pixelMarks,
	int *bigPixel, int *colorMark);


static Window NewWindow( char *progname, double lowX, double lowY,
	double upX, double upY, double asp);
static int DrawWindow(LocalWin *win_info);
static void DelWindow( Window win, LocalWin *win_info);
static void PrintWindow( Window win, LocalWin *win_info);
static int HandlePickPoint( char *progname, XButtonPressedEvent *evt,
	LocalWin *wi, Cursor cur);
static int HandleShowZ ( char *progname, XButtonPressedEvent *evt,
	LocalWin *wi, Cursor cur);
static int HandleZoom( char *progname, XButtonPressedEvent *evt, LocalWin *wi, Cursor cur);
static int TransformCompute (LocalWin *wi);
static void DrawTitle (LocalWin *wi);
static int DrawLegend(LocalWin* wi);
static int DrawGridAndAxis (LocalWin *wi);
static int DrawData (LocalWin *wi);
static double initGrid( double low, double step, int logFlag);
static double stepGrid();
static double RoundUp(double val);
static int WriteValue(char*, double, int, int);

static int dovz (int argc, char** argv);

int
main(int argc, char** argv) {
	char* dbg = getenv("DEBUG");
	if (dbg != nullptr) Trace::debug(std::stoi(dbg));
	Trace trc(1,"main");
	Ftmpdir ftmp;   // create a temp directory if not already

	try {
		return dovz(argc, argv);
	} catch (runtime_error& s) {
		flaps::error(s.what());
	} catch (std::bad_alloc& b) {
		ostringstream os;
		flaps::error("out of memory");
	} catch (std::exception& e) {
		flaps::error("caught exception: ",e.what());
	}
	return flaps::nerrors();
}

static int
dovz (int argc, char** argv) {
// This sets up the hard-wired defaults and reads the X defaults.
// The command line format is: xgraph [host:display].
	Trace trc(2,"dovz");
	Window primary;
	XEvent theEvent;
	LocalWin *win_info;
	Cursor zoomCursor;
	XColor fg_color, bg_color;
	char keys[MAXKEYS];
	char *display = getenv("DISPLAY");
	string disp_name(":0");
	int nbytes, idx, flags;
	Specs& sp = specs();

	Prog_Name = argv[0];

	// Check the DISPLAY variable for validity if the local host -
	// if so set it to ":0" to force X to use unix domain sockets
	string remotehost = get_uname("nodename");
	string disphost = remotehost;
	if (display) {
		disp_name = display;
		string::size_type idx = disp_name.find(':');
		if (idx == string::npos) {
			if (disp_name == "novz" || disp_name == "off" ||
				disp_name == "none")
				return 0;
			flaps::warning("bad DISPLAY variable (", display, "): no colon");
			disp_name = ":0";
			disphost = remotehost;
		} else {
			disphost = disp_name.substr(0,idx);
			if (disphost == remotehost)
				disp_name = ":0";
		}
	}

	trc.dprint("opening display on ",disphost);

	disp = XOpenDisplay(0);
	if (!disp) {
		if (!disp_name.empty()) {
			if (disphost.size() < 1) {
				flaps::error("bad DISPLAY environment variable: \"", disp_name,
					"\" - it should have the form \"host:0\"");
				return (1);
			}
			flaps::error("cannot open display \"", disp_name,
				"\": did you forget to type \"xhost ", remotehost,
				"\" on \"", disphost, "\"?");
		} else {
			flaps::error("cannot open display on local workstation");
		}
		return(1);
	} else {
		trc.dprint("got display ",disp);
	}

	XSetErrorHandler(XErrHandler);

	// Set up hard-wired defaults and allocate spaces
	InitSets();

	// Read X defaults and override hard-coded defaults
	ReadDefaults();

	// Convert command-line options to Flaps program style
	// options and feed them to parse_specs()
	string options = convertOptions (argc, argv);
	if (!parse_specs(options, specs()))
		return(1);

	// Read the data into the data sets
	// Do this before forking so we don't get clobbered
	// by subsequent vz's
	llx = lly = std::numeric_limits<int>::max();
	urx = ury = -std::numeric_limits<int>::max();
	Curves = Plotcurve::getcurves();
	if (Curves.empty()) {
		trc.dprint("quitting: Plotcurve::getcurves returned empty");
		return(1);
	}
	// take vzid from the first curve
	string vzid = Curves[0]->vzid();		// XXX not needed? rm from Curve?
	// cout << "got vzid = " << vzid << endl;
	if (!vzid.empty())
		Settings::defaults.set("vzid", vzid);

	// default title: first curve aid: xname-yname
	if (sp.title.empty()) {
		ostringstream os;
		Plotcurve* cp = Curves[0];
		assert(cp != nullptr);
		os << cp->aid() << ": ";
		if (cp->xparam != nullptr)
			os << cp->xparam->name << '-';
		os << cp->yparam->name;
		sp.title = os.str();
	}

	if (sp.gnuplot) {
		gnuplot();
	}

	/* Reverse Video Hack */
	if (PM_BOOL("ReverseVideo")) ReverseIt();
	hard_init();
	if (PM_BOOL("Debug")) {
		(void) XSynchronize(disp, 1);
		param_dump();
	}

	/* Logarithmic and bounding box computation */
	flags = 0;
#ifdef NEVER // use specs
	if (PM_BOOL("LogX")) flags |= LOG_X;
	if (PM_BOOL("LogY")) flags |= LOG_Y;
#else // NEVER // use specs
	if (sp.xlog) flags |= LOG_X;
	if (sp.ylog) flags |= LOG_Y;
#endif // NEVER // use specs
	Traverse(flags);

	/* Nasty hack here for bar graphs */
	if (PM_BOOL("BarGraph")) {
		double base;

		llx -= PM_DBL("BarWidth");
		urx += PM_DBL("BarWidth");
		base = PM_DBL("BarBase");
		if (base < lly) lly = base;
		if (base > ury) ury = base;
	}

	/* Create initial window */
	xtb_init(disp, screen, PM_PIXEL("Foreground"), PM_PIXEL("Background"),
		 PM_FONT("LabelFont"));
	double aspect_ratio{1.2};
	primary = NewWindow(Prog_Name, sp.xmin, sp.ymin,
		sp.xmax, sp.ymax, aspect_ratio);
	if (!primary)
		throw runtime_error("main window would not open");

	zoomCursor = XCreateFontCursor(disp, XC_sizing);
	fg_color = PM_COLOR("Foreground");
	bg_color = PM_COLOR("Background");
	XRecolorCursor(disp, zoomCursor, &fg_color, &bg_color);

	Num_Windows = 1;
	while (Num_Windows > 0) {
		XNextEvent(disp, &theEvent);
		if (xtb_dispatch(&theEvent) != XTB_NOTDEF) continue;
		if (XFindContext(theEvent.xany.display,
			 theEvent.xany.window,
			 win_context, (caddr_t *) &win_info)) {
			/* Nothing found */
			continue;
		}
		switch (theEvent.type) {
		case Expose:
			if (theEvent.xexpose.count <= 0) {
				XWindowAttributes win_attr;

				XGetWindowAttributes(disp, theEvent.xany.window, &win_attr);
				win_info->dev_info.area_w = win_attr.width;
				win_info->dev_info.area_h = win_attr.height;
				init_X(win_info->dev_info.user_state);
				DrawWindow(win_info);
			}
			break;
		case KeyPress:
			nbytes = XLookupString(&theEvent.xkey, keys, MAXKEYS,
				 (KeySym *) 0, (XComposeStatus *) 0);
			for (idx = 0;  idx < nbytes;  idx++) {
				// pressing 'q' or ctrl-d kills this window
				if (keys[idx] == CONTROL_D) {
					/* Delete this window */
					DelWindow(theEvent.xkey.window, win_info);
				} else if (keys[idx] == 'q') {
					/* Delete this window */
					DelWindow(theEvent.xkey.window, win_info);
				} else if (keys[idx] == CONTROL_C) {
					/* Exit program */
					Num_Windows = 0;
				} else if (keys[idx] == 'h') {
					PrintWindow(theEvent.xany.window, win_info);
				}
			}
			break;
		case ButtonPress:
			/* Handle creating a new window */
			{
				XButtonEvent *ep = (XButtonEvent*)&theEvent;
				int button = ep->button;
				int state = ep->state;
					// right-click or left & ? state 4?
				if (button == 3 ||
					(button == 1 && state == 4)) {
					HandlePickPoint(Prog_Name,
						 &theEvent.xbutton,
						 win_info, zoomCursor);
					// middle-click or shift-left-click: show z-values
				} else if (button == 2 ||
					(button == 1 && (state == 1 || state == 17))) {
					HandleShowZ (Prog_Name, &theEvent.xbutton, win_info, zoomCursor);
				} else {
					// left-click & drag: zoom
					Num_Windows += HandleZoom(Prog_Name,
						 &theEvent.xbutton,
						 win_info, zoomCursor);
				}
			}
			break;
		default:
			(void) fprintf(stderr, "Unknown event type: %x\n", theEvent.type);
			break;
		}
	}
	return 0;
}  // dovz

static void
gnuplot() {
	Trace trc(1,"gnuplot");
	string ftmp{getenv("FTMP")};
	ostringstream os;

	// create a directory under ftmp for gnuplot's data
	string gnudir{vastr(ftmp,"/gnuplot.",getpid())};
	if (mkdir(gnudir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO) == -1) {
		throw std::runtime_error("cannot create temp directory \"" +
					gnudir + "\": " + strerror(errno));
	}

	// move to the gnuplot directory
	Chdir tmp(gnudir);

	// create a file to hold gnuplot commands
	string cmdfilename(".gnuplot");
	ofstream cmdfile(cmdfilename);

	// If there is a z parameter use "splot" (3D)
	bool splot = false;
	if (Curves[0]->zparam == nullptr) {
		splot = false;
		cmdfile << "plot ";
	} else {
		splot = true;
		cmdfile << "splot ";
	}

	// create a file for the plot data, one per aid:cid
	bool first = true;
	for (size_t i=0; i<Curves.size(); i++) {
		Par* xp = Curves[i]->xparam;
		if (xp == nullptr)
			throw runtime_error("must include x param for gnuplot");
		string aid = Curves[i]->aid();
		string cid = Curves[i]->cid();
		string datafilename = vastr(aid,':',cid);
		ofstream datafile(datafilename);
		if (first) {
			first = false;
		} else {
			cmdfile << ", ";
		}
		cmdfile << "\"" << datafilename << "\" with lines";

		if (xp == nullptr)
			continue;
		Par* yp = Curves[i]->yparam;
		Par* zp = splot? Curves[i]->zparam : nullptr;
		for (size_t k = 0; k<xp->nsolns(); k++) {
			if (splot) {
				datafile << xp->solns[k] << ' ' << yp->solns[k]
						<< ' ' << zp->solns[k] << endl;
			} else {
				datafile << xp->solns[k] << ' ' << yp->solns[k] << endl;
			}
		}
		// if 3D add a blank line so that files can be concatenated
		// to form a gridded data file
		if (splot)
			datafile << endl;
		datafile.close();
	}

	cmdfile << "\npause -1 \"press q to quit\"\n";

	cmdfile.close();

	trc.dprint("running gnuplot with (",gnudir,'/',cmdfilename);

	execlp ("gnuplot", "gnuplot", ".gnuplot", nullptr);

	// should not return from execlp:
	os << "cannot start gnuplot: " << strerror(errno);
	throw runtime_error(os.str());
}

#define BLACK_THRES	30000

static void ReversePix(
const char *param_name)		/* Name of color parameter */
/*
 * Looks up `param_name' in the parameters database.  If found, the
 * color is examined and judged to be either black or white based
 * upon its red, green, and blue intensities.  The sense of the
 * color is then reversed and reset to its opposite.
 */
{
	params val;

	if (param_get(param_name, &val)) {
		if ((val.pixv.value.red < BLACK_THRES) &&
		    (val.pixv.value.green < BLACK_THRES) &&
		    (val.pixv.value.blue < BLACK_THRES)) {
			/* Color is black */
			param_reset(param_name, "white");
		} else {
			/* Color is white */
			param_reset(param_name, "black");
		}
	} else {
		(void) fprintf(stderr, "Cannot reverse color `%s'\n", param_name);
	}
}

static void
ReverseIt()
/*
 * This routine attempts to implement reverse video.  It steps through
 * all of the important colors in the parameters database and makes
 * black white (and vice versa).
 */
{
	int i;
	char buf[1024];

	for (i = 0;  i < MAXATTR;  i++) {
		(void) sprintf(buf, "%d.Color", i);
		ReversePix(buf);
	}
	ReversePix("Foreground");
	ReversePix("Border");
	ReversePix("ZeroColor");
	ReversePix("Background");
}


static void
Traverse(int flags) {
/*
 * Traverses through all of the data applying certain options to the
 * data and computing the overall bounding box.  The flags are:
 *   LOG_X	Take the log of the X axis
 *   LOG_Y	Take the log of the Y axis
 */
	double eps = 1.e-10;
	for (size_t i = 0;  i < Curves.size();  i++) {
		Par* xp = Curves[i]->xparam;
		if (!xp) continue;
		Par* yp = Curves[i]->yparam;
		for (size_t k = 0; k<yp->nsolns(); k++) {
			if (flags & LOG_Y) {
				if (yp->solns[k] > 0.0) {
					yp->solns[k] = log10(yp->solns[k]);
				} else {
					flaps::warning("Cannot plot non-positive Y values "
						"when the logarithmic option is selected. Point ",
						k+1," of data set ",i+1," is ",yp->solns[k]);
					yp->solns[k] = eps;
				}
			}
			if (flags & LOG_X) {
				if (xp->solns[k] > 0.0) {
					xp->solns[k] = log10(xp->solns[k]);
				} else {
					flaps::warning("Cannot plot non-positive X values "
						"when the logarithmic option is selected. Point ",
						k+1," of data set ",i+1," is ",xp->solns[k]);
					xp->solns[k] = 0.0;
				}
			}
		}
		// Update global bounding box - include min/max of points
		// adjacent to ensure we use the range necessary to
		// capture all the visible parts of the curve
		for (size_t k = 0; k<yp->nsolns(); k++) {
			double xmin = xp->solns[k];
			double xmax = xp->solns[k];
			double ymin = yp->solns[k];
			double ymax = yp->solns[k];
			if (k > 0) {
				xmin = std::min(xmin, xp->solns[k-1]);
				xmax = std::max(xmax, xp->solns[k-1]);
				ymin = std::min(ymin, yp->solns[k-1]);
				ymax = std::max(ymax, yp->solns[k-1]);
			}
			if (k < yp->nsolns()-1) {
				xmin = std::min(xmin, xp->solns[k+1]);
				xmax = std::max(xmax, xp->solns[k+1]);
				ymin = std::min(ymin, yp->solns[k+1]);
				ymax = std::max(ymax, yp->solns[k+1]);
			}
			updateBoundingBox(xp->solns[k], xmin, xmax,
					yp->solns[k], ymin, ymax, llx, urx, lly, ury);
			// if (xp->values[k] < llx) llx = xp->values[k];
			// if (xp->values[k] > urx) urx = xp->values[k];
			// if (yp->values[k] < lly) lly = yp->values[k];
			// if (yp->values[k] > ury) ury = yp->values[k];
		}
	}
}



/*
 * Button handling functions
 */

/*ARGSUSED*/
xtb_hret del_func(
Window win,			/* Button window    */
int bval,			/* Button value     */
xtb_data info)			/* User information */
/*
 * This routine is called when the `Close' button is pressed in
 * an xgraph window.  It causes the window to go away.
 */
{
	Window the_win = (Window) info;
	LocalWin *win_info;

	xtb_bt_set(win, 1, (char *) 0, 0);
	if (!XFindContext(disp, the_win, win_context, (caddr_t *) &win_info)) {
		if (win_info->flags & HARDCOPY_IN_PROGRESS) {
			do_error("Can't close window while\nhardcopy dialog is posted.\n");
			xtb_bt_set(win, 0, (char *) 0, 0);
		} else {
			DelWindow(the_win, win_info);
		}
	}
	return XTB_HANDLED;
}

/*ARGSUSED*/
xtb_hret hcpy_func(
Window win,			/* Button Window    */
int bval,			/* Button value     */
xtb_data info)			/* User Information */
/*
 * This routine is called when the hardcopy button is pressed
 * in an xgraph window.  It causes the output dialog to be
 * posted.
 */
{
	Trace trc(1,"hcpy_func");
	Window the_win = (Window) info;
	LocalWin *win_info;

	xtb_bt_set(win, 1, (char *) 0, 0);
	if (!XFindContext(disp, the_win, win_context, (caddr_t *) &win_info)) {
		win_info->flags |= HARDCOPY_IN_PROGRESS;
		PrintWindow(the_win, win_info);
		win_info->flags &= (~HARDCOPY_IN_PROGRESS);
	}
	xtb_bt_set(win, 0, (char *) 0, 0);
	return XTB_HANDLED;
}

#define NORMSIZE	600
#define MINDIM		100

static Window
NewWindow(
char *progname,			/* Name of program    */
double lowX, double lowY,		/* Lower left corner  */
double upX, double upY,		/* Upper right corner */
double asp)			/* Aspect ratio       */
/*
 * Creates and maps a new window.  This includes allocating its
 * local structure and associating it with the XId for the window.
 * The aspect ratio is specified as the ratio of width over height.
 */
{
	Window new_window;
	LocalWin *new_info;
	static Cursor theCursor = (Cursor) 0;
	XSizeHints sizehints;
	XSetWindowAttributes wattr;
	XWMHints wmhints;
	XColor fg_color, bg_color;
	int geo_mask;
	int width, height;
	unsigned long wamask;
	char defSpec[120];
	double pad;
	ostringstream os;
	Specs& sp = specs();

	new_info = (LocalWin *) calloc(1,sizeof(LocalWin));

	// watch out for urx-llx small: add a tiny amount
	double eps = 100.0*std::numeric_limits<double>::epsilon()*(abs(llx)+abs(urx));
	if (abs(urx-llx) < eps) {
		llx -= 1000.0*eps;
		urx += 1000.0*eps;
	}
	eps = 100.0*std::numeric_limits<double>::epsilon()*(abs(lly)+abs(ury));
	if (abs(ury-lly) < eps) {
		lly -= 1000.0*eps;
		ury += 1000.0*eps;
	}
	/*
	 * if a lower limit is less than the upper => take
	 * limits from the data
	 */
	if (upX > lowX) {
		new_info->loX = lowX;
		new_info->hiX = upX;
	} else {
		new_info->loX = llx;
		new_info->hiX = urx;
	}
	if (upY > lowY) {
		new_info->loY = lowY;
		new_info->hiY = upY;
	} else {
		new_info->loY = lly;
		new_info->hiY = ury;
	}

	/* Increase the padding for aesthetics */
	if (new_info->hiX - new_info->loX == 0.0) {
		pad = MAX(0.5, fabs(new_info->hiX/2.0));
		new_info->hiX += pad;
		new_info->loX -= pad;
	}
	if (new_info->hiY - new_info->loY == 0) {
		pad = MAX(0.5, fabs(ury/2.0));
		new_info->hiY += pad;
		new_info->loY -= pad;
	}

	/* Add 10% padding to bounding box (div by 20 yeilds 5%) */
	pad = (new_info->hiX - new_info->loX) / 20.0;
	new_info->loX -= pad;  
	new_info->hiX += pad;
	pad = (new_info->hiY - new_info->loY) / 20.0;
	new_info->loY -= pad;  
	new_info->hiY += pad;

	/* Aspect ratio computation */
	if (asp < 1.0) {
		height = NORMSIZE;
		width = ((int) (((double) NORMSIZE) * asp));
	} else {
		width = NORMSIZE;
		height = ((int) (((double) NORMSIZE) / asp));
	}
	height = MAX(MINDIM, height);
	width = MAX(MINDIM, width);
	(void) sprintf(defSpec, "%dx%d+100+100", width, height);

	wamask = CWBackPixel | CWBorderPixel | CWColormap;
	wattr.background_pixel = PM_PIXEL("Background");
	wattr.border_pixel = PM_PIXEL("Border");
	wattr.colormap = cmap;

	sizehints.flags = PPosition|PSize;
	sizehints.x = sizehints.y = 100;
	sizehints.width = width;
	sizehints.height = height;

	new_window = XCreateWindow(disp, RootWindow(disp, screen),
	    sizehints.x, sizehints.y,
	    (unsigned int) sizehints.width,
	    (unsigned int) sizehints.height,
	    (unsigned int) PM_INT("BorderSize"),
	    depth, InputOutput, visual,
	    wamask, &wattr);


	if (new_window == 0)
		return new_window;

	xtb_frame cl_frame, hd_frame;
	struct utsname sysName;
	const char *machine;
	static string title;

	if (sp.title.empty()) {
		if (uname(&sysName) == -1) {
			machine = "unknown machine";
		} else {
			machine = sysName.nodename;
		}
		os.str("");
		os << "Flaps on " << machine;
		title = os.str();
	}

	XStoreName(disp, new_window, sp.title.c_str());
	XSetIconName(disp, new_window, sp.title.c_str());

	wmhints.flags = InputHint | StateHint;
	wmhints.input = True;
	wmhints.initial_state = NormalState;
	XSetWMHints(disp, new_window, &wmhints);

	geo_mask = XParseGeometry(pm_str("Geometry"), &sizehints.x, &sizehints.y,
		 (unsigned int *) &sizehints.width,
		 (unsigned int *) &sizehints.height);
	if (geo_mask & (XValue | YValue)) {
		sizehints.flags = (sizehints.flags & ~PPosition) | USPosition;
	}
	if (geo_mask & (WidthValue | HeightValue)) {
		sizehints.flags = (sizehints.flags & ~PSize) | USSize;
	}
	XSetNormalHints(disp, new_window, &sizehints);

	/* Set device info */
	set_X(new_window, &(new_info->dev_info));

	/* Make buttons */
	xtb_bt_new(new_window, "Close", del_func,
		 (xtb_data) new_window, &cl_frame);
	new_info->close = cl_frame.win;
	XMoveWindow(disp, new_info->close, (int) BTNPAD, (int) BTNPAD);
	xtb_bt_new(new_window, "Hardcopy", hcpy_func,
		 (xtb_data) new_window, &hd_frame);
	new_info->hardcopy = hd_frame.win;
	XMoveWindow(disp, new_info->hardcopy,
		 (int) (BTNPAD + cl_frame.width + BTNINTER),
		 BTNPAD);

	new_info->flags = 0;
	XSelectInput(disp, new_window,
		 ExposureMask|KeyPressMask|ButtonPressMask);
	if (!theCursor) {
		theCursor = XCreateFontCursor(disp, XC_top_left_arrow);
		fg_color = PM_COLOR("Foreground");
		bg_color = PM_COLOR("Background");
		XRecolorCursor(disp, theCursor, &fg_color, &bg_color);
	}
	XDefineCursor(disp, new_window, theCursor);
	if (!win_context) {
		win_context = XUniqueContext();
	}
	XSaveContext(disp, new_window, win_context, (caddr_t) new_info);
	XMapWindow(disp, new_window);
	return new_window;
}


static void
DelWindow( Window win, LocalWin *win_info) {
/*
 * This routine actually deletes the specified window and
 * decrements the window count.
 */
	xtb_data info;

	XDeleteContext(disp, win, win_context);
	xtb_bt_del(win_info->close, &info);
	xtb_bt_del(win_info->hardcopy, &info);
	free((char *) win_info);
	XDestroyWindow(disp, win);
	Num_Windows -= 1;
}

static void
PrintWindow( Window win, LocalWin *win_info) {
/*
 * This routine posts a dialog asking about the hardcopy
 * options desired.  If the user hits `OK',  the hard
 * copy is performed.
 */
	ho_dialog(win, Prog_Name, (char *) win_info);
}


static XRectangle boxEcho;
// GC is X11 graphics context, not generalized coordinate
static GC echoGC = (GC) 0;

#define DRAWBOX \
if (startX < curX) { \
   boxEcho.x = startX; \
   boxEcho.width = curX - startX; \
} else { \
   boxEcho.x = curX; \
   boxEcho.width = startX - curX; \
} \
if (startY < curY) { \
   boxEcho.y = startY; \
   boxEcho.height = curY - startY; \
} else { \
   boxEcho.y = curY; \
   boxEcho.height = startY - curY; \
} \
XDrawRectangles(disp, win, echoGC, &boxEcho, 1);

#define TRANX(xval) \
(((double) ((xval) - wi->XOrgX)) * wi->XUnitsPerPixel + wi->UsrOrgX)

#define TRANY(yval) \
(wi->UsrOppY - (((double) ((yval) - wi->XOrgY)) * wi->YUnitsPerPixel))

int
closestPoint (double screenx, double screeny, LocalWin* wi, size_t& solnidx,
		double& actualx, double& actualy) {
/*
 * Find the closest point...
 * Returns:
 *    index (zero-based) of the curve number which is closest
 *    solnidx - index (zero-based) of the point on the curve, i.e. the
 *             index in the solns arrays of the parameters
 */
	Trace trc(1,"closestPoint");
	int curveno = -1;
	double mindist = std::numeric_limits<double>::max();

	for (size_t idx=0; idx<Curves.size(); idx++) {
		Par* xp = Curves[idx]->xparam;
		if (xp == nullptr) continue;
		trc.dprint("searching ",Curves[idx]->cid());
		Par* yp = Curves[idx]->yparam;
		for (size_t j = 0; j < yp->nsolns(); j++) {
			if (VISIBLE(wi,xp->solns[j],yp->solns[j])) {
				double dx = SCREENX(wi,xp->solns[j]) - screenx;
				double dy = SCREENY(wi,yp->solns[j]) - screeny;
				double dist = dx*dx + dy*dy;
				// XXX inefficient - find mindist first, then save solns
				if (dist < mindist) {
					mindist = dist;
					actualx = xp->solns[j];
					actualy = yp->solns[j];
					solnidx = j;  // zero-based
					curveno = idx;
				}
			}
		}
	}
	return curveno;
}

bool
get_mids(string& nodes, string& coords, string& conn, string& gct) {
	Trace trc(1,"get_mids");
	Specs& sp = specs();
	// first check to see if they are in Matrix::collection()
	Matrix* nodes_p = Matrix::find_desc("nodes");
	Matrix* coords_p = Matrix::find_desc("coords");
	Matrix* conn_p = Matrix::find_desc("conn");
	Matrix* gct_p = Matrix::find_desc("gct");
	if (nodes_p != nullptr) {
		nodes = nodes_p->mid();
		if (coords_p == nullptr || conn_p == nullptr || gct_p == nullptr)
			throw runtime_error("visualization data is incomplete");
		coords = coords_p->mid();
		conn = conn_p->mid();
		gct = gct_p->mid();
	} else {
		// ... then try fetching a list of mids in this analysis (aid)
		vector<string> aids = sp.aids;
		if (aids.empty())
			return false;
			//!! throw runtime_error("cannot process right-click: no aids specified");
		vector<pair<string,string>> mids = Matrix::fetch_mids(aids[0]);
		for (auto& pp : mids) {
			if (pp.first == "nodes")
				nodes = pp.second;
			else if (pp.first == "coords")
				coords = pp.second;
			else if (pp.first == "conn")
				conn = pp.second;
			else if (pp.first == "gct")
				gct = pp.second;
		}
	}
	trc.dprint("returning nodes=\"",nodes,"\"");
	trc.dprint("returning coords=\"",coords,"\"");
	trc.dprint("returning conn=\"",conn,"\"");
	trc.dprint("returning gct=\"",gct,"\"");
	return true;
}

static int
HandlePickPoint( char *progname, XButtonPressedEvent *evt,
	LocalWin *wi, Cursor cur) {
// Respond to a right-click by multiplying the gc transformation matrix
// (gctransform) by the picked gc, add the result to the global vector Disp.
// Create optional universal and output4 files (.op4)
	Trace trc(1,"HandlePickPoint");
	double actualx, actualy;
	int ax, ay;
	int markFlag, pixelMarks, bigPixel, colorMark;
	int startX, startY;
	ostringstream os;
	static bool noamvz{false};  // is amvz available?

	if (noamvz) {
		cerr << "amvz is not available\n";
		return 1;
	}

	// NOTE: GC here means Graphics Context, not Generalized Coordinate :-(
	if (echoGC == (GC) 0) {
		unsigned long gcmask;
		XGCValues gcvals;

		gcmask = GCForeground | GCFunction;
		gcvals.foreground = PM_PIXEL("ZeroColor") ^ PM_PIXEL("Background");
		gcvals.function = GXxor;
		echoGC = XCreateGC(disp, evt->window, gcmask, &gcvals);
	}
	set_mark_flags(&markFlag, &pixelMarks, &bigPixel, &colorMark);
	startX = evt->x;  
	startY = evt->y;
	trc.dprint("x ",startX,", y ",startY);
	wi->dev_info.xg_dot(wi->dev_info.user_state,
	    startX,startY,
	    P_MARK,MARKSTYLE(0),PIXVALUE(0));
	XFlush (disp);

	// Find the closest point: Curve (index) & solns (solnidx)
	size_t solnidx;
	int index = closestPoint (startX, startY, wi, solnidx, actualx, actualy);

	if (index == -1) {
		cout << "no visible points\n";
		trc.dprint("returning 0: no visible points");
		return 0;
	}
	// no x parameter?
	if (Curves[index]->xparam == nullptr)
		return 0;

	Curves[index]->tagged[solnidx] = true;  // for redraws
	trc.dprint("subindex ",solnidx," tagged");

	trc.dprint("index ",index,", cid ",Curves[index]->cid());
	trc.dprint("actual x ",actualx,", y ",actualy);
	ax = SCREENX(wi, actualx);
	ay = SCREENY(wi, actualy);
	trc.dprint("coord of actual pt (screen units): ",actualx, " (",ax,") ",actualy," (",ay,") solnidx ",solnidx);
	/* wi->dev_info.xg_dot(wi->dev_info.user_state,
							ax, ay,
							P_MARK,BWMARK(0),PIXVALUE(0)); */

	// ... and color it
	wi->dev_info.xg_dot(wi->dev_info.user_state, ax, ay,
	    P_TAG,MARKSTYLE(index),PIXVALUE(index));
	/* XFillRectangle (disp, evt->window, echoGC, ax-4, ay-4, 8, 8); */
	// put a number next to the spot
	// XXX this should use PickedPoint instead of instance
	static int instance{1};
	string picked_name{vastr(instance++)};
	wi->dev_info.xg_text(wi->dev_info.user_state, ax, ay,
		 picked_name.c_str(), T_LOWERRIGHT, T_TITLE);

	XFlush (disp);

	trc.dprint("index ",index,",  title ",os.str());

	ModeName = Curves[index]->cid();
	XParValue = Curves[index]->xparam->solns[solnidx];
	YParValue = Curves[index]->yparam->solns[solnidx];

	// calling instance() the first time creates the only Amdata
	static std::thread* amvz_thread{nullptr};
	if (amvz_thread == nullptr) {
		string ufname;
		Settings::defaults.get("ufname", ufname);
		if (!ufname.empty()) {
			Amdata::initialize(ufname);
		} else {
			// defaults:
			string nodes{"nodes"};
			string coords{"coords"};
			string conn{"conn"};
			string gct{"gct"};
			// check to see if they are in Matrix::collection()
			if (get_mids(nodes, coords, conn, gct))
				Amdata::initialize(nodes, coords, conn, gct);
			else {
				cerr << "amvz data is not available\n";
				noamvz = true;
				return 1;
			}
		}
	}
	// add gctransform*ev to the list of selected modes
	// if this is the first time here instance() will create an Amdata
	// which will fetch coord, conn, and gctransform
	Amdata::instance()->add(Curves[index]->params.get_eigenvector(solnidx));
	// add the picked instance to picked_names XXX add name to add()?
	Amdata::instance()->add_picked_name(picked_name);
	// start amvz if not already running
	if (amvz_thread == nullptr) {
		amvz_thread = new std::thread(amvz);
		amvz_thread->detach();
		trc.dprint("started amvz thread: ", amvz_thread->get_id());
	}

	return 0;
} // HandlePickPoint

static int
HandleShowZ ( char *progname, XButtonPressedEvent *evt,
	LocalWin *wi, Cursor cur) {
	Trace trc(1,"HandleShowZ");
// Respond to a middle or shift-left button click by displaying
// all parameter values in a separate window.
// Also write these values to a file, ".vz" which gets
// overwritten each time the middle button is clicked
	double actualx, actualy;
	int markFlag, pixelMarks, bigPixel, colorMark;
	int startX, startY;
	ostringstream os;
	Specs& sp = specs();

	if (echoGC == (GC) 0) {
		unsigned long gcmask;
		XGCValues gcvals;

		gcmask = GCForeground | GCFunction;
		gcvals.foreground = PM_PIXEL("ZeroColor") ^ PM_PIXEL("Background");
		gcvals.function = GXxor;
		echoGC = XCreateGC(disp, evt->window, gcmask, &gcvals);
	}
	set_mark_flags(&markFlag, &pixelMarks, &bigPixel, &colorMark);
	startX = evt->x;  
	startY = evt->y;
	trc.dprint("HandleShowZ: x ",startX,", y ",startY);
	wi->dev_info.xg_dot(wi->dev_info.user_state,
	    startX,startY,
	    P_MARK,MARKSTYLE(0),PIXVALUE(0));
	XFlush (disp);

	// Find the closest point...
	size_t solnidx;
	int index = closestPoint (startX, startY, wi, solnidx, actualx, actualy);
	if (index == -1) {
		printf ("no visible points\n");
		return 0;
	}

	Plotcurve* curve = Curves[index];

	// Display info about the picked point in a message box
	// If there was no z parameter specified, just display
	// the name of the closest set so the user can figure
	// out which curve he is looking at
	int perline = 4;
	int items = 0;
	size_t width = 0;
	bool hasnewline = false;
	string vzfile{".vz"};
	string saveclicks = sp.saveclicksfile;
	std::ios_base::openmode om = std::ios::trunc;
	if (!saveclicks.empty()) {
		om = std::ios::app;
		vzfile = saveclicks;
	}
	ofstream file(vzfile, om);
	for (auto par : curve->params.pmap()) {
		Par* pp = par.second;
		if (!pp->is_nostate())
			width = std::max(width, pp->name.size());
	}
	os << "point " << solnidx+1 << std::endl;
	if (saveclicks.empty())
		file << "point " << solnidx+1 << std::endl;
	// update curve->params with solns[solnidx]
	curve->params.update_from_solns(solnidx);
	curve->params.freshen();
	// write each value twice: to os, and to the .vz file
	// Only display parameters that have a State
	for (auto par : curve->params.pmap()) {
		Par* pp = par.second;
		if (!(pp->is_nostate() || pp->is_eigv())) {
			os << setw(width) << pp->name << " = " << setw(12) << pp->value();
			if (++items >= perline) {
				os << std::endl;
				items = 0;
				hasnewline = true;
			} else {
				os << "       ";
				hasnewline = false;
			}
		}
		if (!saveclicks.empty())
			file << " " << pp->value();
		else
			file << setw(width) << pp->name << setw(24) << pp->value() << endl;
	}
	if (!saveclicks.empty())
		file << endl;
	if (!hasnewline)
		os << endl;
	if (saveclicks.empty())
		msg_box (curve->legend.c_str(), os.str().c_str());
	return 0;
}



static int
HandleZoom( char *progname, XButtonPressedEvent *evt, LocalWin *wi, Cursor cur) {
	Window win, new_win;
	Window root_rtn, child_rtn;
	XEvent theEvent;
	int startX, startY, curX, curY, newX, newY, stopFlag, numwin;
	int root_x, root_y;
	unsigned int mask_rtn;
	double loX, loY, hiX, hiY;

	win = evt->window;
	if (XGrabPointer(disp, win, True,
	    (unsigned int) (ButtonPressMask|ButtonReleaseMask|
	    PointerMotionMask|PointerMotionHintMask),
	    GrabModeAsync, GrabModeAsync,
	    win, cur, CurrentTime) != GrabSuccess) {
		XBell(disp, 0);
		return 0;
	}
	if (echoGC == (GC) 0) {
		unsigned long gcmask;
		XGCValues gcvals;

		gcmask = GCForeground | GCFunction;
		/* gcvals.foreground = PM_PIXEL("ZeroColor") ^ PM_PIXEL("Background"); */
		gcvals.foreground = PM_PIXEL("Foreground") ^ PM_PIXEL("Background");
		gcvals.function = GXxor;
		echoGC = XCreateGC(disp, win, gcmask, &gcvals);
	}
	startX = evt->x;  
	startY = evt->y;
	XQueryPointer(disp, win, &root_rtn, &child_rtn, &root_x, &root_y,
	    &curX, &curY, &mask_rtn);
	/* Draw first box */
	DRAWBOX;
	stopFlag = 0;
	while (!stopFlag) {
		XNextEvent(disp, &theEvent);
		switch (theEvent.xany.type) {
		case MotionNotify:
			XQueryPointer(disp, win, &root_rtn, &child_rtn, &root_x, &root_y,
			    &newX, &newY, &mask_rtn);
			/* Undraw the old one */
			DRAWBOX;
			/* Draw the new one */
			curX = newX;  
			curY = newY;
			DRAWBOX;
			break;
		case ButtonRelease:
			DRAWBOX;
			XUngrabPointer(disp, CurrentTime);
			stopFlag = 1;
			if ((startX-curX != 0) && (startY-curY != 0)) {
				/* Figure out relative bounding box */
				loX = TRANX(startX);   
				loY = TRANY(startY);
				hiX = TRANX(curX);     
				hiY = TRANY(curY);
				if (loX > hiX) {
					double temp;

					temp = hiX;
					hiX = loX;
					loX = temp;
				}
				if (loY > hiY) {
					double temp;

					temp = hiY;
					hiY = loY;
					loY = temp;
				}
				// physical aspect ratio XXX use the orig ratio (1.2) to avoid resizing
				new_win = NewWindow(progname, loX, loY, hiX, hiY, 1.2);
				if (new_win) {
					numwin = 1;
				} else {
					numwin = 0;
				}
			} else {
				numwin = 0;
			}
			break;
		default:
			printf("unknown event: %d\n", theEvent.xany.type);
			break;
		}
	}
	return numwin;
}


static int
InitSets() {
/*
 * Initializes the data sets with default information.  Sets up
 * original values for parameters in parameters package.
 */
	int idx;
	char buf[1024];

	/*
     * Used to do all kinds of searching through visuals, etc.
     * Got complaints -- so back to the simple version.
     */
	visual = DefaultVisual(disp, DefaultScreen(disp));
	cmap = DefaultColormap(disp, DefaultScreen(disp));
	screen = DefaultScreen(disp);
	depth = DefaultDepth(disp, DefaultScreen(disp));

	param_init(disp, cmap);

	param_set("Debug", BOOL, "false");
	param_set("Geometry", STR, DEF_GEOMETRY);
	param_set("ReverseVideo", BOOL, DEF_REVERSE);

	param_set("BorderSize", INT, DEF_BORDER_WIDTH);
	param_set("TitleText", STR, DEF_TITLE_TEXT);
	param_set("XUnitText", STR, DEF_XUNIT_TEXT);
	param_set("YUnitText", STR, DEF_YUNIT_TEXT); /* YUnits */
	param_set("Ticks", BOOL, DEF_TICK_FLAG);

	param_set("Markers", BOOL, DEF_MARKER_FLAG); /* markFlag (-m) */
	param_set("StyleMarkers", BOOL, DEF_DIFFMARK_FLAG); /* colorMark (-M) */
	param_set("PixelMarkers", BOOL, DEF_PIXMARK_FLAG); /* pixelMarks  (-p) */
	param_set("LargePixels", BOOL, DEF_LARGEPIX_FLAG); /* bigPixel (-P) */

	param_set("BoundBox", BOOL, DEF_BB_FLAG);
	param_set("NoLines", BOOL, DEF_NOLINE_FLAG);
	param_set("LogX", BOOL, DEF_LOGX_FLAG);
	param_set("LogY", BOOL, DEF_LOGY_FLAG); /* logYFlag */
	param_set("BarGraph", BOOL, DEF_BAR_FLAG);
	param_set("BarBase", DBL, DEF_BAR_BASE);
	param_set("BarWidth", DBL, DEF_BAR_WIDTH);
	param_set("LineWidth", INT, DEF_LINE_WIDTH);
	param_set("GridSize", INT, DEF_GRID_SIZE);
	param_set("GridStyle", STYLE, DEF_GRID_STYLE);

	param_set("Device", STR, DEF_DEVICE);
	param_set("Disposition", STR, DEF_DISPOSITION);
	param_set("FileOrDev", STR, DEF_FILEORDEV);

	/* Set the user bounding box */
	param_set("XLowLimit", DBL, DEF_LOW_LIMIT);
	param_set("YLowLimit", DBL, DEF_LOW_LIMIT);
	param_set("XHighLimit", DBL, DEF_HIGH_LIMIT);
	param_set("YHighLimit", DBL, DEF_HIGH_LIMIT);

	/* Depends critically on whether the display has color */
	if (depth < 4) {
		/* Its black and white */
		param_set("Background", PIXEL, DEF_BW_BACKGROUND);
		param_set("Border", PIXEL, DEF_BW_BORDER);
		param_set("ZeroColor", PIXEL, DEF_BW_ZEROCOLOR);
		param_set("ZeroWidth", INT, DEF_BW_ZEROWIDTH);
		param_set("ZeroStyle", STYLE, DEF_BW_ZEROSTYLE);
		param_set("Foreground", PIXEL, DEF_BW_FOREGROUND);
		/* Initialize set defaults */
		for (idx = 0;  idx < MAXATTR;  idx++) {
			(void) sprintf(buf, "%d.Style", idx);
			param_set(buf, STYLE, defStyle[idx]);
			(void) sprintf(buf, "%d.Color", idx);
			param_set(buf, PIXEL, DEF_BW_FOREGROUND);
		}
	} else {
		/* Its color */
		param_set("Background", PIXEL, DEF_COL_BACKGROUND);
		param_set("Border", PIXEL, DEF_COL_BORDER);
		param_set("ZeroColor", PIXEL, DEF_COL_ZEROCOLOR);
		param_set("ZeroWidth", INT, DEF_COL_ZEROWIDTH);
		param_set("ZeroStyle", STYLE, DEF_COL_ZEROSTYLE);
		param_set("Foreground", PIXEL, DEF_COL_FOREGROUND);
		/* Initalize attribute colors defaults */
		for (idx = 0;  idx < MAXATTR;  idx++) {
			(void) sprintf(buf, "%d.Style", idx);
			param_set(buf, STYLE, defStyle[idx]);
			(void) sprintf(buf, "%d.Color", idx);
			param_set(buf, PIXEL, defColors[idx]);
		}
	}

	param_set("LabelFont", FONT, DEF_LABEL_FONT);
	param_set("TitleFont", FONT, DEF_TITLE_FONT);

	return 0;
}



static char *def_str;

#define DEF(name, type) \
if ((def_str = XGetDefault(disp, Prog_Name, name))) { \
    param_set(name, type, def_str); \
}

static int
ReadDefaults() {
/*
 * Reads X default values which override the hard-coded defaults
 * set up by InitSets.
 */
	char newname[100];
	int idx;

	DEF("Debug", BOOL);
	DEF("Geometry", STR);
	DEF("Background", PIXEL);
	DEF("BorderSize", INT);
	DEF("Border", PIXEL);
	DEF("GridSize", INT);
	DEF("GridStyle", STYLE);
	DEF("Foreground", PIXEL);
	DEF("ZeroColor", PIXEL);
	DEF("ZeroStyle", STYLE);
	DEF("ZeroWidth", INT);
	DEF("LabelFont", FONT);
	DEF("TitleFont", FONT);
	DEF("Ticks", BOOL);
	DEF("Device", STR);
	DEF("Disposition", STR);
	DEF("FileOrDev", STR);
	DEF("PixelMarkers", BOOL);
	DEF("LargePixels", BOOL);
	DEF("Markers", BOOL);
	DEF("StyleMarkers", BOOL);
	DEF("BoundBox", BOOL);
	DEF("NoLines", BOOL);
	DEF("LineWidth", INT);

	/* Read device specific parameters */
	for (idx = 0;  idx < hard_count;  idx++) {
		sprintf(newname, "%s.Dimension", hard_devices[idx].dev_name);
		DEF(newname, DBL);	/* hard_devices[idx].dev_max_dim */
		sprintf(newname, "%s.OutputTitleFont", hard_devices[idx].dev_name);
		DEF(newname, STR);	/* hard_devices[idx].dev_title_font */
		sprintf(newname, "%s.OutputTitleSize", hard_devices[idx].dev_name);
		DEF(newname, DBL);	/* hard_devices[idx].dev_title_size */
		sprintf(newname, "%s.OutputAxisFont", hard_devices[idx].dev_name);
		DEF(newname, STR);	/* hard_devices[idx].dev_axis_font */
		sprintf(newname, "%s.OutputAxisSize", hard_devices[idx].dev_name);
		DEF(newname, DBL);	/* hard_devices[idx].dev_axis_size */
	}


	/* Read the default line and color attributes */
	for (idx = 0;  idx < MAXATTR;  idx++) {
		(void) sprintf(newname, "%d.Style", idx);
		DEF(newname, STYLE);	/* AllAttrs[idx].lineStyleLen */
		(void) sprintf(newname, "%d.Color", idx);
		DEF(newname, PIXEL);	/* AllAttrs[idx].pixelValue */
	}

	DEF("ReverseVideo", BOOL);
	return 0;
}


#define FS(str)	(void) fprintf(stderr, str)

int argerror(char* err, char* val) {
	(void) fprintf(stderr, "Error: %s: %s\n\n", val, err);

	FS("Usage: vz [-<digit> set_name] [-bar] [-bb] [-bd border_color]\n");
	FS("         [-bg background_color] [-brb bar_base] [-brw bar_width]\n");
	FS("         [-bw bdr_width] [-db] [-fg foreground_color] [-gw grid_size]\n");
	FS("         [-gs grid_style] [-lf label_font] [-lnx] [-lny] [-lw line_width]\n");
	FS("         [-lx x1,x2] [-ly y1,y2] [-m] [-M] [-nl] [-p] [-P] [-rv]\n");
	FS("         [-t title] [-tf title_font] [-tk] [-x x_unit_name]\n");
	FS("         [-y y_unit_name] [-zg zero_color] [-zw zero_size]\n");
	FS("         [=WxH+X+Y] [-display <host>:<disp>.<screen>] file...\n\n");
	FS("-bar   Draw bar graph with base -brb and width -brw\n");
	FS("-bb    Draw bounding box around data\n");
	FS("-db    Turn on debugging\n");
	FS("-lnx   Logarithmic scale for X axis\n");
	FS("-lny   Logarithmic scale for Y axis\n");
	FS("-m -M  Mark points distinctively (M varies with color)\n");
	FS("-nl    Don't draw lines (scatter plot)\n");
	FS("-p -P  Mark points with dot (P means big dot)\n");
	FS("-rv    Reverse video on black and white displays\n");
	FS("-tk    Draw tick marks instead of full grid\n");

	exit(1);
}


static int
DrawWindow(LocalWin *win_info) {
/*
 * Draws the data in the window.  Does not clear the window.
 * The data is scaled so that all of the data will fit.
 * Grid lines are drawn at the nearest power of 10 in engineering
 * notation.  Draws axis numbers along bottom and left hand edges.
 * Centers title at top of window.
 */
	/* Figure out the transformation constants */
	if (TransformCompute(win_info)) {

		if (!Hardcopy) {
			/* Draw the title */
			DrawTitle(win_info);

			/* Draw the legend */
			DrawLegend(win_info);

		}
			/* Draw the axis unit labels,  grid lines,  and grid labels */
		// move this line to the !Hardcopy block if grid not wanted in xfig
			DrawGridAndAxis(win_info);

		/* Draw the data sets themselves */
		DrawData(win_info);
	}
	return 0;
}




static void
DrawTitle (LocalWin *wi) {
/*
 * This routine draws the title of the graph centered in
 * the window.  It is spaced down from the top by an amount
 * specified by the constant PADDING.  The font must be
 * fixed width.  The routine returns the height of the
 * title in pixels.
 */
	Specs& sp = specs();

	wi->dev_info.xg_text(wi->dev_info.user_state,
	    wi->dev_info.area_w/2,
	    wi->dev_info.axis_pad,
	    sp.title.c_str(), T_TOP, T_TITLE);
}




static int
TransformCompute (LocalWin *wi)
/*
 * This routine figures out how to draw the axis labels and grid lines.
 * Both linear and logarithmic axes are supported.  Axis labels are
 * drawn in engineering notation.  The power of the axes are labeled
 * in the normal axis labeling spots.  The routine also figures
 * out the necessary transformation information for the display
 * of the points (it touches XOrgX, XOrgY, UsrOrgX, UsrOrgY, and
 * UnitsPerPixel).
 */
{
	double bbCenX, bbCenY, bbHalfWidth, bbHalfHeight;
	int tempSize, maxName, leftWidth = 0;
	/*
     * First,  we figure out the origin in the X window.  Above
     * the space we have the title and the Y axis unit label.
     * To the left of the space we have the Y axis grid labels.
     */

	wi->XOrgX = wi->dev_info.bdr_pad + (7 * wi->dev_info.axis_width)
	    + wi->dev_info.bdr_pad;
	wi->XOrgY = wi->dev_info.bdr_pad + wi->dev_info.title_height
	    + wi->dev_info.bdr_pad + wi->dev_info.axis_height
	    + wi->dev_info.axis_height/2 + wi->dev_info.bdr_pad;

	/*
     * Now we find the lower right corner.  Below the space we
     * have the X axis grid labels.  To the right of the space we
     * have the X axis unit label and the legend.  We assume the 
     * worst case size for the unit label.
     */

	maxName = 0;

	for (size_t idx = 0;  idx < Curves.size();  idx++) {
		if (Curves[idx]->xparam == nullptr) continue;
		tempSize = (int)Curves[idx]->legend.size();
		maxName = MAX(tempSize, maxName);
	}

	tempSize = (maxName+1)*wi->dev_info.axis_width + 2 + wi->dev_info.bdr_pad;
	if (tempSize > leftWidth)
		leftWidth = tempSize;

	wi->XOppX = wi->dev_info.area_w - wi->dev_info.bdr_pad - leftWidth;
	wi->XOppY = wi->dev_info.area_h - wi->dev_info.bdr_pad
	    - 3*(wi->dev_info.axis_height) - wi->dev_info.bdr_pad;

	if ((wi->XOrgX >= wi->XOppX) || (wi->XOrgY >= wi->XOppY)) {
		do_error("Drawing area is too small\n");
		return 0;
	}

	/* 
     * We now have a bounding box for the drawing region.
     * Figure out the units per pixel using the data set bounding box.
     */
	wi->XUnitsPerPixel = (wi->hiX - wi->loX)/((double) (wi->XOppX - wi->XOrgX));
	wi->YUnitsPerPixel = (wi->hiY - wi->loY)/((double) (wi->XOppY - wi->XOrgY));

	/*
     * Find origin in user coordinate space.  We keep the center of
     * the original bounding box in the same place.
     */
	bbCenX = (wi->loX + wi->hiX) / 2.0;
	bbCenY = (wi->loY + wi->hiY) / 2.0;
	bbHalfWidth = ((double) (wi->XOppX - wi->XOrgX))/2.0 * wi->XUnitsPerPixel;
	bbHalfHeight = ((double) (wi->XOppY - wi->XOrgY))/2.0 * wi->YUnitsPerPixel;
	wi->UsrOrgX = bbCenX - bbHalfWidth;
	wi->UsrOrgY = bbCenY - bbHalfHeight;
	wi->UsrOppX = bbCenX + bbHalfWidth;
	wi->UsrOppY = bbCenY + bbHalfHeight;

	/*
     * Everything is defined so we can now use the SCREENX and SCREENY
     * transformations.
     */
	return 1;
}

static int
DrawGridAndAxis(LocalWin *wi) {
/*
 * This routine draws grid line labels in engineering notation,
 * the grid lines themselves,  and unit labels on the axes.
 */
	Trace trc(1,"DrawGridAndAxis");
	Specs& sp = specs();
	int expX, expY;		/* Engineering powers */
	int startX, startY;
	int Yspot, Xspot;
	char power[10], value[10];
	double Xincr, Yincr, Xstart, Ystart, Yindex, Xindex, larger;
	XSegment segs[2];
	int tickFlag = PM_BOOL("Ticks");
#ifdef NEVER // use specs
	int logXFlag = PM_BOOL("LogX");
	int logYFlag = PM_BOOL("LogY");
#else // NEVER // use specs
	int logXFlag = sp.xlog ? 1 : 0;
	int logYFlag = sp.ylog ? 1 : 0;
#endif // NEVER // use specs
	// char *XUnitText = pm_str("XUnitText");
	// char *YUnitText = pm_str("YUnitText");
	string xunittext;
	Settings::defaults.get("xunittext", xunittext);
	string yunittext;
	Settings::defaults.get("yunittext", yunittext);

	/*
     * Grid display powers are computed by taking the log of
     * the largest numbers and rounding down to the nearest
     * multiple of 3.
     */
	if (logXFlag) {
		expX = 0;
	} else {
		if (fabs(wi->UsrOrgX) > fabs(wi->UsrOppX)) {
			larger = fabs(wi->UsrOrgX);
		} else {
			larger = fabs(wi->UsrOppX);
		}
		expX = ((int) floor(nlog10(larger)/3.0)) * 3;
	}
	if (logYFlag) {
		expY = 0;
	} else {
		if (fabs(wi->UsrOrgY) > fabs(wi->UsrOppY)) {
			larger = fabs(wi->UsrOrgY);
		} else {
			larger = fabs(wi->UsrOppY);
		}
		expY = ((int) floor(nlog10(larger)/3.0)) * 3;
	}

	/*
     * With the powers computed,  we can draw the axis labels.
     */
	if (expY != 0) {
		string final = vastr(yunittext, " x 10");
		Xspot = wi->dev_info.bdr_pad +
		    ((yunittext.size()+5) * wi->dev_info.axis_width);
		Yspot = wi->dev_info.bdr_pad * 2 + wi->dev_info.title_height +
		    wi->dev_info.axis_height;
		wi->dev_info.xg_text(wi->dev_info.user_state,
		    Xspot, Yspot, final.c_str(), T_RIGHT, T_AXIS);
		(void) sprintf(power, "%d", expY);
		wi->dev_info.xg_text(wi->dev_info.user_state,
		    Xspot, Yspot, power, T_LOWERLEFT, T_AXIS);
	} else {
		Yspot = wi->dev_info.bdr_pad * 2 + wi->dev_info.title_height +
		    wi->dev_info.axis_height/2;
		wi->dev_info.xg_text(wi->dev_info.user_state,
		    wi->dev_info.bdr_pad, Yspot, yunittext.c_str(),
		    T_UPPERLEFT, T_AXIS);
	}
/*
 * X-axis unit label
 */
	/* startX = (wi->dev_info.area_w)/2 + strlen(XUnitText); */
	startX = (wi->XOppY + wi->XOrgY +
		xunittext.size()*wi->dev_info.axis_width)/2;
	startY = wi->XOppY + 3*(wi->dev_info.axis_height) - 2*wi->dev_info.bdr_pad;
	if (expX != 0) {
		(void) sprintf(power, "%d", expX);
		startX -= (strlen(power) * wi->dev_info.axis_width);
		wi->dev_info.xg_text(wi->dev_info.user_state,
		    startX, startY, power, T_LOWERLEFT, T_AXIS);
		string final = vastr(xunittext, " x 10");
		wi->dev_info.xg_text(wi->dev_info.user_state,
		    startX, startY, final.c_str(), T_RIGHT, T_AXIS);
	} else {
		wi->dev_info.xg_text(wi->dev_info.user_state,
		    startX, startY, xunittext.c_str(), T_RIGHT, T_AXIS);
	}

	/* 
     * First,  the grid line labels
     */
	Yincr = (wi->dev_info.axis_pad + wi->dev_info.axis_height) * wi->YUnitsPerPixel;
	Ystart = initGrid(wi->UsrOrgY, Yincr, logYFlag);
	for (Yindex = Ystart;  Yindex < wi->UsrOppY;  Yindex = stepGrid()) {
		Yspot = SCREENY(wi, Yindex);
		/* Write the axis label */
		WriteValue(value, Yindex, expY, logYFlag);
		wi->dev_info.xg_text(wi->dev_info.user_state,
		    wi->dev_info.bdr_pad +
		    (7 * wi->dev_info.axis_width),
		    Yspot, value, T_RIGHT, T_AXIS);
	}

	Xincr = (wi->dev_info.axis_pad + (wi->dev_info.axis_width * 7)) * wi->XUnitsPerPixel;
	Xstart = initGrid(wi->UsrOrgX, Xincr, logXFlag);
	Yspot = wi->XOppY + (wi->dev_info.axis_height);
	for (Xindex = Xstart;  Xindex < wi->UsrOppX;  Xindex = stepGrid()) {
		Xspot = SCREENX(wi, Xindex);
		/* Write the axis label */
		WriteValue(value, Xindex, expX, logXFlag);
		wi->dev_info.xg_text(wi->dev_info.user_state,
		    Xspot, Yspot,
		    value, T_BOTTOM, T_AXIS);
	}

	/*
     * Now,  the grid lines or tick marks
     */
	Yincr = (wi->dev_info.axis_pad + wi->dev_info.axis_height) * wi->YUnitsPerPixel;
	Ystart = initGrid(wi->UsrOrgY, Yincr, logYFlag);
	for (Yindex = Ystart;  Yindex < wi->UsrOppY;  Yindex = stepGrid()) {
		Yspot = SCREENY(wi, Yindex);
		/* Draw the grid line or tick mark */
		if (tickFlag) {
			segs[0].x1 = wi->XOrgX;
			segs[0].x2 = wi->XOrgX + wi->dev_info.tick_len;
			segs[1].x1 = wi->XOppX - wi->dev_info.tick_len;
			segs[1].x2 = wi->XOppX;
			segs[0].y1 = segs[0].y2 = segs[1].y1 = segs[1].y2 = Yspot;
		} else {
			segs[0].x1 = wi->XOrgX;  
			segs[0].x2 = wi->XOppX;
			segs[0].y1 = segs[0].y2 = Yspot;
		}
		if ((ABS(Yindex) < ZERO_THRES) && !logYFlag) {
			wi->dev_info.xg_seg(wi->dev_info.user_state,
			    1, segs, PM_INT("ZeroWidth"),
			    L_ZERO, 0, 0);
			if (tickFlag) {
				wi->dev_info.xg_seg(wi->dev_info.user_state,
				    1, &(segs[1]), PM_INT("ZeroWidth"),
				    L_ZERO, 0, 0);
			}
		} else {
			wi->dev_info.xg_seg(wi->dev_info.user_state,
			    1, segs, PM_INT("GridSize"),
			    L_AXIS, 0, 0);
			if (tickFlag) {
				wi->dev_info.xg_seg(wi->dev_info.user_state,
				    1, &(segs[1]), PM_INT("GridSize"),
				    L_AXIS, 0, 0);
			}
		}
	}

	Xincr = (wi->dev_info.axis_pad + (wi->dev_info.axis_width * 7)) * wi->XUnitsPerPixel;
	Xstart = initGrid(wi->UsrOrgX, Xincr, logXFlag);
	for (Xindex = Xstart;  Xindex < wi->UsrOppX;  Xindex = stepGrid()) {
		Xspot = SCREENX(wi, Xindex);
		/* Draw the grid line or tick marks */
		if (tickFlag) {
			segs[0].x1 = segs[0].x2 = segs[1].x1 = segs[1].x2 = Xspot;
			segs[0].y1 = wi->XOrgY;
			segs[0].y2 = wi->XOrgY + wi->dev_info.tick_len;
			segs[1].y1 = wi->XOppY - wi->dev_info.tick_len;
			segs[1].y2 = wi->XOppY;
		} else {
			segs[0].x1 = segs[0].x2 = Xspot;
			segs[0].y1 = wi->XOrgY; 
			segs[0].y2 = wi->XOppY;
		}
		if ((ABS(Xindex) < ZERO_THRES) && !logXFlag) {
			wi->dev_info.xg_seg(wi->dev_info.user_state,
			    1, segs, PM_INT("ZeroWidth"), L_ZERO, 0, 0);
			if (tickFlag) {
				wi->dev_info.xg_seg(wi->dev_info.user_state,
				    1, &(segs[1]), PM_INT("ZeroWidth"),
				    L_ZERO, 0, 0);

			}
		} else {
			wi->dev_info.xg_seg(wi->dev_info.user_state,
			    1, segs, PM_INT("GridSize"), L_AXIS, 0, 0);
			if (tickFlag) {
				wi->dev_info.xg_seg(wi->dev_info.user_state,
				    1, &(segs[1]), PM_INT("GridSize"), L_AXIS, 0, 0);
			}
		}
	}
	/* Check to see if he wants a bounding box */
	if (PM_BOOL("BoundBox")) {
		XSegment bb[4];

		/* Draw bounding box */
		bb[0].x1 = bb[0].x2 = bb[1].x1 = bb[3].x2 = wi->XOrgX;
		bb[0].y1 = bb[2].y2 = bb[3].y1 = bb[3].y2 = wi->XOrgY;
		bb[1].x2 = bb[2].x1 = bb[2].x2 = bb[3].x1 = wi->XOppX;
		bb[0].y2 = bb[1].y1 = bb[1].y2 = bb[2].y1 = wi->XOppY;
		wi->dev_info.xg_seg(wi->dev_info.user_state,
		    4, bb, PM_INT("GridSize"), L_AXIS, 0, 0);
	}
	return 0;
}

static double gridBase, gridStep, gridJuke[101];
static int gridNJuke, gridCurJuke;

#define ADD_GRID(val)	(gridJuke[gridNJuke++] = log10(val))

static double
initGrid(
double low,			/* desired low value          */
double step,			/* desired step (user coords) */
int logFlag)			/* is axis logarithmic?       */
{
	double ratio, x;

	gridNJuke = gridCurJuke = 0;
	gridJuke[gridNJuke++] = 0.0;

	if (logFlag) {
		ratio = pow(10.0, step);
		gridBase = floor(low);
		gridStep = ceil(step);
		if (ratio <= 3.0) {
			if (ratio > 2.0) {
				ADD_GRID(3.0);
			} else if (ratio > 1.333) {
				ADD_GRID(2.0);	
				ADD_GRID(5.0);
			} else if (ratio > 1.25) {
				ADD_GRID(1.5);	
				ADD_GRID(2.0);	
				ADD_GRID(3.0);
				ADD_GRID(5.0);	
				ADD_GRID(7.0);
			} else {
				for (x = 1.0; x < 10.0 && (x+.5)/(x+.4) >= ratio; x += .5) {
					ADD_GRID(x + .1);	
					ADD_GRID(x + .2);
					ADD_GRID(x + .3);	
					ADD_GRID(x + .4);
					ADD_GRID(x + .5);
				}
				if (floor(x) != x) ADD_GRID(x += .5);
				for ( ; x < 10.0 && (x+1.0)/(x+.5) >= ratio; x += 1.0) {
					ADD_GRID(x + .5);	
					ADD_GRID(x + 1.0);
				}
				for ( ; x < 10.0 && (x+1.0)/x >= ratio; x += 1.0) {
					ADD_GRID(x + 1.0);
				}
				if (x == 7.0) {
					gridNJuke--;
					x = 6.0;
				}
				if (x < 7.0) {
					ADD_GRID(x + 2.0);
				}
				if (x == 10.0) gridNJuke--;
			}
			x = low - gridBase;
			//for (gridCurJuke = -1; x >= gridJuke[gridCurJuke+1]; gridCurJuke++){ }
			for (int i=0; i<gridNJuke; i++){
				if (x < gridJuke[i]) {
					gridCurJuke = i-1;
					break;
				}
			}
		}
	} else {
		gridStep = RoundUp(step);
		gridBase = floor(low / gridStep) * gridStep;
	}
	return(stepGrid());
}

static double
stepGrid() {
	if (++gridCurJuke >= gridNJuke) {
		gridCurJuke = 0;
		gridBase += gridStep;
	}
	return(gridBase + gridJuke[gridCurJuke]);
}

static double
RoundUp(double val)
/*
 * This routine rounds up the given positive number such that
 * it is some power of ten times either 1, 2, or 5.  It is
 * used to find increments for grid lines.
 */
{
	int exponent, idx;
	double eps = FLT_MIN;

	if (val < eps)
		return 1.0;
	if (val > 1.0/eps)
		return 1.0;

	exponent = (int) floor(nlog10(val));
	if (exponent < 0) {
		for (idx = exponent;  idx < 0; idx++) {
			val *= 10.0;
		}
	} else {
		for (idx = 0;  idx < exponent; idx++) {
			val /= 10.0;
		}
	}
	if (val > 5.0) val = 10.0;
	else if (val > 2.0) val = 5.0;
	else if (val > 1.0) val = 2.0;
	else val = 1.0;
	if (exponent < 0) {
		for (idx = exponent;  idx < 0;  idx++) {
			val /= 10.0;
		}
	} else {
		for (idx = 0;  idx < exponent;  idx++) {
			val *= 10.0;
		}
	}
	return val;
}

static int
WriteValue(
	char *str,			/* String to write into */
	double val,			/* Value to print       */
	int expv,			/* Exponent             */
	int logFlag) {			/* Is this a log axis?  */
/*
 * Writes the value provided into the string in a fixed format
 * consisting of seven characters.  The format is:
 *   -ddd.dd
 */
	int idx;

	if (logFlag) {
		if (val == floor(val)) {
			(void) sprintf(str, "%.0e", pow(10.0, val));
		} else {
			(void) sprintf(str, "%.2g", pow(10.0, val - floor(val)));
		}
	} else {
		if (expv < 0) {
			for (idx = expv;  idx < 0;  idx++) {
				val *= 10.0;
			}
		} else {
			for (idx = 0;  idx < expv;  idx++) {
				val /= 10.0;
			}
		}
		(void) sprintf(str, "%.2f", val);
	}
	return 0;
}

bool
inRange(vector<Par*> const& parlist, int a, int b) {
// are any parameters out of range at both points a and b?
	for (size_t i=0; i<parlist.size(); i++) {
		if (parlist[i]->is_fixed())
			continue;
		parlist[i]->update_from_solns(a);
		bool rval = parlist[i]->inRange();
		if (rval)
			continue;
		parlist[i]->update_from_solns(b);
		if (!parlist[i]->inRange()) {
			return false;
		}
	}
	return true;
}

#define LEFT_CODE	0x01
#define RIGHT_CODE	0x02
#define BOTTOM_CODE	0x04
#define TOP_CODE	0x08

/* Clipping algorithm from Neumann and Sproull by Cohen and Sutherland */
#define C_CODE(xval, yval, rtn) \
rtn = 0; \
if ((xval) < wi->UsrOrgX) rtn = LEFT_CODE; \
else if ((xval) > wi->UsrOppX) rtn = RIGHT_CODE; \
if ((yval) < wi->UsrOrgY) rtn |= BOTTOM_CODE; \
else if ((yval) > wi->UsrOppY) rtn |= TOP_CODE

static int
DrawData(LocalWin* wi) {
// This routine draws the data sets themselves using the macros
// for translating coordinates.
	Trace trc(1,"DrawData");
	double sx1{0.0}, sy1{0.0}, sx2{0.0}, sy2{0.0};
	double tx{0.0}, ty{0.0};
	double wval = 0.0;
	int code1, code2, cd;
	int markFlag, pixelMarks, bigPixel, colorMark;
	int lineWidth = PM_INT("LineWidth");
	bool isTagged;
	size_t i;

	set_mark_flags(&markFlag, &pixelMarks, &bigPixel, &colorMark);
	// idx: loop over each curve
	for (size_t idx = 0;  idx < Curves.size();  idx++) {
		trc.dprint("Curves[",idx,"]: ",Curves[idx]->summary());
		Par* xp = Curves[idx]->xparam;
		Par* wp = Curves[idx]->wparam;

		vector<XSegment> xsegs;		/* Point space for X */
		Par* yp = Curves[idx]->yparam;
		int color = Curves[idx]->color;
		int wcolor = color+1;
		int currentColor = color;
		int style = LINESTYLE(Curves[idx]->style);
		int nseg = (int)yp->nsolns()-1;
		trc.dprint("y ",yp->name,", color ",color,", style ",style, ", number of segments: ",nseg);

		// subindex: loop over each x value, adding to xsegs until
		// either all x's processed or w value changes sign
		for (int subindex = 0; subindex<nseg; subindex++) {
			if (wp) {
				// if color has changed dump xsegs
				wval = wp->solns[subindex];
				if (subindex > 1 && !xsegs.empty() &&
						wval*wp->solns[subindex-1] < 0.0) {
					trc.dprint("w sign change: dump ",xsegs.size()," segs");
					wi->dev_info.xg_seg(wi->dev_info.user_state,
						 xsegs.size(), &xsegs[0], lineWidth, L_VAR,
						 style, currentColor);
					// mark data points XXX move to xg_seg's?
					int mark_style = P_MARK;
					int type = MARKSTYLE(idx);
					if (wval > 0.0) {
						//mark_style = P_TAG;
						type = MARKSTYLE(idx+2);
					}
						/* Distinctive markers */
					for (i=0; i<xsegs.size(); i++)
						wi->dev_info.xg_dot(wi->dev_info.user_state,
							 xsegs[i].x1, xsegs[i].y1,
							 mark_style, type, currentColor);
					// and delete accumulated segments
					xsegs.clear();
				}
				// set w color, no change if wval=0
				if (wval < 0.0) {
					currentColor = wcolor;
				} else if (wval > 0.0) {
					currentColor = color;
				}
				trc.dprint("w value ",wval,", using color ",currentColor);
			}
			lineWidth = PM_INT("LineWidth");

			/* Put segment in (sx1,sy1) (sx2,sy2) */
			if (xp == nullptr) {
				sx1 = subindex;
				sx2 = subindex+1;
			} else {
				sx1 = xp->solns[subindex];
				sx2 = xp->solns[subindex+1];
			}
			sy1 = yp->solns[subindex];
			sy2 = yp->solns[subindex+1];
			trc.dprint("seg ",subindex,": x(",sx1,',',sx2,"), y(",sy1,',',sy2,")");
			isTagged = Curves[idx]->tagged[subindex];
			if (isTagged)
				trc.dprint("idx ",idx,", subindex ",subindex," is tagged");
			/* Now clip to current window boundary */
			C_CODE(sx1, sy1, code1);
			C_CODE(sx2, sy2, code2);
			while (code1 || code2) {
				if (code1 & code2) break;
				cd = (code1 ? code1 : code2);
				if (cd & LEFT_CODE) {	/* Crosses left edge */
					ty = sy1 + (sy2 - sy1) * (wi->UsrOrgX - sx1) / (sx2 - sx1);
					tx = wi->UsrOrgX;
				} else if (cd & RIGHT_CODE) { /* Crosses right edge */
					ty = sy1 + (sy2 - sy1) * (wi->UsrOppX - sx1) / (sx2 - sx1);
					tx = wi->UsrOppX;
				} else if (cd & BOTTOM_CODE) { /* Crosses bottom edge */
					tx = sx1 + (sx2 - sx1) * (wi->UsrOrgY - sy1) / (sy2 - sy1);
					ty = wi->UsrOrgY;
				} else if (cd & TOP_CODE) { /* Crosses top edge */
					tx = sx1 + (sx2 - sx1) * (wi->UsrOppY - sy1) / (sy2 - sy1);
					ty = wi->UsrOppY;
				}
				if (cd == code1) {
					sx1 = tx;  
					sy1 = ty;
					C_CODE(sx1, sy1, code1);
				} else {
					sx2 = tx;  
					sy2 = ty;
					C_CODE(sx2, sy2, code2);
				}
			}
			if (!code1 && !code2) {
				/* Add segment to list */
				XSegment xs;
				xs.x1 = SCREENX(wi, sx1);
				xs.y1 = SCREENY(wi, sy1);
				xs.x2 = SCREENX(wi, sx2);
				xs.y2 = SCREENY(wi, sy2);
				xsegs.push_back(xs);
			}
		}

		// Draw last segments
		trc.dprint("finishing with ",xsegs.size()," segments");
		wi->dev_info.xg_seg(wi->dev_info.user_state, xsegs.size(),
				&xsegs[0], lineWidth, L_VAR, style, currentColor);

		// ... and their markers XXX move to xg_seg's?
		int mark_style = P_MARK;
		int type = MARKSTYLE(idx);
		if (wval > 0.0) {
			// mark_style = P_TAG;
			type = MARKSTYLE(idx+2);
		}
		// Distinctive markers
		for (i=0; i<xsegs.size(); i++)
			wi->dev_info.xg_dot(wi->dev_info.user_state,
				 xsegs[i].x1, xsegs[i].y1,
				 mark_style, type, currentColor);
	}
	return 0;
}



static int
DrawLegend(LocalWin* wi) {
// This draws a legend of the data sets displayed.  Only those that
// will fit are drawn.
	Trace trc(1,"DrawLegend");
	int spot;
	size_t idx;
	XSegment leg_line;
	int markFlag, pixelMarks, bigPixel, colorMark;
	string file;
	string lastfile;

	set_mark_flags(&markFlag, &pixelMarks, &bigPixel, &colorMark);
	/*
	 * Compute how long to make the lines
	 */
	int maxleg{0};
	for (auto curve : Curves) {
		maxleg = std::max(maxleg, (int)curve->legend.size());
	}
	int linelen = maxleg * wi->dev_info.axis_width;
	trc.dprint("max legend ",maxleg,", linelen ",linelen);

	leg_line.x1 = wi->XOppX + wi->dev_info.bdr_pad;
	leg_line.x2 = leg_line.x1 + linelen;
	spot = wi->XOrgY;

	for (idx = 0;  idx < Curves.size();  idx++) {
		if (!Curves[idx]->xparam) continue;
		/*
		 * the line and mark...
		 */
		int color = Curves[idx]->color;
		int style = Curves[idx]->style;
		if ( (spot + wi->dev_info.axis_height + 2 < wi->XOppY)) {
			/* Meets the criteria */
			leg_line.y1 = leg_line.y2 = spot - wi->dev_info.legend_pad;
			wi->dev_info.xg_seg(wi->dev_info.user_state,
			    1, &leg_line, 1, L_VAR, style, color);
			if (markFlag && !pixelMarks) {
				wi->dev_info.xg_dot(wi->dev_info.user_state,
				    leg_line.x1, leg_line.y1,
				    P_MARK, MARKSTYLE(idx), color);

			}
			/* spot += 2 + wi->dev_info.axis_height + wi->dev_info.bdr_pad; */
			/*
			 * ...then the legend
			 */
			wi->dev_info.xg_text(wi->dev_info.user_state,
			    wi->XOppX + wi->dev_info.bdr_pad,
			    spot+2,
				 Curves[idx]->legend.c_str(),
			    // Curves[idx]->cid().c_str(),
			    T_UPPERLEFT, T_AXIS);
			spot += 2 + wi->dev_info.axis_height + wi->dev_info.bdr_pad;
		}
	}
	return 0;
}



static void
set_mark_flags( int *markFlag, int *pixelMarks, int *bigPixel, int *colorMark)
/*
 * Determines the values of the old boolean flags based on the
 * new values in the parameters database.
 */
{
	*markFlag = 0;  
	*pixelMarks = 0;  
	*colorMark = 0;  
	*bigPixel = 0;
	if (PM_BOOL("Markers")) {
		*markFlag = 1;  
		*pixelMarks = 0;  
		*colorMark = 0;
	}
	if (PM_BOOL("PixelMarkers")) {
		*markFlag = 1;  
		*pixelMarks = 1;  
		*bigPixel = 0;
	}
	if (PM_BOOL("LargePixels")) {
		*markFlag = 1;  
		*pixelMarks = 1;  
		*bigPixel = 1;
	}
	if (PM_BOOL("StyleMarkers")) {
		*markFlag = 1;  
		*pixelMarks = 0;  
		*colorMark = 1;
	}
}



#define RND(val)	((int) ((val) + 0.5))

/*ARGSUSED*/
void
do_hardcopy(
char *prog,			/* Program name for Xdefaults    */
xtb_data info,			/* Some state information        */
int (*init_fun)(FILE*,int,int,char*,double,char*,double,int,xgOut*,char*),
const char *dev_spec,			/* Device specification (if any) */
char *file_or_dev,		/* Filename or device spec       */
double maxdim,			/* Maximum dimension in cm       */
char *ti_fam, double ti_size, char* ax_fam, double ax_size,
int doc_p)			/* Documentation predicate       */
/*
 * This routine resets the function pointers to those specified
 * by `init_fun' and causes a screen redisplay.  If `dev_spec'
 * is non-zero,  it will be considered a sprintf string with
 * one %s which will be filled in with `file_or_dev' and fed
 * to popen(3) to obtain a stream.  Otherwise,  `file_or_dev'
 * is considered to be a file and is opened for writing.  The
 * resulting stream is fed to the initialization routine for
 * the device.
 */
{
	LocalWin *curWin = (LocalWin *) info;
	LocalWin thisWin;
	FILE *out_stream;
	char buf[MAXBUFSIZE];
	// char err[MAXBUFSIZE];
	char ierr[ERRBUFSIZE];
	char tilde[MAXBUFSIZE*10];
	int final_w, final_h, flags;
	double ratio;

	if (dev_spec) {
		(void) sprintf(buf, dev_spec, file_or_dev);
		out_stream = popen(buf, "w");
		if (!out_stream) {
			// sprintf(err, "Unable to issue command:\n  %s\n", buf);
			do_error(vastr("Unable to issue command \"",buf,"\""));
			return;
		}
	} else {
		tildeExpand(tilde, file_or_dev);
		out_stream = fopen(tilde, "w");
		if (!out_stream) {
			// sprintf(err, "Unable to open file `%s'\n", tilde);
			do_error(vastr("Unable to open file ",tilde));
			return;
		}
	}
	thisWin = *curWin;
	ratio = ((double) thisWin.dev_info.area_w) /
	    ((double) thisWin.dev_info.area_h);
	if (thisWin.dev_info.area_w > thisWin.dev_info.area_h) {
		final_w = RND(maxdim * 10000.0);
		final_h = RND(maxdim/ratio * 10000.0);
	} else {
		final_w = RND(maxdim * ratio * 10000.0);
		final_h = RND(maxdim * 10000.0);
	}
	ierr[0] = '\0';
	flags = 0;
	if (doc_p) flags |= D_DOCU;
	if ((*init_fun)(out_stream, final_w, final_h, ti_fam, ti_size,
	    ax_fam, ax_size, flags, &(thisWin.dev_info), ierr)) {
		Hardcopy = true;
		DrawWindow(&thisWin);
		Hardcopy = false;
		if (thisWin.dev_info.xg_end) {
			thisWin.dev_info.xg_end(thisWin.dev_info.user_state);
		}
	} else {
		do_error(ierr);
	}
	if (dev_spec) {
		(void) pclose(out_stream);
	} else {
		(void) fclose(out_stream);
	}
}


static char *
tildeExpand(
char *out,			/* Output space for expanded file name */
char *in)			/* Filename with tilde                 */
/*
 * This routine expands out a file name passed in `in' and places
 * the expanded version in `out'.  It returns `out'.
 */
{
	char username[50], *userPntr;
	struct passwd *userRecord;

	out[0] = '\0';

	/* Skip over the white space in the initial path */
	while ((*in == ' ') || (*in == '\t')) in++;

	/* Tilde? */
	if (in[0] == TILDE) {
		/* Copy user name into 'username' */
		in++;  
		userPntr = &(username[0]);
		while ((*in != '\0') && (*in != '/')) {
			*(userPntr++) = *(in++);
		}
		*(userPntr) = '\0';
		/* See if we have to fill in the user name ourselves */
		if (strlen(username) == 0) {
			userRecord = getpwuid(getuid());
		} else {
			userRecord = getpwnam(username);
		}
		if (userRecord) {
			/* Found user in passwd file.  Concatenate user directory */
			strcat(out, userRecord->pw_dir);
		}
	}

	/* Concantenate remaining portion of file name */
	strcat(out, in);
	return out;
}



#define ERR_MSG_SIZE	2048

/*ARGSUSED*/
static int
XErrHandler(Display *disp_ptr, XErrorEvent *evt) {
/*
 * Displays a nicely formatted message and core dumps.
 */
	char err_buf[ERR_MSG_SIZE], mesg[ERR_MSG_SIZE], number[ERR_MSG_SIZE];
	const char *mtype = "XlibMessage";

	XGetErrorText(disp_ptr, evt->error_code, err_buf, ERR_MSG_SIZE);
	(void) fprintf(stderr, "X Error: %s\n", err_buf);
	XGetErrorDatabaseText(disp_ptr, mtype, "MajorCode",
	    "Request Major code %d", mesg, ERR_MSG_SIZE);
	(void) fprintf(stderr, mesg, evt->request_code);
	(void) sprintf(number, "%d", evt->request_code);
	XGetErrorDatabaseText(disp_ptr, "XRequest", number, "", err_buf,
	    ERR_MSG_SIZE);
	(void) fprintf(stderr, " (%s)\n", err_buf);

	/* abort(); */
	return 0;
}

