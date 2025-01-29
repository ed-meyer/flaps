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
 * Dialog Boxes
 *
 * This file constructs the hardcopy and error dialog
 * boxes.  It uses the mini-toolbox given in xtb.c.
 */

#include <string>

#include "config.h"
#include "hard_devices.h"
#include "params.h"

using namespace std;

static int getline(char** tptr, char* lptr);

#define MAXCHBUF	1024

#ifdef SHADOW
#define gray_width 16
#define gray_height 16
static short gray_bits[] = {
	0x5555, 0xaaaa, 0x5555, 0xaaaa,
	0x5555, 0xaaaa, 0x5555, 0xaaaa,
	0x5555, 0xaaaa, 0x5555, 0xaaaa,
	0x5555, 0xaaaa, 0x5555, 0xaaaa};
#endif

static void make_msg_box(char const* text, char const* title, xtb_frame *frame);
static void del_msg_box(Window msg);


#define D_VPAD	2
#define D_HPAD	2
#define D_INT	4
#define D_BRDR	2
#define D_INP	35
#define D_DSP	10
#define D_FS	10

struct d_info {
	char *prog;			/* Program name              */
	xtb_data cookie;		/* Info used by do_harcopy   */
	Window choices;		/* Output device choices     */
	Window fod;			/* File or device flag       */
	Window fodspec;		/* File or device spec       */
	Window docu_p;		/* Document predicate        */
	Window dimspec;		/* Maximum dimension spec    */
	Window tf_family;		/* Title font family spec    */
	Window tf_size;		/* Title font size spec      */
	Window af_family;		/* Axis font family spec     */
	Window af_size;		/* Axis font size spec       */
};

#define	BACKSPACE	0010
#define DELETE		0177
#define CONTROL_U	0025
#define CONTROL_X	0030

/*ARGSUSED*/
static xtb_hret df_fun(
Window win,			/* Widget window   */
int ch,				/* Typed character */
char* text,			/* Copy of text    */
xtb_data val			/* User info       */
)
/*
 * This is the handler function for the text widget for
 * specifing the file or device name.  It supports simple
 * line editing operations.
 */
{
	if ((ch == BACKSPACE) || (ch == DELETE)) {
		if (!xtb_ti_dch(win)) XBell(disp, 0);
	} else if ((ch == CONTROL_U) || (ch == CONTROL_X)) {
		(void) xtb_ti_set(win, "", (xtb_data) 0);
	} else {
		/* Insert if printable - ascii dependent */
		if ((ch < ' ') || (ch >= DELETE) || !xtb_ti_ins(win, ch)) {
			XBell(disp, 0);
		}
	}
	return XTB_HANDLED;
}

/*ARGSUSED*/
static xtb_hret ok_fun(
Window win,			/* Button window     */
int bval,			/* Button value      */
xtb_data info) 			/* Local button info */
/*
 * This is the handler function for when the `Ok' button
 * is hit.  It sets the button,  does the hardcopy output,
 * and turns the button off.  It returns a status which
 * deactivates the dialog.
 */
{
	struct d_info *real_info = (struct d_info *) info;
	int val, dev_p, doc_p;
	char file_or_dev[MAXCHBUF], dim_spec[MAXCHBUF];
	const char *dev_spec;
	char tfam[MAXCHBUF], afam[MAXCHBUF];
	char tsizstr[MAXCHBUF], asizstr[MAXCHBUF];
	double centimeters, tsize, asize;
	xtb_hret rtn;

	(void) xtb_bt_set(win, 1, (xtb_data) 0, 0);
	val = xtb_br_get(real_info->choices);
	if ((val >= 0) && (val < hard_count)) {
		dev_p = xtb_br_get(real_info->fod);
		if ((dev_p == 0) || (dev_p == 1)) {
			xtb_ti_get(real_info->fodspec, file_or_dev, (xtb_data *) 0);
			doc_p = xtb_bt_get(real_info->docu_p, (xtb_data *) 0, (int *) 0);
			xtb_ti_get(real_info->dimspec, dim_spec, (xtb_data *) 0);
			if (sscanf(dim_spec, "%lf", &centimeters) == 1) {
				xtb_ti_get(real_info->tf_family, tfam, (xtb_data *) 0);
				xtb_ti_get(real_info->af_family, afam, (xtb_data *) 0);
				xtb_ti_get(real_info->tf_size, tsizstr, (xtb_data *) 0);
				if (sscanf(tsizstr, "%lf", &tsize) == 1) {
					xtb_ti_get(real_info->af_size, asizstr, (xtb_data *) 0);
					if (sscanf(asizstr, "%lf", &asize) == 1) {
						/* Got all the info */
						if (dev_p) {
							dev_spec = (char *) 0;
						} else {
							dev_spec = hard_devices[val].dev_spec;
						}
						do_hardcopy(real_info->prog, real_info->cookie,
						    hard_devices[val].dev_init, dev_spec,
						    file_or_dev, centimeters,
						    tfam, tsize, afam, asize, doc_p);
						rtn = XTB_STOP;
					} else {
						/* Bad axis size */
						do_error("Bad axis font size\n");
						rtn = XTB_HANDLED;
					}
				} else {
					/* Bad title size */
					do_error("Bad title font size\n");
					rtn = XTB_HANDLED;
				}
			} else {
				/* Bad max dimension */
				do_error("Bad maximum dimension\n");
				rtn = XTB_HANDLED;
			}
		} else {
			/* Bad device spec */
			do_error("Must specify `To File' or `To Device'\n");
			rtn = XTB_HANDLED;
		}
	} else {
		/* Bad value spec */
		do_error("Must specify an output device\n");
		rtn = XTB_HANDLED;
	}
	(void) xtb_bt_set(win, 0, (xtb_data) 0, 0);
	return rtn;
}

/*ARGSUSED*/
static xtb_hret can_fun(
Window win,			/* Button window     */
int val,			/* Button value      */
xtb_data info)			/* Local button info */
/*
 * This is the handler function for the cancel button.  It
 * turns itself on and off and then exits with a status 
 * which kills the dialog.
 */
{
	(void) xtb_bt_set(win, 1, (xtb_data) 0, 0);
	(void) xtb_bt_set(win, 0, (xtb_data) 0, 0);
	return XTB_STOP;
}


/*ARGSUSED*/
static xtb_hret docu_fun(
Window win,			/* Button window     */
int val,			/* Button value      */
xtb_data info)			/* Local button info */
/*
 * This is the handler function for the document button.  It
 * toggles it's state.
 */
{
	int state;

	state = xtb_bt_get(win, (xtb_data *) 0, (int *) 0);
	xtb_bt_set(win, (state == 0), (xtb_data) 0, 0);
	return XTB_HANDLED;
}

/*ARGSUSED*/
static xtb_hret dev_fun(
Window win,			/* Button row window */
int old,			/* Previous button   */
int current,			/* Current button    */
xtb_data info)			/* User data         */
/*
 * This routine swaps the device specific information
 * in the dialog based on what device is selected.  The
 * information passed is the information for the whole
 * dialog.
 */
{
	struct d_info *data = (struct d_info *) info;
	char text[MAXCHBUF];
	int fod_spot, inactive;

	fod_spot = xtb_br_get(data->fod);
	if ((old >= 0) && (old < hard_count)) {
		/* Save old info */
		xtb_ti_get(data->fodspec, text, (xtb_data *) 0);
		if (fod_spot == 1) {
			strncpy(hard_devices[old].dev_file, text, MFNAME-1);
		} else if (fod_spot == 0) {
			strncpy(hard_devices[old].dev_printer, text, MFNAME-1);
		}
		if (xtb_bt_get(data->docu_p, (xtb_data *) 0, &inactive)) {
			if (inactive) hard_devices[old].dev_docu = NONE;
			else hard_devices[old].dev_docu = YES;
		} else if (inactive) hard_devices[old].dev_docu = NONE;
		else hard_devices[old].dev_docu = NO;
		xtb_ti_get(data->dimspec, text, (xtb_data *) 0);
		if (sscanf(text, "%lf", &hard_devices[old].dev_max_dim) != 1) {
			do_error("Warning: can't read maximum dimension");
		}
		xtb_ti_get(data->tf_family, text, (xtb_data *) 0);
		strncpy(hard_devices[old].dev_title_font, text, MFNAME-1);
		xtb_ti_get(data->tf_size, text, (xtb_data *) 0);
		if (sscanf(text, "%lf", &hard_devices[old].dev_title_size) != 1) {
			do_error("Warning: can't read title font size");
		}
		xtb_ti_get(data->af_family, text, (xtb_data *) 0);
		strncpy(hard_devices[old].dev_axis_font, text, MFNAME-1);
		xtb_ti_get(data->af_size, text, (xtb_data *) 0);
		if (sscanf(text, "%lf", &hard_devices[old].dev_axis_size) != 1) {
			do_error("Warning: can't read axis font size");
		}
	}
	/* Insert new info */
	if ((current >= 0) && (current < hard_count)) {
		if (fod_spot == 1) {
			xtb_ti_set(data->fodspec, hard_devices[current].dev_file, (xtb_data) 0);
		} else if (fod_spot == 0) {
			xtb_ti_set(data->fodspec, hard_devices[current].dev_printer,
			    (xtb_data) 0);
		} else {
			xtb_ti_set(data->fodspec, "", (xtb_data) 0);
		}
		switch (hard_devices[current].dev_docu) {
		case NONE:
			(void) xtb_bt_set(data->docu_p, 0, (xtb_data) 0, 1);
			break;
		case NO:
			(void) xtb_bt_set(data->docu_p, 0, (xtb_data) 0, 0);
			break;
		case YES:
			(void) xtb_bt_set(data->docu_p, 1, (xtb_data) 0, 0);
			break;
		}
		(void) sprintf(text, "%lg", hard_devices[current].dev_max_dim);
		xtb_ti_set(data->dimspec, text, (xtb_data) 0);
		xtb_ti_set(data->tf_family, hard_devices[current].dev_title_font,
		    (xtb_data) 0);
		(void) sprintf(text, "%lg", hard_devices[current].dev_title_size);
		xtb_ti_set(data->tf_size, text, (xtb_data) 0);
		xtb_ti_set(data->af_family, hard_devices[current].dev_axis_font,
		    (xtb_data) 0);
		(void) sprintf(text, "%lg", hard_devices[current].dev_axis_size);
		xtb_ti_set(data->af_size, text, (xtb_data) 0);
	}
	return XTB_HANDLED;
}

/*ARGSUSED*/
static xtb_hret fd_fun(
Window win,			/* Button row window */
int old,			/* Previous button   */
int current,			/* Current button    */
xtb_data info)			/* User data         */
/*
 * This routine swaps the default file or device names
 * based on the state of the file or device buttons.
 * The information passed is the information for the whole
 * dialog.
 */
{
	struct d_info *data = (struct d_info *) info;
	char text[MAXCHBUF];
	int which_one;

	which_one = xtb_br_get(data->choices);
	if ((which_one >= 0) && (which_one < hard_count)) {
		if (old == 0) {
			/* Save into device space */
			xtb_ti_get(data->fodspec, text, (xtb_data *) 0);
			strncpy(hard_devices[which_one].dev_printer, text, MFNAME-1);
		} else if (old == 1) {
			/* Save into file space */
			xtb_ti_get(data->fodspec, text, (xtb_data *) 0);
			which_one = xtb_br_get(data->choices);
			strncpy(hard_devices[which_one].dev_file, text, MFNAME-1);
		}
		if (current == 0) {
			/* Restore into device */
			xtb_ti_set(data->fodspec, hard_devices[which_one].dev_printer,
			    (xtb_data *) 0);
		} else if (current == 1) {
			xtb_ti_set(data->fodspec, hard_devices[which_one].dev_file,
			    (xtb_data *) 0);
		}
	}
	return XTB_HANDLED;
}


/* Indices for frames made in make_dialog */
enum d_frames_defn {
	TITLE_F, ODEVLBL_F, ODEVROW_F, DISPLBL_F, DISPROW_F, FDLBL_F,
	FDINP_F, OPTLBL_F, DOCU_F, MDIMLBL_F, MDIMI_F, TFLBL_F, TFFAMLBL_F, TFFAM_F,
	TFSIZLBL_F, TFSIZ_F, AFLBL_F, AFFAMLBL_F, AFFAM_F, AFSIZLBL_F, AFSIZ_F,
	OK_F, CAN_F, BAR_F, LAST_F
} 
d_frames;

#define AF(ix)	af[(int) (ix)]

#define BAR_SLACK	10

static void make_dialog(
Window win,			/* Parent window          */
Window spawned,			/* Spawned from window    */
char *prog,			/* Program name           */
xtb_data cookie,		/* Info for do_hardcopy  */
xtb_frame *okbtn,		/* Frame for OK button    */
xtb_frame *frame)		/* Returned window/size   */
/*
 * This routine constructs a new dialog for asking the user about
 * hardcopy devices.  The dialog and its size is returned in 
 * `frame'.  The window of the `ok' button is returned in `btnwin'.  
 * This can be used to reset some of the button state to reuse the dialog.
 */
{
	Window overall;
	xtb_frame AF(LAST_F);
	xtb_fmt *def, *cntrl, *mindim, *tfarea, *afarea;
	Cursor diag_cursor;
	XColor fg_color, bg_color;
	XSizeHints hints;
	unsigned long wamask;
	XSetWindowAttributes wattr;
	struct d_info *info;
	int i, found, max_width, which_one, fodi;
	const char **names;
	static const char *fodstrs[] = { 
		"To Device", "To File" 	};
	XFontStruct *bigFont = PM_FONT("TitleFont");
	XFontStruct *medFont = PM_FONT("LabelFont");
	char *Odevice = pm_str("Device");
	char *Odisp = pm_str("Disposition");
	char *OfileDev = pm_str("FileOrDev");

	wamask = CWBackPixel | CWBorderPixel | CWOverrideRedirect |
	    CWSaveUnder| CWColormap;
	wattr.background_pixel = PM_PIXEL("Background");
	wattr.border_pixel = PM_PIXEL("Border");
	wattr.override_redirect = True;
	wattr.save_under = True;
	wattr.colormap = cmap;
	overall = XCreateWindow(disp, win, 0, 0, 1, 1, D_BRDR,
	    depth, InputOutput, visual,
	    wamask, &wattr);
	frame->win = overall;
	frame->width = frame->height = frame->x_loc = frame->y_loc = 0;
	XStoreName(disp, overall, "Hardcopy Dialog");
	XSetTransientForHint(disp, spawned, overall);
	info = (struct d_info *) calloc(1,sizeof(struct d_info));
	info->prog = prog;
	info->cookie = cookie;

	/* Make all frames */
	xtb_to_new(overall, "Hardcopy Options", bigFont, &AF(TITLE_F));
	xtb_to_new(overall, "Output device:", medFont, &AF(ODEVLBL_F));
	found = -1;
	names = (const char **) calloc(1,(unsigned) (sizeof(char *) * hard_count));
	// names = new char*[hard_count];
	for (i = 0;  i < hard_count;  i++) {
		names[i] = hard_devices[i].dev_name;
		if (strcmp(Odevice, names[i]) == 0) {
			found = i;
		}
	}
	xtb_br_new(overall, hard_count, names, found,
	    dev_fun, (xtb_data) info, &AF(ODEVROW_F));
	info->choices = AF(ODEVROW_F).win;
	xtb_to_new(overall, "Disposition:", medFont, &AF(DISPLBL_F));
	found = -1;
	for (i = 0;  i < 2;  i++) {
		if (strcmp(Odisp, fodstrs[i]) == 0) {
			found = i;
		}
	}
	xtb_br_new(overall, 2, fodstrs, found,
	    fd_fun, (xtb_data) info, &AF(DISPROW_F));
	info->fod = AF(DISPROW_F).win;
	xtb_to_new(overall, "File or Device Name:", medFont, &AF(FDLBL_F));
	xtb_ti_new(overall, "", D_INP, df_fun, (xtb_data) 0, &AF(FDINP_F));
	if (OfileDev && strlen(OfileDev)) {
		which_one = xtb_br_get(info->choices);
		if ((which_one >= 0) && (which_one < hard_count)) {
			fodi = xtb_br_get(info->fod);
			if (fodi == 0) {
				strncpy(hard_devices[which_one].dev_printer, OfileDev, MFNAME-1);
			} else if (fodi == 1) {
				strncpy(hard_devices[which_one].dev_file, OfileDev, MFNAME-1);
			}
		}
	}
	info->fodspec = AF(FDINP_F).win;
	xtb_to_new(overall, "Optional Parameters", bigFont, &AF(OPTLBL_F));
	xtb_bt_new(overall, "Include in Document", docu_fun, (xtb_data) 0, &AF(DOCU_F));
	info->docu_p = AF(DOCU_F).win;
	xtb_to_new(overall, "Maximum Dimension (cm):", medFont, &AF(MDIMLBL_F));
	xtb_ti_new(overall, "", D_DSP, df_fun, (xtb_data) 0, &AF(MDIMI_F));
	info->dimspec = AF(MDIMI_F).win;
	xtb_to_new(overall, "Title Font", medFont, &AF(TFLBL_F));
	xtb_to_new(overall, "Family:", medFont, &AF(TFFAMLBL_F));
	xtb_ti_new(overall, "", MFNAME, df_fun, (xtb_data) 0, &AF(TFFAM_F));
	info->tf_family = AF(TFFAM_F).win;
	xtb_to_new(overall, "Size (pnts):", medFont, &AF(TFSIZLBL_F));
	xtb_ti_new(overall, "", D_FS, df_fun, (xtb_data) 0, &AF(TFSIZ_F));
	info->tf_size = AF(TFSIZ_F).win;
	xtb_to_new(overall, "Axis Font", medFont, &AF(AFLBL_F));
	xtb_to_new(overall, "Family:", medFont, &AF(AFFAMLBL_F));
	xtb_ti_new(overall, "", MFNAME, df_fun, (xtb_data) 0, &AF(AFFAM_F));
	info->af_family = AF(AFFAM_F).win;
	xtb_to_new(overall, "Size (pnts):", medFont, &AF(AFSIZLBL_F));
	xtb_ti_new(overall, "", D_FS, df_fun, (xtb_data) 0, &AF(AFSIZ_F));
	info->af_size = AF(AFSIZ_F).win;
	xtb_bt_new(overall, "  Ok  ", ok_fun, (xtb_data) info, &AF(OK_F));
	xtb_bt_new(overall, "Cancel", can_fun, (xtb_data) 0, &AF(CAN_F));
	/* Dividing bar */
	max_width = 0;
	for (i = 0;  i < ((int) BAR_F);  i++) {
		if ((int)AF(i).width > max_width) max_width = AF(i).width;
	}
	xtb_bk_new(overall, max_width - BAR_SLACK, 1, &AF(BAR_F));

	/* Set device specific info */
	(void) dev_fun(info->choices, -1,xtb_br_get(info->choices),(xtb_data) info);
	(void) fd_fun(info->fod, -1, xtb_br_get(info->fod), (xtb_data) info);

	/* 
     * Now place elements - could make this one expression but pcc
     * is too wimpy.
     */
	cntrl = xtb_vert(XTB_LEFT, D_VPAD, D_INT,
	    xtb_hort(XTB_CENTER, D_HPAD, D_INT,
	    xtb_w(&AF(ODEVLBL_F)),
	    xtb_w(&AF(ODEVROW_F)), (xtb_fmt*)0),
	    xtb_hort(XTB_CENTER, D_HPAD, D_INT,
	    xtb_w(&AF(DISPLBL_F)),
	    xtb_w(&AF(DISPROW_F)), (xtb_fmt*)0),
	    xtb_hort(XTB_CENTER, D_HPAD, D_INT,
	    xtb_w(&AF(FDLBL_F)),
	    xtb_w(&AF(FDINP_F)), (xtb_fmt*)0),
	    (xtb_fmt*)0);

	mindim = xtb_vert(XTB_LEFT, D_VPAD, D_INT,
	    xtb_hort(XTB_CENTER, D_HPAD, D_INT,
	    xtb_w(&AF(MDIMLBL_F)),
	    xtb_w(&AF(MDIMI_F)),
	    (xtb_fmt*)0),
	    (xtb_fmt*)0);

	tfarea = xtb_vert(XTB_LEFT, D_VPAD, D_INT,
	    xtb_hort(XTB_CENTER, D_HPAD, D_INT,
	    xtb_w(&AF(TFFAMLBL_F)),
	    xtb_w(&AF(TFFAM_F)),
	    xtb_w(&AF(TFSIZLBL_F)),
	    xtb_w(&AF(TFSIZ_F)), (xtb_fmt*)0),
	    (xtb_fmt*)0);

	afarea = xtb_vert(XTB_LEFT, D_VPAD, D_INT,
	    xtb_hort(XTB_CENTER, D_HPAD, D_INT,
	    xtb_w(&AF(AFFAMLBL_F)),
	    xtb_w(&AF(AFFAM_F)),
	    xtb_w(&AF(AFSIZLBL_F)),
	    xtb_w(&AF(AFSIZ_F)), (xtb_fmt*)0),
	    (xtb_fmt*)0);

	def = xtb_fmt_do
	    (xtb_vert(XTB_CENTER, D_VPAD, D_INT,
	    xtb_w(&AF(TITLE_F)),
	    cntrl,
	    xtb_w(&AF(BAR_F)),
	    xtb_w(&AF(OPTLBL_F)),
	    mindim,
	    xtb_w(&AF(DOCU_F)),
	    xtb_w(&AF(TFLBL_F)),
	    tfarea,
	    xtb_w(&AF(AFLBL_F)),
	    afarea,
	    xtb_hort(XTB_CENTER, D_HPAD, D_INT,
	    xtb_w(&AF(OK_F)), xtb_w(&AF(CAN_F)), (xtb_fmt*)0),
	    (xtb_fmt*)0),
	    &frame->width, &frame->height);
	xtb_mv_frames(LAST_F, af);
	xtb_fmt_free(def);

	/* Make window large enough to contain the info */
	XResizeWindow(disp, overall, frame->width, frame->height);
	hints.flags = PSize;
	hints.width = frame->width;
	hints.height = frame->height;
	XSetNormalHints(disp, overall, &hints);
	diag_cursor = XCreateFontCursor(disp, XC_dotbox);
	fg_color.pixel = PM_PIXEL("Foreground");
	XQueryColor(disp, cmap, &fg_color);
	bg_color.pixel = PM_PIXEL("Background");
	XQueryColor(disp, cmap, &bg_color);
	XRecolorCursor(disp, diag_cursor, &fg_color, &bg_color);
	XDefineCursor(disp, overall, diag_cursor);
	frame->width += (2 * D_BRDR);
	frame->height += (2 * D_BRDR);
	*okbtn = AF(OK_F);
}



#ifdef SHADOW
Window make_shadow(int w, int h)
/*
 * Makes a shadow window for a pop-up window of the specified size.
 * Needs hint work as well.  Try no background window.
 */
{
	Window shadow;
	Bitmap gray_bm;
	Pixmap gray_pm;

	gray_bm = XStoreBitmap(gray_width, gray_height, gray_bits);
	gray_pm = XMakePixmap(gray_bm, PM_PIXEL("Foreground"), PM_PIXEL("Background"));
	shadow = XCreateWindow(RootWindow, 0, 0, w, h, 0, BlackPixmap, gray_pm);
	XFreePixmap(gray_pm);
	XFreeBitmap(gray_bm);
	return shadow;
}
#endif



#define SH_W	5
#define SH_H	5

void
ho_dialog(
Window parent,			/* Parent window              */
char *prog,			/* Program name               */
xtb_data cookie)		/* Info passed to do_hardcopy */
/*
 * Asks the user about hardcopy devices.  A table of hardcopy
 * device names and function pointers to their initialization
 * functions is assumed to be in the global `hard_devices' and
 * `hard_count'.  Returns a non-zero status if the dialog was
 * sucessfully posted.  Calls do_hardcopy in xgraph to actually
 * output information.
 */
{
#ifdef SHADOW
	static Window shadow;
#endif
	static Window dummy;
	static xtb_frame overall = { 
		(Window) 0, 0, 0, 0, 0 	};
	static xtb_frame okbtn;
	XEvent evt;
	XWindowAttributes winInfo;
	XSizeHints hints;
	struct d_info *info;

	if (!overall.win) {
		make_dialog(RootWindow(disp, screen), parent, prog, cookie,
		    &okbtn, &overall);
#ifdef SHADOW
		shadow = make_shadow(d_w, d_h);
#endif
	} else {
		/* Change the button information */
		(void) xtb_bt_get(okbtn.win, (xtb_data *) &info, (int *) 0);
		info->prog = prog;
		info->cookie = cookie;
	}
	XGetWindowAttributes(disp, parent, &winInfo);
	XTranslateCoordinates(disp, parent, RootWindow(disp, screen),
	    0, 0, &winInfo.x, &winInfo.y, &dummy);
	XMoveWindow(disp, overall.win,
	    (int) (winInfo.x + winInfo.width/2 - overall.width/2),
	    (int) (winInfo.y + winInfo.height/2 - overall.height/2));
	hints.flags = PPosition;
	hints.x = winInfo.x + winInfo.width/2 - overall.width/2;
	hints.y = winInfo.y + winInfo.height/2 - overall.height/2;
	XSetNormalHints(disp, overall.win, &hints);
#ifdef SHADOW
	XMoveWindow(disp, shadow, winInfo.x + winInfo.width/2 - d_w/2 + SH_W,
	    winInfo.y + winInfo.height/2 - d_h/2 + SH_H);
	hints.flags = PPosition;
	hints.x = winInfo.x + winInfo.width/2 - d_w/2 + SH_W;
	hints.y = winInfo.y + winInfo.height/2 - d_h/2 + SH_H;
	XSetNormalHints(disp, shadow, &hints);
	XRaiseWindow(disp, shadow);
	XMapWindow(disp, shadow);
#endif
	XRaiseWindow(disp, overall.win);
	XMapWindow(disp, overall.win);
	do {
		XNextEvent(disp, &evt);
	} while (xtb_dispatch(&evt) != XTB_STOP);
	XUnmapWindow(disp, overall.win);
#ifdef SHADOW
	XUnmapWindow(disp, shadow);
#endif
}



/*ARGSUSED*/
static xtb_hret
err_func(
Window win,			/* Button window     */
int bval,			/* Button value      */
xtb_data info)			/* Local button info */
/*
 * Handler function for button in error box.  Simply stops dialog.
 */
{
	(void) xtb_bt_set(win, 1, (xtb_data) 0, 0);
	(void) xtb_bt_set(win, 0, (xtb_data) 0, 0);
	return XTB_STOP;
}


struct err_info {
	Window title;
	Window contbtn;
	int num_lines;
	int alloc_lines;
	Window *lines;
};

#define E_LINES	2
#define E_VPAD	3
#define E_HPAD	3
#define E_INTER	1
#define ML	256

static void
make_msg_box(
char const* txt,			/* Error text    */
char const* ti,			/* Title text    */
xtb_frame *frame)		/* Returned frame */
// Makes an message box with a title.
{
	XSizeHints hints;
	struct err_info *new_info;
	xtb_frame tf, cf, lf;
	char *lineptr, line[ML];
	int y, i;
	unsigned long wamask;
	char* text = strdup(txt);
	char* title = strdup(ti);

	XSetWindowAttributes wattr;
	XFontStruct *bigFont = PM_FONT("TitleFont");
	XFontStruct *medFont = PM_FONT("LabelFont");

	wamask = CWBackPixel | CWBorderPixel | CWOverrideRedirect |
	    CWSaveUnder | CWColormap;
	wattr.background_pixel = PM_PIXEL("Background");
	wattr.border_pixel = PM_PIXEL("Border");
	wattr.override_redirect = True;
	wattr.save_under = True;
	wattr.colormap = cmap;
	frame->win = XCreateWindow(disp, RootWindow(disp, screen),
	    0, 0, 1, 1, D_BRDR,
	    depth, InputOutput, visual,
	    wamask, &wattr);
	frame->x_loc = frame->y_loc = frame->width = frame->height = 0;
	XStoreName(disp, frame->win, "Error Dialog");
	XSetTransientForHint(disp, RootWindow(disp, screen), frame->win);
	new_info = (struct err_info *) calloc(1,(unsigned) sizeof(struct err_info));
	xtb_to_new(frame->win, title, bigFont, &tf);
	new_info->title = tf.win;
	if (tf.width > frame->width) frame->width = tf.width;

	xtb_bt_new(frame->win, "Dismiss", err_func, (xtb_data) 0, &cf);
	new_info->contbtn = cf.win;
	if (cf.width > frame->width) frame->width = cf.width;

	new_info->alloc_lines = E_LINES;
	new_info->num_lines = 0;
	new_info->lines = (Window *) calloc(1,(unsigned) (sizeof(Window) * E_LINES));

	lineptr = text;
	while (getline(&lineptr, line)) {
		if (new_info->num_lines >= new_info->alloc_lines) {
			new_info->alloc_lines *= 2;
			new_info->lines = (Window *) realloc((char *) new_info->lines,
			    (unsigned)
			    (new_info->alloc_lines *
			    sizeof(Window)));
		}
		xtb_to_new(frame->win, line, medFont, &lf);
		new_info->lines[new_info->num_lines] = lf.win;
		new_info->num_lines += 1;
		if (lf.width > frame->width) frame->width = lf.width;
	}


	/* Placement */
	frame->width += (2 * E_HPAD);
	y = E_VPAD;
	/* Title */
	XMoveWindow(disp, new_info->title, (int) (frame->width/2 - tf.width/2), y);
	y += (tf.height + E_INTER);
	/* All lines */
	for (i = 0;  i < new_info->num_lines;  i++) {
		XMoveWindow(disp, new_info->lines[i], E_HPAD, y);
		y += (lf.height + E_INTER);
	}
	/* Button */
	XMoveWindow(disp, new_info->contbtn, (int) (frame->width/2 - cf.width/2), y);
	y += (cf.height + E_INTER);

	/* Make dialog the right size */
	y += (E_VPAD - E_INTER);
	XResizeWindow(disp, frame->win, frame->width, (unsigned int) y);
	hints.flags = PSize;
	hints.width = frame->width;
	hints.height = (unsigned int) y;
	XSetNormalHints(disp, frame->win, &hints);
	frame->width += (2 * D_BRDR);
	frame->height = y + (2 * D_BRDR);
	/* xtb_register(frame->win, (xtb_hret (*)()) 0, (xtb_data) new_info); */
	xtb_register(frame->win,  0, (xtb_data) new_info);
}


void
msg_box(char const* title, char const* text) {
// This posts a dialog that contains lines of text and a continue
// button.  The text may be multiple lines.  The dialog is remade
// each time.
#ifdef SHADOW
	Window shadow;
#endif
	XWindowAttributes info;
	XEvent evt;
	XSizeHints hints;
	xtb_frame text_frame;

	// text = formatText (text);

	make_msg_box(text, title, &text_frame);
#ifdef SHADOW
	shadow = make_shadow(w, h);
#endif
	XGetWindowAttributes(disp, RootWindow(disp, screen), &info);
	XMoveWindow(disp, text_frame.win, (int) (info.width/2 - text_frame.width/2),
	     (int) (info.height/2 - text_frame.height/2));
	hints.flags = PPosition;
	hints.x = info.width/2 - text_frame.width/2;
	hints.y = info.height/2 - text_frame.height/2;
	XSetNormalHints(disp, text_frame.win, &hints);
#ifdef SHADOW
	XMoveWindow(disp, shadow, info.width/2 - w/2 + SH_W,
	    info.height/2 - h/2 + SH_H);
	hints.flags = PPosition;
	hints.x = info.width/2 - w/2 + SH_W;
	hints.y = info.height/2 - h/2 + SH_H;
	XSetNormalHints(disp, text_frame.win, &hints);
	XRaiseWindow(disp, shadow);
	XMapWindow(disp, shadow);
#endif
	XRaiseWindow(disp, text_frame.win);
	XMapWindow(disp, text_frame.win);
	do {
		XNextEvent(disp, &evt);
		// a shortcut for deleting the message box: press x
		if (evt.type == KeyPress) {
			char keys[2];
			XLookupString(&evt.xkey, keys, 2, nullptr, nullptr);
			if (keys[0] == 'x')
				break;
		}
	} while (xtb_dispatch(&evt) != XTB_STOP);
#ifdef SHADOW
	XDestroyWindow(disp, shadow);
#endif
	del_msg_box(text_frame.win);
}

void
do_error(string const& err_text) {
//do_error(char const* err_text) {
	msg_box("Error:", err_text.c_str());
}



static int
getline(char** tptr, char* lptr)
/*
 * Returns next line from tptr.  Assumes `lptr' is large enough to
 * accept the line.
 */
{
	char *start;

	start = *tptr;
	while (*tptr && **tptr && (**tptr != '\n')) {
		(*tptr)++;
	}
	if (*tptr > start) {
		(void) strncpy(lptr, start, (*tptr - start));
		lptr[*tptr-start] = '\0';
		if (**tptr == '\n') (*tptr)++;
		return 1;
	} else {
		return 0;
	}
}



static void
del_msg_box(Window msg)
/*
 * Deletes all components of an msg box
 */
{
	struct err_info *info;
	char *dummy;
	int i;

	if (xtb_unregister(msg, (xtb_data *) &info)) {
		xtb_to_del(info->title);
		xtb_bt_del(info->contbtn, (xtb_data *) &dummy);
		for (i = 0;  i < info->num_lines;  i++) {
			xtb_to_del(info->lines[i]);
		}
		free((char*)info->lines);
		free((char*)info);
		XDestroyWindow(disp, msg);
	}
}
