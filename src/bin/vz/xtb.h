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
 * xtb - a mini-toolbox for X11
 *
 * David Harrison
 * University of California, Berkeley
 * 1988
 */

#ifndef _XTB_
#define _XTB_

#include "copyright.h"

typedef void* xtb_data;

/* Handler function return codes */
typedef enum xtb_hret_defn { XTB_NOTDEF, XTB_HANDLED, XTB_STOP } xtb_hret;

/* Basic return value */
typedef struct xtb_frame_defn {
    Window win;
    int x_loc, y_loc;
    unsigned int width, height;
} xtb_frame;

void xtb_init(Display *disp, int scrn,
			 unsigned long foreground,
			 unsigned long background,
			 XFontStruct *font);
   /* Initializes mini-toolbox */

/*
 * Basic event handling
 */

void xtb_register(Window win,
			     xtb_hret (*func)(XEvent *evt, xtb_data info),
			     xtb_data info);
   /* Registers call-back function */

xtb_data xtb_lookup(Window win);
   /* Returns data associated with window */

xtb_hret xtb_dispatch(XEvent *evt);
   /* Dispatches events for mini-toolbox */

int xtb_unregister(Window win, xtb_data *info);
   /* Unregisters a call-back function */

/*
 * Command button frame
 */

void xtb_bt_new(Window win, char const* text,
		xtb_hret (*func)(Window win, int state, xtb_data val),
		xtb_data val, xtb_frame *frame);
   /* Creates new button  */

int xtb_bt_get(Window win, xtb_data *stuff, int *na);
   /* Returns state of button */
int xtb_bt_set(Window win, int val, xtb_data stuff, int na);
   /* Sets state of button */
void xtb_bt_del(Window win, xtb_data *info);
   /* Deletes a button */

/*
 * Button row frame - built on top of buttons
 */

void xtb_br_new(Window win, int cnt, const char *lbls[], int init,
			xtb_hret (*func)(Window win, int prev, int thisn, xtb_data val),
			xtb_data val, xtb_frame *frame);
   /* Creates a new button row frame */

int xtb_br_get(Window win);
   /* Returns currently selected button */
void xtb_br_del(Window win);
   /* Deletes a button row */

/*
 * Text output (label) frames
 */

void xtb_to_new(Window win, char const* text,
			   XFontStruct *ft, xtb_frame *frame);
   /* Create new text output frame */
void xtb_to_del(Window win);

/*
 * Text input (editable text) frames
 */

#define MAXCHBUF	1024

void xtb_ti_new(Window win, const char *text, int maxchar,
		xtb_hret (*func)(Window win, int ch, char *textcopy, xtb_data val),
		xtb_data val, xtb_frame *frame);
   /* Creates a new text input frame */		   

void xtb_ti_get(Window win, char text[MAXCHBUF], xtb_data *val);
   /* Returns state of text input frame */
int xtb_ti_set(Window win, char const* text, xtb_data val);
   /* Sets the state of text input frame */
int xtb_ti_ins(Window win, int ch);
   /* Inserts character onto end of text input frame */
int xtb_ti_dch(Window win);
   /* Deletes character from end of text input frame */
void xtb_ti_del(Window win, xtb_data *info);
   /* Deletes an text input frame */

/*
 * Block frame
 */

void xtb_bk_new(Window win, unsigned width, unsigned height, xtb_frame *frame);
   /* Makes a new block frame */
void xtb_bk_del(Window win);
   /* Deletes a block frame */


/*
 * Formatting support
 */

#define MAX_BRANCH	50

typedef enum xtb_fmt_types_defn { W_TYPE, A_TYPE } xtb_fmt_types;
typedef enum xtb_fmt_dir_defn { HORIZONTAL, VERTICAL } xtb_fmt_dir;
typedef enum xtb_just_defn {
    XTB_CENTER=0, XTB_LEFT, XTB_RIGHT, XTB_TOP, XTB_BOTTOM
} xtb_just;

typedef struct xtb_fmt_widget_defn {
    xtb_fmt_types type;		/* W_TYPE */
    xtb_frame *w;
} xtb_fmt_widget;

typedef struct xtb_fmt_align_defn {
    xtb_fmt_types type;		/* A_TYPE */
    xtb_fmt_dir dir;		/* HORIZONTAL or VERTICAL */
    int padding;		/* Outside padding        */
    int interspace;		/* Internal padding       */
    xtb_just just;		/* Justification          */
    int ni;			/* Number of items */
    union xtb_fmt_defn *items[MAX_BRANCH]; /* Branches themselves */
} xtb_fmt_align;

typedef union xtb_fmt_defn {
    xtb_fmt_types type;		/* W_TYPE or A_TYPE */
    xtb_fmt_widget wid;
    xtb_fmt_align align;
} xtb_fmt;

#define NE	0

xtb_fmt* xtb_w(xtb_frame *w);
   /* Returns formatting structure for frame */
xtb_fmt* xtb_hort(xtb_just just, int padding, int interspace, ...);
   /* Varargs routine for horizontal formatting */
xtb_fmt* xtb_vert(xtb_just just, int padding, int interspace, ...);
   /* Varargs routine for vertical formatting */
xtb_fmt* xtb_fmt_do(xtb_fmt *def, unsigned *w, unsigned *h);
   /* Carries out formatting */
void xtb_mv_frames(int nf, xtb_frame frames[]);
   /* Actually moves widgets */
void xtb_fmt_free(xtb_fmt *def);
   /* Frees resources claimed by xtb_w, xtb_hort, and xtb_vert */


#endif /* _XTB_ */
