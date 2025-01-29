//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef XGOUT_H
#define XGOUT_H
/*
 * Output Device Information
 *
 * This file contains definitions for output device interfaces
 * to the graphing program xgraph.
 */

typedef struct xg_out xgOut;

#include "plotcurve.h"


/* Passed device option flags */
#define D_DOCU		0x01

/* Returned device capability flags */
#define D_COLOR		0x01

/* Text justifications */
#define T_CENTER	0
#define T_LEFT		1
#define T_UPPERLEFT	2
#define T_TOP		3
#define T_UPPERRIGHT	4
#define T_RIGHT		5
#define T_LOWERRIGHT	6
#define T_BOTTOM	7
#define T_LOWERLEFT	8

/* Text styles */
#define T_AXIS		0
#define T_TITLE		1

/* Line Styles */
#define L_AXIS		0
#define L_ZERO		1
#define L_VAR		2

/* Marker Styles */
#define P_PIXEL		0
#define P_DOT		1
#define P_MARK		2
#define P_TAG		3

/* Output device information returned by initialization routine */

struct xg_out {
    int dev_flags;		/* Device characteristic flags           */
    int area_w, area_h;		/* Width and height in pixels            */
    int bdr_pad;		/* Padding from border                   */
    int axis_pad;		/* Extra space around axis labels        */
    int tick_len;		/* Length of tick mark on axis           */
    int legend_pad;		/* Top of legend text to legend line     */
    int axis_width;		/* Width of big character of axis font   */
    int axis_height;		/* Height of big character of axis font  */
    int title_width;		/* Width of big character of title font  */
    int title_height;		/* Height of big character of title font */
    int max_segs;		/* Maximum number of segments in group   */

    void (*xg_text)(char*,int,int,char const*,int,int);	/* Draws text at a location */
    void (*xg_seg)(char*,int,XSegment*,int,int,int,int);		/* Draws a series of segments            */
    void (*xg_dot)(char*,int,int,int,int,int);		/* Draws a dot or marker at a location   */
    void (*xg_end)(char*);		/* Stops the drawing sequence            */

    char *user_state;		/* User supplied data                    */
};

void init_X(char* user_state);

#define ERRBUFSIZE	2048

#endif /* XGOUT_H */
