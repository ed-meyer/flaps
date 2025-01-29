
/* $Id: matview_glob.c,v 1.52 1999/03/11 22:09:19 kohl Exp $ */


/*
 *                      MatView Version 1.0:
 *                 Simple Scalable Matrix Viewer
 *          Oak Ridge National Laboratory, Oak Ridge TN.
 *                   Author:  James Arthur Kohl
 *                  (C) 1998 All Rights Reserved
 *
 *                              NOTICE
 *
 * Permission to use, copy, modify, and distribute this software and
 * its documentation for any purpose and without fee is hereby granted
 * provided that the above copyright notice appear in all copies and
 * that both the copyright notice and this permission notice appear
 * in supporting documentation.
 *
 * Neither the Institution, Oak Ridge National Laboratory, nor the
 * Authors make any representations about the suitability of this
 * software for any purpose.  This software is provided ``as is''
 * without express or implied warranty.
 *
 * MatView was funded by the U.S. Department of Energy.
 */


/* MatView Globals Declarations */


#include "matview.h"


/* MatView Globals */

REGION MAX_REGION;

int PLANE_CHOICE;


/* TCL / TK Globals */

Tcl_Interp		*interp = 0;

Tk_Window		Top;

Display			*Disp;

char			**Argv;
int 			Argc;

char			*MENUVAR_AUTO_COLOR;
char			*MENUVAR_GREYSCALE;
char			*MENUVAR_FILTER_VALUES;
char			*MENUVAR_MIN_MAX;
char			*MENUVAR_CUMULATIVE;
char			*MENUVAR_ABSVALS;
char			*MENUVAR_PATTERN;
char			*MENUVAR_TRANSPOSE;
char			*MENUVAR_NOLOADVALUES;

int				BORDER_SPACE;
int				FRAME_BORDER;
int				ROW_HEIGHT;
int				COL_WIDTH;

int				PTR_VALID;

int				PTR_X;
int				PTR_Y;

/* MatView Globals */

REGION_STACK	REGIONS_SAVE;
REGION_STACK	REGIONS;

VIEW			PRINT_VIEW;
VIEW			THUMB_VIEW;
VIEW			MAIN_VIEW;
VIEW			SEP_VIEW;

double			TMP_VALUE;

void			*TMP_PTR;

char			TMP_COORD1[1024];
char			TMP_COORD2[1024];
char			TMP_COORD3[1024];
char			TMP_COORD4[1024];

char			TMP_CMD[2048];

char			*MATVIEW_MATFILE;

char			*MATVIEW_DIR;

char			*MV;
char			*MV_M;
char			*MV_G;

char			*PROGRESS_HEIGHT;
char			*PROGRESS_WIDTH;

MATRIX			MAT;

double			*SAVE_FIELD;
int				SAVE_FIELD_NELEMS;

double			COLOR_MIN_VALUE;
double			COLOR_MAX_VALUE;

int				AUTO_COLOR_MAP;

int				GREYSCALE;

int				FILTER_VALUES;

int				COMPLEX_VIEW;

int				PROGRESS_TICKS;

char			*ZERO_COLOR;

double			COLOR_DIST;

int				*COLORIDS;
int				NCOLORIDS;

int				GRID_STATUS;
int				GRID_ROWS;
int				GRID_ROFF;
int				GRID_COLS;
int				GRID_COFF;

int				VIEW_ABS_VALS;

int				VIEW_PATTERN;

int				TRANSPOSE_MATRIX;

int				SHOW_CUMULATIVE;

int				SHOW_MIN_MAX;

double			CUMULATIVE_MIN;
double			CUMULATIVE_MAX;
int				CUMULATIVE_SET;

int				NOLOAD_VALUES;

int				REDUCTION_FUNCTION;

int				NAVIGATE_SIZE;

int				THUMB_BORDER;

int				USE_SEP_VIEW;

int				REGION_ID;

int				DATA_MODE;

int				VIS_MODE;

int				NO_SHELL;

int 			vflag;

#ifdef FILE_BUFFERING

int				FILE_BUF_SIZE;

#endif

