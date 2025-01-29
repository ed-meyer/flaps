/* $Id: matview_glob.h,v 1.53 1999/03/11 22:09:20 kohl Exp $ */


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


/* Externals */

extern	int 	errno;


/* MatView Globals */

extern	REGION MAX_REGION;

extern	int PLANE_CHOICE;


/* TCL / TK Globals */

extern	Tcl_Interp		*interp;

extern	Tk_Window		Top;

extern	Display			*Disp;

extern	char			**Argv;
extern	int 			Argc;

extern	char			*MENUVAR_AUTO_COLOR;
extern	char			*MENUVAR_GREYSCALE;
extern	char			*MENUVAR_FILTER_VALUES;
extern	char			*MENUVAR_MIN_MAX;
extern	char			*MENUVAR_CUMULATIVE;
extern	char			*MENUVAR_ABSVALS;
extern	char			*MENUVAR_PATTERN;
extern	char			*MENUVAR_TRANSPOSE;
extern	char			*MENUVAR_NOLOADVALUES;

extern	int				BORDER_SPACE;
extern	int				FRAME_BORDER;
extern	int				ROW_HEIGHT;
extern	int				COL_WIDTH;

extern	int				PTR_VALID;

extern	int				PTR_X;
extern	int				PTR_Y;

/* MatView Globals */

extern	REGION_STACK	REGIONS_SAVE;
extern	REGION_STACK	REGIONS;

extern	VIEW			PRINT_VIEW;
extern	VIEW			THUMB_VIEW;
extern	VIEW			MAIN_VIEW;
extern	VIEW			SEP_VIEW;

extern	double			TMP_VALUE;

extern	void			*TMP_PTR;

extern	char			TMP_COORD1[255];
extern	char			TMP_COORD2[255];
extern	char			TMP_COORD3[255];
extern	char			TMP_COORD4[255];

extern	char			TMP_CMD[2048];

extern	char			*MATVIEW_MATFILE;

extern	char			*MATVIEW_DIR;

extern	char			*MV;
extern	char			*MV_M;
extern	char			*MV_G;

extern	char			*PROGRESS_HEIGHT;
extern	char			*PROGRESS_WIDTH;

extern	MATRIX			MAT;

extern	double			*SAVE_FIELD;
extern	int				SAVE_FIELD_NELEMS;

extern	double			COLOR_MIN_VALUE;
extern	double			COLOR_MAX_VALUE;

extern	int				AUTO_COLOR_MAP;

extern	int				GREYSCALE;

extern	int				FILTER_VALUES;

extern	int				COMPLEX_VIEW;

extern	int				PROGRESS_TICKS;

extern	char			*ZERO_COLOR;

extern	double			COLOR_DIST;

extern	int				*COLORIDS;
extern	int				NCOLORIDS;

extern	int				GRID_STATUS;
extern	int				GRID_ROWS;
extern	int				GRID_ROFF;
extern	int				GRID_COLS;
extern	int				GRID_COFF;

extern	int				VIEW_ABS_VALS;

extern	int				VIEW_PATTERN;

extern	int				TRANSPOSE_MATRIX;

extern	int				SHOW_CUMULATIVE;

extern	int				SHOW_MIN_MAX;

extern	double			CUMULATIVE_MIN;
extern	double			CUMULATIVE_MAX;
extern	int				CUMULATIVE_SET;

extern	int				NOLOAD_VALUES;

extern	int				REDUCTION_FUNCTION;

extern	int				NAVIGATE_SIZE;

extern	int				THUMB_BORDER;

extern	int				USE_SEP_VIEW;

extern	int				REGION_ID;

extern	int				DATA_MODE;

extern	int				VIS_MODE;

extern	int				NO_SHELL;

extern	int 			vflag;

#ifdef FILE_BUFFERING

extern	int 			FILE_BUF_SIZE;

#endif

