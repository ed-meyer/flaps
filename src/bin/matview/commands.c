/* $Id: commands.c,v 1.127 2001/07/12 18:38:32 kohl Exp $ */


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

#include <stdlib.h>
#include <unistd.h>

#include "matview.h"
#include "matview_glob.h"


/* ARGSUSED */
int set_mat_name_proc(ClientData clientData,Tcl_Interp* itp,int argc,char** argv )
{
	if ( argc != 2 )
	{
		Tcl_SetResult( itp, "usage: set_mat_name <name>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	if ( strcmp( argv[1], "" ) &&
		( MATVIEW_MATFILE == NULL
			|| strcmp( argv[1], MATVIEW_MATFILE ) ) )
	{
		if ( MATVIEW_MATFILE != NULL )
			FREE( MATVIEW_MATFILE );

		MATVIEW_MATFILE = copy_str( argv[1] );

		load_mat( FALSE );
	}

	return( TCL_OK );
}


/* ARGSUSED */
int
get_mat_name_proc(ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	Tcl_SetResult( itp, copy_str( MATVIEW_MATFILE
		? MATVIEW_MATFILE : "" ), TCL_DYNAMIC );

	return( TCL_OK );
}


/* ARGSUSED */
int
reload_mat_proc(ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	load_mat( TRUE );

	return( TCL_OK );
}


/* ARGSUSED */
int
set_plane_choice_proc (ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	if ( argc != 2 )
	{
		Tcl_SetResult( itp, "usage: set_plane_choice < XY | YZ | XZ >",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	if ( !strcmp( argv[1], "XY" ) )
		PLANE_CHOICE = PLANE_XY;

	else if ( !strcmp( argv[1], "YZ" ) )
		PLANE_CHOICE = PLANE_YZ;

	else if ( !strcmp( argv[1], "XZ" ) )
		PLANE_CHOICE = PLANE_XZ;

	else
	{
		Tcl_SetResult( itp, "usage: set_plane_choice < XY | YZ | XZ >",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	return( TCL_OK );
}


/* ARGSUSED */
int
set_color_range_proc(ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	void **args;

	char *minstr;
	char *maxstr;

	double min, max;

	if ( argc != 3 )
	{
		Tcl_SetResult( itp, "usage: set_color_range <min> <max>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	minstr = argv[1];
	maxstr = argv[2];

	min = atof( minstr );
	max = atof( maxstr );

	if ( !EQUIV( min, COLOR_MIN_VALUE, 0.00001 )
		|| !EQUIV( max, COLOR_MAX_VALUE, 0.00001 ) )
	{
		COLOR_MIN_VALUE = min;
		COLOR_MAX_VALUE = max;

		if ( COLOR_MAX_VALUE < COLOR_MIN_VALUE )
			COLOR_MAX_VALUE = COLOR_MIN_VALUE;

		update_color_range();

		args = make_args( "Color Range", MAIN_VIEW );
		Tk_DoWhenIdle( recolor_matrix, (ClientData) args );

		args = make_args( "Color Range", THUMB_VIEW );
		Tk_DoWhenIdle( recolor_matrix, (ClientData) args );

		Tcl_SetResult( itp, "1", TCL_STATIC );
	}

	else
		Tcl_SetResult( itp, "0", TCL_STATIC );

	return( TCL_OK );
}


/* ARGSUSED */
int
set_color_function_proc(ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	static double last_t = -1.0;

	void **args;

	char *tstr;

	double t;

	if ( argc != 2 )
	{
		Tcl_SetResult( itp, "usage: set_color_function <value>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	tstr = argv[1];

	t = atof( tstr );

	t = ( t < 0.0 ) ? 0.0 : t;
	t = ( t > 1.0 ) ? 1.0 : t;

	if ( !EQUIV( t, last_t, 0.00001 ) )
	{
		last_t = t;

		COLOR_DIST = tan( PI * t / 2.0 );

		/* Draw Color Function Graph */

		draw_color_graph( itp );

		/* Redraw Matrix */

		args = make_args( "Color Function", MAIN_VIEW );
		Tk_DoWhenIdle( recolor_matrix, (ClientData) args );

		args = make_args( "Color Function", THUMB_VIEW );
		Tk_DoWhenIdle( recolor_matrix, (ClientData) args );

		Tcl_SetResult( itp, "1", TCL_STATIC );
	}

	else
		Tcl_SetResult( itp, "0", TCL_STATIC );

	return( TCL_OK );
}


/* ARGSUSED */
int
set_grid_proc(ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	char *onoff;
	char *tmp;

	int nrows, ncols;
	int roff, coff;
	int changed;

	if ( argc != 6 )
	{
		Tcl_SetResult( itp,
			"usage: set_grid <onoff> <nrows> <roff> <ncols> <coff>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	onoff = argv[1];

	tmp = argv[2];
	nrows = atoi( tmp );

	tmp = argv[3];
	roff = atoi( tmp );

	tmp = argv[4];
	ncols = atoi( tmp );

	tmp = argv[5];
	coff = atoi( tmp );

	changed = 0;

	if ( !strcmp( onoff, "on" ) )
	{
		if ( nrows < 1 )
		{
			nrows = 1;

			UPDATE_IENTRY( "grid_row_entry", nrows );
		}

		if ( ncols < 1 )
		{
			ncols = 1;

			UPDATE_IENTRY( "grid_col_entry", ncols );
		}

		if ( !GRID_STATUS )
		{
			GRID_STATUS = TRUE;

			changed++;
		}

		if ( nrows != GRID_ROWS || ncols != GRID_COLS
			|| roff != GRID_ROFF || coff != GRID_COFF )
		{
			GRID_ROWS = nrows;
			GRID_ROFF = roff;
			GRID_COLS = ncols;
			GRID_COFF = coff;

			changed++;
		}
	}

	else if ( GRID_STATUS )
	{
		GRID_STATUS = FALSE;

		changed++;
	}

	if ( changed )
	{
		handle_grid( THUMB_VIEW );
		handle_grid( MAIN_VIEW );
	}

	return( TCL_OK );
}


/* ARGSUSED */
int
set_complex_view_proc(ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	char *cv;

	int changed;
	int tmp;

	if ( argc != 2 )
	{
		Tcl_SetResult( itp,
			"usage: set_complex_view [ re | im | magn | phase ]",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	cv = argv[1];

	if ( !strcmp( cv, "re" ) )
		tmp = COMPLEX_RE;

	else if ( !strcmp( cv, "im" ) )
		tmp = COMPLEX_IM;

	else if ( !strcmp( cv, "magn" ) )
		tmp = COMPLEX_MAGN;

	else if ( !strcmp( cv, "phase" ) )
		tmp = COMPLEX_PHASE;

	changed = 0;

	if ( COMPLEX_VIEW != tmp )
	{
		COMPLEX_VIEW = tmp;

		changed++;
	}

	if ( changed )
		handle_frame( MAIN_VIEW, THUMB_VIEW );

	return( TCL_OK );
}


/* ARGSUSED */
int get_complex_view_proc(ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	char *tmp;

	if ( COMPLEX_VIEW == COMPLEX_RE )
		tmp = "re";

	else if ( COMPLEX_VIEW == COMPLEX_IM )
		tmp = "im";

	else if ( COMPLEX_VIEW == COMPLEX_MAGN )
		tmp = "magn";

	else if ( COMPLEX_VIEW == COMPLEX_PHASE )
		tmp = "phase";

	/* default */
	else
		tmp = "re";

	Tcl_SetResult( itp, tmp, TCL_STATIC );

	return( TCL_OK );
}


/* ARGSUSED */
int
set_zero_color_proc(ClientData clientData, Tcl_Interp* itp, int argc, char** argv ) {
	void **args;

	char *color;

	int changed;

	if ( argc != 2 )
	{
		Tcl_SetResult( itp, "usage: set_zero_color <color>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	color = argv[1];

	changed = 0;

	if ( !strcmp( color, "none" ) )
	{
		if ( ZERO_COLOR != NULL )
		{
			FREE( ZERO_COLOR );
			ZERO_COLOR = (char *) NULL;

			changed++;
		}
	}

	else
	{
		if ( ZERO_COLOR != NULL )
		{
			if ( strcmp( ZERO_COLOR, color ) )
			{
				FREE( ZERO_COLOR );

				ZERO_COLOR = copy_str( color );

				changed++;
			}
		}

		else
		{
			ZERO_COLOR = copy_str( color );

			changed++;
		}
	}

	if ( changed )
	{
		args = make_args( "Zero Color", MAIN_VIEW );
		Tk_DoWhenIdle( recolor_matrix, (ClientData) args );

		args = make_args( "Zero Color", THUMB_VIEW );
		Tk_DoWhenIdle( recolor_matrix, (ClientData) args );
	}

	return( TCL_OK );
}


void
draw_color_graph( Tcl_Interp* itp ) {
	char cmd[255];

	double f;
	double x;

	int x1, y1;
	int x2, y2;
	int ns;
	int fb;
	int i;

	/* Manually Get Color Map Canvas & Size, In Case Startup */

	const char* mv_g = GET_TCL_GLOBAL( itp, "MV_G" );

	const char* tmp = GET_TCL_GLOBAL( itp, "frame_border" );
	fb = atoi( tmp );

	tmp = GET_TCL_GLOBAL( itp, "navigate_size" );
	ns = atoi( tmp );

	ns -= 4 * fb;

	/* Delete Any Existing Lines */

	sprintf( cmd, "%s delete all", mv_g );
	Tcl_Eval( itp, cmd );

	/* Draw New Graph */

	x2 = y2 = -1;

	for ( i=0 ; i <= ns ; i++ )
	{
		x = (double) i / (double) ns;

		f = pow( x, COLOR_DIST );

		x1 = ( 2 * fb ) + (int) ( x * (double) ns );
		y1 = ( 2 * fb ) + ns - (int) ( f * (double) ns );

		if ( x2 != -1 && y2 != -1 )
		{
			sprintf( cmd, "%s create line %d %d %d %d -fill %s",
				mv_g, x1, y1, x2, y2, "purple" );
			Tcl_Eval( itp, cmd );
		}

		x2 = x1;
		y2 = y1;
	}
}


void **
make_args( char* str, VIEW V ) {
	void **tmp;

	tmp = (void **) MALLOC( (unsigned) 2 * sizeof( void * ) );
	memcheck( tmp, "Make Args Pointer" );

	tmp[0] = (void *) str;
	tmp[1] = (void *) V;

	return( tmp );
}


void
recolor_matrix( clientData )
ClientData clientData;
{
	VIEW V;

	RECT R;

	GRID G;

	void **args;

	char color[16];
	char tmp[255];

	char *grid_color;
	char *label;

	int progress;
	int inst;
	int seq;
	int i;

	/* Get Label & View from Client Data */

	args = (void **) clientData;

	label = (char *) args[0];

	V = (VIEW) args[1];

	FREE( args );

	/* Make Sure There's Something to Draw */

	if ( MAT->loaded == ML_UNLOADED || V->rectids == NULL )
		return;

	/* Grab Frame Instance Number & Bump/Grab Drawing Sequence Number */
	/* (Double Trouble!) */

	inst = V->instance;

	seq = ++(V->sequence);

	/* Prep for Recolor */

	sprintf( tmp, "setMsg {Redrawing for New %s...}", label );

	Tcl_Eval( interp, tmp );
	INST_CHECK( V, inst, V->instance, FALSE );
	INST_CHECK( V, seq, V->sequence, FALSE );

	update_progress( V, 0 );
	INST_CHECK( V, inst, V->instance, FALSE );
	INST_CHECK( V, seq, V->sequence, FALSE );

	hide_image( V );
	INST_CHECK( V, inst, V->instance, TRUE );
	INST_CHECK( V, seq, V->sequence, TRUE );

	progress = V->nrectids / PROGRESS_TICKS;
	progress = ( progress ) ? progress : 1;

	/* Draw Zero Rect Background */

	sprintf( tmp, "clear_image %s %d %d", V->image, V->wt, V->ht );
	Tcl_Eval( interp, tmp );

	if ( ZERO_COLOR != NULL )
		strcpy( color, ZERO_COLOR );

	else
		color_value( color, 0.0 );

	DRAW_IMAGE_RECT( interp, V->image,
		V->win->xoff, V->win->yoff,
		( V->iscale * V->nx ) + V->win->xoff,
		( V->iscale * V->ny ) + V->win->yoff,
		color );

	DRAW_IMAGE_OUTLINE( interp, V->image,
		V->win->xoff - 1, V->win->yoff - 1,
		( V->iscale * V->nx ) + V->win->xoff,
		( V->iscale * V->ny ) + V->win->yoff,
		"black" );

	/* Draw Rects */

	R = &(V->rectids[0]);

	for ( i=0 ; i < V->nrectids ; i++ )
	{
		if ( !( i % progress ) )
		{
			update_progress( V, i / progress );
			INST_CHECK( V, inst, V->instance, TRUE );
			INST_CHECK( V, seq, V->sequence, TRUE );
		}

		if ( R->drawn == V->instance && R->value != 0.0
			&& ( !FILTER_VALUES
				|| ( R->value >= COLOR_MIN_VALUE
					&& R->value <= COLOR_MAX_VALUE ) ) )
		{
			color_value( color, R->value );

			if ( R->x1 != R->x2 && R->y1 != R->y2 )
			{
				DRAW_IMAGE_RECT( interp, V->image,
					R->x1, R->y1, R->x2, R->y2,
					color );
			}

			else
			{
				DRAW_IMAGE_POINT( interp, V->image,
					R->x1, R->y1, color );
			}
		}

		R++;
	}

	/* Draw Grid Lines */

	if ( V->grids != NULL )
	{
		grid_color =
			( ZERO_COLOR != NULL && !strcmp( ZERO_COLOR, "black" ) )
			? "white" : "black";

		G = &(V->grids[0]);

		for ( i=0 ; i < V->ngrids ; i++ )
		{
			if ( G->id != -1 )
				LINE_COLOR( interp, V->win->canvas, G->id, grid_color );
			
			G++;
		}
	}

	/* Done, Show Image */

	show_image( V );
	INST_CHECK( V, inst, V->instance, FALSE );
	INST_CHECK( V, seq, V->sequence, FALSE );

	Tcl_Eval( interp,
		"setMsg {New Color Range Redraw Complete.}" );
	INST_CHECK( V, inst, V->instance, FALSE );
	INST_CHECK( V, seq, V->sequence, FALSE );

	update_progress( V, 0 );
	INST_CHECK( V, inst, V->instance, FALSE );
	INST_CHECK( V, seq, V->sequence, FALSE );
}


void
update_color_range()
{
	UPDATE_FENTRY( "color_min", COLOR_MIN_VALUE );
	UPDATE_FENTRY( "color_max", COLOR_MAX_VALUE );
}


/* ARGSUSED */
int draw_color_map_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	static double color_dist_save = -1.0;
	static int equiv_save = -1;
	static int csize_save = -1;
	static int grey_save = -1;

	char color[16];

	double cstart;
	double cstep;
	double c;

	int *C, *C2;

	int *newids;

	int ncolorids;
	int csize;
	int equiv;
	int wt;
	int x;
	int i;

	/* Manually Get Color Map Canvas, In Case Startup */

	const char* mv_m = GET_TCL_GLOBAL( itp, "MV_M" );

	/* Get New Map Size */

	const char* color_size = GET_TCL_GLOBAL( itp, "color_size" );

	csize = atoi( color_size );

	csize -= 4 * FRAME_BORDER;

	if ( !csize )
	{
		Tcl_SetResult( itp, "Error:  Color Map Width of Zero",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	/* Verify Map Needs Redraw */

	if ( EQUIV( COLOR_MIN_VALUE, COLOR_MAX_VALUE, 0.000001 ) )
	{
		cstart = COLOR_MIN_VALUE - 1.0;
		cstep = 1.0;

		ncolorids = 3;

		equiv = 1;
	}

	else
	{
		cstart = COLOR_MIN_VALUE;
		cstep = ( COLOR_MAX_VALUE - COLOR_MIN_VALUE )
			 / (double) csize;

		ncolorids = csize;

		equiv = 0;
	}

	if ( csize == csize_save && equiv == equiv_save
		&& color_dist_save == COLOR_DIST
		&& grey_save == GREYSCALE )
	{
		return( TCL_OK );
	}

	/* Verify Allocation of Color Ids */

	if ( ncolorids > NCOLORIDS || COLORIDS == NULL )
	{
		newids = (int *) MALLOC( (unsigned) ncolorids * sizeof(int) );
		memcheck( newids, "Color Map ID Array" );

		C = &(newids[0]);

		if ( COLORIDS != NULL )
		{
			C2 = &(COLORIDS[0]);

			for ( i=0 ; i < NCOLORIDS ; i++ )
				*C++ = *C2++;

			for ( i=NCOLORIDS ; i < ncolorids ; i++ )
				*C++ = -1;

			FREE( COLORIDS );
		}

		else
		{
			for ( i=0 ; i < ncolorids ; i++ )
				*C++ = -1;
		}

		COLORIDS = newids;

		NCOLORIDS = ncolorids;
	}

	/* Draw New Map */

	x = FRAME_BORDER + 1;

	c = cstart;

	for ( i=0 ; i < ncolorids ; i++ )
	{
		/* Set Color Bar Width - Different for Binary Case... */

		if ( equiv )
		{
			switch ( i )
			{
				case 0:
				case 2:
					wt = ( csize - BORDER_SPACE ) / 2;
					break;
				
				case 1:
					wt = BORDER_SPACE;
					break;
				
				default:
					wt = 1;
					break;
			}
		}

		else
			wt = 1;

		/* Draw Color Bar */

		color_value( color, c );

		if ( COLORIDS[i] != -1 )
		{
			if ( csize != csize_save )
			{
				RECT_COORDS( itp, mv_m, COLORIDS[i],
					x, 0, x + wt, 100 );
			}

			RECT_COLOR( itp, mv_m, COLORIDS[i], color, color );
		}

		else
		{
			COLORIDS[i] = DRAW_RECT( itp, mv_m, x, 0, x + wt, 100,
				color, color, 1, "color_bar" );
		}

		c += cstep;

		x += wt;
	}

	/* Free Any Leftover Color Rects */

	for ( i=ncolorids ; i < NCOLORIDS ; i++ )
	{
		if ( COLORIDS[i] != -1 )
		{
			DELETE_ITEM( itp, mv_m, COLORIDS[i] );

			COLORIDS[i] = -1;
		}
	}

	color_dist_save = COLOR_DIST;
	grey_save = GREYSCALE;
	csize_save = csize;
	equiv_save = equiv;

	return( TCL_OK );
}


/* ARGSUSED */
int set_vis_region_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	int xmin, xmax, xcell;
	int ymin, ymax, ycell;
	int xi, yi, zi;
	int zval;

	if ( MAT->loaded == ML_UNLOADED )
		return( TCL_OK );

	if ( argc != 8 )
	{
		Tcl_SetResult( itp,
		"usage: set_vis_region xmin xmax xcell ymin ymax ycell zval",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	xmin =	atoi( argv[1] );
	xmax =	atoi( argv[2] );
	xcell =	atoi( argv[3] );

	ymin =	atoi( argv[4] );
	ymax =	atoi( argv[5] );
	ycell =	atoi( argv[6] );

	zval =	atoi( argv[7] );

	if ( !get_plane_indices( &xi, &yi, &zi ) )
	{
		Tcl_SetResult( itp,
			"Error:  Invalid Plane Selection in set_vis_region{}",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	if ( REGION_MIN( MAIN_VIEW->visregion, xi ) != xmin
		|| REGION_MAX( MAIN_VIEW->visregion, xi ) != xmax
		|| CELL_OF_DIM( MAIN_VIEW->viscell, xi ) != xcell
		|| REGION_MIN( MAIN_VIEW->visregion, yi ) != ymin
		|| REGION_MAX( MAIN_VIEW->visregion, yi ) != ymax
		|| CELL_OF_DIM( MAIN_VIEW->viscell, yi ) != ycell
		|| REGION_MIN( MAIN_VIEW->visregion, zi ) != zval
		|| REGION_MAX( MAIN_VIEW->visregion, zi ) != zval
		|| CELL_OF_DIM( MAIN_VIEW->viscell, zi ) != 1 )
	{
		REGION_MIN( MAIN_VIEW->visregion, xi ) = xmin;
		REGION_MAX( MAIN_VIEW->visregion, xi ) = xmax;

		CELL_OF_DIM( MAIN_VIEW->viscell, xi ) = xcell;

		REGION_MIN( MAIN_VIEW->visregion, yi ) = ymin;
		REGION_MAX( MAIN_VIEW->visregion, yi ) = ymax;

		CELL_OF_DIM( MAIN_VIEW->viscell, yi ) = ycell;

		REGION_MIN( MAIN_VIEW->visregion, zi ) =
			REGION_MAX( MAIN_VIEW->visregion, zi ) = zval;

		CELL_OF_DIM( MAIN_VIEW->viscell, zi ) = 1;

		if ( check_region( MAIN_VIEW->visregion, MAIN_VIEW->viscell,
				MAT->glb, MAT->gub, MAT->dim ) )
			update_vis_region();

		if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
			draw_thumb_border();
	}

	return( TCL_OK );
}


/* ARGSUSED */
int set_sep_view_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *onoff;

	char cmd[255];

	int tmp;

	if ( argc != 2 )
		goto err;

	onoff = argv[1];

	if ( !strcmp( onoff, "OFF" ) )
		tmp = FALSE;

	else if ( !strcmp( onoff, "ON" ) )
		tmp = TRUE;

	else
	{
		printf( "set_sep_view{}:  Unknown Option \"%s\".\n", onoff );

		goto err;
	}

	if ( tmp != USE_SEP_VIEW )
	{
		USE_SEP_VIEW = tmp;

		sprintf( cmd, "set_view_size %s -1 -1 -1 -1",
			MAIN_VIEW->win->toplevel );
		Tcl_Eval( interp, cmd );

		update_progress( MAIN_VIEW, 0 );

		swap_views( &SEP_VIEW, &MAIN_VIEW );

		if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
			draw_thumb_border();
	}

	return( TCL_OK );

err:

	Tcl_SetResult( itp, "usage: set_sep_view < ON | OFF >",
		TCL_STATIC );

	return( TCL_ERROR );
}


/* ARGSUSED */
int set_data_mode_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *mode;

	int tmp;

	if ( argc != 2 )
		goto err;

	mode = argv[1];

	if ( !strcmp( mode, "matrix" ) )
		tmp = DATA_MODE_MATRIX;

	else if ( !strcmp( mode, "field" ) )
		tmp = DATA_MODE_FIELD;

	else
	{
		printf( "set_data_mode{}:  Unknown Data Mode \"%s\".\n", mode );

		goto err;
	}

	if ( tmp != DATA_MODE )
	{
		DATA_MODE = tmp;

		Tcl_Eval( interp, "layout_main_panel" );

		handle_frame( MAIN_VIEW, THUMB_VIEW );
	}

	return( TCL_OK );

err:

	Tcl_SetResult( itp, "usage: set_data_mode [ matrix | field ]",
		TCL_STATIC );

	return( TCL_ERROR );
}


/* ARGSUSED */
int get_data_mode_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *tmp;

	if ( DATA_MODE == DATA_MODE_MATRIX )
		tmp = "matrix";

	else if ( DATA_MODE == DATA_MODE_FIELD )
		tmp = "field";

	/* default */
	else
		tmp = "matrix";

	Tcl_SetResult( itp, tmp, TCL_STATIC );

	return( TCL_OK );
}


/* ARGSUSED */
int set_vis_mode_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *mode;

	int tmp;

	if ( argc != 2 )
		goto err;

	mode = argv[1];

	if ( !strcmp( mode, "normal" ) )
		tmp = VIS_MODE_NORMAL;

	else if ( !strcmp( mode, "diff" ) )
	{
		if ( SAVE_FIELD != NULL )
			tmp = VIS_MODE_DIFF;
		
		else
		{
			printf( "set_vis_mode{}:  No Saved Field - " );
			printf( "Reverting to Normal Mode\n" );

			tmp = VIS_MODE_NORMAL;

			SET_TCL_GLOBAL( itp, "vis_mode", "normal" );
		}
	}

	else
	{
		printf( "set_vis_mode{}:  Unknown Vis Mode \"%s\".\n", mode );

		goto err;
	}

	if ( tmp != VIS_MODE )
	{
		VIS_MODE = tmp;

		handle_frame( MAIN_VIEW, THUMB_VIEW );
	}

	return( TCL_OK );

err:

	Tcl_SetResult( itp, "usage: set_vis_mode [ normal | diff ]",
		TCL_STATIC );

	return( TCL_ERROR );
}


/* ARGSUSED */
int get_vis_mode_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *tmp;

	if ( VIS_MODE == VIS_MODE_NORMAL )
		tmp = "normal";

	else if ( VIS_MODE == VIS_MODE_DIFF )
		tmp = "diff";

	/* default */
	else
		tmp = "normal";

	Tcl_SetResult( itp, tmp, TCL_STATIC );

	return( TCL_OK );
}


/* ARGSUSED */
int set_opt_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char usage[255];

	void **args;

	char *onoff;
	char *opt;

	double min, max;

	int thumb_redraw;
	int recolor;
	int redraw;
	int tmp;

	if ( argc != 3 )
		goto err;

	opt = argv[1];
	onoff = argv[2];

	/* Decode On vs. Off */

	if ( !strcmp( onoff, "ON" ) )
		tmp = TRUE;

	else if ( !strcmp( onoff, "OFF" ) )
		tmp = FALSE;

	else if ( strcmp( onoff, "none" )
			|| strcmp( opt, "reset_cumulative" ) )
		goto err;

	/* Set Option */

	thumb_redraw = 0;
	recolor = 0;
	redraw = 0;

	if ( !strcmp( opt, "auto_color" ) )
	{
		AUTO_COLOR_MAP = tmp;

		if ( tmp && MAIN_VIEW->field != NULL )
		{
			compute_view_min_max( MAIN_VIEW, &min, &max );

			handle_min_max( min, max );

			recolor++;
		}
	}

	else if ( !strcmp( opt, "greyscale" ) )
	{
		GREYSCALE = tmp;

		recolor++;
	}

	else if ( !strcmp( opt, "filter_values" ) )
	{
		FILTER_VALUES = tmp;

		recolor++;
	}

	else if ( !strcmp( opt, "absvals" ) )
	{
		VIEW_ABS_VALS = tmp;

		thumb_redraw++;
	}

	else if ( !strcmp( opt, "pattern" ) )
	{
		VIEW_PATTERN = tmp;

		thumb_redraw++;
	}

	else if ( !strcmp( opt, "transpose" ) )
	{
		TRANSPOSE_MATRIX = tmp;

		if ( !transpose_matrix() )
			return( TCL_OK );

		else
			thumb_redraw++;
	}

	else if ( !strcmp( opt, "min_max" ) )
	{
		SHOW_MIN_MAX = tmp;

		if ( MAIN_VIEW->field != NULL )
		{
			compute_view_min_max( MAIN_VIEW, &min, &max );

			handle_min_max( min, max );
		}
	}

	else if ( !strcmp( opt, "cumulative" ) )
	{
		SHOW_CUMULATIVE = tmp;

		if ( MAIN_VIEW->field != NULL )
		{
			compute_view_min_max( MAIN_VIEW, &min, &max );

			handle_min_max( min, max );
		}
	}

	else if ( !strcmp( opt, "reset_cumulative" ) )
	{
		CUMULATIVE_SET = 0;

		if ( MAIN_VIEW->field != NULL )
		{
			compute_view_min_max( MAIN_VIEW, &min, &max );

			handle_min_max( min, max );
		}
	}

	else if ( !strcmp( opt, "noload_values" ) )
	{
		NOLOAD_VALUES = tmp;

		return( TCL_OK );
	}

	else
		goto err;

	if ( thumb_redraw )
		handle_frame( MAIN_VIEW, THUMB_VIEW );

	else if ( redraw )
	{
		if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
		{
			args = make_args( "Option Change", THUMB_VIEW );
			Tk_DoWhenIdle( recolor_matrix, (ClientData) args );
		}
	}

	else if ( recolor )
	{
		args = make_args( "Option Change", MAIN_VIEW );
		Tk_DoWhenIdle( recolor_matrix, (ClientData) args );

		args = make_args( "Option Change", THUMB_VIEW );
		Tk_DoWhenIdle( recolor_matrix, (ClientData) args );
	}

	return( TCL_OK );

err:

	sprintf( usage,
		"usage: set_opt [%s|%s|%s|%s|%s] [ON|OFF]",
		"auto_color", "transpose", "min_max", "cumulative",
		"reset_cumulative" );

	Tcl_SetResult( itp, copy_str( usage ), TCL_DYNAMIC );

	return( TCL_ERROR );
}


/* ARGSUSED */
int set_reduction_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *rf;

	if ( argc != 2 )
		goto err;

	rf = argv[1];

	if ( !strcmp( rf, "avg" ) )
		REDUCTION_FUNCTION = RF_AVG;

	else if ( !strcmp( rf, "min" ) )
		REDUCTION_FUNCTION = RF_MIN;

	else if ( !strcmp( rf, "max" ) )
		REDUCTION_FUNCTION = RF_MAX;

	else if ( !strcmp( rf, "sum" ) )
		REDUCTION_FUNCTION = RF_SUM;

	else if ( !strcmp( rf, "midpoint" ) )
		REDUCTION_FUNCTION = RF_MIDPOINT;

	else if ( !strcmp( rf, "count" ) )
		REDUCTION_FUNCTION = RF_COUNT;

	else if ( !strcmp( rf, "std" ) )
		REDUCTION_FUNCTION = RF_STD;

	else if ( !strcmp( rf, "var" ) )
		REDUCTION_FUNCTION = RF_VAR;

	else
		goto err;

	handle_frame( MAIN_VIEW, THUMB_VIEW );

	return( TCL_OK );

err:

	Tcl_SetResult( itp,
		"usage: set_reduction [Avg|Min|Max|Sum|Mean|Std|Var]",
		TCL_STATIC );

	return( TCL_ERROR );
}


/* ARGSUSED */
int get_reduction_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *tmp;

	if ( REDUCTION_FUNCTION == RF_AVG )
		tmp = "avg";

	else if ( REDUCTION_FUNCTION == RF_MIN )
		tmp = "min";

	else if ( REDUCTION_FUNCTION == RF_MAX )
		tmp = "max";

	else if ( REDUCTION_FUNCTION == RF_SUM )
		tmp = "sum";

	else if ( REDUCTION_FUNCTION == RF_MIDPOINT )
		tmp = "midpoint";

	else if ( REDUCTION_FUNCTION == RF_COUNT )
		tmp = "count";

	else if ( REDUCTION_FUNCTION == RF_STD )
		tmp = "std";

	else if ( REDUCTION_FUNCTION == RF_VAR )
		tmp = "var";

	/* default */
	else
		tmp = "sum";

	Tcl_SetResult( itp, tmp, TCL_STATIC );

	return( TCL_OK );
}


void
update_vis_region()
{
	int xi, yi, zi;

	if ( !get_plane_indices( &xi, &yi, &zi ) )
	{
		printf( "\nError:  Invalid Plane Selection in %s\n\n",
			"update_vis_reg()" );

		return;
	}

	UPDATE_IENTRY( "xmin_entry",
		REGION_MIN( MAIN_VIEW->visregion, xi ) );
	UPDATE_IENTRY( "xmax_entry",
		REGION_MAX( MAIN_VIEW->visregion, xi ) );
	UPDATE_IENTRY( "xcell_entry",
		CELL_OF_DIM( MAIN_VIEW->viscell, xi ) );

	UPDATE_IENTRY( "ymin_entry",
		REGION_MIN( MAIN_VIEW->visregion, yi ) );
	UPDATE_IENTRY( "ymax_entry",
		REGION_MAX( MAIN_VIEW->visregion, yi ) );
	UPDATE_IENTRY( "ycell_entry",
		CELL_OF_DIM( MAIN_VIEW->viscell, yi ) );

	UPDATE_IENTRY( "zval_entry",
		REGION_MIN( MAIN_VIEW->visregion, zi ) );
}


int
get_plane_indices( xi, yi, zi )
int *xi;
int *yi;
int *zi;
{
	switch ( PLANE_CHOICE )
	{
		case PLANE_XY:
			*xi = 0;
			*yi = 1;
			*zi = 2;
			break;

		case PLANE_YZ:
			*xi = 1;
			*yi = 2;
			*zi = 0;
			break;

		case PLANE_XZ:
			*xi = 0;
			*yi = 2;
			*zi = 1;
			break;

		default:
		{
			printf( "\nError in %s:  Invalid Plane Selection (%d)\n\n",
				"get_plane_indices()", PLANE_CHOICE );

			return( FALSE );
		}
	}

	return( TRUE );
}


/* ARGSUSED */
int resize_mat_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
		draw_thumb_border();

	return( TCL_OK );
}


/* ARGSUSED */
int save_field_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *cmd;

	int do_usage = 0;

	if ( argc < 2 || argc > 3 )
		do_usage++;

	else
	{
		cmd = argv[1];

		if ( !strcmp( cmd, "set" ) )
			set_current_field();

		else if ( !strcmp( cmd, "free" ) )
			free_current_field();

		else if ( !strcmp( cmd, "save" ) )
		{
			if ( argc != 3 || !strcmp( argv[2], "" ) )
				do_usage++;
			
			else
				save_current_field( argv[2] );
		}

		else if ( !strcmp( cmd, "load" ) )
		{
			if ( argc != 3 || !strcmp( argv[2], "" ) )
				do_usage++;
			
			else
				load_current_field( argv[2] );
		}

		else
			do_usage++;
	}

	if ( do_usage )
	{
		Tcl_SetResult( itp,
		"usage: save_field [set] [free] [save <file>] [load <file>]",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	return( TCL_OK );
}


/* ARGSUSED */
int canvas_handle_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	static int SAVEID = -1;

	static int markx = -1;
	static int marky = -1;

	RECT R;

	char cmd[1024];

	char *color;
	char *op;

	int xi, yi, zi;
	int coords[3];

	int xmin, xmax;
	int ymin, ymax;

	int offset;
	int dx, dy;
	int rx, ry;
	int r, c;
	int x, y;
	int flag;

	if ( MAIN_VIEW->win == NULL )
		return( TCL_OK );

	if ( argc != 4 )
	{
		Tcl_SetResult( itp, "usage: canvas_handle <op> <x> <y>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	op = argv[1];

	x = atoi( argv[2] );
	y = atoi( argv[3] );

	/* Munge coords a little, so pointer's not in the way... */

	x -= 3;
	y -= 3;

	if ( vflag )
	{
		printf( "(%d, %d) -> (%d, %d)\n", x, y,
			SCREEN_TO_REGION( x, MAIN_VIEW->win->xoff,
				REGION_MIN( MAIN_VIEW->visregion, 0 ),
				REGION_MAX( MAIN_VIEW->visregion, 0 ),
				MAIN_VIEW->nx, MAIN_VIEW->iscale ),
			SCREEN_TO_REGION( y, MAIN_VIEW->win->yoff,
				REGION_MIN( MAIN_VIEW->visregion, 1 ),
				REGION_MAX( MAIN_VIEW->visregion, 1 ),
				MAIN_VIEW->ny, MAIN_VIEW->iscale ) );
	}

	/* Handle Command */

	if ( !strcmp( op, "query_press" ) )
	{
		PTR_X = x;
		PTR_Y = y;

		R = find_rect( MAIN_VIEW, x, y );

		if ( R != NULL )
		{
			if ( get_plane_indices( &xi, &yi, &zi ) )
			{
				color = ( ZERO_COLOR != NULL
						&& !strcmp( ZERO_COLOR, "black" ) )
					? "white" : "black";

				SAVEID = DRAW_RECT( itp, MAIN_VIEW->win->canvas,
					R->x1, R->y1, R->x2, R->y2, "", color, 1,
					"query" );

				coords[xi] = SCREEN_TO_REGION( x, MAIN_VIEW->win->xoff,
					REGION_MIN( MAIN_VIEW->visregion, xi ),
					REGION_MAX( MAIN_VIEW->visregion, xi ),
					MAIN_VIEW->nx, MAIN_VIEW->iscale );
				coords[yi] = SCREEN_TO_REGION( y, MAIN_VIEW->win->yoff,
					REGION_MIN( MAIN_VIEW->visregion, yi ),
					REGION_MAX( MAIN_VIEW->visregion, yi ),
					MAIN_VIEW->ny, MAIN_VIEW->iscale );
				coords[zi] = SCREEN_TO_REGION( 0, 0,
					REGION_MIN( MAIN_VIEW->visregion, zi ),
					REGION_MAX( MAIN_VIEW->visregion, zi ),
					1, MAIN_VIEW->iscale );

				sprintf( cmd,
					"query_handle %s %s %d %d \"(%d,%d,%d) = %lg\"",
					MAIN_VIEW->win->toplevel, MAIN_VIEW->win->canvas,
					x, y, coords[2], coords[1], coords[0], R->value );

				Tcl_Eval( itp, cmd );
			}
		}
	}

	else if ( !strcmp( op, "query_slide" ) )
	{
		PTR_X = x;
		PTR_Y = y;

		R = find_rect( MAIN_VIEW, x, y );

		if ( R != NULL )
		{
			if ( get_plane_indices( &xi, &yi, &zi ) )
			{
				if ( SAVEID != -1 )
				{
					RECT_COORDS( itp, MAIN_VIEW->win->canvas, SAVEID,
						R->x1, R->y1, R->x2, R->y2 );
				}

				else
				{
					color = ( ZERO_COLOR != NULL
							&& !strcmp( ZERO_COLOR, "black" ) )
						? "white" : "black";

					SAVEID = DRAW_RECT( itp, MAIN_VIEW->win->canvas,
						R->x1, R->y1, R->x2, R->y2, "", color, 1,
						"query" );
				}

				coords[xi] = SCREEN_TO_REGION( x, MAIN_VIEW->win->xoff,
					REGION_MIN( MAIN_VIEW->visregion, xi ),
					REGION_MAX( MAIN_VIEW->visregion, xi ),
					MAIN_VIEW->nx, MAIN_VIEW->iscale );
				coords[yi] = SCREEN_TO_REGION( y, MAIN_VIEW->win->yoff,
					REGION_MIN( MAIN_VIEW->visregion, yi ),
					REGION_MAX( MAIN_VIEW->visregion, yi ),
					MAIN_VIEW->ny, MAIN_VIEW->iscale );
				coords[zi] = SCREEN_TO_REGION( 0, 0,
					REGION_MIN( MAIN_VIEW->visregion, zi ),
					REGION_MAX( MAIN_VIEW->visregion, zi ),
					1, MAIN_VIEW->iscale );

				sprintf( cmd,
					"query_handle %s %s %d %d \"(%d,%d,%d) = %lg\"",
					MAIN_VIEW->win->toplevel, MAIN_VIEW->win->canvas,
					x, y, coords[2], coords[1], coords[0], R->value );

				Tcl_Eval( itp, cmd );
			}
		}

		else
		{
			if ( SAVEID != -1 )
			{
				DELETE_ITEM( itp, MAIN_VIEW->win->canvas, SAVEID );

				SAVEID = -1;
			}

			sprintf( cmd, "query_handle %s %s -1 -1 \"\"",
				MAIN_VIEW->win->toplevel, MAIN_VIEW->win->canvas );

			Tcl_Eval( itp, cmd );
		}
	}

	else if ( !strcmp( op, "query_release" ) )
	{
		PTR_X = x;
		PTR_Y = y;

		if ( SAVEID != -1 )
		{
			DELETE_ITEM( itp, MAIN_VIEW->win->canvas, SAVEID );

			SAVEID = -1;
		}

		sprintf( cmd, "query_handle %s %s -1 -1 \"\"",
			MAIN_VIEW->win->toplevel, MAIN_VIEW->win->canvas );

		Tcl_Eval( itp, cmd );
	}

	else if ( !strcmp( op, "zoom_press" ) )
	{
		PTR_X = x;
		PTR_Y = y;

		/* Create Zoom Rect */

		if ( REGION_ID != -1 )
			DELETE_ITEM( itp, MAIN_VIEW->win->canvas, REGION_ID );
	
		color = ( ZERO_COLOR != NULL && !strcmp( ZERO_COLOR, "black" ) )
			? "white" : "black";

		REGION_ID = DRAW_RECT( itp, MAIN_VIEW->win->canvas,
			x, y, x, y, "", color, 1, "zoom" );

		markx = x;
		marky = y;
	}

	else if ( !strcmp( op, "zoom_slide" ) )
	{
		PTR_X = x;
		PTR_Y = y;

		/* Adjust Zoom Rect */

		if ( REGION_ID != -1 )
		{
			RECT_COORDS( itp, MAIN_VIEW->win->canvas, REGION_ID,
				markx, marky, x, y );
		}
	}

	else if ( !strcmp( op, "zoom_release" ) )
	{
		PTR_X = x;
		PTR_Y = y;

		/* Set Desired Vis Region */

		DELETE_ITEM( itp, MAIN_VIEW->win->canvas, REGION_ID );

		if ( MAT->loaded == ML_UNLOADED )
			return( TCL_OK );

		if ( markx != x && marky != y )
		{
			/* Determine Min/Max Coords */

			if ( markx < x )
			{
				xmin = markx;
				xmax = x;
			}

			else
			{
				xmin = x;
				xmax = markx;
			}

			if ( marky < y )
			{
				ymin = marky;
				ymax = y;
			}

			else
			{
				ymin = y;
				ymax = marky;
			}

			/* Set New Region Bounds */

			xmin = SCREEN_TO_REGION( xmin, MAIN_VIEW->win->xoff,
				REGION_MIN( MAIN_VIEW->visregion, 0 ),
				REGION_MAX( MAIN_VIEW->visregion, 0 ),
				MAIN_VIEW->nx, MAIN_VIEW->iscale );

			ymin = SCREEN_TO_REGION( ymin, MAIN_VIEW->win->yoff,
				REGION_MIN( MAIN_VIEW->visregion, 1 ),
				REGION_MAX( MAIN_VIEW->visregion, 1 ),
				MAIN_VIEW->ny, MAIN_VIEW->iscale );

			xmax = SCREEN_TO_REGION( xmax, MAIN_VIEW->win->xoff,
				REGION_MIN( MAIN_VIEW->visregion, 0 ),
				REGION_MAX( MAIN_VIEW->visregion, 0 ),
				MAIN_VIEW->nx, MAIN_VIEW->iscale );

			ymax = SCREEN_TO_REGION( ymax, MAIN_VIEW->win->yoff,
				REGION_MIN( MAIN_VIEW->visregion, 1 ),
				REGION_MAX( MAIN_VIEW->visregion, 1 ),
				MAIN_VIEW->ny, MAIN_VIEW->iscale );

			if ( REGION_MIN( MAIN_VIEW->visregion, 0 ) != xmin
				|| REGION_MIN( MAIN_VIEW->visregion, 1 ) != ymin
				|| REGION_MAX( MAIN_VIEW->visregion, 0 ) != xmax
				|| REGION_MAX( MAIN_VIEW->visregion, 1 ) != ymax )
			{
				/* Save Old Region */

				push_region( MAIN_VIEW->visregion );

				/* Check New Region */

				REGION_MIN( MAIN_VIEW->visregion, 0 ) = xmin;
				REGION_MIN( MAIN_VIEW->visregion, 1 ) = ymin;

				REGION_MAX( MAIN_VIEW->visregion, 0 ) = xmax;
				REGION_MAX( MAIN_VIEW->visregion, 1 ) = ymax;

				check_region( MAIN_VIEW->visregion, MAIN_VIEW->viscell,
					MAT->glb, MAT->gub, MAT->dim  );

				/* Send New Vis Region */

				update_vis_region();

				if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
					draw_thumb_border();
			}
		}

		else if ( GRID_STATUS && MAIN_VIEW->do_grid )
		{
			/* Undo Coords Munge... */

			x += 3;
			y += 3;

			/* Determine Region Coords */

			rx = SCREEN_TO_REGION( x, MAIN_VIEW->win->xoff,
				REGION_MIN( MAIN_VIEW->visregion, 0 ),
				REGION_MAX( MAIN_VIEW->visregion, 0 ),
				MAIN_VIEW->nx, MAIN_VIEW->iscale );

			ry = SCREEN_TO_REGION( y, MAIN_VIEW->win->yoff,
				REGION_MIN( MAIN_VIEW->visregion, 1 ),
				REGION_MAX( MAIN_VIEW->visregion, 1 ),
				MAIN_VIEW->ny, MAIN_VIEW->iscale );

			/* Determine Current Grid Cell */

			if ( rx < GRID_COFF )
			{
				xmin = REGION_MIN( MAX_REGION, 0 );
				xmax = xmin + GRID_COFF - 1;
			}

			else
			{
				offset = REGION_MIN( MAX_REGION, 0 ) + GRID_COFF;

				c = ( rx - offset ) / GRID_COLS;

				xmin = ( c * GRID_COLS ) + offset;
				xmax = xmin + GRID_COLS - 1;
			}

			if ( ry < GRID_ROFF )
			{
				ymin = REGION_MIN( MAX_REGION, 1 );
				ymax = ymin + GRID_ROFF - 1;
			}

			else
			{
				offset = REGION_MIN( MAX_REGION, 1 ) + GRID_ROFF;

				r = ( ry - offset ) / GRID_ROWS;

				ymin = ( r * GRID_ROWS ) + offset;
				ymax = ymin + GRID_ROWS - 1;
			}

			if ( REGION_MIN( MAIN_VIEW->visregion, 0 ) != xmin
				|| REGION_MIN( MAIN_VIEW->visregion, 1 ) != ymin
				|| REGION_MAX( MAIN_VIEW->visregion, 0 ) != xmax
				|| REGION_MAX( MAIN_VIEW->visregion, 1 ) != ymax )
			{
				/* Save Old Region */

				push_region( MAIN_VIEW->visregion );

				/* Check New Vis Region */

				REGION_MIN( MAIN_VIEW->visregion, 0 ) = xmin;
				REGION_MIN( MAIN_VIEW->visregion, 1 ) = ymin;

				REGION_MAX( MAIN_VIEW->visregion, 0 ) = xmax;
				REGION_MAX( MAIN_VIEW->visregion, 1 ) = ymax;

				check_region( MAIN_VIEW->visregion, MAIN_VIEW->viscell,
					MAT->glb, MAT->gub, MAT->dim );

				/* Set New Vis Region */

				update_vis_region();

				if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
					draw_thumb_border();
			}
		}
	}

	else if ( !strcmp( op, "nav_press" ) )
	{
		if ( MAT->loaded != ML_UNLOADED )
		{
			dx = x + 3 - ( ( NAVIGATE_SIZE / 2 ) - FRAME_BORDER );
			dy = y + 3 - ( ( NAVIGATE_SIZE / 2 ) - FRAME_BORDER );

			dx = ( abs( dx ) < NAVIGATE_SIZE / 10 )
				? 0 : ( dx / abs( dx ) );
			dy = ( abs( dy ) < NAVIGATE_SIZE / 10 )
				? 0 : ( dy / abs( dy ) );

			if ( dx != 0 || dy != 0 )
				nudge_view( dx, dy );
		}
	}

	else if ( !strcmp( op, "Up" ) || !strcmp( op, "KP_Up" ) )
		nudge_view( 0, -1 );

	else if ( !strcmp( op, "Down" ) || !strcmp( op, "KP_Down" ) )
		nudge_view( 0, 1 );

	else if ( !strcmp( op, "Left" ) || !strcmp( op, "KP_Left" ) )
		nudge_view( -1, 0 );

	else if ( !strcmp( op, "Right" ) || !strcmp( op, "KP_Right" ) )
		nudge_view( 1, 0 );

	else if ( !strcmp( op, "pop" ) )
	{
		if ( double_click( "pop" ) )
		{
			flag = 0;

			while ( pop_region( MAIN_VIEW->visregion ) )
				flag++;

			if ( flag )
			{
				update_vis_region();

				if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
					draw_thumb_border();
			}
		}

		else if ( pop_region( MAIN_VIEW->visregion ) )
		{
			update_vis_region();

			if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
				draw_thumb_border();
		}
	}

	else if ( !strcmp( op, "Prior" ) || !strcmp( op, "Page_Up" )
		|| !strcmp( op, "KP_Prior" ) || !strcmp( op, "KP_Page_Up" ) )
	{
		if ( pop_region( MAIN_VIEW->visregion ) )
		{
			update_vis_region();

			if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
				draw_thumb_border();
		}
	}

	else if ( !strcmp( op, "POP" )
		|| !strcmp( op, "Home" ) || !strcmp( op, "KP_Home" ) )
	{
		flag = 0;

		while ( pop_region( MAIN_VIEW->visregion ) )
			flag++;

		if ( flag )
		{
			update_vis_region();

			if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
				draw_thumb_border();
		}
	}

	else if ( !strcmp( op, "Next" ) || !strcmp( op, "Page_Down" )
		|| !strcmp( op, "KP_Next" ) || !strcmp( op, "KP_Page_Down" ) )
	{
		if ( repush_region( MAIN_VIEW->visregion ) )
		{
			update_vis_region();

			if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
				draw_thumb_border();
		}
	}

	else if ( !strcmp( op, "End" ) || !strcmp( op, "KP_End" ) )
	{
		flag = 0;

		while ( repush_region( MAIN_VIEW->visregion ) )
			flag++;
		
		if ( flag )
		{
			update_vis_region();

			if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
				draw_thumb_border();
		}
	}

	else if ( !strcmp( op, "enter" ) )
	{
		PTR_VALID = TRUE;

		PTR_X = x;
		PTR_Y = y;
	}

	else if ( !strcmp( op, "move" ) )
	{
		PTR_X = x;
		PTR_Y = y;
	}

	else if ( !strcmp( op, "leave" ) )
	{
		PTR_VALID = FALSE;

		PTR_X = -1;
		PTR_Y = -1;
	}

	else if ( strlen( op ) == 1 )
	{
		switch ( op[0] )
		{
			case 'q':
			case 'Q':
			{
				Tcl_Eval( itp, "quit" );

				break;
			}

			default:
			{
				if ( vflag )
				{
					printf(
						"Error in %s: Unknown Key '%c' (%d)\n",
						"canvas_handle{}", op[0], op[0] );
				}
				
				break;
			}
		}
	}

	else if ( vflag )
		printf( "Error in canvas_handle{}: Unknown Op = \"%s\"\n", op );

	return( TCL_OK );
}


void
nudge_view( dx, dy )
int dx, dy;
{
	REGION tmpregion;

	copy_region( tmpregion, MAIN_VIEW->visregion );

	REGION_MIN( tmpregion, 0 ) += dx;
	REGION_MAX( tmpregion, 0 ) += dx;

	REGION_MIN( tmpregion, 1 ) += dy;
	REGION_MAX( tmpregion, 1 ) += dy;

	check_in_bounds( tmpregion );

	if ( !match_region( tmpregion, MAIN_VIEW->visregion,
		MAT->dim ) )
	{
		copy_region( MAIN_VIEW->visregion, tmpregion );

		update_vis_region();

		if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
			draw_thumb_border();
	}
}


/* ARGSUSED */
int thumb_nav_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	static int markx = -1;
	static int marky = -1;

	REGION tmpregion;

	char *color;
	char *op;

	int xmin, xmax;
	int ymin, ymax;

	int vx, vy;
	int x, y;

	if ( THUMB_VIEW->win == NULL )
		return( TCL_OK );

	if ( argc != 4 )
	{
		Tcl_SetResult( itp, "usage: thumb_nav <op> <x> <y>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	if ( MAT->loaded == ML_UNLOADED )
		return( TCL_OK );

	op = argv[1];

	x = atoi( argv[2] );
	y = atoi( argv[3] );

	/* Munge coords a little, so pointer's not in the way... */

	x -= 3;
	y -= 3;

	/* Handle Release Command */

	if ( !strcmp( op, "press" )
		|| !strcmp( op, "slide" )
		|| !strcmp( op, "release" ) )
	{
		vx = SCREEN_TO_REGION( x, THUMB_VIEW->win->xoff,
			REGION_MIN( THUMB_VIEW->visregion, 0 ),
			REGION_MAX( THUMB_VIEW->visregion, 0 ),
			THUMB_VIEW->nx, THUMB_VIEW->iscale );

		vy = SCREEN_TO_REGION( y, THUMB_VIEW->win->yoff,
			REGION_MIN( THUMB_VIEW->visregion, 1 ),
			REGION_MAX( THUMB_VIEW->visregion, 1 ),
			THUMB_VIEW->ny, THUMB_VIEW->iscale );

		copy_region( tmpregion, MAIN_VIEW->visregion );

		REGION_MIN( tmpregion, 0 ) = vx;
		REGION_MIN( tmpregion, 1 ) = vy;

		REGION_MAX( tmpregion, 0 ) = vx +
			( REGION_MAX( MAIN_VIEW->visregion, 0 )
				- REGION_MIN( MAIN_VIEW->visregion, 0 ) );

		REGION_MAX( tmpregion, 1 ) = vy +
			( REGION_MAX( MAIN_VIEW->visregion, 1 )
				- REGION_MIN( MAIN_VIEW->visregion, 1 ) );

		check_in_bounds( tmpregion );

		if ( !strcmp( op, "press" ) || !strcmp( op, "slide" ) )
			draw_thumb_border_region( tmpregion );

		else if ( !strcmp( op, "release" ) )
		{
			if ( !match_region( tmpregion, MAIN_VIEW->visregion,
				MAT->dim ) )
			{
				copy_region( MAIN_VIEW->visregion, tmpregion );

				update_vis_region();

				if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
					draw_thumb_border();
			}

			else
				draw_thumb_border_region( tmpregion );
		}
	}

	else if ( !strcmp( op, "region_press" ) )
	{
		if ( THUMB_BORDER == -1 )
		{
			color = ( ZERO_COLOR != NULL
					&& !strcmp( ZERO_COLOR, "black" ) )
				? "white" : "black";

			THUMB_BORDER = DRAW_RECT( itp, THUMB_VIEW->win->canvas,
				x, y, x, y, "", color, 1, "region" );
		}

		else
		{
			RECT_COORDS( itp, THUMB_VIEW->win->canvas, THUMB_BORDER,
				x, y, x, y );
		}

		markx = x;
		marky = y;
	}

	else if ( !strcmp( op, "region_slide" ) )
	{
		if ( THUMB_BORDER == -1 )
		{
			color = ( ZERO_COLOR != NULL
					&& !strcmp( ZERO_COLOR, "black" ) )
				? "white" : "black";

			THUMB_BORDER = DRAW_RECT( itp, THUMB_VIEW->win->canvas,
				markx, marky, x, y, "", color, 1, "region" );
		}

		else
		{
			RECT_COORDS( itp, THUMB_VIEW->win->canvas, THUMB_BORDER,
				markx, marky, x, y );
		}
	}

	else if ( !strcmp( op, "region_release" ) )
	{
		if ( markx != x && marky != y )
		{
			/* Save Old Region */

			push_region( MAIN_VIEW->visregion );

			/* Determine Min/Max Coords */

			if ( markx < x )
			{
				xmin = markx;
				xmax = x;
			}

			else
			{
				xmin = x;
				xmax = markx;
			}

			if ( marky < y )
			{
				ymin = marky;
				ymax = y;
			}

			else
			{
				ymin = y;
				ymax = marky;
			}

			/* Set New Region Bounds */

			xmin = SCREEN_TO_REGION( xmin, THUMB_VIEW->win->xoff,
				REGION_MIN( THUMB_VIEW->visregion, 0 ),
				REGION_MAX( THUMB_VIEW->visregion, 0 ),
				THUMB_VIEW->nx, THUMB_VIEW->iscale );

			ymin = SCREEN_TO_REGION( ymin, THUMB_VIEW->win->yoff,
				REGION_MIN( THUMB_VIEW->visregion, 1 ),
				REGION_MAX( THUMB_VIEW->visregion, 1 ),
				THUMB_VIEW->ny, THUMB_VIEW->iscale );

			xmax = SCREEN_TO_REGION( xmax, THUMB_VIEW->win->xoff,
				REGION_MIN( THUMB_VIEW->visregion, 0 ),
				REGION_MAX( THUMB_VIEW->visregion, 0 ),
				THUMB_VIEW->nx, THUMB_VIEW->iscale );

			ymax = SCREEN_TO_REGION( ymax, THUMB_VIEW->win->yoff,
				REGION_MIN( THUMB_VIEW->visregion, 1 ),
				REGION_MAX( THUMB_VIEW->visregion, 1 ),
				THUMB_VIEW->ny, THUMB_VIEW->iscale );

			REGION_MIN( MAIN_VIEW->visregion, 0 ) = xmin;
			REGION_MIN( MAIN_VIEW->visregion, 1 ) = ymin;

			REGION_MAX( MAIN_VIEW->visregion, 0 ) = xmax;
			REGION_MAX( MAIN_VIEW->visregion, 1 ) = ymax;

			check_region( MAIN_VIEW->visregion, MAIN_VIEW->viscell,
				MAT->glb, MAT->gub, MAT->dim );

			/* Send New Vis Region */

			update_vis_region();

			if ( handle_frame( MAIN_VIEW, (VIEW) NULL ) )
				draw_thumb_border();
		}

		else
			draw_thumb_border();
	}

	else if ( vflag )
		printf( "Error in thumb_nav{}: Unknown Op = \"%s\"\n", op );

	return( TCL_OK );
}


int
check_in_bounds( visregion )
REGION visregion;
{
	int changed = 0;
	int tmp;

	if ( REGION_MIN( visregion, 0 ) < MAT->glb[0] )
	{
		tmp = MAT->glb[0] - REGION_MIN( visregion, 0 );

		REGION_MIN( visregion, 0 ) += tmp;
		REGION_MAX( visregion, 0 ) += tmp;
	}

	else if ( REGION_MAX( visregion, 0 ) > MAT->gub[0] )
	{
		tmp = REGION_MAX( visregion, 0 ) - MAT->gub[0];

		REGION_MIN( visregion, 0 ) -= tmp;
		REGION_MAX( visregion, 0 ) -= tmp;
	}

	if ( REGION_MIN( visregion, 1 ) < MAT->glb[1] )
	{
		tmp = MAT->glb[1] - REGION_MIN( visregion, 1 );

		REGION_MIN( visregion, 1 ) += tmp;
		REGION_MAX( visregion, 1 ) += tmp;
	}

	else if ( REGION_MAX( visregion, 1 ) > MAT->gub[1] )
	{
		tmp = REGION_MAX( visregion, 1 ) - MAT->gub[1];

		REGION_MIN( visregion, 1 ) -= tmp;
		REGION_MAX( visregion, 1 ) -= tmp;
	}

	return( changed );
}


RECT
find_rect( V, x, y )
VIEW V;
int x;
int y;
{
	RECT R;

	int xx, yy;
	int index;

	/* Make Sure We're Loaded... */

	if ( MAT->loaded == ML_UNLOADED || V->rectids == NULL )
		return( (RECT) NULL );

	/* Screw the Linear Search, Calculate It! */

	xx = ( x - V->win->xoff ) / V->iscale;
	yy = ( y - V->win->yoff ) / V->iscale;

	if ( xx >= 0 && xx < V->nx && yy >= 0 && yy < V->ny )
	{
		index = ( yy * V->nx ) + xx;

		/* One Last Sanity Check */

		if ( index < V->nrectids )
		{
			R = &(V->rectids[index]);

			return( R );
		}
	}

	return( (RECT) NULL );
}


/* ARGSUSED */
int do_print_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char fname[256];
	char usage[256];
	char cmd[512];

	char *format;
	char *orient;
	char *media;
	char *dest;
	char *file;
	char *pcmd;
	char *pp;

	int seq;

	if ( argc < 5 )
	{
		snprintf( usage, 255,
			"usage: do_print %s %s <orient> <media> <args...>",
			"[ print | preview | preview_print | cancel ]",
			"[ printer | file ]" );

		Tcl_SetResult( itp, usage, TCL_STATIC );

		return( TCL_ERROR );
	}

	pp = argv[1];
	dest = argv[2];
	orient = argv[3];
	media = argv[4];

	/* Check for Cancel... */

	if ( !strcmp( pp, "cancel" ) )
	{
		(PRINT_VIEW->sequence)++;

		return( TCL_OK );
	}

	/* Get Starting Instance Number */

	seq = PRINT_VIEW->sequence;

	/* Verify Canvas Image ID */

	if ( PRINT_VIEW->image_id == -1 )
	{
		const char* image_id = GET_TCL_GLOBAL( interp, "pv_image_id" );
		PRINT_VIEW->image_id = atoi( image_id );
	}

	/* Execute Print / Save Function */

	if ( !strcmp( dest, "printer" ) )
	{
		if ( argc != 6 )
		{
			Tcl_SetResult( itp,
		"usage: do_print <cmd> printer <orient> <media> <print_cmd>",
				TCL_STATIC );

			return( TCL_ERROR );
		}

		pcmd = argv[5];

		/* Draw Image on Print View */

		if ( strcmp( pp, "preview_print" ) )
		{
			prepare_print( orient, media, "ps" );

			INST_CHECK_ALT_RETURN( seq, PRINT_VIEW->sequence,
				( Tcl_Eval( itp, "setMsg {Print Cancelled.}" ),
					update_progress( PRINT_VIEW, 0 ) ), TCL_OK );
		}

		/* Just Return if Only Preview */

		if ( !strcmp( pp, "preview" ) )
			return( TCL_OK );

		/* Send to Printer */

		snprintf( fname, 255, ".mv_print.%d.ps", getpid() );

		save_image_to_file( fname, "ps" );

		INST_CHECK_ALT_RETURN( seq, PRINT_VIEW->sequence,
			( Tcl_Eval( itp, "setMsg {Print Cancelled.}" ),
				update_progress( PRINT_VIEW, 0 ) ), TCL_OK );

		snprintf( cmd, 500, "setMsg {Executing Print Command \"%s\"...}",
			pcmd );
		Tcl_Eval( itp, cmd );

		INST_CHECK_ALT_RETURN( seq, PRINT_VIEW->sequence,
			( Tcl_Eval( itp, "setMsg {Print Cancelled.}" ),
				update_progress( PRINT_VIEW, 0 ) ), TCL_OK );

		snprintf( cmd, 500, "%s %s", pcmd, fname );

		system( cmd );

		sleep( 3 );

		unlink( fname );

		Tcl_Eval( itp, "setMsg {Print Command Completed.}" );
	}

	else if ( !strcmp( dest, "file" ) )
	{
		if ( argc != 7 )
		{
			Tcl_SetResult( itp,
	"usage: do_print <cmd> file <orient> <media> <file_name> <format>",
				TCL_STATIC );

			return( TCL_ERROR );
		}

		file = argv[5];
		format = argv[6];

		/* Draw Image on Print View */

		if ( strcmp( pp, "preview_print" ) )
		{
			prepare_print( orient, media, format );

			INST_CHECK_ALT_RETURN( seq, PRINT_VIEW->sequence,
				( Tcl_Eval( itp, "setMsg {Print Cancelled.}" ),
					update_progress( PRINT_VIEW, 0 ) ), TCL_OK );
		}

		/* Just Return if Only Preview */

		if ( !strcmp( pp, "preview" ) )
			return( TCL_OK );

		/* Save Image to File */

		save_image_to_file( file, format );
	}

	return( TCL_OK );
}


void
prepare_print( char* orient, char* media, char* format ) {
	char cmd[256];

	int save_grey;
	int wt, ht;
	int tmp;

	/* Calculate Preview Width & Height (Portrait) */

	wt = ( MAIN_VIEW->iscale * MAIN_VIEW->nx )
		+ ( 2 * MAIN_VIEW->win->xoff );

	ht = ( MAIN_VIEW->iscale * MAIN_VIEW->ny )
		+ ( 2 * MAIN_VIEW->win->yoff );

	/* Handle Landscape Mode */

	if ( !strcmp( orient, "landscape" ) )
	{
		PRINT_VIEW->draw_type = DRAW_LANDSCAPE;

		tmp = wt;
		wt = ht;
		ht = tmp;
	}

	else
		PRINT_VIEW->draw_type = DRAW_PORTRAIT;

	/* Handle Media Type */

	save_grey = GREYSCALE;

	if ( !strcmp( media, "greyscale" ) )
		GREYSCALE = TRUE;

	else
		GREYSCALE = FALSE;

	/* Copy Current Main View to Print View */

	copy_region( PRINT_VIEW->visregion, MAIN_VIEW->visregion );

	copy_cell( PRINT_VIEW->viscell, MAIN_VIEW->viscell );

	snprintf( cmd, 255, "%d", ht );
	SET_TCL_GLOBAL( interp, "pview_field_height", cmd );

	snprintf( cmd, 255, "%d", wt );
	SET_TCL_GLOBAL( interp, "pview_field_width", cmd );

	Tcl_Eval( interp, "layout_pv_panel" );

	/* Set Draw Type */

	if ( !strcmp( format, "ps" ) )
		PRINT_VIEW->draw_type |= DRAW_CANVAS;
	
	else if ( !strcmp( format, "ppm" ) )
		PRINT_VIEW->draw_type |= DRAW_IMAGE;

	/* Draw Print View */

	handle_frame( PRINT_VIEW, (VIEW) NULL );

	/* Reset Greyscale */

	GREYSCALE = save_grey;
}


/* ARGSUSED */
int interface_param_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char tmp[256];

	char *value;
	char *what;

	if ( argc != 3 )
	{
		Tcl_SetResult( itp,
			"usage: interfaceParamHandle <what> <value>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	what	= argv[1];
	value	= argv[2];

	if ( !strcmp( what, "progress_ticks" ) )
		PROGRESS_TICKS = atoi( value );

	else
	{
		snprintf( tmp, 255,
			"Error in interfaceParamHandle{}:  unknown param \"%s\".",
			what );

		Tcl_SetResult( itp, copy_str( tmp ), TCL_DYNAMIC );

		return( TCL_ERROR );
	}

	return( TCL_OK );
}


/* ARGSUSED */
int register_menuvar_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *menuvar;
	char *index;
	char *menu;

	if ( argc != 4 )
	{
		Tcl_SetResult( itp,
			"usage: register_menuvar <menu> <index> <menuvar>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	menu	= argv[1];
	index	= argv[2];
	menuvar	= argv[3];

	if ( !strcmp( menu, ".matview.options_menu" ) )
	{
		if ( !strcmp( index, "auto_color" ) )
		{
			if ( MENUVAR_AUTO_COLOR != NULL )
				FREE( MENUVAR_AUTO_COLOR );

			MENUVAR_AUTO_COLOR = copy_str( menuvar );
		}

		else if ( !strcmp( index, "greyscale" ) )
		{
			if ( MENUVAR_GREYSCALE != NULL )
				FREE( MENUVAR_GREYSCALE );

			MENUVAR_GREYSCALE = copy_str( menuvar );
		}

		else if ( !strcmp( index, "filter_values" ) )
		{
			if ( MENUVAR_FILTER_VALUES != NULL )
				FREE( MENUVAR_FILTER_VALUES );

			MENUVAR_FILTER_VALUES = copy_str( menuvar );
		}

		else if ( !strcmp( index, "min_max" ) )
		{
			if ( MENUVAR_MIN_MAX != NULL )
				FREE( MENUVAR_MIN_MAX );

			MENUVAR_MIN_MAX = copy_str( menuvar );
		}

		else if ( !strcmp( index, "cumulative" ) )
		{
			if ( MENUVAR_CUMULATIVE != NULL )
				FREE( MENUVAR_CUMULATIVE );

			MENUVAR_CUMULATIVE = copy_str( menuvar );
		}

		else if ( !strcmp( index, "absvals" ) )
		{
			if ( MENUVAR_ABSVALS != NULL )
				FREE( MENUVAR_ABSVALS );

			MENUVAR_ABSVALS = copy_str( menuvar );
		}

		else if ( !strcmp( index, "pattern" ) )
		{
			if ( MENUVAR_PATTERN != NULL )
				FREE( MENUVAR_PATTERN );

			MENUVAR_PATTERN = copy_str( menuvar );
		}

		else if ( !strcmp( index, "transpose" ) )
		{
			if ( MENUVAR_TRANSPOSE != NULL )
				FREE( MENUVAR_TRANSPOSE );

			MENUVAR_TRANSPOSE = copy_str( menuvar );
		}
	}

	else if ( !strcmp( menu, ".matview.file_menu" ) )
	{
		if ( !strcmp( index, "noload_values" ) )
		{
			if ( MENUVAR_NOLOADVALUES != NULL )
				FREE( MENUVAR_NOLOADVALUES );

			MENUVAR_NOLOADVALUES = copy_str( menuvar );
		}
	}

	return( TCL_OK );
}


/* ARGSUSED */
int load_bitmap_file_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	Pixmap bitmap;

	char tmp[1024];

	if ( argc != 2 )
	{
		Tcl_SetResult( itp, "usage: load_bitmap_file <file>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	snprintf( tmp, 1023, "@%s", argv[1] );

	bitmap = Tk_GetBitmap( itp, Top, Tk_GetUid( tmp ) );

	if ( bitmap <= 0 )
	{
		snprintf( tmp, 1023, "Error Opening Bitmap %s", argv[1] );

		Tcl_SetResult( itp, copy_str( tmp ), TCL_DYNAMIC );

		return( TCL_ERROR );
	}

	return( TCL_OK );
}


/* ARGSUSED */
int fix_help_line_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char tmp[1024];

	char *result;
	char *ptr;
	char *str;

	if ( argc != 2 )
	{
		Tcl_SetResult( itp, "usage: fix_help_line <line>", TCL_STATIC );

		return( TCL_ERROR );
	}

	str = argv[1];

	ptr = tmp;

	while ( *str != '\0' )
	{
		if ( *str == '' )
			*ptr++ = '{';

		else if ( *str == '' )
			*ptr++ = '}';

		else
			*ptr++ = *str;

		str++;
	}

	*ptr = '\0';

	result = copy_str( tmp );

	Tcl_SetResult( itp, result, TCL_DYNAMIC );

	return( TCL_OK );
}


/* ARGSUSED */
int double_click_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *id;

	if ( argc != 2 )
	{
		Tcl_SetResult( itp, "usage: double_click <id>", TCL_STATIC );

		return( TCL_ERROR );
	}

	id = argv[1];

	if ( double_click( id ) )
		Tcl_SetResult( itp, "1", TCL_STATIC );

	else
		Tcl_SetResult( itp, "0", TCL_STATIC );

	return( TCL_OK );
}


int
double_click( id )
char *id;
{
	static struct timeval last = { -1, -1 };
	static char *last_id = (char *) NULL;

	struct timeval tmp;

	float diff;

	int dc;

	gettimeofday( &tmp, (struct timezone *) NULL );

	if ( last_id == NULL || strcmp( id, last_id ) )
	{
		if ( last_id != NULL )
			FREE( last_id );

		last_id = copy_str( id );

		dc = FALSE;
	}

	else
	{
		diff = ( (float) ( tmp.tv_sec - last.tv_sec ) * 1000000.0 )
			+ (float) ( tmp.tv_usec - last.tv_usec );

		if ( diff < 500000.0 )
			dc = TRUE;

		else
			dc = FALSE;
	}

	last.tv_sec = tmp.tv_sec;
	last.tv_usec = tmp.tv_usec;

	return( dc );
}


/* ARGSUSED */
int strip_label_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char tmp[1024];

	char *result;
	char *ptr;
	char *str;

	if ( argc != 2 )
	{
		Tcl_SetResult( itp, "usage: strip_label <str>", TCL_STATIC );

		return( TCL_ERROR );
	}

	str = argv[1];

	ptr = tmp;

	while ( *str != '\0' )
	{
		if ( *str == '-' || *str == ' ' || *str == '.' )
			*ptr++ = '_';

		else
			*ptr++ = *str;

		str++;
	}

	*ptr = '\0';

	result = copy_str( tmp );

	Tcl_SetResult( itp, result, TCL_DYNAMIC );

	return( TCL_OK );
}


/* ARGSUSED */
int title_info_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char tmp[256];

	snprintf( tmp, 255, "%s", MATVIEW_VERSION );

	Tcl_SetResult( itp, copy_str( tmp ), TCL_DYNAMIC );

	return( TCL_OK );
}


/* ARGSUSED */
int pad_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char *result;
	char *ptr;
	char *str;

	int len;
	int nd;
	int nz;
	int i;

	if ( argc < 3 )
	{
		Tcl_SetResult( itp,
			"usage: pad <number_str> <desired_length>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	str = argv[1];

	len = atoi( argv[2] );

	nd = strlen( str );

	nz = len - nd;

	result = (char *) MALLOC( (unsigned) (len + 1) * sizeof(char) );
	memcheck( result, "Numerical Pad String" );

	ptr = result;

	for ( i=0 ; i < nz ; i++ )
		*ptr++ = '0';

	snprintf( (char *) (result + nz), len-nz, "%s", str );

	Tcl_SetResult( itp, result, TCL_DYNAMIC );

	return( TCL_OK );
}


static char *about_text[] =
{
	"                     MatView Version 1.0:",
	"                Simple Scalable Matrix Viewer",
	"         Oak Ridge National Laboratory, Oak Ridge TN.",
	"                  Author:  James Arthur Kohl",
	"                 (C) 1998 All Rights Reserved",
	"",
	"                             NOTICE",
	"",
	"Permission to use, copy, modify, and distribute this software",
	"and its documentation for any purpose and without fee is hereby",
	"granted provided that the above copyright notice appear in all",
	"copies and that both the copyright notice and this permission",
	"notice appear in supporting documentation.",
	"",
	"Neither the Institution, Oak Ridge National Laboratory, nor the",
	"Authors make any representations about the suitability of this",
	"software for any purpose.  This software is provided ``as is''",
	"without express or implied warranty.",
 	"",
	"MatView was funded by the U.S. Department of Energy.",
	"",
	"Acknowledgements:",
	"",
	"Many Thanks to Tammy Kolda, Chuck Romine and Noel",
	"Nachtigal for their suggestions in helping to develop",
	"this tool!",
	"",
	"Any questions, problems, or suggestions regarding",
	"MatView should be directed to:",
	"",
	"    kohlja@ornl.gov",
	"",
	"    James Arthur Kohl, Ph.D.",
	"    P.O. Box 2008, Bldg 6012, MS 6367",
	"    Oak Ridge National Laboratory",
	"    Oak Ridge, TN 37831-6367",
	"",
	"I usually provide timely responses unless I'm on travel.",
	"",
	"Thanks much, and Good Luck!  :-)"
};

static int about_text_size
	= sizeof( about_text ) / sizeof( about_text[0] );

/* ARGSUSED */
int fill_about_text_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	char tmp[256];

	char *result[2];

	char *canvas;
	char *color;
	char *font;
	char *dims;

	int len;
	int wt;
	int y;
	int i;

	if ( argc != 4 )
	{
		Tcl_SetResult( itp,
			"usage: fill_about_text <canvas> <color> <font>",
			TCL_STATIC );

		return( TCL_ERROR );
	}

	canvas = argv[1];
	color = argv[2];
	font = argv[3];

	y = BORDER_SPACE;

	wt = -1;

	for ( i=0 ; i < about_text_size ; i++ )
	{
		DRAW_TEXT( itp, canvas, BORDER_SPACE, y, about_text[i],
			color, font );

		len = strlen( about_text[i] );

		wt = ( len > wt ) ? len : wt;

		y += ROW_HEIGHT;
	}

	/* Pass Back Text Height & Width */

	wt *= COL_WIDTH;

	snprintf( tmp, 255, "%d", y );
	result[0] = copy_str( tmp );

	snprintf( tmp, 255, "%d", wt );
	result[1] = copy_str( tmp );

	/* XXX gcc says result is an incompatible pointer */
	dims = Tcl_Merge( 2, (const char* const*)result );

	Tcl_SetResult( itp, dims, TCL_DYNAMIC );

	FREE( result[0] );
	FREE( result[1] );

	return( TCL_OK );
}


/* ARGSUSED */
void
quit_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	// printf( "Quitting MatView.\n" );

	exit( 0 );
}

