
/* $Id: view.c,v 1.114 2005/03/02 15:14:30 kohl Exp $ */


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


#include "matview.h"
#include "matview_glob.h"


void
do_load_mat( clientData )
ClientData clientData;
{
	long lreload;
	int reload;

	lreload = (long) clientData;
	reload = (int) lreload;

	load_mat( reload );
}


int
load_mat( reload )
int reload;
{
	char tmp[1024];

	int i;

	if ( MATVIEW_MATFILE != NULL )
	{
		/* Bump Frame Instance Number */

		(MAIN_VIEW->instance)++;
		(THUMB_VIEW->instance)++;

		/* Read Matrix From File */

		MAT->loaded = ML_UNLOADED;

		if ( !read_mat_file( MATVIEW_MATFILE ) )
		{
			clear_view( MAIN_VIEW );
			clear_view( THUMB_VIEW );

			return( FALSE );
		}

		if ( !reload )
		{
			while ( pop_region( MAIN_VIEW->visregion ) );

			/* Reset Plane Choice */

			SET_TCL_GLOBAL( interp, "plane_choice", "XY" );

			Tcl_Eval( interp, "layout_main_panel" );

			Tcl_Eval( interp, "set_plane_choice XY" );
		}

		/* Set Max Region */

		for ( i=0 ; i < MAX_DIM ; i++ )
		{
			if ( i < MAT->dim )
			{
				REGION_MIN( MAX_REGION, i ) = MAT->glb[i];
				REGION_MAX( MAX_REGION, i ) = MAT->gub[i];
			}

			else
			{
				REGION_MIN( MAX_REGION, i ) =
					REGION_MAX( MAX_REGION, i ) =
						0;
			}
		}

		/* Set Initial Vis Region */

		for ( i=0 ; i < MAT->dim && i < 2 ; i++ )
		{
			REGION_MIN( THUMB_VIEW->visregion, i ) = MAT->glb[i];
			REGION_MAX( THUMB_VIEW->visregion, i ) = MAT->gub[i];

			CELL_OF_DIM( THUMB_VIEW->viscell, i ) = 1;

			if ( !reload )
			{
				REGION_MIN( MAIN_VIEW->visregion, i ) = MAT->glb[i];
				REGION_MAX( MAIN_VIEW->visregion, i ) = MAT->gub[i];

				CELL_OF_DIM( MAIN_VIEW->viscell, i ) = 1;
			}
		}

		for ( i=((MAT->dim < 2) ? MAT->dim : 2) ; i < MAX_DIM ; i++ )
		{
			if ( i < MAT->dim )
			{
				REGION_MIN( THUMB_VIEW->visregion, i ) =
					REGION_MAX( THUMB_VIEW->visregion, i ) =
						MAT->glb[i];

				if ( !reload )
				{
					REGION_MIN( MAIN_VIEW->visregion, i ) =
						REGION_MAX( MAIN_VIEW->visregion, i ) =
							MAT->glb[i];
				}
			}

			else
			{
				REGION_MIN( THUMB_VIEW->visregion, i ) =
					REGION_MAX( THUMB_VIEW->visregion, i ) =
						0;

				if ( !reload )
				{
					REGION_MIN( MAIN_VIEW->visregion, i ) =
						REGION_MAX( MAIN_VIEW->visregion, i ) =
							0;
				}
			}

			CELL_OF_DIM( THUMB_VIEW->viscell, i ) = 1;

			if ( !reload )
				CELL_OF_DIM( MAIN_VIEW->viscell, i ) = 1;
		}

		/* Make Sure Things Are O.K. After the Reload... */

		if ( reload )
		{
			check_region( MAIN_VIEW->visregion, MAIN_VIEW->viscell,
				MAT->glb, MAT->gub, MAT->dim );
		}

		if ( MAT->alloc_size == MAT->nelems )
			MAT->loaded = ML_LOADED;
		else
			MAT->loaded = ML_PARTIAL;

		update_vis_region();

		CUMULATIVE_SET = 0;

		if ( MAT->complex )
			Tcl_Eval( interp, "set_complex ON" );
		else
			Tcl_Eval( interp, "set_complex OFF" );

		handle_frame( MAIN_VIEW, THUMB_VIEW );

		sprintf( tmp,
			"setMsg {Matrix File \"%s\" Successfully Loaded.}",
			MATVIEW_MATFILE );
		Tcl_Eval( interp, tmp );
	}

	return( TRUE );
}


int
handle_frame( V, thumb )
VIEW V;
VIEW thumb;
{
	RECT R, R2;

	RECT newrects;

	char color[16];

	const char *tmp;

	double *dptr;
	double *ptr;
	int nelems;

	double min, max;
	double xscale;
	double yscale;
	double fval;

	int max_indices[MAX_DIM];
	int indices[MAX_DIM];

	int x1, y1, x2, y2;
	int xi, yi, zi;
	int progress;
	int vis_mode;
	int newsize;
	int nrects;
	int nx, ny;
	int height;
	int width;
	int inst;
	int i;

	/* Check to See if a Matrix is Loaded */

	if ( MAT->loaded == ML_UNLOADED )
		return( FALSE );

	/* Set View Redraw Flags */

	V->need_redraw = TRUE;

	if ( thumb )
		thumb->need_redraw = TRUE;

	/* Get Proper Plane Indices */

	if ( !get_plane_indices( &xi, &yi, &zi ) )
		return( FALSE );

	/* Calc Total Num of Rects Needed */

	nx = ( ( REGION_MAX( V->visregion, xi )
				- REGION_MIN( V->visregion, xi ) + 1 )
			+ ( CELL_OF_DIM( V->viscell, xi ) - 1 ) )
		/ CELL_OF_DIM( V->viscell, xi );

	ny = ( ( REGION_MAX( V->visregion, yi )
				- REGION_MIN( V->visregion, yi ) + 1 )
			+ ( CELL_OF_DIM( V->viscell, yi ) - 1 ) )
		/ CELL_OF_DIM( V->viscell, yi );

	/* Get Canvas Dimensions */

	tmp = GET_TCL_GLOBAL( interp, V->win->htvar );
	height = atoi( tmp );
	height -= 2 * V->win->yoff;

	tmp = GET_TCL_GLOBAL( interp, V->win->wtvar );
	width = atoi( tmp );
	width -= 2 * V->win->xoff;

	/* Swap for Landscape Scaling... */

	if ( V->draw_type & DRAW_LANDSCAPE )
	{
		i = height;
		height = width;
		width = i;
	}

	/* Determine Scale Factor */

	xscale = (double) width / (double) nx;
	yscale = (double) height / (double) ny;

	V->scale = ( xscale < yscale ) ? xscale : yscale;

	if ( V->scale < 1.0 )
	{
		nx = (int) ceil( V->scale * (double) nx );
		ny = (int) ceil( V->scale * (double) ny );
	}

	nrects = nx * ny;

	/* Determine Integer Scale for Drawing... */

	V->iscale = (int) floor( V->scale );
	V->iscale = ( V->iscale ) ? V->iscale : 1;

	/* Set Reduce Status Indicator */

	handle_reduce_status( V );

	/* Allocate Rect Array */

	if ( V->rectids == NULL || V->rectidsize < nrects )
	{
		newsize = V->rectidsize ? V->rectidsize : 1;

		while ( newsize < nrects )
			newsize *= 2;

		newrects = (RECT) MALLOC( (unsigned) newsize
	 		* sizeof( struct rect_struct ) );
		memcheck( newrects, "Rectangle IDs List" );

		R = &(newrects[0]);

		if ( V->rectids != NULL )
		{
			R2 = &(V->rectids[0]);

			for ( i=0 ; i < V->rectidsize ; i++ )
			{
				R->x1 = R2->x1;
				R->y1 = R2->y1;
				R->x2 = R2->x2;
				R->y2 = R2->y2;

				R->drawn = R2->drawn;

				R++;  R2++;
			}

			for ( i=(V->rectidsize) ; i < newsize ; i++ )
			{
				R->drawn = -1;

				R++;
			}

			FREE( V->rectids );
		}

		else
		{
			for ( i=0 ; i < newsize ; i++ )
			{
				R->drawn = -1;

				R++;
			}
		}

		V->rectids = newrects;

		V->rectidsize = newsize;
	}

	/* Reset Drawn Flags for (What Will Be) Unused Rects */

	else if ( nrects < V->nrectids )
	{
		R = &(V->rectids[nrects]);

		for ( i=nrects ; i < V->nrectids ; R++, i++ )
			R->drawn = -1;
	}

	V->nrectids = nrects;

	/* Grab & Increment Frame Instance Number */

	inst = ++(V->instance);

	/* Generate Viewing Field */

	if ( !generate_view_field( V, nx, ny, xi, yi ) )
		return( FALSE );
	INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

	/* Draw It */

	switch ( VIS_MODE )
	{
		case VIS_MODE_NORMAL:
			tmp = "setMsg {Drawing Data Frame...}";
			break;

		case VIS_MODE_DIFF:
			tmp = "setMsg {Drawing Differential Data Frame...}";
			break;

		default:
			tmp = "setMsg {Drawing Data Frame...}";
			break;
	}

	vis_mode = ( VIS_MODE == VIS_MODE_DIFF && !(V->do_diff_mode) )
		? VIS_MODE_NORMAL : VIS_MODE;

	Tcl_Eval( interp, tmp );
	INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

	update_progress( V, 0 );
	INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

	/* Automatically Determine Color Map Range */

	if ( V->do_auto_color &&
		( AUTO_COLOR_MAP || SHOW_MIN_MAX || SHOW_CUMULATIVE ) )
	{
		compute_view_min_max( V, &min, &max );

		handle_min_max( min, max );
	}

	/* Initialize Local & Max Index */

	for ( i=0 ; i < MAT->dim ; i++ )
		max_indices[i] = indices[i] = 0;

	max_indices[xi] = V->nx;
	max_indices[yi] = V->ny;

	/* Set Data Pointers */

	ptr = V->field;

	switch ( vis_mode )
	{
		case VIS_MODE_NORMAL:
			dptr = (double *) NULL;
			break;

		case VIS_MODE_DIFF:
			dptr = SAVE_FIELD;
			nelems = SAVE_FIELD_NELEMS - 1;
			break;
	}

	/* Determine Progress Factor */

	progress = V->nrectids / PROGRESS_TICKS;
	progress = ( progress ) ? progress : 1;

	/* Hide Viewing Image */

	hide_image( V );
	INST_CHECK_RETURN( V, inst, V->instance, TRUE, FALSE );

	/* Clear View */

	clear_view( V );

	/* Draw Zero Rect Background */

	if ( ZERO_COLOR != NULL )
		strcpy( color, ZERO_COLOR );

	else
		color_value( color, 0.0 );

	x1 = V->win->xoff;
	y1 = V->win->yoff;

	x2 = ( V->iscale * V->nx ) + V->win->xoff;
	y2 = ( V->iscale * V->ny ) + V->win->yoff;

	if ( V->draw_type & DRAW_LANDSCAPE )
	{
		landscape_coords( &x1, &y1, &x2, &y2,
			( V->iscale * V->ny ) + ( 2 * V->win->yoff ) );
	}

	if ( V->draw_type & DRAW_CANVAS )
	{
		DRAW_RECT( interp, V->win->canvas,
			x1 - 1, y1 - 1, x2, y2, color, "black", 1,
			"view" );
	}

	else
	{
		DRAW_IMAGE_RECT( interp, V->image,
			x1, y1, x2, y2, color );

		DRAW_IMAGE_OUTLINE( interp, V->image,
			x1 - 1, y1 - 1, x2, y2, "black" );
	}

	/* Draw Rects */

	R = &(V->rectids[0]);

	i = 0;

	do
	{
		/* Update Progress Meter */

		if ( !( i % progress ) )
		{
			update_progress( V, i / progress );
			INST_CHECK_RETURN( V, inst, V->instance, TRUE, FALSE );
		}

		/* Get Data Value */

		if ( dptr != NULL )
		{
			fval = (*ptr) - (*dptr++);

			if ( i >= nelems )
				dptr = (double *) NULL;
		}

		else
			fval = *ptr;

		R->value = fval;

		/* Get Coords */

		R->x1 = ( V->iscale * indices[xi] ) + V->win->xoff;
		R->y1 = ( V->iscale * indices[yi] ) + V->win->yoff;

		/* Note: *not* - 1, make sure rects show up in PostScript! */
		/* -> appears to be identical to - 1 in PPM format (whoa...) */
		R->x2 = R->x1 + V->iscale;
		R->y2 = R->y1 + V->iscale;

		/* Update Rectangle Image */

		if ( R->value != 0.0
			&& ( !FILTER_VALUES
				|| ( R->value >= COLOR_MIN_VALUE
					&& R->value <= COLOR_MAX_VALUE ) ) )
		{
			color_value( color, fval );

			if ( V->draw_type & DRAW_LANDSCAPE )
			{
				landscape_coords(
					&(R->x1), &(R->y1), &(R->x2), &(R->y2),
					( V->iscale * V->ny ) + ( 2 * V->win->yoff ) );
			}

			if ( V->draw_type & DRAW_CANVAS )
			{
				DRAW_RECT( interp, V->win->canvas,
					R->x1, R->y1, R->x2, R->y2, color, "", 1,
					"view" );
			}

			else
			{
				if ( R->x1 != R->x2 && R->y1 != R->y2 )
				{
					DRAW_IMAGE_RECT( interp, V->image,
						R->x1, R->y1, R->x2, R->y2, color );
				}

				else
				{
					DRAW_IMAGE_POINT( interp, V->image,
						R->x1, R->y1, color );
				}
			}
		}

		R->drawn = V->instance;

		/* Increment Ptrs */

		R++;  i++;

		ptr++;
	}
	while ( i < V->nrectids
		&& incr_field_index( MAT->dim, indices, max_indices,
			ROW_MAJOR ) );

	/* Make Sure All Old Rects Reset */

	while ( i < V->nrectids )
	{
		R->drawn = -1;

		R++;  i++;
	}

	/* Check On Any Grid... */

	handle_grid( V );

	/* Show Viewing Image */

	show_image( V );
	INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

	/* Reset Progress Bar */

	update_progress( V, 0 );
	INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

	/* Done */

	switch ( VIS_MODE )
	{
		case VIS_MODE_NORMAL:
			tmp = "setMsg {Data Frame Done.}";
			break;

		case VIS_MODE_DIFF:
			tmp = "setMsg {Differential Data Frame Done.}";
			break;

		default:
			tmp = "setMsg {Data Frame Done.}";
			break;
	}

	Tcl_Eval( interp, tmp );
	INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

	/* View Has Been Completely Drawn. */

	V->need_redraw = FALSE;

	/* Handle Additional Thumb View If Present... */

	if ( thumb != NULL )
	{
		if ( !handle_frame( thumb, (VIEW) NULL ) )
			return( FALSE );

		draw_thumb_border();
	}

	/* Check for Pre-empted Views */

	check_need_redraw();

	return( TRUE );
}


void
check_need_redraw()
{
	if ( THUMB_VIEW->need_redraw )
		handle_frame( THUMB_VIEW, (VIEW) NULL );

	else if ( MAIN_VIEW->need_redraw )
		handle_frame( MAIN_VIEW, (VIEW) NULL );

	else if ( PRINT_VIEW->need_redraw )
		handle_frame( PRINT_VIEW, (VIEW) NULL );
}


void
clear_view( V )
VIEW V;
{
	char cmd[255];

	/* Clear the Heck Out of the Viewing Area */

	Tcl_Eval( interp, "set_cursor busy" );

	sprintf( cmd, "%s delete view", V->win->canvas );
	Tcl_Eval( interp, cmd );

	sprintf( cmd, "clear_image %s %d %d", V->image, V->wt, V->ht );
	Tcl_Eval( interp, cmd );

	Tcl_Eval( interp, "set_cursor reset" );
}


void
handle_reduce_status( V )
VIEW V;
{
	char cmd[255];

	char *rf;

	int i;

	/* Determine if Reduction is Needed */

	if ( V->scale < 1.0 )
		V->need_reduction = TRUE;

	else
	{
		V->need_reduction = FALSE;

		for ( i=0 ; i < MAT->dim && !(V->need_reduction) ; i++ )
		{
			if ( CELL_OF_DIM( V->viscell, i ) != 1 )
				V->need_reduction = TRUE;
		}
	}

	/* Set Interface Status */

	if ( V == MAIN_VIEW )
	{
		if ( V->need_reduction )
		{
			switch ( REDUCTION_FUNCTION )
			{
				case RF_AVG:		rf = "Avg";			break;
				case RF_MIN:		rf = "Min";			break;
				case RF_MAX:		rf = "Max";			break;
				case RF_SUM:		rf = "Sum";			break;
				case RF_MIDPOINT:	rf = "Midpoint";	break;
				case RF_COUNT:		rf = "Count";		break;
				case RF_STD:		rf = "Std";			break;
				case RF_VAR:		rf = "Var";			break;

				/* defer to avg... */
				default:		rf = "Avg";		break;
			}
		}

		else
			rf = "none";

		sprintf( cmd, "reduceStatus %s", rf );

		Tcl_Eval( interp, cmd );
	}
}


int
generate_view_field( V, nx, ny, xi, yi )
VIEW V;
int nx, ny;
int xi, yi;
{
	TEMPLATE T;

	COMPLEX C;

	DATA D;

	char cmd[255];

	double *F1;
	double *F2;
	double *F;

	double eff_scale;
	double xfactor;
	double yfactor;
	double val;

	int x1, y1, x2, y2;
	int multi_pass;
	int progress;
	int symm_cnt;
	int use_tmps;
	int xx, yy;
	int index;
	int x, y;
	int i, j;
	int inst;
	int tmp;

	/* Make Sure Any Partial Matrix is Reset to Beginning... */

	if ( MAT->loaded == ML_PARTIAL && MAT->first_loaded != 0 )
	{
		if ( !read_mat_file( MATVIEW_MATFILE ) )
			return( FALSE );
	}

	inst = V->instance;

	Tcl_Eval( interp, "setMsg {Generating View Field...}" );
	INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

	update_progress( V, 0 );
	INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

	/* Determine Stat Requirements */

	multi_pass = 0;

	if ( V->need_reduction
		&& REDUCTION_FUNCTION != RF_MIN
		&& REDUCTION_FUNCTION != RF_MAX
		&& REDUCTION_FUNCTION != RF_SUM )
	{
		use_tmps = 1;

		if ( REDUCTION_FUNCTION == RF_STD
				|| REDUCTION_FUNCTION == RF_VAR )
			multi_pass++;
	}

	else
	{
		if ( V->ftmp1 )
		{
			FREE( V->ftmp1 );
			V->ftmp1 = (double *) NULL;
		}

		if ( V->ftmp2 )
		{
			FREE( V->ftmp2 );
			V->ftmp2 = (double *) NULL;
		}

		use_tmps = 0;
	}

	/* (Re-)Allocate Viewing Field */

	if ( nx != V->nx || ny != V->ny
		|| REDUCTION_FUNCTION != V->last_rf )
	{
		V->last_rf = REDUCTION_FUNCTION;

		if ( V->field != NULL )
			FREE( V->field );

		V->field = (double *) MALLOC( (unsigned) (nx * ny)
			* sizeof(double) );
		memcheck( V->field, "Viewing Field Array" );

		if ( use_tmps )
		{
			if ( V->ftmp1 != NULL )
				FREE( V->ftmp1 );

			V->ftmp1 = (double *) MALLOC( (unsigned) (nx * ny)
				* sizeof(double) );
			memcheck( V->ftmp1, "Viewing Field Array" );

			if ( V->ftmp2 != NULL )
				FREE( V->ftmp2 );

			V->ftmp2 = (double *) MALLOC( (unsigned) (nx * ny)
				* sizeof(double) );
			memcheck( V->ftmp2, "Viewing Field Array" );
		}

		V->nx = nx;
		V->ny = ny;

		V->field_nelems = V->nx * V->ny;

		if ( V->do_auto_color )
		{
			sprintf( cmd, "set_view_size %s %d %d %d %d",
				V->win->toplevel,
				REGION_MAX( V->visregion, yi )
					- REGION_MIN( V->visregion, yi ) + 1,
				REGION_MAX( V->visregion, xi )
					- REGION_MIN( V->visregion, xi ) + 1,
				V->ny, V->nx );
			Tcl_Eval( interp, cmd );
		}
	}

	/* Clear Viewing Field */

#ifdef NO_BZERO

	F = V->field;

	if ( use_tmps )
	{
		F1 = V->ftmp1;
		F2 = V->ftmp2;

		for ( i=0 ; i < V->field_nelems ; i++ )
		{
			*F1++ = 0.0;
			*F2++ = 0.0;
			*F++ = 0.0;
		}
	}

	else
	{
		for ( i=0 ; i < V->field_nelems ; i++ )
			*F++ = 0.0;
	}

#else

	if ( use_tmps )
	{
		bzero( (void *) V->ftmp1,
			(size_t) V->field_nelems * sizeof( double ) );
		bzero( (void *) V->ftmp2,
			(size_t) V->field_nelems * sizeof( double ) );
		bzero( (void *) V->field,
			(size_t) V->field_nelems * sizeof( double ) );
	}

	else
	{
		bzero( (void *) V->field,
			(size_t) V->field_nelems * sizeof( double ) );
	}

#endif

	/* Calculate Progress Factor */

	if ( multi_pass )
		tmp = ( 2 * MAT->nelems ) + V->field_nelems;

	else
		tmp = MAT->nelems;

	if ( VIEW_ABS_VALS )
		tmp += V->field_nelems;

	progress = tmp / PROGRESS_TICKS;
	progress = ( progress ) ? progress : 1;

	/* Determine Effective Scale */

	eff_scale = ( V->scale < 1.0 ) ? V->scale : 1.0;

	xfactor = eff_scale / ( (double) CELL_OF_DIM( V->viscell, xi ) );
	yfactor = eff_scale / ( (double) CELL_OF_DIM( V->viscell, yi ) );

	/* Fill In Viewing Field */

	if ( MAT->complex )
		C = MAT->complex;

	else if ( MAT->data )
		D = MAT->data;

	else
		T = MAT->template;

	/* Set Symmetry Loop Count */

	symm_cnt = ( MAT->symm == MS_UNSYMMETRIC ) ? 1 : 2;

	for ( i=0 ; i < MAT->nelems ; i++ )
	{
		/* Check for Partial Load */

		if ( i > MAT->last_loaded )
		{
			/* Get Next Chunk of Matrix Data */

			if ( !load_rest_of_mat( V ) )
				return( FALSE );

			/* Re-Display Generating View Field Message */

			Tcl_Eval( interp, "setMsg {Generating View Field...}" );
			INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

			/* Reset Data Pointers */

			if ( MAT->complex )
				C = MAT->complex;

			else if ( MAT->data )
				D = MAT->data;

			else
				T = MAT->template;
		}

		if ( !( i % progress ) )
		{
			update_progress( V, i / progress );
			INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );
		}

		/* Get Coords & Value */

		if ( MAT->complex )
		{
			xx = C->x;
			yy = C->y;

			if ( VIEW_PATTERN )
				val = 1.0;
			
			else
			{
				switch ( COMPLEX_VIEW )
				{
					case COMPLEX_RE:
						val = C->val;
						break;

					case COMPLEX_IM:
						val = C->jval;
						break;

					case COMPLEX_MAGN:
						val = hypot( C->val, C->jval );
						break;

					case COMPLEX_PHASE:
						val = atan2( C->jval, C->val );
						break;

					default:
						val = C->val;
						break;
				}
			}

			C++;
		}

		else if ( MAT->data )
		{
			xx = D->x;
			yy = D->y;

			val = ( VIEW_PATTERN ) ? 1.0 : D->val;

			D++;
		}

		else
		{
			xx = T->x;
			yy = T->y;

			val = 1.0;

			T++;
		}

		for ( j=0 ; j < symm_cnt ; j++ )
		{
			if ( j )
			{
				/* Don't Duplicate Diagonal Data Dummy! (CHR :-) */
				if ( xx == yy )
					break;

				switch ( MAT->symm )
				{
					case MS_SYMMETRIC:
					{
						tmp = xx;
						xx = yy;
						yy = tmp;

						break;
					}

					case MS_SKEW_SYMMETRIC:
					{
						tmp = xx;
						xx = yy;
						yy = tmp;

						if ( !VIEW_PATTERN && MAT->template == NULL )
							val *= -1.0;

						break;
					}

					case MS_HERMITIAN:
					{
						tmp = xx;
						xx = yy;
						yy = tmp;

						if ( !VIEW_PATTERN && MAT->complex )
						{
							switch ( COMPLEX_VIEW )
							{
								case COMPLEX_IM:
								case COMPLEX_PHASE:
									val *= -1.0;
									break;

								default:
									break;
							}
						}

						break;
					}
				}
			}

			/* In Region */

			if ( xx >= REGION_MIN( V->visregion, xi )
				&& xx <= REGION_MAX( V->visregion, xi )
				&& yy >= REGION_MIN( V->visregion, yi )
				&& yy <= REGION_MAX( V->visregion, yi ) )
			{
				/* Calculate View Field Coords */

				x = (int) ( (double) ( xx -
						REGION_MIN( V->visregion, xi ) )
					* xfactor );
				y = (int) ( (double) ( yy -
						REGION_MIN( V->visregion, yi ) )
					* yfactor );

				index = ( y * V->nx ) + x;

				if ( index > V->field_nelems )
					printf( "Whoa!  Index %d > V->field_nelems %d...\n",
						index, V->field_nelems );

				/* Process Reduction Function */

				F = &(V->field[ index ]);

				if ( V->need_reduction )
				{
					if ( use_tmps )
					{
						F1 = &(V->ftmp1[ index ]);
						F2 = &(V->ftmp2[ index ]);
					}

					switch ( REDUCTION_FUNCTION )
					{
						case RF_AVG:
						case RF_STD:
						case RF_VAR:
							if ( (*F2) == 0.0 )
							{
								x1 = (int) ceil( (double) x / xfactor );
								x2 = (int) ceil( (double) (x+1)
									/ xfactor );

								if ( x2 >
										REGION_MAX( V->visregion, xi ) )
									x2 = REGION_MAX( V->visregion, xi );

								y1 = (int) ceil( (double) y / yfactor );
								y2 = (int) ceil( (double) (y+1)
									/ yfactor );

								if ( y2 >
										REGION_MAX( V->visregion, yi ) )
									y2 = REGION_MAX( V->visregion, yi );

								*F2 = (double) ( ( x2 - x1 )
									* ( y2 - y1 ) );
							}
							*F1 += val;
							*F = (*F1) / (*F2);
							break;

						case RF_MIN:
							if ( (*F) == 0.0 )
								*F = val;
							else
								*F = ( val < (*F) ) ? val : *F;
							break;

						case RF_MAX:
							if ( (*F) == 0.0 )
								*F = val;
							else
								*F = ( val > (*F) ) ? val : *F;
							break;

						case RF_SUM:
							*F += val;
							break;

						case RF_MIDPOINT:
							if ( (*F1) == 0.0 && (*F2) == 0.0 )
								*F1 = *F2 = val;
							else
							{
								*F1 = ( val < (*F1) ) ? val : *F1;
								*F2 = ( val > (*F2) ) ? val : *F2;
							}
							*F = (*F1) + ( ( (*F2) - (*F1) ) / 2.0 );
							break;

						case RF_COUNT:
							if ( val != 0.0 )
								(*F)++;
							break;

						/* defer to avg, & skip the zillion warnings */
						default:
							if ( (*F2) == 0.0 )
							{
								*F2 =
								(ceil( (double) x / xfactor )
									- ceil( (double) (x+1) / xfactor ))
								* (ceil( (double) y / yfactor )
									- ceil( (double) (y+1) / yfactor ));
							}
							*F1 += val;
							*F = (*F1) / (*F2);
							break;
					}
				}

				else
					*F = val;
			}
		}
	}

	if ( V->need_reduction &&
		( REDUCTION_FUNCTION == RF_STD
			|| REDUCTION_FUNCTION == RF_VAR ) )
	{
		F1 = V->ftmp1;

		for ( i=0 ; i < V->field_nelems ; i++ )
			*F1++ = 0.0;

		if ( MAT->loaded == ML_PARTIAL )
		{
			if ( !read_mat_file( MATVIEW_MATFILE ) )
				return( FALSE );
		}

		if ( MAT->complex )
			C = MAT->complex;

		else if ( MAT->data )
			D = MAT->data;

		else
			T = MAT->template;

		for ( i=0 ; i < MAT->nelems ; i++ )
		{
			/* Check for Partial Load */

			if ( i > MAT->last_loaded )
			{
				/* Get Next Chunk of Matrix Data */

				load_rest_of_mat();

				/* Reset Data Pointers */

				if ( MAT->complex )
					C = MAT->complex;

				else if ( MAT->data )
					D = MAT->data;

				else
					T = MAT->template;
			}

			if ( !( ( MAT->nelems + i ) % progress ) )
			{
				update_progress( V, ( MAT->nelems + i ) / progress );
				INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );
			}

			/* Get Coords & Value */

			if ( MAT->complex )
			{
				xx = C->x;
				yy = C->y;

				val = ( VIEW_PATTERN ) ? 1.0 : C->val;

				C++;
			}

			else if ( MAT->data )
			{
				xx = D->x;
				yy = D->y;

				val = ( VIEW_PATTERN ) ? 1.0 : D->val;

				D++;
			}

			else
			{
				xx = T->x;
				yy = T->y;

				val = 1.0;

				T++;
			}

			for ( j=0 ; j < symm_cnt ; j++ )
			{
				if ( j )
				{
					switch ( MAT->symm )
					{
						case MS_SYMMETRIC:
						{
							tmp = xx;
							xx = yy;
							yy = tmp;

							break;
						}

						case MS_SKEW_SYMMETRIC:
						{
							tmp = xx;
							xx = yy;
							yy = tmp;

							if ( !VIEW_PATTERN
									&& MAT->template == NULL )
								val *= -1.0;

							break;
						}

						case MS_HERMITIAN:
						{
							tmp = xx;
							xx = yy;
							yy = tmp;

							if ( MAT->complex )
							{
								/* Im part only: val *= -1.0; */
							}

							break;
						}
					}
				}

				/* In Region */

				if ( xx >= REGION_MIN( V->visregion, xi )
					&& xx <= REGION_MAX( V->visregion, xi )
					&& yy >= REGION_MIN( V->visregion, yi )
					&& yy <= REGION_MAX( V->visregion, yi ) )
				{
					/* Calculate View Field Coords */

					x = (int) ( (double)
							( xx - REGION_MIN( V->visregion, xi ) )
						* eff_scale
							/ (double) CELL_OF_DIM( V->viscell, xi ) );
					y = (int) ( (double)
							( yy - REGION_MIN( V->visregion, yi ) )
						* eff_scale
							/ (double) CELL_OF_DIM( V->viscell, yi ) );

					index = ( y * V->nx ) + x;

					if ( index > V->field_nelems )
						printf(
							"Whoa!  Index %d > V->field_nelems %d...\n",
							index, V->field_nelems );

					/* Process Reduction Function */

					F1 = &(V->ftmp1[ index ]);

					*F1 += val * val;
				}
			}
		}

		F2 = &(V->ftmp2[0]);
		F1 = &(V->ftmp1[0]);
		F = &(V->field[0]);

		for ( i=0 ; i < V->field_nelems ; i++ )
		{
			if ( !( ( (2 * MAT->nelems) + i ) % progress ) )
			{
				update_progress( V,
					( (2 * MAT->nelems) + i ) / progress );
				INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );
			}

			if ( *F2 )
			{
				if ( (*F2) > 1 )
				{
					*F = ( (*F1) - ( (*F2) * (*F) * (*F) ) )
						/ ( (*F2) - 1 );

					if ( REDUCTION_FUNCTION == RF_STD )
						*F = sqrt( *F );
				}

				else
					*F = 0;
			}

			F1++;  F2++;  F++;
		}
	}

	else if ( VIEW_ABS_VALS )
	{
		F = &(V->field[0]);

		for ( i=0 ; i < V->field_nelems ; i++ )
		{
			if ( !( ( MAT->nelems + i ) % progress ) )
			{
				update_progress( V, ( MAT->nelems + i ) / progress );
				INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );
			}

			if ( *F < 0.0 )
				*F = fabs( *F );

			F++;
		}
	}

	return( TRUE );
}


void
compute_view_min_max( V, minptr, maxptr )
VIEW V;
double *minptr;
double *maxptr;
{
	double *dptr;
	double *ptr;

	double min, max;
	double fval;

	int vis_mode;
	int nelems;
	int i;

	/* Determine Vis Mode */

	vis_mode = ( VIS_MODE == VIS_MODE_DIFF && !(V->do_diff_mode) )
		? VIS_MODE_NORMAL : VIS_MODE;

	/* Set Data Pointers */

	ptr = V->field;

	switch ( vis_mode )
	{
		case VIS_MODE_NORMAL:
			dptr = (double *) NULL;
			break;

		case VIS_MODE_DIFF:
			dptr = SAVE_FIELD;
			nelems = SAVE_FIELD_NELEMS - 1;
			break;
	}

	/* Determine Min & Max Values */

	i = 0;

	do
	{
		if ( dptr != NULL )
		{
			fval = (*ptr) - (*dptr++);

			if ( i >= nelems )
				dptr = (double *) NULL;
		}

		else
			fval = *ptr;

		if ( !i )
			min = max = fval;

		else
		{
			min = ( fval < min ) ? fval : min;
			max = ( fval > max ) ? fval : max;
		}

		ptr++;

		i++;
	}
	while ( i < V->nrectids );

	if ( minptr )
		*minptr = min;
	
	if ( maxptr )
		*maxptr = max;
}


void
handle_min_max( min, max )
double min, max;
{
	static double minsave = 1.0;
	static double maxsave = -1.0;

	static int showsave = -1;

	char minstr[512];
	char maxstr[512];
	char cmd[1024];

	int changed = 0;
	int show;

	/* Handle Auto Color Map */

	if ( AUTO_COLOR_MAP )
	{
		if ( COLOR_MIN_VALUE != min || COLOR_MAX_VALUE != max )
		{
			COLOR_MIN_VALUE = min;
			COLOR_MAX_VALUE = max;

			update_color_range();

			Tcl_Eval( interp, "draw_color_map" );
		}
	}

	/* Force Redraw of Min & Max Display If State Has Changed... */

	show = ( ( SHOW_MIN_MAX ) ? 1 : 0 )
		+ ( ( ( SHOW_CUMULATIVE ) ? 1 : 0 ) * 10 );

	if ( show != showsave )
	{
		showsave = show;
		changed++;
	}

	/* Handle Min & Max, Cumulative Displays */

	if ( SHOW_CUMULATIVE )
	{
		if ( !CUMULATIVE_SET )
		{
			CUMULATIVE_MIN = min;
			CUMULATIVE_MAX = max;

			CUMULATIVE_SET++;

			changed++;
		}

		else
		{
			if ( min < CUMULATIVE_MIN )
			{
				CUMULATIVE_MIN = min;

				changed++;
			}

			if ( max > CUMULATIVE_MAX )
			{
				CUMULATIVE_MAX = max;

				changed++;
			}
		}
	}

	if ( SHOW_MIN_MAX || SHOW_CUMULATIVE )
	{
		if ( !SHOW_CUMULATIVE )
		{
			if ( min != minsave || max != maxsave )
			{
				sprintf( minstr, "Min: %lg", min );
				sprintf( maxstr, "Max: %lg", max );

				changed++;
			}
		}

		else if ( !SHOW_MIN_MAX )
		{
			if ( changed )
			{
				sprintf( minstr, "Cumulative Min: %lg",
					CUMULATIVE_MIN );
				sprintf( maxstr, "Cumulative Max: %lg",
					CUMULATIVE_MAX );
			}
		}

		else
		{
			if ( changed || min != minsave || max != maxsave )
			{
				sprintf( minstr, "Min: %lg / %lg",
					min, CUMULATIVE_MIN );
				sprintf( maxstr, "Max: %lg / %lg",
					max, CUMULATIVE_MAX );

				changed++;
			}
		}

		if ( changed )
		{
			sprintf( cmd, "%s.minval configure -text {%s}",
				MV, minstr );
			Tcl_Eval( interp, cmd );

			sprintf( cmd, "%s.maxval configure -text {%s}",
				MV, maxstr );
			Tcl_Eval( interp, cmd );
		}
	}
}


void
handle_grid( V )
VIEW V;
{
	GRID G, G2;

	GRID newgrids;

	char *color;

	int nrows, ncols;
	int newsize;
	int ngrids;
	int x1, x2;
	int y1, y2;
	int index;
	int start;
	int x, y;
	int mod;
	int i;

	if ( GRID_STATUS && V->do_grid )
	{
		/* Determine # of Grid Lines (Horiz & Vert) */

		nrows = (int) ceil( (double) ( REGION_MAX( V->visregion, 0 )
				- REGION_MIN( V->visregion, 0 ) + 1 )
			/ (double) GRID_ROWS ) + 1;

		ncols = (int) ceil( (double) ( REGION_MAX( V->visregion, 1 )
				- REGION_MIN( V->visregion, 1 ) + 1 )
			/ (double) GRID_COLS ) + 1;

		ngrids = nrows + ncols;

		/* Allocate Grid Structs */

		if ( !(V->draw_type & DRAW_IMAGE) )
		{
			if ( ngrids > V->gridsize || V->grids == NULL )
			{
				newsize = V->gridsize ? V->gridsize : 1;

				while ( newsize < ngrids )
					newsize *= 2;

				newgrids = (GRID) MALLOC( (unsigned) newsize
					* sizeof( struct grid_struct ) );
				memcheck( newgrids, "Grid Line Array" );

				G = &(newgrids[0]);

				if ( V->grids != NULL )
				{
					G2 = &(V->grids[0]);

					for ( i=0 ; i < V->ngrids ; i++ )
					{
						G->id = G2->id;

						G->x1 = G2->x1;
						G->y1 = G2->y1;
						G->x2 = G2->x2;
						G->y2 = G2->y2;

						G++;  G2++;
					}

					for ( i=(V->ngrids) ; i < newsize ; i++ )
					{
						G->id = -1;

						G->x1 = G->y1 = -1;
						G->x2 = G->y2 = -1;

						G++;
					}

					FREE( V->grids );
				}

				else
				{
					for ( i=0 ; i < newsize ; i++ )
					{
						G->id = -1;

						G->x1 = G->y1 = -1;
						G->x2 = G->y2 = -1;

						G++;
					}
				}

				V->grids = newgrids;

				V->gridsize = newsize;
			}

			V->ngrids = ngrids;
		}

		/* Set Grid Line Color */

		color = ( ZERO_COLOR != NULL && !strcmp( ZERO_COLOR, "black" ) )
			? "white" : "black";

		/* Initialize Grid Line Ptr */

		G = ( V->grids != NULL ) ? &(V->grids[0]) : NULL;

		/* Draw Column Grid Lines */

		start = REGION_MIN( V->visregion, 0 );

		index = start - REGION_MIN( MAX_REGION, 0 );

		mod = ( index - GRID_COFF ) % GRID_COLS;

		if ( mod )
			start -= mod;

		while ( start < GRID_COFF
			|| start < REGION_MIN( V->visregion, 0 ) )
		{
			start += GRID_COLS;
		}

		x = start;

		while ( x <= REGION_MAX( V->visregion, 0 ) )
		{
			y1 = V->win->yoff;
			y2 = ( V->iscale * V->ny ) + V->win->yoff;

			x1 = x2 = REGION_TO_SCREEN( x, V->win->xoff,
				REGION_MIN( V->visregion, 0 ),
				REGION_MAX( V->visregion, 0 ),
				V->nx, V->iscale );

			if ( V->draw_type & DRAW_LANDSCAPE )
			{
				landscape_coords( &x1, &y1, &x2, &y2,
					( V->iscale * V->ny ) + ( 2 * V->win->yoff ) );
			}

			if ( V->draw_type & DRAW_IMAGE )
			{
				if ( V->draw_type & DRAW_LANDSCAPE )
					y2++;
				else
					x2++;

				DRAW_IMAGE_RECT( interp, V->image,
					x1, y1, x2, y2, color );
			}

			else
			{
				if ( G->id != -1 )
				{
					LINE_COORDS( interp, V->win->canvas, G->id,
						x1, y1, x2, y2 );

					LINE_COLOR( interp, V->win->canvas, G->id, color );
				}

				else
				{
					G->id = DRAW_LINE( interp, V->win->canvas,
						x1, y1, x2, y2, color, 1 );
				}

				RAISE_ITEM( interp, V->win->canvas, G->id );

				G++;
			}

			x += GRID_COLS;
		}

		/* Draw Row Grid Lines */

		start = REGION_MIN( V->visregion, 1 );

		index = start - REGION_MIN( MAX_REGION, 1 );

		mod = ( index - GRID_ROFF ) % GRID_ROWS;

		if ( mod )
			start -= mod;

		while ( start < GRID_ROFF
			|| start < REGION_MIN( V->visregion, 1 ) )
		{
			start += GRID_ROWS;
		}

		y = start;

		while ( y <= REGION_MAX( V->visregion, 1 ) )
		{
			x1 = V->win->xoff;
			x2 = ( V->iscale * V->nx ) + V->win->xoff;

			y1 = y2 = REGION_TO_SCREEN( y, V->win->yoff,
				REGION_MIN( V->visregion, 1 ),
				REGION_MAX( V->visregion, 1 ),
				V->ny, V->iscale );

			if ( V->draw_type & DRAW_LANDSCAPE )
			{
				landscape_coords( &x1, &y1, &x2, &y2,
					( V->iscale * V->ny ) + ( 2 * V->win->yoff ) );
			}

			if ( V->draw_type & DRAW_IMAGE )
			{
				if ( V->draw_type & DRAW_LANDSCAPE )
					x2++;
				else
					y2++;

				DRAW_IMAGE_RECT( interp, V->image,
					x1, y1, x2, y2, color );
			}

			else
			{
				if ( G->id != -1 )
				{
					LINE_COORDS( interp, V->win->canvas, G->id,
						x1, y1, x2, y2 );

					LINE_COLOR( interp, V->win->canvas, G->id, color );
				}

				else
				{
					G->id = DRAW_LINE( interp, V->win->canvas,
						x1, y1, x2, y2, color, 1 );
				}

				RAISE_ITEM( interp, V->win->canvas, G->id );

				G++;
			}

			y += GRID_ROWS;
		}

		/* Free Any Leftover Grid Lines */

		if ( G )
		{
			while ( G->id != -1 )
			{
				DELETE_ITEM( interp, V->win->canvas, G->id );

				G->id = -1;

				G++;
			}
		}
	}

	else if ( V->grids != NULL )
	{
		G = &(V->grids[0]);

		for ( i=0 ; i < V->ngrids ; i++ )
		{
			if ( G->id != -1 )
			{
				DELETE_ITEM( interp, V->win->canvas, G->id );

				G->id = -1;
			}

			G++;
		}
	}
}


int
transpose_matrix()
{
	TEMPLATE T;

	COMPLEX C;

	DATA D;

	int tmp;
	int i;

	/* Check for Partial Load */

	if ( MAT->loaded == ML_PARTIAL )
	{
		load_mat( TRUE );

		return( FALSE );
	}

	/* Transpose Data / Template */

	if ( MAT->complex != NULL )
	{
		C = MAT->complex;

		for ( i=0 ; i < MAT->nelems ; i++ )
		{
			tmp = C->x;
			C->x = C->y;
			C->y = tmp;

			C++;
		}
	}

	else if ( MAT->data != NULL )
	{
		D = MAT->data;

		for ( i=0 ; i < MAT->nelems ; i++ )
		{
			tmp = D->x;
			D->x = D->y;
			D->y = tmp;

			D++;
		}
	}

	else if ( MAT->template != NULL )
	{
		T = MAT->template;

		for ( i=0 ; i < MAT->nelems ; i++ )
		{
			tmp = T->x;
			T->x = T->y;
			T->y = tmp;

			T++;
		}
	}

	/* Transpose GLB & GUB */

	tmp = MAT->glb[0];
	MAT->glb[0] = MAT->glb[1];
	MAT->glb[1] = tmp;

	tmp = MAT->gub[0];
	MAT->gub[0] = MAT->gub[1];
	MAT->gub[1] = tmp;

	/* Transpose Regions */

	/* Max */

	tmp = REGION_MIN( MAX_REGION, 0 );
	REGION_MIN( MAX_REGION, 0 ) = REGION_MIN( MAX_REGION, 1 );
	REGION_MIN( MAX_REGION, 1 ) = tmp;

	tmp = REGION_MAX( MAX_REGION, 0 );
	REGION_MAX( MAX_REGION, 0 ) = REGION_MAX( MAX_REGION, 1 );
	REGION_MAX( MAX_REGION, 1 ) = tmp;

	/* Thumb */

	tmp = REGION_MIN( THUMB_VIEW->visregion, 0 );
	REGION_MIN( THUMB_VIEW->visregion, 0 ) =
		REGION_MIN( THUMB_VIEW->visregion, 1 );
	REGION_MIN( THUMB_VIEW->visregion, 1 ) = tmp;

	tmp = REGION_MAX( THUMB_VIEW->visregion, 0 );
	REGION_MAX( THUMB_VIEW->visregion, 0 ) =
		REGION_MAX( THUMB_VIEW->visregion, 1 );
	REGION_MAX( THUMB_VIEW->visregion, 1 ) = tmp;

	tmp = CELL_OF_DIM( THUMB_VIEW->viscell, 0 );
	CELL_OF_DIM( THUMB_VIEW->viscell, 0 ) =
		CELL_OF_DIM( THUMB_VIEW->viscell, 1 );
	CELL_OF_DIM( THUMB_VIEW->viscell, 1 ) = tmp;

	/* Main */

	tmp = REGION_MIN( MAIN_VIEW->visregion, 0 );
	REGION_MIN( MAIN_VIEW->visregion, 0 ) =
		REGION_MIN( MAIN_VIEW->visregion, 1 );
	REGION_MIN( MAIN_VIEW->visregion, 1 ) = tmp;

	tmp = REGION_MAX( MAIN_VIEW->visregion, 0 );
	REGION_MAX( MAIN_VIEW->visregion, 0 ) =
		REGION_MAX( MAIN_VIEW->visregion, 1 );
	REGION_MAX( MAIN_VIEW->visregion, 1 ) = tmp;

	tmp = CELL_OF_DIM( MAIN_VIEW->viscell, 0 );
	CELL_OF_DIM( MAIN_VIEW->viscell, 0 ) =
		CELL_OF_DIM( MAIN_VIEW->viscell, 1 );
	CELL_OF_DIM( MAIN_VIEW->viscell, 1 ) = tmp;

	/* Update Interface */

	update_vis_region();

	/* Return Redraw Status */

	return( TRUE );
}


void
draw_thumb_border()
{
	draw_thumb_border_region( MAIN_VIEW->visregion );
}


void
draw_thumb_border_region( visregion )
REGION visregion;
{
	char *color;

	int x1, y1, x2, y2;

	if ( MAT->loaded == ML_UNLOADED )
		return;

	x1 = REGION_TO_SCREEN( REGION_MIN( visregion, 0 ),
		THUMB_VIEW->win->xoff,
		REGION_MIN( THUMB_VIEW->visregion, 0 ),
		REGION_MAX( THUMB_VIEW->visregion, 0 ),
		THUMB_VIEW->nx, THUMB_VIEW->iscale ) - 1;

	y1 = REGION_TO_SCREEN( REGION_MIN( visregion, 1 ),
		THUMB_VIEW->win->yoff,
		REGION_MIN( THUMB_VIEW->visregion, 1 ),
		REGION_MAX( THUMB_VIEW->visregion, 1 ),
		THUMB_VIEW->ny, THUMB_VIEW->iscale ) - 1;

	x2 = REGION_TO_SCREEN( REGION_MAX( visregion, 0 ) + 1,
		THUMB_VIEW->win->xoff,
		REGION_MIN( THUMB_VIEW->visregion, 0 ),
		REGION_MAX( THUMB_VIEW->visregion, 0 ),
		THUMB_VIEW->nx, THUMB_VIEW->iscale ) + 1;

	y2 = REGION_TO_SCREEN( REGION_MAX( visregion, 1 ) + 1,
		THUMB_VIEW->win->yoff,
		REGION_MIN( THUMB_VIEW->visregion, 1 ),
		REGION_MAX( THUMB_VIEW->visregion, 1 ),
		THUMB_VIEW->ny, THUMB_VIEW->iscale ) + 1;

	if ( THUMB_BORDER == -1 )
	{
		color = ( ZERO_COLOR != NULL && !strcmp( ZERO_COLOR, "black" ) )
			? "white" : "black";

		THUMB_BORDER = DRAW_RECT( interp, THUMB_VIEW->win->canvas,
			x1, y1, x2, y2, "", color, 1, "thumb_border" );
	}

	else
	{
		RECT_COORDS( interp, THUMB_VIEW->win->canvas, THUMB_BORDER,
			x1, y1, x2, y2 );
	}

	RAISE_ITEM( interp, THUMB_VIEW->win->canvas, THUMB_BORDER );
}


void
update_progress( V, percent )
VIEW V;
int percent;
{
	char cmd[255];

	int width;

	width = atoi( PROGRESS_WIDTH );

	sprintf( cmd, "%s coords %d 0 0 %d %s",
		V->progress_canvas, V->progress_bar,
		percent * width / PROGRESS_TICKS, PROGRESS_HEIGHT );

	Tcl_Eval( interp, cmd );

	Tcl_Eval( interp, "update" );
}


void
hide_image( V )
VIEW V;
{
	char cmd[255];

	if ( !(V->hide_image)++ )
	{
		sprintf( cmd, "%s raise %d", V->win->canvas, V->rect_id );

		Tcl_Eval( interp, cmd );

		Tcl_Eval( interp, "update" );
	}
}


void
show_image( V )
VIEW V;
{
	char cmd[255];

	if ( !(--(V->hide_image)) )
	{
		sprintf( cmd, "%s lower %d", V->win->canvas, V->rect_id );

		Tcl_Eval( interp, cmd );

		Tcl_Eval( interp, "update" );
	}
}


void
set_current_field()
{
	double *F, *S;

	int i;

	/* Allocate Saved Field Storage & Copy Data */

	if ( SAVE_FIELD != NULL )
		FREE( SAVE_FIELD );

	SAVE_FIELD = (double *) MALLOC( (unsigned) MAIN_VIEW->field_nelems
		* sizeof(double) );

	if ( SAVE_FIELD == NULL )
	{
		printf( "\nError: Insufficient Memory to Save Field!\n\n" );
	}

	else
	{
		S = SAVE_FIELD;
		F = MAIN_VIEW->field;

		for ( i=0 ; i < MAIN_VIEW->field_nelems ; i++ )
			*S++ = *F++;

		SAVE_FIELD_NELEMS = MAIN_VIEW->field_nelems;
	}
}


void
free_current_field()
{
	int i;

	if ( SAVE_FIELD != NULL )
	{
		for ( i=0 ; i < SAVE_FIELD_NELEMS ; i++ )
			SAVE_FIELD[i] = -1.0;

		FREE( SAVE_FIELD );

		SAVE_FIELD = (double *) NULL;
	}

	SAVE_FIELD_NELEMS = -1;

	if ( VIS_MODE == VIS_MODE_DIFF )
		Tcl_Eval( interp, "set_vis_mode diff" );
}


void
save_current_field( file )
char *file;
{
	FILE *fp;

	double *fptr;

	int i;

	if ( SAVE_FIELD == NULL )
		set_current_field();

	fp = fopen( file, "w" );

	if ( !filecheck( fp, file ) )
		return;

	fprintf( fp, "%d\n", SAVE_FIELD_NELEMS );

	fptr = SAVE_FIELD;

	for ( i=1 ; i <= SAVE_FIELD_NELEMS ; i++ )
		fprintf( fp, "%lf%c", *fptr++, ( !(i % 8) ? '\n' : ' ' ) );

	fclose( fp );
}


void
load_current_field( file )
char *file;
{
	FILE *fp;

	double *fptr;

	int nelems = -1;
	int i;

	fp = fopen( file, "r" );

	if ( !filecheck( fp, file ) )
		return;

	if ( fscanf( fp, "%d", &nelems ) != 1 || nelems <= 0 )
	{
		printf( "%s:  Error Parsing Saved Field File - nelems=%d...\n",
			"load_current_field()", nelems );

		return;
	}

	if ( SAVE_FIELD != NULL )
		FREE( SAVE_FIELD );

	SAVE_FIELD = (double *) MALLOC( (unsigned) nelems
		* sizeof(double) );

	if ( SAVE_FIELD == NULL )
		printf( "\nError: Insufficient Memory to Load Field!\n\n" );

	else
	{
		/* Initialize Saved Field Storage */

		fptr = SAVE_FIELD;

		for ( i=0 ; i < nelems ; i++ )
			*fptr++ = 0.0;

		SAVE_FIELD_NELEMS = nelems;

		/* Read in Field */

		fptr = SAVE_FIELD;

		for ( i=0 ; i < nelems ; i++ )
		{
			if ( fscanf( fp, "%lf", fptr++ ) != 1 )
			{
				printf(
				"%s:  Error Parsing Saved Field File (%d of %d)...\n",
					"load_current_field()", i, nelems );

				break;
			}
		}
	}

	fclose( fp );
}


void
save_image_to_file( file, format )
char *file;
char *format;
{
	char cmd[255];

	int x1, y1, x2, y2;

	/* Set Busy Cursor */

	Tcl_Eval( interp, "set_cursor busy" );

	/* Dump Postscript to File from Canvas */

	if ( !strcmp( format, "ps" ) )
	{
		Tcl_Eval( interp,
			"setMsg {Saving View Image to PostScript File...}" );

		sprintf( cmd, "%s postscript -file %s",
			PRINT_VIEW->win->canvas, file );

		Tcl_Eval( interp, cmd );

		Tcl_Eval( interp, "setMsg {PostScript File Written.}" );
	}

	/* Dump PPM to File from Image */

	else if ( !strcmp( format, "ppm" ) )
	{
		Tcl_Eval( interp,
			"setMsg {Saving View Image to PPM File...}" );

		/* Calculate Bounds */

		x1 = PRINT_VIEW->win->xoff;
		y1 = PRINT_VIEW->win->yoff;

		x2 = ( PRINT_VIEW->iscale * PRINT_VIEW->nx )
			+ PRINT_VIEW->win->xoff;
		y2 = ( PRINT_VIEW->iscale * PRINT_VIEW->ny )
			+ PRINT_VIEW->win->yoff;

		/* Swap Bounds for Landscape */

		if ( PRINT_VIEW->draw_type & DRAW_LANDSCAPE )
		{
			landscape_coords( &x1, &y1, &x2, &y2,
				( PRINT_VIEW->iscale * PRINT_VIEW->ny )
					+ ( 2 * PRINT_VIEW->win->yoff ) );
		}

		/* Dump It */

		sprintf( cmd, "%s write %s -format ppm -from %d %d %d %d",
			PRINT_VIEW->image, file,
			x1 - 1, y1 - 1, x2, y2 );

		Tcl_Eval( interp, cmd );

		Tcl_Eval( interp, "setMsg {PPM File Written.}" );
	}

	else
	{
		sprintf( cmd, "setMsg {Error:  Unkown Print Format \"%s\"}",
			format );

		Tcl_Eval( interp, cmd );
	}

	/* Reset Cursor */

	Tcl_Eval( interp, "set_cursor reset" );
}


void
color_value( color, value )
char color[16];
double value;
{
	char Rstr[4];
	char Gstr[4];
	char Bstr[4];

	double ratio;
	double f;

	double R;
	double G;
	double B;

	if ( EQUIV( COLOR_MIN_VALUE, COLOR_MAX_VALUE, 0.000001 ) )
	{
		if ( value > COLOR_MIN_VALUE )
		{
			if ( GREYSCALE )
			{
				R = 0.0;
				G = 0.0;
				B = 0.0;
			}

			else
			{
				R = 255.0;
				G = 0.0;
				B = 0.0;
			}
		}

		else if ( value < COLOR_MIN_VALUE )
		{
			if ( GREYSCALE )
			{
				R = 255.0;
				G = 255.0;
				B = 255.0;
			}

			else
			{
				R = 0.0;
				G = 0.0;
				B = 255.0;
			}
		}

		else
		{
			if ( GREYSCALE )
			{
				R = 127.0;
				G = 127.0;
				B = 127.0;
			}

			else
			{
				R = 0.0;
				G = 255.0;
				B = 0.0;
			}
		}
	}

	else
	{
		ratio = ( value - COLOR_MIN_VALUE )
			 / ( COLOR_MAX_VALUE - COLOR_MIN_VALUE );

		ratio = ( ratio < 0.0 ) ? 0.0 : ratio;

		f = pow( ratio, COLOR_DIST );

		if ( GREYSCALE )
		{
			B = 255.0 * ( 1.0 - f );

			B = ( B < 0 ) ? 0.0 : B;
			B = ( B > 255.0 ) ? 255.0 : B;

			R = G = B;
		}

		else
		{
			B = 255.0 * ( 1.0 - f );

			B = ( B < 0 ) ? 0.0 : B;
			B = ( B > 255.0 ) ? 255.0 : B;

			if ( f > 0.5 )
				G = 255.0 * ( 1.0 - ( ( f - 0.5 ) * 2.0 ) );

			else
				G = 255.0 * f * 2.0;

			G = ( G < 0 ) ? 0.0 : G;
			G = ( G > 255.0 ) ? 255.0 : G;

			R = 255.0 * f;

			R = ( R < 0 ) ? 0.0 : R;
			R = ( R > 255.0 ) ? 255.0 : R;
		}
	}

	if ( (int) R > 0xf )
		sprintf( Rstr, "%2x", (int) R );
	else
		sprintf( Rstr, "0%1x", (int) R );

	if ( (int) G > 0xf )
		sprintf( Gstr, "%2x", (int) G );
	else
		sprintf( Gstr, "0%1x", (int) G );

	if ( (int) B > 0xf )
		sprintf( Bstr, "%2x", (int) B );
	else
		sprintf( Bstr, "0%1x", (int) B );

	sprintf( color, "#%s%s%s", Rstr, Gstr, Bstr );
}

