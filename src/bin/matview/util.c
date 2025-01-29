
/* $Id: util.c,v 1.39 2000/06/19 18:30:46 kohl Exp $ */


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


MATRIX
create_matrix()
{
	MATRIX tmp;

	int i;

	tmp = (MATRIX) MALLOC( sizeof( struct matrix_struct ) );
	memcheck( tmp, "Matrix Object" );

	tmp->template = (TEMPLATE) NULL;
	tmp->complex = (COMPLEX) NULL;
	tmp->data = (DATA) NULL;

	tmp->alloc_size = -1;
	tmp->nelems = -1;
	tmp->dim = -1;

	for ( i=0 ; i < MAX_DIM ; i++ )
	{
		tmp->glb[i] = -1;
		tmp->gub[i] = -1;
	}

	tmp->symm = -1;

	tmp->fpsave = (FILE *) NULL;
	tmp->is_compressed = -1;
	tmp->read_mat_func = (ifp) NULL;
	tmp->load_inst = -1;
	tmp->first_loaded = -1;
	tmp->last_loaded = -1;
	tmp->loaded = -1;

	return( tmp );
}


void
free_matrix( ptr )
MATRIX *ptr;
{
	MATRIX M;

	int i;

	if ( ptr == NULL || *ptr == NULL )
		return;
	
	M = *ptr;

	if ( M->template != NULL )
	{
		FREE( M->template );
		M->template = (TEMPLATE) NULL;
	}

	if ( M->complex != NULL )
	{
		FREE( M->complex );
		M->complex = (COMPLEX) NULL;
	}

	if ( M->data != NULL )
	{
		FREE( M->data );
		M->data = (DATA) NULL;
	}

	M->alloc_size = -1;
	M->nelems = -1;
	M->dim = -1;

	for ( i=0 ; i < MAX_DIM ; i++ )
	{
		M->glb[i] = -1;
		M->gub[i] = -1;
	}

	M->symm = -1;

	M->fpsave = (FILE *) NULL;
	M->is_compressed = -1;
	M->read_mat_func = (ifp) NULL;
	M->load_inst = -1;
	M->first_loaded = -1;
	M->last_loaded = -1;
	M->loaded = -1;

	FREE( M );

	*ptr = (MATRIX) NULL;
}


VIEW
create_view()
{
	VIEW tmp;

	tmp = (VIEW) MALLOC( sizeof( struct view_struct ) );
	memcheck( tmp, "View Object" );

	tmp->win = (WIN) NULL;

	init_region( tmp->visregion );

	init_cell( tmp->viscell );

	tmp->rectids = (RECT) NULL;
	tmp->rectidsize = 0;
	tmp->nrectids = 0;

	tmp->nx = tmp->ny = -1;

	tmp->grids = (GRID) NULL;
	tmp->gridsize = 0;
	tmp->ngrids = 0;

	tmp->do_auto_color = -1;
	tmp->do_diff_mode = -1;
	tmp->do_grid = -1;

	tmp->field = (double *) NULL;
	tmp->ftmp1 = (double *) NULL;
	tmp->ftmp2 = (double *) NULL;
	tmp->field_nelems = -1;

	tmp->image = (char *) NULL;

	tmp->wt = tmp->ht = -1;

	tmp->scale = -1.0;

	tmp->iscale = -1;

	tmp->progress_canvas = (char *) NULL;
	tmp->progress_bar = -1;

	tmp->need_reduction = -1;
	tmp->need_redraw = -1;
	tmp->hide_image = 0;
	tmp->draw_type = 0;
	tmp->instance = 0;
	tmp->sequence = 0;
	tmp->fileinst = 0;
	tmp->last_rf = -1;

	tmp->image_id = -1;
	tmp->rect_id = -1;

	return( tmp );
}


void
free_view( ptr )
VIEW *ptr;
{
	VIEW V;

	if ( ptr == NULL || *ptr == NULL )
		return;
	
	V = *ptr;

	if ( V->win != NULL )
		free_win( &(V->win) );

	init_region( V->visregion );

	init_cell( V->viscell );

	if ( V->rectids != NULL )
	{
		FREE( V->rectids );
		V->rectids = (RECT) NULL;
	}

	V->rectidsize = 0;
	V->nrectids = 0;

	V->nx = V->ny = -1;

	if ( V->grids != NULL )
	{
		FREE( V->grids );
		V->grids = (GRID) NULL;
	}

	V->gridsize = 0;
	V->ngrids = 0;

	V->do_auto_color = -1;
	V->do_diff_mode = -1;
	V->do_grid = -1;

	if ( V->field != NULL )
	{
		FREE( V->field );
		V->field = (double *) NULL;
	}

	if ( V->ftmp1 != NULL )
	{
		FREE( V->ftmp1 );
		V->ftmp1 = (double *) NULL;
	}

	if ( V->ftmp2 != NULL )
	{
		FREE( V->ftmp2 );
		V->ftmp2 = (double *) NULL;
	}

	V->field_nelems = -1;

	V->image = (char *) NULL;

	V->wt = V->ht = -1;

	V->scale = -1.0;

	V->iscale = -1;

	V->progress_canvas = (char *) NULL;
	V->progress_bar = -1;

	V->need_reduction = -1;
	V->need_redraw = -1;
	V->hide_image = 0;
	V->draw_type = 0;
	V->instance = 0;
	V->sequence = 0;
	V->fileinst = 0;
	V->last_rf = -1;

	V->image_id = -1;
	V->rect_id = -1;

	FREE( V );

	*ptr = (VIEW) NULL;
}


void
swap_views( Vptr1, Vptr2 )
VIEW *Vptr1, *Vptr2;
{
	VIEW V1, V2;

	GRID G;

	int i;

	if ( Vptr1 == NULL || *Vptr1 == NULL
			|| Vptr2 == NULL || *Vptr2 == NULL )
		return;

	V1 = *Vptr1;
	V2 = *Vptr2;

	/* Copy Latest Region & Cell Over */

	copy_region( V1->visregion, V2->visregion );

	copy_cell( V1->viscell, V2->viscell );

	/* Pass Rectid Storage Over */

	V1->rectids = V2->rectids;
	V1->rectidsize = V2->rectidsize;
	V1->nrectids = V2->nrectids;

	V2->rectids = (RECT) NULL;
	V2->rectidsize = 0;
	V2->nrectids = 0;

	/* Wipe Out Any Grid Lines */

	if ( V2->grids != NULL )
	{
		G = &(V2->grids[0]);

		for ( i=0 ; i < V2->ngrids ; i++ )
		{
			if ( G->id != -1 )
			{
				DELETE_ITEM( interp, V2->win->canvas, G->id );

				G->id = -1;
			}

			G++;
		}
	}

	/* Pass Grid Storage Over */

	V1->grids = V2->grids;
	V1->gridsize = V2->gridsize;
	V1->ngrids = V2->ngrids;

	V2->grids = (GRID) NULL;
	V2->gridsize = 0;
	V2->ngrids = 0;

	/* Free Up Any Field Arrays */

	if ( V2->field != NULL )
	{
		FREE( V2->field );
		V2->field = (double *) NULL;
	}

	if ( V2->ftmp1 != NULL )
	{
		FREE( V2->ftmp1 );
		V2->ftmp1 = (double *) NULL;
	}

	if ( V2->ftmp2 != NULL )
	{
		FREE( V2->ftmp2 );
		V2->ftmp2 = (double *) NULL;
	}

	V1->field_nelems = -1;
	V1->nx = V1->ny = -1;
	V1->last_rf = -1;

	V2->field_nelems = -1;
	V2->nx = V2->ny = -1;
	V2->last_rf = -1;

	/* Bump Frame Instance & Drawing Sequence Numbers */
	/* Keep Consistent (Unique) Between Views... */

	V1->instance = V2->instance + 1;
	V1->sequence = V2->sequence + 1;

	/* Pass Along File Instance Number (Do Not Disturb!) */

	V1->fileinst = V2->fileinst;

	/* Swap Ptrs */

	*Vptr1 = V2;
	*Vptr2 = V1;
}


WIN
create_win()
{
	WIN tmp;

	tmp = (WIN) MALLOC( sizeof( struct win_struct ) );
	memcheck( tmp, "Window Display Object" );

	tmp->toplevel = (char *) NULL;

	tmp->canvas = (char *) NULL;

	tmp->htvar = (char *) NULL;
	tmp->wtvar = (char *) NULL;

	tmp->scale = -1.0;

	tmp->csx = tmp->csy = -1;

	init_region( tmp->region );

	tmp->xoff = tmp->yoff = -1;

	tmp->border_id = -1;

	return( tmp );
}


void
free_win( ptr )
WIN *ptr;
{
	WIN W;

	if ( ptr == NULL || *ptr == NULL )
		return;
	
	W = *ptr;

	W->toplevel = (char *) NULL;

	W->canvas = (char *) NULL;

	W->htvar = (char *) NULL;
	W->wtvar = (char *) NULL;

	W->scale = -1.0;

	W->csx = W->csy = -1;

	init_region( W->region );

	W->xoff = W->yoff = -1;

	W->border_id = -1;

	FREE( W );

	*ptr = (WIN) NULL;
}


/* Region & Cell Size Routines */


void
init_region( R )
REGION R;
{
	int i;

	for ( i=0 ; i < MAX_DIM ; i++ )
	{
		REGION_MIN( R, i ) = -1;
		REGION_MAX( R, i ) = -1;
	}
}


void
copy_region( R1, R2 )
REGION R1, R2;
{
	int i;

	for ( i=0 ; i < MAX_DIM ; i++ )
	{
		REGION_MIN( R1, i ) = REGION_MIN( R2, i );
		REGION_MAX( R1, i ) = REGION_MAX( R2, i );
	}
}


int
match_region( R1, R2, dim )
REGION R1, R2;
int dim;
{
	int i;

	for ( i=0 ; i < dim ; i++ )
	{
		if ( REGION_MIN( R1, i ) != REGION_MIN( R2, i ) )
			return( FALSE );

		if ( REGION_MAX( R1, i ) != REGION_MAX( R2, i ) )
			return( FALSE );
	}

	return( TRUE );
}


void
init_cell( C )
CELL C;
{
	int i;

	for ( i=0 ; i < CELL_SIZE( MAX_DIM ) ; i++ )
		CELL_OF_DIM( C, i ) = -1;
}


void
copy_cell( C1, C2 )
CELL C1, C2;
{
	int i;

	for ( i=0 ; i < CELL_SIZE( MAX_DIM ) ; i++ )
		C1[i] = C2[i];
}


int
match_cell( C1, C2, dim )
CELL C1, C2;
int dim;
{
	int i;

	for ( i=0 ; i < dim ; i++ )
	{
		if ( C1[i] != C2[i] )
			return( FALSE );
	}

	return( TRUE );
}


int
check_region( R, C, glb, gub, dim )
REGION R;
CELL C;
int *glb;
int *gub;
int dim;
{
	int changed;
	int diff;
	int i;

	changed = FALSE;

	/* Check Against Decomposition Global Bounds */

	for ( i=0 ; i < dim ; i++ )
	{
		if ( REGION_MIN( R, i ) < glb[i] )
		{
			REGION_MIN( R, i ) = glb[i];

			changed = TRUE;
		}

		if ( REGION_MIN( R, i ) > gub[i] )
		{
			REGION_MIN( R, i ) = gub[i];

			changed = TRUE;
		}

		if ( REGION_MAX( R, i ) > gub[i] )
		{
			REGION_MAX( R, i ) = gub[i];

			changed = TRUE;
		}

		if ( REGION_MAX( R, i ) < glb[i] )
		{
			REGION_MAX( R, i ) = glb[i];

			changed = TRUE;
		}
	}

	/* Check Divisibility with Cell Size */

	for ( i=0 ; i < dim ; i++ )
	{
		/* First, Make Sure Cell Size is Positive...  D-Oh! */

		if ( CELL_OF_DIM( C, i ) <= 0 )
		{
			CELL_OF_DIM( C, i ) = 1;

			changed = TRUE;
		}

		/* Leave Min Where It Is, Use Min as Offset for Divisibility */

		/* Bring ( Max - Min ) Down to Last Divisible By Cell */

		diff = REGION_MAX( R, i ) - REGION_MIN( R, i );

		while ( !DIVISIBLE( diff, CELL_OF_DIM( C, i ) ) && diff > 0 )
		{
			REGION_MAX( R, i )--;

			changed = TRUE;

			diff = REGION_MAX( R, i ) - REGION_MIN( R, i );
		}
	}

	/* Check for Valid Region */

	for ( i=0 ; i < dim ; i++ )
	{
		/* Bogus Region, Adjust Bounds to Fix */

		if ( REGION_MIN( R, i ) > REGION_MAX( R, i ) )
		{
			REGION_MAX( R, i ) = REGION_MIN( R, i );

			changed = TRUE;
		}
	}

	return( changed );
}


REGION_STACK create_region_stack()
{
	REGION_STACK tmp;

	int i;

	tmp = (REGION_STACK) MALLOC( sizeof( struct region_stack_struct ) );
	memcheck( tmp, "Region Stack Entry" );

	for ( i=0 ; i < REGION_SIZE( MAX_DIM ) ; i++ )
		tmp->region[i] = -1;
	
	tmp->next = (REGION_STACK) NULL;

	return( tmp );
}


void
free_region_stack( ptr )
REGION_STACK *ptr;
{
	REGION_STACK RS;

	int i;

	if ( ptr == NULL || *ptr == NULL )
		return;
	
	RS = *ptr;

	for ( i=0 ; i < REGION_SIZE( MAX_DIM ) ; i++ )
		RS->region[i] = -1;
	
	RS->next = (REGION_STACK) NULL;

	FREE( RS );

	*ptr = (REGION_STACK) NULL;
}


void
push_region( R )
REGION R;
{
	REGION_STACK RS;

	/* Create New Region Stack Entry */

	RS = create_region_stack();

	copy_region( RS->region, R );

	/* Add to Regular Region Stack List */

	RS->next = REGIONS;

	REGIONS = RS;

	/* Clean Up Any Saved Re-Push Regions */

	while ( REGIONS_SAVE )
	{
		RS = REGIONS_SAVE;

		REGIONS_SAVE = REGIONS_SAVE->next;

		free_region_stack( &RS );
	}
}


int
pop_region( R )
REGION R;
{
	REGION_STACK RS;

	REGION tmp;

	if ( REGIONS != NULL )
	{
		/* Save Current Region */

		copy_region( tmp, R );

		/* Get Region Values */

		copy_region( R, REGIONS->region );

		/* Remove From Region Stack */

		RS = REGIONS;

		REGIONS = REGIONS->next;

		/* Save on Re-Push Stack */

		RS->next = REGIONS_SAVE;

		REGIONS_SAVE = RS;

		copy_region( RS->region, tmp );

		/* Yep, We Popped */

		return( TRUE );
	}

	return( FALSE );
}


int
repush_region( R )
REGION R;
{
	REGION_STACK RS;

	REGION tmp;

	if ( REGIONS_SAVE != NULL )
	{
		/* Save Current Region */

		copy_region( tmp, R );

		/* Get Region Values */

		copy_region( R, REGIONS_SAVE->region );

		/* Remove From Re-Push Region Stack */

		RS = REGIONS_SAVE;

		REGIONS_SAVE = REGIONS_SAVE->next;

		/* Save Back on Regular Stack */

		RS->next = REGIONS;

		REGIONS = RS;

		copy_region( RS->region, tmp );

		/* Yep, We Re-Pushed */

		return( TRUE );
	}

	return( FALSE );
}


int
incr_field_index( dim, indices, max, storage_order )
int dim;
int *indices;
int *max;
int storage_order;
{
	int index;
	int flag;
	int incr;

	if ( storage_order == ROW_MAJOR )
	{
		index = 0;
		incr = 1;
	}

	else if ( storage_order == COLUMN_MAJOR )
	{
		index = dim - 1;
		incr = -1;
	}

	else
	{
		printf( "%s: Unknown Data Storage Order, ID = %d",
			"incr_field_index()", storage_order );

		return( FALSE );
	}

	do
	{
		flag = 0;

		(indices[index])++;

		if ( indices[index] >= max[index] )
		{
			indices[index] = 0;

			flag++;
		}

		if ( flag )
		{
			index += incr;

			if ( index >= dim || index < 0 )
				return( FALSE );
		}
	}
	while ( flag );

	return( TRUE );
}


void
landscape_coords( x1, y1, x2, y2, ht )
int *x1, *y1, *x2, *y2;
int ht;
{
	int tmp;

	tmp = *x1;
	*x1 = ht - *y1;
	*y1 = tmp;

	tmp = *x2;
	*x2 = ht - *y2;
	*y2 = tmp;
}


/* Other Utilities */


char *
center_str( str, width )
char *str;
int width;
{
	char *tmp;

	int i, j;
	int len;

	tmp = (char *) MALLOC( (unsigned) (width + 1) * sizeof(char) );
	memcheck( tmp, "Center String" );

	len = strlen( str );

	if ( len > width )
	{
		strcpy( tmp, (char *) ( str + ( len - width ) ) );

		tmp[0] = '<';
	}

	else
	{
		j = ( width - len ) / 2;

		for ( i=0 ; i < j ; i++ )
			tmp[i] = ' ';

		strcpy( (char *) ( tmp + j ), str );

		for ( i=(j + len) ; i < width ; i++ )
			tmp[i] = ' ';

		tmp[width] = '\0';
	}

	return( tmp );
}


char *
pad_num( num, max )
long num;
long max;
{
	char tmp[1024];

	char *ptr;
	char *str;

	int nd;
	int nz;
	int i;

	sprintf( tmp, "%ld", num );

	nd = strlen( tmp );

	nz = max - nd;

	str = (char *) MALLOC( (unsigned) (max + 1) * sizeof(char) );
	memcheck( str, "Numerical Pad String" );

	ptr = str;

	for ( i=0 ; i < nz ; i++ )
		*ptr++ = '0';

	sprintf( (char *) (str + nz), "%ld", num );

	return( str );
}


char *
trunc_str( str, len, mark_incl )
char *str;
int len;
int mark_incl;
{
	char tmp[1024];

	char *ptr;

	int slen;

	sprintf( tmp, "%s", str );

	slen = strlen( str );

	if ( slen > len )
	{
		if ( mark_incl )
			ptr = tmp + len - 1;
		
		else
			ptr = tmp + len;

		*ptr++ = '*';

		*ptr = '\0';
	}

	return( copy_str( tmp ) );
}


char *
upper_str( str )
char *str;
{
	char *ptr;
	char *tmp;

	tmp = copy_str( str );

	ptr = tmp;

	while ( *ptr != '\0' )
	{
		if ( *ptr >= 'a' && *ptr <= 'z' )
			*ptr += 'A' - 'a';

		ptr++;
	}

	return( tmp );
}


char *
suffix_of( str )
char *str;
{
	char *ptr;

	int len;

	if ( str == NULL )
		return( copy_str( "" ) );

	len = strlen( str );

	ptr = &(str[len]);

	while ( ptr != str && *ptr != '.' )
		ptr--;

	if ( *ptr == '.' )
		return( copy_str( ptr + 1 ) );
	
	else
		return( copy_str( "" ) );
}


int
append_cmd_str( result, str, length, err, cmd_flag )
char *result;
char *str;
int length;
char *err;
int cmd_flag;
{
	char *strcat();

	char *errstr;

	int elen;
	int rlen;
	int len;
	int rc;
	int i;

	rlen = strlen( result );

	len = rlen + strlen( str ) + 1;

	if ( len <= length )
		strcat( result, str );

	else
	{
		if ( cmd_flag == TRUE )
		{
			if ( err == NULL )
			{
				printf( "\nError in Append Cmd String: Null Cmd\n\n" );

				return( TCL_ERROR );
			}

			elen = len + strlen( err ) + 5;

			errstr = (char *) MALLOC( (unsigned) elen
				* sizeof(char) );
			memcheck( errstr, "Append Cmd String Err String" );

			sprintf( errstr, "%s { %s }", err, result );

			rc = Tcl_Eval( interp, errstr );

			if ( rc == TCL_ERROR )
				printf( "\nError Executing \"%s\"\n\n", errstr );

			FREE( errstr );

			if ( strlen( str ) <= length )
				strcpy( result, str );

			else
			{
				printf( "\nError Append Cmd String \"%s\" Too Long\n\n",
					str );

				strcpy( result, "" );
			}

			return( rc );
		}

		else
		{
			if ( err != NULL )
				printf( "\nError in Append Cmd String: %s\n\n", err );

			i = rlen;

			while ( i < length - 4 )
			{
				result[i] = str[ i - rlen ];

				i++;
			}

			sprintf( (char *) (result + (length - 4)), "..." );
		}
	}

	return( TCL_OK );
}


double
distance( x1, y1, x2, y2 )
double x1, y1, x2, y2;
{
	return( sqrt(((x2 - x1) * (x2 - x1)) + ((y2 - y1) * (y2 - y1))) );
}


double
angle_of( x1, y1, x2, y2 )
double x1, y1, x2, y2;
{
	double slope;
	double at;

	/* North is heading angle 0.0 */
	/* y's grow downward from the top of the screen */

	if ( x2 == x1 )
	{
		if ( y2 > y1 )
			at = -1.0 * PI;

		else
			at = 0.0;
	}

	else
	{
		slope = (y1 - y2) / (x2 - x1);

		at = atan( slope );

		if ( x1 > x2 )
			at += PI / 2.0;

		else
			at -= PI / 2.0;
	}

	return( at );
}


int
split_line( line, av, ac )
char *line;
char ***av;
int *ac;
{
	char *ptr, *ptr2;
	char *tmp;

	char match;

	int len;
	int cnt;
	int eos;

	if ( line == NULL )
	{
		if ( av )
			*av = (char **) NULL;

		if ( ac )
			*ac = 0;

		return( FALSE );
	}

	/* Determine # of tokens */

	tmp = copy_str( line );

	len = strlen( tmp );

	if ( tmp[ len - 1 ] == '\n' )
		tmp[ len - 1 ] = '\0';

	ptr = strtok( tmp, " \t" );

	cnt = 0;

	while ( ptr != NULL )
	{
		cnt++;

		if ( *ptr == '{' || *ptr == '"' )	/* matching } */
		{
			match = ( *ptr == '{' ) ? '}' : '"';

			while ( *ptr != match && *ptr != '\0' )
				ptr++;
			
			if ( *ptr != match )
			{
				ptr++;

				while ( *ptr != match && *ptr != '\0' )
					ptr++;
			}
			
			if ( *ptr != '\0' ) ptr++;

			if ( *ptr == ' ' || *ptr == '\t' ) ptr++;

			ptr = strtok( ptr, " \t" );
		}

		else
			ptr = strtok( (char *) NULL, " \t" );
	}

	FREE( tmp );

	/* Allocate List of Ptrs */

	if ( cnt > 0 )
	{
		if ( av )
		{
			*av = (char **) MALLOC( (unsigned) cnt * sizeof(char *) );
			memcheck( *av, "String List String Pointer Array" );
		}

		if ( ac )
			*ac = cnt;
	}

	else
	{
		if ( av )
			*av = (char **) NULL;

		if ( ac )
			*ac = 0;

		return( FALSE );
	}

	/* Fill List */

	tmp = copy_str( line );

	if ( tmp[ len - 1 ] == '\n' )
		tmp[ len - 1 ] = '\0';

	ptr = strtok( tmp, " \t" );

	cnt = 0;

	while ( ptr != NULL )
	{
		if ( *ptr == '{' || *ptr == '"' )	/* matching } */
		{
			match = ( *ptr == '{' ) ? '}' : '"';

			ptr++;

			ptr2 = ptr;

			while ( *ptr2 != match && *ptr2 != '\0' )
				ptr2++;

			if ( *ptr2 != match )
			{
				*ptr2++ = ' ';

				while ( *ptr2 != match && *ptr2 != '\0' )
					ptr2++;
			}

			eos = 0;

			if ( *ptr2 == match )
				*ptr2 = '\0';
			
			else
				eos++;

			if ( av )
				(*av)[cnt] = copy_str( ptr );

			if ( !eos )
			{
				ptr2++;

				ptr = strtok( ptr2, " \t" );
			}

			else
				ptr = (char *) NULL;
		}

		else
		{
			if ( av )
				(*av)[cnt] = copy_str( ptr );

			ptr = strtok( (char *) NULL, " \t" );
		}

		cnt++;
	}

	FREE( tmp );

	return( TRUE );
}


char *
copy_str( str )
char *str;
{
	char *tmp;

	tmp = (char *) MALLOC( (unsigned) ( strlen( str ) + 1 )
		* sizeof(char) );
	memcheck( tmp, "Copy String" );

	strcpy( tmp, str );

	return( tmp );
}


int
compare( x, y )
char *x, *y;
{
	while ( *x != '\0' )
	{
		if ( *x != *y )
			return( FALSE );

		x++; y++;
	}

	return( TRUE );
}


int
filecheck( fp, name )
FILE *fp;
char *name;
{
	char msg[1024];
	char tmp[1024];

	if ( fp == NULL )
	{
		sprintf( tmp, "Error Opening File \"%s\".", name );

		sprintf( msg, "setMsg {%s}", tmp );

		Tcl_Eval( interp, msg );

		fprintf( stderr, "\n%s\n\n", tmp );

		return( FALSE );
	}

	return( TRUE );
}


void
memcheck( ptr, name )
char *ptr;
char *name;
{
	if ( ptr == NULL )
	{
		fprintf( stderr, "\nError Allocating Memory for \"%s\".\n\n",
			name );

		exit( -1 );
	}
}

