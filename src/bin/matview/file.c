
/* $Id: file.c,v 1.51 2001/03/26 19:44:43 kohl Exp $ */


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

#ifdef READ_MATLAB
#include "mat.h"
#endif

int open_mat_file( file, line, mat, is_compressed );
int read_mm_coordinate_pattern( mat, contflag, inst );
int read_coordinate( mat, contflag, inst );
int parse_format_spec( fmt, label, num, size, factor );
int realloc_mat( nelems, has_values );

/* Pointer to first (header) line of currently open matrix file */
char *firstline;


int
read_mat_file( file )
char *file;
{
#ifdef READ_MATLAB
	MATFile *pmat;
#endif

	FILE *mat = (FILE *) NULL;

#ifdef FILE_BUFFERING
	char *line;
#else
	char line[1024];
#endif

	char cmd[255];

	char **av;
	int ac;

	int is_compressed;
	int inst;
	int ret;

	/* Bump & Grab File Instance Number */

	inst = ++(MAIN_VIEW->fileinst);

	sprintf( cmd, "setMsg {Reading Matrix File \"%s\"...}", file );
	Tcl_Eval( interp, cmd );
	INST_CHECK_RETURN( MAIN_VIEW, inst, MAIN_VIEW->fileinst,
		FALSE, FALSE );

	update_progress( MAIN_VIEW, 0 );
	INST_CHECK_RETURN( MAIN_VIEW, inst, MAIN_VIEW->fileinst,
		FALSE, FALSE );

#ifdef READ_MATLAB

	/* Check for MatLab Files */

	pmat = matOpen( file, "r" );

	if ( pmat != NULL )
	{
		ret = read_matlab_mat( pmat, FALSE, inst );

		matClose( pmat );

		if ( ret )
			return( ret );
	}

#endif

	/* Open Matrix File and Read Header Line */

	if ( !open_mat_file( file, line, &mat, &is_compressed ) )
		return( FALSE );

	/* Parse Header */

	if ( !split_line( line, &av, &ac ) )
	{
		( is_compressed ) ? pclose( mat ) : fclose( mat );

		Tcl_Eval( interp,
			"setMsg {Error Parsing Matrix File Header.}" );

		return( FALSE );
	}

	/* Determine File Type */

	ret = -1;

	/* Matrix Market */

	if ( ac > 1 && !strcmp( av[0], "%%MatrixMarket" ) )
	{
		/* Sparse Coordinate Matrix */

		if ( ac >= 5
			&& !strcmp( av[1], "matrix" )
			&& !strcmp( av[2], "coordinate" ) )
		{
			if ( !strcmp( av[3], "pattern" ) )
				ret = read_mm_coordinate_pattern( mat, FALSE, inst );
		
			else if ( !strcasecmp( av[3], "real" )
				|| !strcmp( av[3], "integer" ) )
			{
				ret = read_mm_coordinate_real( mat, FALSE, inst );
			}

			else if ( !strcasecmp( av[3], "complex" ) )
				ret = read_mm_coordinate_complex( mat, FALSE, inst );

			if ( ret )
			{
				if ( !strcmp( av[4], "general" ) )
					MAT->symm = MS_UNSYMMETRIC;

				else if ( !strcmp( av[4], "symmetric" ) )
					MAT->symm = MS_SYMMETRIC;

				else if ( !strcmp( av[4], "skew-symmetric" ) )
					MAT->symm = MS_SKEW_SYMMETRIC;

				else if ( !strcmp( av[4], "Hermitian" ) )
					MAT->symm = MS_HERMITIAN;
			}
		}
	}

	else
	{
		/* Try Harwell-Boeing */

		ret = read_harwell_boeing( mat, FALSE, inst );

		/* Try Generic Coordinate Format */

		if ( !ret )
		{
			/* Need to re-open file after read_harwell_boeing()... */

			if ( !open_mat_file( file, line, &mat, &is_compressed ) )
				return( FALSE );

			ret = read_coordinate( mat, FALSE, inst );
		}

		/* Tried Everything We Know... */

		if ( !ret )
			ret = -1;
	}

	if ( MAT->nelems == MAT->alloc_size )
		( is_compressed ) ? pclose( mat ) : fclose( mat );

	MAT->is_compressed = is_compressed;

	if ( ret < 0 )
	{
		Tcl_Eval( interp,
			"setMsg {Error:  Unknown or Unsupported File Format...}" );

		ret = FALSE;
	}

	return( ret );
}


int
read_mm_coordinate_pattern( mat, contflag, inst )
FILE *mat;
int contflag;
int inst;
{
	return( read_mm_coordinate( mat, MATRIX_PATTERN, contflag, inst ) );
}


int
read_mm_coordinate_real( mat, contflag, inst )
FILE *mat;
int contflag;
int inst;
{
	return( read_mm_coordinate( mat, MATRIX_REAL, contflag, inst ) );
}


int
read_mm_coordinate_complex( mat, contflag, inst )
FILE *mat;
int contflag;
int inst;
{
	return( read_mm_coordinate( mat, MATRIX_COMPLEX, contflag, inst ) );
}


int
read_mm_coordinate( mat, has_values, contflag, inst )
FILE *mat;
int has_values;
int contflag;
int inst;
{
	static int cnt;

	TEMPLATE T;

	COMPLEX C;

	DATA D;

#ifdef FILE_BUFFERING
	char *line;
#else
	char line[1024];
#endif

	char *ptr;

	double jval;
	double val;

	int xmax, ymax;
	int progress;
	int istart;
	int nelems;
	int x, y;
	int err;
	int num;
	int tmp;
	int ok;
	int i;

	/* Check Whether Re-Entered */

	if ( !contflag )
	{
		/* Read Matrix Bounds & Number of Elements */

		cnt = 1;  /* for header */

		ok = 0;

		do
		{
			cnt++;

#ifdef FILE_BUFFERING
			if ( (line = read_line( mat )) == NULL )
#else
			if ( fgets( line, 1024, mat ) == NULL )
#endif
			{
				printf( "\nError Reading Matrix Market File.\n\n" );

				return( FALSE );
			}

			if ( line[0] != '%' && line[0] != '\0'
				&& split_line( line, (char ***) NULL, (int *) NULL ) )
			{
				if ( sscanf( line, "%d %d %d\n",
					&ymax, &xmax, &nelems ) != 3 )
				{
					printf( "\nError Parsing %s Matrix Dimensions.\n\n",
						"Matrix Market" );

					return( FALSE );
				}

				else
					ok++;
			}
		}
		while ( !ok );

		if ( TRANSPOSE_MATRIX )
		{
			tmp = xmax;
			xmax = ymax;
			ymax = tmp;
		}

		/* Set Up Matrix */

		free_mat();

		if ( !alloc_mat( nelems, has_values ) )
			return( FALSE );

		MAT->glb[0] = 1;
		MAT->glb[1] = 1;

		MAT->gub[0] = xmax;
		MAT->gub[1] = ymax;

		MAT->dim = 2;
	}

	/* Calculate Progress Factor */

	progress = MAT->nelems / PROGRESS_TICKS;
	progress = ( progress ) ? progress : 1;

	/* Read In Data */

	if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
		C = MAT->complex;
	
	else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
		D = MAT->data;
	
	else
		T = MAT->template;

	istart = ( !contflag ) ? 0 : MAT->last_loaded + 1;

	num = 0;

	for ( i=istart ; i < MAT->nelems ; i++ )
	{
		if ( !( i % progress ) )
		{
			update_progress( MAIN_VIEW, i / progress );
			INST_CHECK_RETURN( MAIN_VIEW, inst, MAIN_VIEW->fileinst,
				FALSE, FALSE );
		}

		/* Check for Partial Load Overflow */

		if ( i - istart >= MAT->alloc_size )
		{
			MAT->fpsave = mat;

			if ( has_values == MATRIX_COMPLEX )
				MAT->read_mat_func = read_mm_coordinate_complex;
	
			else if ( has_values == MATRIX_REAL )
				MAT->read_mat_func = read_mm_coordinate_real;
	
			else
				MAT->read_mat_func = read_mm_coordinate_pattern;

			MAT->load_inst = inst;

			MAT->first_loaded = istart;
			MAT->last_loaded = MAT->first_loaded + MAT->alloc_size - 1;

			return( TRUE );
		}

		/* Read Valid Line */

		ok = 0;

		do
		{
			cnt++;

#ifdef FILE_BUFFERING
			if ( (line = read_line( mat )) == NULL )
#else
			if ( fgets( line, 1024, mat ) == NULL )
#endif
			{
				printf(
					"\nError Reading Matrix Market File (line %d).\n\n",
					cnt );

				free_mat();

				return( FALSE );
			}

			if ( line[0] != '%' && line[0] != '\0'
				&& split_line( line, (char ***) NULL, (int *) NULL ) )
			{
				ok++;
			}
		}
		while ( !ok );

		/* Munge Goofy nnnD[+-]nn Fortran Exponent Format...  :-Q */

		ptr = line;

		while ( *ptr != '\0' )
		{
			if ( ( *ptr == 'D' || *ptr == 'd' )
					&& ( *(ptr+1) == '+' || *(ptr+1) == '-' ) )
				*ptr = 'E';
			ptr++;
		}

		/* Parse Out Coords & Any Data */

		err = 0;

		switch ( has_values )
		{
			case MATRIX_COMPLEX:
			{
				if ( sscanf( line, "%d %d %lf %lf\n",
						&y, &x, &val, &jval ) != 4 )
					err++;

				break;
			}

			case MATRIX_REAL:
			{
				if ( sscanf( line, "%d %d %lf\n", &y, &x, &val ) != 3 )
					err++;

				break;
			}

			case MATRIX_PATTERN:
			default:
			{
				if ( sscanf( line, "%d %d\n", &y, &x ) != 2 )
					err++;

				break;
			}
		}

		if ( err )
		{
			printf(
				"\nError Parsing Matrix Market File (line %d).\n\n",
				cnt );

			free_mat();

			return( FALSE );
		}

		/* Save Data in Global Structure */

		if ( TRANSPOSE_MATRIX )
		{
			tmp = x;
			x = y;
			y = tmp;
		}

		if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
		{
			C->x = x;
			C->y = y;

			C->val = val;
			C->jval = jval;

			C++;
		}

		else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
		{
			D->x = x;
			D->y = y;

			D->val = val;

			D++;
		}

		else
		{
			T->x = x;
			T->y = y;

			T++;
		}

		num++;
	}

	MAT->fpsave = mat;

	if ( has_values == MATRIX_COMPLEX )
		MAT->read_mat_func = read_mm_coordinate_complex;
	
	else if ( has_values == MATRIX_REAL )
		MAT->read_mat_func = read_mm_coordinate_real;
	
	else
		MAT->read_mat_func = read_mm_coordinate_pattern;

	MAT->load_inst = inst;

	MAT->first_loaded = istart;
	MAT->last_loaded = MAT->first_loaded + num - 1;

	return( TRUE );
}


int
read_harwell_boeing( mat, contflag, inst )
FILE *mat;
int contflag;
int inst;
{
	static int has_values;

	static char value[80];

	static char valfmt[20];

	static int *colptrs;
	static int ncolptrs;

	static int *rowinds;
	static int nrowinds;

	static int valcrd;

	static int start, stop;
	static int line_offset;
	static int progress;
	static int r, c;
	static int cnt;

	TEMPLATE T;

	COMPLEX C;

	DATA D;

	char line[82];
	char fmt[32];

	char title[73];
	char type[4];
	char key[9];

	char ptrfmt[20];
	char indfmt[20];
	char rhsfmt[20];

	char eolch;

	int *iptr;

	int totcrd;
	int ptrcrd;
	int indcrd;
	int rhscrd;
	int nnzero;
	int neltvl;


	int nrow, ncol;
	int cstart;
	int factor;
	int nelems;
	int nlines;
	int nvals;
	int i, j;
	int x, y;
	int bail;
	int done;
	int size;
	int num;
	int tmp;
	int r_c;

	/* Check Whether Re-Entered */

	if ( !contflag )
	{
		/* Extract Title & Key */

		strncpy( title, firstline, 72 );
		title[72] = '\0';

		strncpy( key, firstline + 72, 8 );
		key[8] = '\0';

		/* Read Numbers of Lines */

		if ( fscanf( mat, "%14d%14d%14d%14d%14d\n",
				&totcrd, &ptrcrd, &indcrd, &valcrd, &rhscrd ) != 5 )
			return( FALSE );

		/* Read Matrix Type & Size */

		if ( fscanf( mat, "%3c%*11c%14d%14d%14d",
				type, &nrow, &ncol, &nnzero ) != 4 )
			return( FALSE );

		nelems = nnzero;

		/* Try to Read NELTVL Value... */
		/* (O.K. if omitted, usual practice.) */

		if ( fscanf( mat, "%14d\n", &neltvl ) != 1 )
			neltvl = 0;

		/* Check Matrix Type */

		if ( type[2] == 'E' )
		{
			printf( "\nError: Parsing of %s Not Yet Supported...\n\n",
				"Elemental (Unassembled) Harwell-Boeing Matrix Files" );

			return( FALSE );
		}

		/* Set Matrix Values */

		switch ( type[0] )
		{
			case 'C':	has_values = MATRIX_COMPLEX;	break;
			case 'R':	has_values = MATRIX_REAL;		break;
			case 'P':	has_values = MATRIX_PATTERN;	break;

			default:
			{
				printf( "\nError: Unknown %s, '%c' != %s.  %s.\n\n",
					"Harwell-Boeing Matrix Values Type", type[0],
					"(R,C,P)", "Parse Failed" );

				return( FALSE );
			}
		}

		/* Set Matrix Symmetry */

		switch ( type[1] )
		{
			case 'S':	MAT->symm = MS_SYMMETRIC;		break;
			case 'U':	MAT->symm = MS_UNSYMMETRIC;		break;
			case 'H':	MAT->symm = MS_HERMITIAN;		break;
			case 'Z':	MAT->symm = MS_SKEW_SYMMETRIC;	break;
			case 'R':	MAT->symm = MS_UNSYMMETRIC;		break;

			default:
			{
				printf( "\nError: Unknown %s, '%c' != %s.  %s.\n\n",
					"Harwell-Boeing Matrix Structure Type", type[0],
					"(S,U,H,Z,R)", "Parse Failed" );

				return( FALSE );
			}
		}

		/* Read Formats */

		if ( fscanf( mat, "%16c%16c", ptrfmt, indfmt ) != 2 )
		{
			printf( "\nError Reading %s.  Parse Failed.\n\n",
				"Harwell-Boeing Format for Pointer or Row Indices" );

			return( FALSE );
		}

		if ( valcrd > 0 && fscanf( mat, "%20c", valfmt ) != 1 )
		{
			printf( "\nError Reading %s.  Parse Failed.\n\n",
			"Harwell-Boeing Format for Values of Coefficient Matrix" );

			return( FALSE );
		}

		if ( rhscrd > 0 && fscanf( mat, "%20c", rhsfmt ) != 1 )
		{
			printf( "\nError Reading %s.  Parse Failed.\n\n",
				"Harwell-Boeing Format for Right-Hand Sides" );

			return( FALSE );
		}

		while ( (eolch = fgetc( mat )) != '\n' && eolch != (char) EOF );

		/* Now Parsing Line... */

		cnt = 5;

		/* Ignore RHS Details */

		if ( rhscrd > 0 )
		{
			if ( fgets( line, 81, mat ) == NULL )
			{
				printf(
					"\nError: Premature End of File %s (line %d).\n\n",
					"Skipping Harwell-Boeing Right-Hand Side Values",
					cnt );

				return( FALSE );
			}

			cnt++;
		}

		/* Calculate Progress Factor */

		if ( has_values != MATRIX_PATTERN )
		{
			progress = ( cnt + ptrcrd + indcrd + valcrd )
				/ PROGRESS_TICKS;
			progress = ( progress ) ? progress : 1;
		}

		else
		{
			progress = ( cnt + ptrcrd + indcrd + ncol )
				/ PROGRESS_TICKS;
			progress = ( progress ) ? progress : 1;
		}

		/* Read In Column Pointers */

		if ( !parse_format_spec( ptrfmt,
				"Harwell-Boeing Column Pointer",
				&num, &size, &factor ) )
			return( FALSE );

		sprintf( fmt, "%%%dc", size );

		ncolptrs = ptrcrd * num;

		if ( !ncolptrs )
		{
			printf(
			"\nError Parsing %s File: %s (PTRCRD=%d num/line=%d).\n\n",
				"Harwell-Boeing", "No Column Pointers?", ptrcrd, num );

			return( FALSE );
		}

		colptrs = (int *) MALLOC( (unsigned) ncolptrs * sizeof( int ) );
		memcheck( colptrs, "Harwell-Boeing Column Pointers" );

		iptr = colptrs;

		done = 0;

		for ( i=0 ; i < ptrcrd ; i++ )
		{
			if ( !( cnt % progress ) )
			{
				update_progress( MAIN_VIEW, cnt / progress );
				INST_CHECK_RETURN( MAIN_VIEW, inst,
					MAIN_VIEW->fileinst, FALSE, FALSE );
			}

			for ( j=0 ; j < num && !done ; j++ )
			{
				if ( fscanf( mat, fmt, value ) != 1 )
				{
					printf( "\nError Reading %s (line %d).\n\n",
						"Harwell-Boeing Column Pointer Field", cnt );

					if ( i || j )
						printf( "Last Value Read:  %s\n\n", value );

					return( FALSE );
				}

				if ( !is_bogus_value( mat, value, size, &bail ) )
				{
					value[ size ] = '\0';

					*iptr++ = factor * atoi( value );
				}

				else
				{
					if ( bail )
					{
						printf(
						"\nError Parsing %s File: %s (line %d).\n\n",
							"Harwell-Boeing",
							"Premature End of Column Pointer Values",
							cnt );

						return( FALSE );
					}

					done++;
				}
			}

			if ( !done )
			{
				while ( (eolch = fgetc( mat )) != '\n'
					&& eolch != (char) EOF );
			
				if ( eolch == (char) EOF )
				{
					printf(
						"\nError Parsing %s File: %s (line %d).\n\n",
						"Harwell-Boeing",
						"Premature End of Column Pointer Values", cnt );

					return( FALSE );
				}
			}

			else if ( i != ptrcrd - 1 )
			{
				printf(
				"\nError Parsing %s File: %s - %d of %d (line %d).\n\n",
					"Harwell-Boeing", "Too Few Column Pointer Lines",
					i + 1, ptrcrd, cnt );

				return( FALSE );
			}

			cnt++;
		}

		/* Read In Row Indices */

		if ( !parse_format_spec( indfmt, "Harwell-Boeing Row Indices",
				&num, &size, &factor ) )
			return( FALSE );

		sprintf( fmt, "%%%dc", size );

		nrowinds = indcrd * num;

		if ( !nrowinds )
		{
			printf(
			"\nError Parsing %s File: %s (INDCRD=%d num/line=%d).\n\n",
				"Harwell-Boeing", "No Row Indices?", indcrd, num );

			return( FALSE );
		}

		rowinds = (int *) MALLOC( (unsigned) nrowinds * sizeof( int ) );
		memcheck( rowinds, "Harwell-Boeing Column Pointers" );

		iptr = rowinds;

		done = 0;

		for ( i=0 ; i < indcrd ; i++ )
		{
			if ( !( cnt % progress ) )
			{
				update_progress( MAIN_VIEW, cnt / progress );
				INST_CHECK_RETURN( MAIN_VIEW, inst,
					MAIN_VIEW->fileinst, FALSE, FALSE );
			}

			for ( j=0 ; j < num && !done ; j++ )
			{
				if ( fscanf( mat, fmt, value ) != 1 )
				{
					if ( has_values != MATRIX_PATTERN
						|| i != indcrd - 1 )
					{
						printf( "\nError Reading %s %s (line %d).\n\n",
							"Harwell-Boeing", "Row Index Field", cnt );

						if ( i || j )
							printf( "Last Value Read:  %s\n\n", value );

						return( FALSE );
					}

					done++;

					break;
				}

				if ( !is_bogus_value( mat, value, size, &bail ) )
				{
					value[ size ] = '\0';

					*iptr++ = factor * atoi( value );
				}

				else
				{
					if ( bail &&
							( has_values != MATRIX_PATTERN
								|| i != indcrd - 1 ) )
					{
						printf(
						"\nError Parsing %s File: %s (line %d).\n\n",
							"Harwell-Boeing",
							"Premature End of Row Index Values", cnt );

						return( FALSE );
					}

					done++;
				}
			}

			if ( !done )
			{
				while ( (eolch = fgetc( mat )) != '\n'
					&& eolch != (char) EOF );
			
				if ( eolch == (char) EOF )
				{
					printf(
						"\nError Parsing %s File: %s (line %d).\n\n",
						"Harwell-Boeing",
						"Premature End of Row Index Values", cnt );

					return( FALSE );
				}
			}

			else if ( i != indcrd - 1 )
			{
				printf(
				"\nError Parsing %s File: %s - %d of %d (line %d).\n\n",
					"Harwell-Boeing", "Too Few Row Index Lines",
					i + 1, indcrd, cnt );

				return( FALSE );
			}

			cnt++;
		}

		/* Set Up Matrix */

		free_mat();

		if ( !alloc_mat( nelems, has_values ) )
			return( FALSE );

		MAT->glb[0] = 1;
		MAT->glb[1] = 1;

		MAT->gub[0] = ncol;
		MAT->gub[1] = nrow;

		MAT->dim = 2;
	}

	/* Set Non-Zero Values */

	if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
		C = MAT->complex;

	else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
		D = MAT->data;

	else
		T = MAT->template;

	/* Handle Template Versus Real Data */

	if ( has_values != MATRIX_PATTERN )
	{
		/* Matrix Data Format */

		if ( !parse_format_spec( valfmt, "Harwell-Boeing Data Values",
				&num, &size, &factor ) )
			return( FALSE );

		sprintf( fmt, "%%%dc", size );

		if ( !contflag )
		{
			/* Initialize Row & Column Indices */

			c = -1;

			do
			{
				c++;

				if ( c < ncolptrs )
				{
					start = colptrs[ c ];
					stop = colptrs[ c + 1 ];
				}
			}
			while ( start == stop && c < ncolptrs );

			r = start;

			/* Initialize Line Offset to Zero */

			line_offset = 0;
		}

		/* Crunch Number of Matrix Data Lines to Read */

		if ( MAT->alloc_size == MAT->nelems )
			nlines = valcrd;

		else
		{
			if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
				nlines = ( 2 * MAT->alloc_size ) / num;

			else
				nlines = MAT->alloc_size / num;
		}

		/* Read in Data */

		done = 0;

		r_c = 0;

		for ( i=0 ; i < nlines && !done ; i++ )
		{
			if ( !( cnt % progress ) )
			{
				update_progress( MAIN_VIEW, cnt / progress );
				INST_CHECK_RETURN( MAIN_VIEW, inst,
					MAIN_VIEW->fileinst, FALSE, FALSE );
			}

			for ( j=0 ; j < num && !done ; j++ )
			{
				if ( fscanf( mat, fmt, value ) != 1 )
				{
					if ( line_offset + i != valcrd - 1 )
					{
						printf( "\nError Reading %s (line %d).\n\n",
							"Harwell-Boeing Data Value Field", cnt );

						if ( line_offset + i || j )
							printf( "Last Value Read:  %s\n\n", value );

						return( FALSE );
					}

					done++;

					break;
				}

				if ( !is_bogus_value( mat, value, size, &bail ) )
				{
					if ( !r_c )
					{
						x = c + 1;
						y = rowinds[ r - 1 ];

						if ( TRANSPOSE_MATRIX )
						{
							tmp = x;
							x = y;
							y = tmp;
						}

						if ( NOLOAD_VALUES )
						{
							T->x = x;
							T->y = y;

							T++;
						}

						else if ( has_values == MATRIX_REAL )
						{
							value[ size ] = '\0';

							D->x = x;
							D->y = y;

							D->val = (double) factor * atof( value );

							D++;
						}

						else if ( has_values == MATRIX_COMPLEX )
						{
							value[ size ] = '\0';

							C->x = x;
							C->y = y;

							C->val = (double) factor * atof( value );

							/* Don't Increment C, Wait for Im Part */
						}

						/* Bump Up Row & Col Indices */

						if ( ++r == stop )
						{
							do
							{
								c++;

								if ( c < ncolptrs )
								{
									start = colptrs[ c ];
									stop = colptrs[ c + 1 ];
								}
							}
							while ( start == stop && c < ncolptrs );

							r = start;
						}
					}

					else if ( has_values == MATRIX_COMPLEX )
					{
						value[ size ] = '\0';

						C->jval = (double) factor * atof( value );

						C++;
					}

					/* Flip Real vs. Complex Flag */

					if ( has_values == MATRIX_COMPLEX )
						r_c = 1 - r_c;
				}

				else
				{
					if ( bail && line_offset + i != valcrd - 1 )
					{
						printf(
						"\nError Parsing %s File: %s (line %d).\n\n",
							"Harwell-Boeing",
							"Premature End of Data Values", cnt );

						return( FALSE );
					}

					done++;
				}
			}

			if ( !done )
			{
				while ( (eolch = fgetc( mat )) != '\n'
					&& eolch != (char) EOF );
			
				if ( eolch == (char) EOF )
				{
					printf(
						"\nError Parsing %s File: %s (line %d).\n\n",
						"Harwell-Boeing",
						"Premature End of Data Values", cnt );

					return( FALSE );
				}
			}

			else if ( line_offset + i != valcrd - 1 )
			{
				printf(
				"\nError Parsing %s File: %s - %d of %d (line %d).\n\n",
					"Harwell-Boeing", "Too Few Data Value Lines",
					line_offset + i + 1, valcrd, cnt );

				return( FALSE );
			}

			cnt++;
		}

		line_offset = nlines;

		nvals = num * nlines;
	}

	else
	{
		cstart = ( !contflag ) ? 0 : c;

		num = 0;

		for ( c=cstart ; c < ncol ; c++ )
		{
			if ( !( cnt % progress ) )
			{
				update_progress( MAIN_VIEW, cnt / progress );
				INST_CHECK_RETURN( MAIN_VIEW, inst,
					MAIN_VIEW->fileinst, FALSE, FALSE );
			}

			start = colptrs[ c ];
			stop = colptrs[ c + 1 ];

			/* Check for Partial Load Overflow */

			if ( num + ( stop - start ) > MAT->alloc_size )
			{
				MAT->fpsave = mat;

				MAT->read_mat_func = read_harwell_boeing;

				MAT->load_inst = inst;

				MAT->first_loaded = ( MAT->last_loaded )
					? ( MAT->last_loaded + 1 ) : 0;
				MAT->last_loaded = MAT->first_loaded + num - 1;

				return( TRUE );
			}

			if ( start != stop )
			{
				for ( r=start ; r < stop ; r++ )
				{
					x = c + 1;
					y = rowinds[ r - 1 ];

					if ( TRANSPOSE_MATRIX )
					{
						tmp = x;
						x = y;
						y = tmp;
					}

					T->x = x;
					T->y = y;

					num++;

					T++;
				}
			}

			cnt++;
		}

		nvals = num;
	}

	/* Done */

	MAT->fpsave = mat;

	MAT->read_mat_func = read_harwell_boeing;

	MAT->load_inst = inst;

	MAT->first_loaded = ( MAT->last_loaded )
		? ( MAT->last_loaded + 1 ) : 0;
	MAT->last_loaded = MAT->first_loaded + nvals - 1;

	return( TRUE );
}


int
parse_format_spec( fmt, label, num, size, factor )
char *fmt;
char *label;
int *num;
int *size;
int *factor;
{
	char *ptr;

	char type;

	int val;

	/* Skip White Space */

	ptr = fmt;

	while ( *ptr == ' ' || *ptr == '\t' )
		ptr++;

	/* Skip Opening Parentheses */

	if ( *ptr == '(' /* matching ) */ )
		ptr++;

	/* Default Value Factor */

	*factor = 1;

	/* Parse Number of Elements and Type */

	if ( sscanf( ptr, "%d%1c", &val, &type ) != 2 )
	{
		printf( "\nError Parsing %s Format:  %s  (%s).\n\n",
			label, "Number of Elements & Type", fmt );

		return( FALSE );
	}

	/* Check for Goofy Fortran 'P' Scaling Factor */

	if ( type == 'P' || type == 'p' )
	{
		/* Calculate Data Value Factor */

		*factor = (int) pow( 10.0, (double) val );

		/* Advance Past Type */

		while ( *ptr != type )
			ptr++;

		ptr++;

		if ( *ptr == ',' )
			ptr++;

		/* Read Actual Number of Elements and Type */

		if ( sscanf( ptr, "%d%1c", &val, &type ) != 2 )
		{
			printf( "\nError Parsing %s Format:  %s  (%s).\n\n",
				label, "Number of Elements & Type", fmt );

			return( FALSE );
		}
	}

	*num = val;

	/* Advance Past Type */

	while ( *ptr != type )
		ptr++;

	ptr++;

	/* Parse Format String */

	if ( sscanf( ptr, "%d", size ) != 1 )
	{
		printf( "\nError Parsing %s Format:  Field Size  (%s).\n\n",
			label, fmt );

		return( FALSE );
	}

	return( TRUE );
}


int
is_bogus_value( mat, value, size, really_bogus )
FILE *mat;
char *value;
int size;
int *really_bogus;
{
	char eolch;

	int non_white;
	int i, j;

	/* Initialize "Really Bogus" Flag */

	if ( really_bogus != NULL )
		*really_bogus = FALSE;

	/* Check Through String */

	non_white = 0;

	for ( i=0 ; i < size ; i++ )
	{
		/* Premature End of Line */

		if ( value[i] == '\n' )
		{
			/* Push Rest of Last Read Back Onto Stream */

			for ( j=(size - 1) ; j > i ; j-- )
			{
				if ( ungetc( value[j], mat ) == EOF )
				{
					if ( really_bogus != NULL )
						*really_bogus = TRUE;

					return( TRUE );
				}
			}

			return( TRUE );
		}

		/* O.K., Non White Space Character */

		if ( value[i] != ' ' && value[i] != '\t' )
			non_white++;

		/* Munge Goofy nnnD[+-]nn Fortran Exponent Format...  :-Q */

		if ( ( value[i] == 'D' || value[i] == 'd' )
				&& ( value[i+1] == '+' || value[i+1] == '-' ) )
			value[i] = 'E';
	}

	/* Return If Something There... */

	if ( non_white )
		return( FALSE );

	/* If Nothing There, Clear Off Rest of Line */

	while ( (eolch = fgetc( mat )) != '\n' && eolch != (char) EOF );

	if ( eolch == (char) EOF && really_bogus != NULL )
		*really_bogus = TRUE;

	return( TRUE );
}


int
read_coordinate( mat, contflag, inst )
FILE *mat;
int contflag;
int inst;
{
	static int has_values;
	static int nelems;
	static int cnt;

	TEMPLATE T;

	COMPLEX C;

	DATA D;

#ifdef FILE_BUFFERING
	char *line;
#else
	char line[1024];
#endif

	char *lptr;

	double a, b, c, d;

	double jval;
	double val;

	int xmax, ymax;
	int progress;
	int is_first;
	int istart;
	int done;
	int x, y;
	int err;
	int num;
	int tmp;
	int ok;
	int i;

	/* Check Whether Re-Entered */

	is_first = 0;

	if ( !contflag )
	{
		/* Initialize Matrix Info  to Defaults */

		xmax = ymax = -1;
		has_values = -1;
		nelems = -1;

		/* Parse Lines Until You Hit Data */

		cnt = 1;	/* first line */

		lptr = firstline;

		do
		{
			num = sscanf( lptr, "%lf %lf %lf %lf", &a, &b, &c, &d );
		
			switch ( num )
			{
				case 0:
					break;

				case 1:
					nelems = (int) a;
					break;

				case 2:
					has_values = MATRIX_PATTERN;
					break;

				case 3:
					has_values = MATRIX_REAL;
					break;

				case 4:
					has_values = MATRIX_COMPLEX;
					break;
			}

			/* Try Next Line */

			if ( has_values < 0 )
			{
				/* Read Valid Line */

				ok = 0;

				do
				{
					cnt++;

#ifdef FILE_BUFFERING
					if ( (line = read_line( mat )) == NULL )
#else
					if ( fgets( line, 1024, mat ) == NULL )
#endif
					{
						printf(
						"\nError Reading %s Matrix File (line %d).\n\n",
							"Coordinate", cnt );

						return( FALSE );
					}

					if ( line[0] != '%' && line[0] != '\0'
						&& split_line( line,
							(char ***) NULL, (int *) NULL ) )
					{
						ok++;
					}
				}
				while ( !ok );

				lptr = line;
			}
		}
		while ( has_values < 0 );

		/* Verify # of Elements (i.e. Guess) */

		tmp = ( nelems < 0 ) ? 10000 : nelems;

		/* Set Up Matrix */

		free_mat();

		if ( !alloc_mat( tmp, has_values ) )
			return( FALSE );

		/* Set Up First Element */

		if ( TRANSPOSE_MATRIX )
		{
			x = (int) b;
			y = (int) a;
		}

		else
		{
			x = (int) a;
			y = (int) b;
		}

		xmax = x;
		ymax = y;

		is_first++;
	}

	else
	{
		xmax = MAT->gub[0];
		ymax = MAT->gub[1];
	}

	/* Calculate Progress Factor */

	progress = MAT->nelems / PROGRESS_TICKS;
	progress = ( progress ) ? progress : 1;

	/* Read In Data */

	if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
	{
		C = MAT->complex;

		if ( is_first )
		{
			C->x = x;
			C->y = y;

			C->val = c;
			C->jval = d;

			C++;
		}
	}
	
	else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
	{
		D = MAT->data;

		if ( is_first )
		{
			D->x = x;
			D->y = y;

			D->val = c;

			D++;
		}
	}
	
	else
	{
		T = MAT->template;

		if ( is_first )
		{
			T->x = x;
			T->y = y;

			T++;
		}
	}

	istart = ( !contflag ) ? 0 : MAT->last_loaded + 1;

	i = istart;

	num = 0;

	if ( is_first )
	{
		i++;

		num++;
	}

	done = 0;

	do
	{
		if ( !( i % progress ) )
		{
			update_progress( MAIN_VIEW, i / progress );
			INST_CHECK_RETURN( MAIN_VIEW, inst, MAIN_VIEW->fileinst,
				FALSE, FALSE );
		}

		/* Check for Partial Load Overflow */

		if ( i - istart >= MAT->alloc_size )
		{
			/* Make Sure Actual Allocation Problem... */

			if ( MAT->alloc_size != MAT->nelems )
			{
				/* Set Temporary Matrix Info */

				MAT->glb[0] = 1;
				MAT->glb[1] = 1;

				MAT->gub[0] = xmax;
				MAT->gub[1] = ymax;

				MAT->dim = 2;

				/* Set Partial Read Info */

				MAT->fpsave = mat;

				MAT->read_mat_func = read_coordinate;
	
				MAT->load_inst = inst;

				MAT->first_loaded = istart;
				MAT->last_loaded =
					MAT->first_loaded + MAT->alloc_size - 1;

				return( TRUE );
			}

			/* We Guessed Wrong - Re-Allocate...  :-Q */

			else
			{
				/* Double Current Allocation */

				tmp = 2 * MAT->nelems;

				if ( !realloc_mat( tmp, has_values ) )
					return( FALSE );

				/* Reset Data Pointer (D-Ohhhh! :-Q) */

				if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
					C = MAT->complex + ( i - istart );
	
				else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
					D = MAT->data + ( i - istart );
	
				else
					T = MAT->template + ( i - istart );

				/* Re-Calculate Progress Factor */

				progress = MAT->nelems / PROGRESS_TICKS;
				progress = ( progress ) ? progress : 1;
			}
		}

		/* Read Valid Line */

		ok = 0;

		do
		{
#ifdef FILE_BUFFERING
			if ( (line = read_line( mat )) == NULL )
#else
			if ( fgets( line, 1024, mat ) == NULL )
#endif
			{
				/* Premature EOF */

				if ( cnt < nelems )
				{
					printf(
						"\nError Reading %s Matrix File (line %d).\n\n",
						"Coordinate", cnt );

					free_mat();

					return( FALSE );
				}

				/* Finally Found EOF :-) */

				else
				{
					nelems = i;

					done++;
				}
			}

			if ( done
				|| ( line[0] != '%' && line[0] != '\0'
					&& split_line( line,
						(char ***) NULL, (int *) NULL ) ) )
			{
				ok++;
			}

			cnt++;
		}
		while ( !ok );

		/* Parse Out Coords & Any Data */

		if ( !done )
		{
			err = 0;

			switch ( has_values )
			{
				case MATRIX_COMPLEX:
				{
					if ( sscanf( line, "%d %d %lf %lf\n",
							&y, &x, &val, &jval ) != 4 )
						err++;

					break;
				}

				case MATRIX_REAL:
				{
					if ( sscanf( line, "%d %d %lf\n",
							&y, &x, &val ) != 3 )
						err++;

					break;
				}

				case MATRIX_PATTERN:
				default:
				{
					if ( sscanf( line, "%d %d\n", &y, &x ) != 2 )
						err++;

					break;
				}
			}

			if ( err )
			{
				printf( "\nError Parsing %s Matrix File (line %d).\n\n",
					"Coordinate", cnt );

				free_mat();

				return( FALSE );
			}

			/* Save Data in Global Structure */

			if ( TRANSPOSE_MATRIX )
			{
				tmp = x;
				x = y;
				y = tmp;
			}

			if ( x > xmax ) xmax = x;
			if ( y > ymax ) ymax = y;

			if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
			{
				C->x = x;
				C->y = y;

				C->val = val;
				C->jval = jval;

				C++;
			}

			else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
			{
				D->x = x;
				D->y = y;

				D->val = val;

				D++;
			}

			else
			{
				T->x = x;
				T->y = y;

				T++;
			}

			if ( ++i >= nelems && nelems > 0 )
				done++;

			num++;
		}
	}
	while ( !done );

	/* Set Actual # of Elements */

	MAT->nelems = nelems;

	/* Set Matrix Info */

	MAT->glb[0] = 1;
	MAT->glb[1] = 1;

	MAT->gub[0] = xmax;
	MAT->gub[1] = ymax;

	MAT->dim = 2;

	/* Set Partial Read Info */

	MAT->fpsave = mat;

	MAT->read_mat_func = read_coordinate;
	
	MAT->load_inst = inst;

	MAT->first_loaded = istart;
	MAT->last_loaded = MAT->first_loaded + num - 1;

	return( TRUE );
}


#ifdef READ_MATLAB

int
read_matlab_mat( pmat, contflag, inst )
MATFile *pmat;
int contflag;
int inst;
{
	static mxArray *pa;

	static double *idata;
	static double *data;

	static int *ir, *jc;

	static int ncol;
	static int nrow;

	static int has_values;
	static int c, r;
	static int cnt;

	TEMPLATE T;

	COMPLEX C;

	DATA D;

	int start, stop;
	int progress;
	int cstart;
	int nelems;
	int x, y;
	int num;
	int tmp;
	int i;

	const int *dims;
	int ndim;

	if ( !contflag )
	{
		/* Get Matrix */

		pa = matGetNextArray( pmat );

		if ( pa == NULL )
			return( FALSE );

		/* Get Dimensions */

		ndim = mxGetNumberOfDimensions( pa );

		dims = mxGetDimensions( pa );

		for ( i=0 ; i < ndim ; i++ )
		{
			MAT->glb[i] = 1;
			MAT->gub[i] = dims[i];
		}

		MAT->dim = ndim;

		/* Check for Complex vs. Real Data */

		has_values = ( mxIsComplex( pa ) )
			? MATRIX_COMPLEX : MATRIX_REAL;
	}
	
	/* Read In Sparse Matrix */

	if ( mxIsSparse( pa ) )
	{
		if ( !contflag )
		{
			/* Get Matrix Parameters */

			ncol = mxGetN( pa );

			ir = mxGetIr( pa );
			jc = mxGetJc( pa );

			data = mxGetPr( pa );

			if ( has_values == MATRIX_COMPLEX )
				idata = mxGetPi( pa );

			/* Determine # of Non-Zeros */

			nelems = jc[ ncol ] - jc[0];

			/* Set Up Matrix */

			free_mat();

			if ( !alloc_mat( nelems, has_values ) )
				return( FALSE );

			MAT->symm = MS_UNSYMMETRIC;
		}

		/* Calculate Progress Factor */

		progress = MAT->nelems / PROGRESS_TICKS;
		progress = ( progress ) ? progress : 1;

		/* Set Non-Zero Values */

		if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
			C = MAT->complex;
	
		else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
			D = MAT->data;
	
		else
			T = MAT->template;

		if ( !contflag )
		{
			cstart = 0;

			cnt = 0;
		}

		else
			cstart = c;

		num = 0;

		for ( c=cstart ; c < ncol ; c++ )
		{
			start = jc[ c ];
			stop = jc[ c + 1 ];

			if ( start != stop )
			{
				/* Check for Partial Load Overflow */

				if ( num + ( stop - start ) > MAT->alloc_size )
				{
					MAT->fpsave = (FILE *) pmat;

					MAT->read_mat_func = read_matlab_mat;

					MAT->load_inst = inst;

					MAT->first_loaded = ( MAT->last_loaded )
						? ( MAT->last_loaded + 1 ) : 0;
					MAT->last_loaded = MAT->first_loaded + num - 1;

					return( TRUE );
				}

				for ( r=start ; r < stop ; r++ )
				{
					if ( !( cnt % progress ) )
					{
						update_progress( MAIN_VIEW, cnt / progress );
						INST_CHECK_RETURN( MAIN_VIEW, inst,
							MAIN_VIEW->fileinst, FALSE, FALSE );
					}

					x = c + 1;
					y = ir[ r ] + 1;

					if ( TRANSPOSE_MATRIX )
					{
						tmp = x;
						x = y;
						y = tmp;
					}

					if ( NOLOAD_VALUES )
					{
						T->x = x;
						T->y = y;

						cnt++;  /* for progress */

						T++;
					}

					else if ( has_values == MATRIX_REAL )
					{
						D->x = x;
						D->y = y;

						D->val = data[ cnt++ ];

						D++;
					}

					else if ( has_values == MATRIX_COMPLEX )
					{
						C->x = x;
						C->y = y;

						C->val = data[ cnt ];
						C->jval = idata[ cnt ];

						cnt++;

						C++;
					}

					num++;
				}
			}
		}

		if ( pa )
			mxDestroyArray( pa );
	}

	/* Read In Full Matrix (huhuhuhuh...) */

	else
	{
		if ( !contflag )
		{
			/* Get Matrix Parameters */

			ncol = mxGetN( pa );
			nrow = mxGetM( pa );

			data = mxGetPr( pa );

			if ( has_values == MATRIX_COMPLEX )
				idata = mxGetPi( pa );

			/* Determine # of Non-Zeros */

			nelems = ncol * nrow;

			/* Set Up Matrix */

			free_mat();

			if ( !alloc_mat( nelems, has_values ) )
				return( FALSE );

			MAT->symm = MS_UNSYMMETRIC;
		}

		/* Calculate Progress Factor */

		progress = MAT->nelems / PROGRESS_TICKS;
		progress = ( progress ) ? progress : 1;

		/* Set Non-Zero Values */

		if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
			C = MAT->complex;
	
		else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
			D = MAT->data;
	
		else
			T = MAT->template;

		if ( !contflag )
		{
			cstart = 0;

			cnt = 0;
		}

		else
			cstart = c;

		num = 0;

		for ( c=cstart ; c < ncol ; c++ )
		{
			/* Check for Partial Load Overflow */

			if ( num + nrow > MAT->alloc_size )
			{
				MAT->fpsave = (FILE *) pmat;

				MAT->read_mat_func = read_matlab_mat;

				MAT->load_inst = inst;

				MAT->first_loaded = ( MAT->last_loaded )
					? ( MAT->last_loaded + 1 ) : 0;
				MAT->last_loaded = MAT->first_loaded + num - 1;

				return( TRUE );
			}

			for ( r=0 ; r < nrow ; r++ )
			{
				if ( !( cnt % progress ) )
				{
					update_progress( MAIN_VIEW, cnt / progress );
					INST_CHECK_RETURN( MAIN_VIEW, inst,
						MAIN_VIEW->fileinst, FALSE, FALSE );
				}

				x = c + 1;
				y = r + 1;

				if ( TRANSPOSE_MATRIX )
				{
					tmp = x;
					x = y;
					y = tmp;
				}

				if ( NOLOAD_VALUES )
				{
					T->x = x;
					T->y = y;

					cnt++;  /* for progress */

					T++;
				}

				else if ( has_values == MATRIX_REAL )
				{
					D->x = x;
					D->y = y;

					D->val = data[ cnt ];

					cnt++;

					D++;
				}

				else if ( has_values == MATRIX_COMPLEX )
				{
					C->x = x;
					C->y = y;

					C->val = data[ cnt ];
					C->jval = idata[ cnt ];

					cnt++;

					C++;
				}

				num++;
			}
		}

		if ( pa )
			mxDestroyArray( pa );
	}

	MAT->fpsave = (FILE *) pmat;

	MAT->read_mat_func = read_matlab_mat;

	MAT->load_inst = inst;

	MAT->first_loaded = ( MAT->last_loaded )
		? ( MAT->last_loaded + 1 ) : 0;
	MAT->last_loaded = MAT->first_loaded + MAT->alloc_size - 1;

	return( TRUE );
}

#endif


int
open_mat_file( file, line, mat, is_compressed )
char *file;
char *line;
FILE **mat;
int *is_compressed;
{
	char cmd[255];

	char *uncompress;
	char *suffix;

	/* Close Any Open File */

	if ( *mat != NULL )
		( *is_compressed ) ? pclose( *mat ) : fclose( *mat );

	/* Open File - For Compressed, Make Sure It Exists & Is Readable */

	*mat = fopen( file, "r" );

	if ( !filecheck( *mat, file ) )
		return( FALSE );

	/* Check for Compressed or Gzipped Files */

	suffix = suffix_of( file );

	*is_compressed = 0;

	if ( !strcmp( suffix, "Z" ) )
	{
		uncompress = "zcat";

		(*is_compressed)++;
	}

	else if ( !strcmp( suffix, "gz" ) )
	{
		uncompress = "gzcat";

		(*is_compressed)++;
	}

	if ( *is_compressed )
	{
		fclose( *mat );

		sprintf( cmd, "%s %s", uncompress, file );

		*mat = popen( cmd, "r" );

		if ( !filecheck( *mat, cmd ) )
			return( FALSE );
	}

	/* Read Header Line */

#ifdef FILE_BUFFERING
	if ( (line = read_line( *mat )) == NULL )
#else
	if ( fgets( line, 1024, *mat ) == NULL )
#endif
	{
		( *is_compressed ) ? pclose( *mat ) : fclose( *mat );

		if ( *is_compressed )
		{
			sprintf( cmd,
				"setMsg {%s - \"%s\" Not Found?}",
				"Error Uncompressing Matrix File",
				uncompress );

			Tcl_Eval( interp, cmd );
		}

		else
		{
			Tcl_Eval( interp,
				"setMsg {Error Reading Matrix File Header.}" );
		}

		return( FALSE );
	}

	/* Save Pointer to First Header Line */

	firstline = line;

	return( TRUE );
}


int
alloc_mat( nelems, has_values )
int nelems;
int has_values;
{
	char tmp[255];

	unsigned size;

	int alloc_size;
	int valid;

	/* Initial Defaults */

	alloc_size = nelems;

	valid = TRUE;

	/* Verify Base # of Elements is Under Hard-coded Limits */

	if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
		size = sizeof( struct complex_struct );
	else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
		size = sizeof( struct data_struct );
	else
		size = sizeof( struct template_struct );

	if ( (unsigned) alloc_size * size > MAX_MAT_SIZE * 1000000 )
		alloc_size = ( 1000000 * MAX_MAT_SIZE ) / size;

	/* Try to Allocate Matrix Memory */

	do
	{
		if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
		{
			MAT->complex = (COMPLEX) MALLOC( (unsigned) alloc_size
				* sizeof( struct complex_struct ) );

			if ( MAT->complex == NULL )
				valid = FALSE;
		}

		else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
		{
			MAT->data = (DATA) MALLOC( (unsigned) alloc_size
				* sizeof( struct data_struct ) );

			if ( MAT->data == NULL )
				valid = FALSE;
		}

		else
		{
			MAT->template = (TEMPLATE) MALLOC( (unsigned) alloc_size
				* sizeof( struct template_struct ) );

			if ( MAT->template == NULL )
				valid = FALSE;
		}

		if ( !valid )
			alloc_size /= 2;
	}
	while ( !valid && alloc_size > 0 );

	if ( !valid )
	{
		sprintf( tmp,
			"setMsg {Error Allocating Memory for Matrix (\"%s\").}",
			MATVIEW_MATFILE );
		Tcl_Eval( interp, tmp );
	}

	else
	{
		MAT->alloc_size = alloc_size;

		MAT->nelems = nelems;

		MAT->first_loaded = -1;
		MAT->last_loaded = 0;
	}

	return( valid );
}


int
realloc_mat( nelems, has_values )
int nelems;
int has_values;
{
	TEMPLATE TSAVE;

	COMPLEX CSAVE;

	DATA DSAVE;

	int old_nelems;

#ifdef NO_BCOPY

	TEMPLATE Tsrc, Tdst;

	COMPLEX Csrc, Cdst;

	DATA Dsrc, Ddst;

	int i;

#endif

	/* Already Maxed Out Allocation... */

	if ( MAT->alloc_size != MAT->nelems )
	{
		MAT->nelems = nelems;

		return( TRUE );
	}

	/* Save Existing Data */

	if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
	{
		CSAVE = MAT->complex;

		MAT->complex = (COMPLEX) NULL;
	}

	else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
	{
		DSAVE = MAT->data;

		MAT->data = (DATA) NULL;
	}

	else
	{
		TSAVE = MAT->template;

		MAT->template = (TEMPLATE) NULL;
	}

	old_nelems = MAT->nelems;

	/* Re-Allocate Matrix */

	if ( !alloc_mat( nelems, has_values ) )
		return( FALSE );

	/* Copy Old Data to New Memory */

	if ( has_values == MATRIX_COMPLEX && !NOLOAD_VALUES )
	{
#ifdef NO_BCOPY
		Csrc = CSAVE;
		Cdst = MAT->complex;

		for ( i=0 ; i < old_nelems ; i++ )
		{
			Cdst->x = Csrc->x;
			Cdst->y = Csrc->y;

			Cdst->val = Csrc->val;
			Cdst->jval = Csrc->jval;

			Cdst++;  Csrc++;
		}
#else
		bcopy( (void *) CSAVE, (void *) MAT->complex,
			(size_t) old_nelems * sizeof( struct complex_struct ) );
#endif

		FREE( CSAVE );
	}

	else if ( has_values == MATRIX_REAL && !NOLOAD_VALUES )
	{
#ifdef NO_BCOPY
		Dsrc = DSAVE;
		Ddst = MAT->data;

		for ( i=0 ; i < old_nelems ; i++ )
		{
			Ddst->x = Dsrc->x;
			Ddst->y = Dsrc->y;

			Ddst->val = Dsrc->val;

			Ddst++;  Dsrc++;
		}
#else
		bcopy( (void *) DSAVE, (void *) MAT->data,
			(size_t) old_nelems * sizeof( struct data_struct ) );
#endif

		FREE( DSAVE );
	}

	else
	{
#ifdef NO_BCOPY
		Tsrc = TSAVE;
		Tdst = MAT->template;

		for ( i=0 ; i < old_nelems ; i++ )
		{
			Tdst->x = Tsrc->x;
			Tdst->y = Tsrc->y;

			Tdst++;  Tsrc++;
		}
#else
		bcopy( (void *) TSAVE, (void *) MAT->template,
			(size_t) old_nelems * sizeof( struct template_struct ) );
#endif

		FREE( TSAVE );
	}

	return( TRUE );
}


void
free_mat()
{
	if ( MAT->template != NULL )
	{
		FREE( MAT->template );
		MAT->template = (TEMPLATE) NULL;
	}

	if ( MAT->complex != NULL )
	{
		FREE( MAT->complex );
		MAT->complex = (COMPLEX) NULL;
	}

	if ( MAT->data != NULL )
	{
		FREE( MAT->data );
		MAT->data = (DATA) NULL;
	}
}


int
load_rest_of_mat( V )
VIEW V;
{
	int inst;
	int ret;

	ret = FALSE;

	if ( MAT->read_mat_func != NULL )
	{
		inst = V->instance;
		Tcl_Eval( interp, "setMsg {Reading More of Matrix File...}" );
		INST_CHECK_RETURN( V, inst, V->instance, FALSE, FALSE );

		MAT->first_loaded = -1;  /* reset in case interrupted */

		inst = MAT->load_inst;

		ret = (MAT->read_mat_func)( MAT->fpsave, TRUE, inst );

		if ( MAT->last_loaded == MAT->nelems )
		{
			( MAT->is_compressed )
				? pclose( MAT->fpsave ) : fclose( MAT->fpsave );

			MAT->fpsave = (FILE *) NULL;
		}
	}

	else
	{
		Tcl_Eval( interp,
			"setMsg {Error Reading Matrix File...  Abort Read.}" );
	}

	return( ret );
}


#ifdef FILE_BUFFERING

/* File Buffering Stuff */

char *
read_line( fp )
FILE *fp;
{
	static FILE *fp_save = NULL;

	static char *bufptr = NULL;
	static char *endptr = NULL;
	static char *buf = NULL;

	static int buffering = TRUE;
	static int data_size = -1;
	static int size = -1;

	static char line[1024];

	static int cnt = 0;

	int fd;

	/* Check for Valid File. */

	if ( fp == NULL )
		return( (char *) NULL );
	
	else if ( fp != fp_save )
	{
		if ( size > 0 )
			buffering = TRUE;

		bufptr = NULL;
		endptr = NULL;

		fp_save = fp;
	}

	/* Kick-Start Buffer */

	if ( size < 0 )
	{
		size = FILE_BUF_SIZE * 1024;

		if ( size > 0 )
		{
			buf = (char *) MALLOC( (unsigned) ( size + 1 )
				* sizeof( char ) );

			if ( buf == NULL )
			{
				fprintf( stderr,
					"\nError Allocating Internal File Buffer.\n\n" );
				fprintf( stderr,
					"File Buffering Disabled.\n\n" );

				buffering = FALSE;

				size = 0;
			}
		}

		else
		{
			buffering = FALSE;

			size = 0;
		}
	}

	/* Handle Non-Buffered Case */

	if ( !buffering )
	{
		if ( fgets( line, 1024, fp ) == NULL )
			return( (char *) NULL );

		return( line );
	}

	/* Do Buffering */

	else
	{
		/* Load Buffer With Initial Data */

		if ( bufptr == NULL )
		{
			fd = fileno( fp );

			data_size = read( fd, buf, size );

			if ( !data_size )
				return( (char *) NULL );

			endptr = bufptr = buf;

			cnt = 0;

			while ( *endptr != '\n' && cnt < data_size )
			{
				endptr++;  cnt++;
			}

			if ( cnt >= data_size )
				return( (char *) NULL );

			*endptr = '\0';

			return( bufptr );
		}

		/* Continue Parsing Buffered Data */

		else
		{
			bufptr = endptr + 1;

			endptr++;

			cnt++;

			while ( *endptr != '\n' && cnt < data_size )
			{
				endptr++;  cnt++;
			}

			/* Ran Out of Data - Read More */

			if ( cnt >= data_size )
			{
				*endptr = '\0';

				strcpy( buf, bufptr );

				cnt = ( endptr - bufptr ) / sizeof( char );

				endptr = buf + cnt;

				bufptr = buf;

				/* Read More Data */

				fd = fileno( fp );

				data_size = read( fd, endptr, size - cnt );

				if ( !data_size )
				{
					if ( bufptr != endptr )
						return( bufptr );
					
					else
						return( (char *) NULL );
				}

				while ( *endptr != '\n' && cnt < data_size )
				{
					endptr++;  cnt++;
				}

				if ( cnt >= data_size )
					return( (char *) NULL );
				
				else
					return( bufptr );
			}

			/* Return Next Line */

			*endptr = '\0';

			return( bufptr );
		}
	}
}

#endif

