/* $Id: matview.c,v 1.94 2000/06/19 18:30:45 kohl Exp $ */


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
#include "time.h"


/* TCL Initialization Routine */

int Tcl_AppInit();


/* TCL C Command Routines */

int get_matview_default_dir_proc();
int set_color_function_proc();
int get_complex_view_proc();
int set_complex_view_proc();
int register_menuvar_proc();
int load_bitmap_file_proc();
int set_plane_choice_proc();
int fill_about_text_proc();
int interface_param_proc();
int set_color_range_proc();
int draw_color_map_proc();
int set_vis_region_proc();
int set_zero_color_proc();
int canvas_handle_proc();
int fix_help_line_proc();
int get_data_mode_proc();
int set_data_mode_proc();
int get_reduction_proc();
int set_reduction_proc();
int double_click_proc();
int get_mat_name_proc();
int set_mat_name_proc();
int set_sep_view_proc();
int get_vis_mode_proc();
int set_vis_mode_proc();
int strip_label_proc();
int reload_mat_proc();
int resize_mat_proc();
int save_field_proc();
int title_info_proc();
int thumb_nav_proc();
int do_print_proc();
int set_grid_proc();
int set_opt_proc();
int pad_proc();


void usage();
void quit_proc();

int tcl_init();
int read_args();
int window_init();


/* Private Storage */

static int complex_was_set = 0;


/* MAIN */

int
main( argc, argv )
int argc;
char **argv;
{
	char **mung_argv;
	int mung_argc;
	int i;

	/* Mung Argv to Circumvent TK Stupidity */

	mung_argc = argc + 1;
	mung_argv = (char **) MALLOC( (unsigned) ( mung_argc + 1 )
		* sizeof( char * ) );
	memcheck( mung_argv, "Munged Argv List" );

	mung_argv[0] = argv[0];

	mung_argv[1] = copy_str( "--" );

	for ( i=2 ; i < mung_argc ; i++ )
		mung_argv[i] = argv[ i - 1 ];

	mung_argv[ mung_argc ] = (char *) NULL;

	/* Pass to TK */

	Tk_Main( mung_argc, mung_argv, Tcl_AppInit );

	return( 0 );
}


int Tcl_AppInit( itp )
Tcl_Interp *itp;
{
	/* Save Interpreter */

	interp = itp;

	/* Initialize TCL & TK */

	if ( tcl_init() == TCL_ERROR )
		return( TCL_ERROR );

	/* Get Main Window */

	Top = Tk_MainWindow( interp );

	if ( Top == NULL )
		return( TCL_ERROR );

	Disp = Tk_Display( Top );

	/* Initialize Program Constants & Structs */

	if ( program_init() == TCL_ERROR )
		return( TCL_ERROR );

	/* Handle Command Line Args */

	if ( read_args() == TCL_ERROR )
		return( TCL_ERROR );

	/* Initialize Modules */

	if ( window_init() == TCL_ERROR )
		return( TCL_ERROR );

	/* If There's a Matrix File, Schedule to Load It... */

	if ( MATVIEW_MATFILE != NULL )
		Tk_DoWhenIdle( do_load_mat, (ClientData) FALSE );

	/* Set Up User-Specific Startup File */

#if ( TCL_MAJOR_VERSION > 7 ) || ( TCL_MINOR_VERSION >= 5 )

	SET_TCL_GLOBAL( interp, "tcl_rcFileName", MATVIEW_RC_FILE );

#else

	tcl_RcFileName = MATVIEW_RC_FILE;

#endif

	/* If No Tk Application Shell Desired, Don't Return... */

	if ( NO_SHELL )
	{
		/* Source User-Specific Startup File */
		Tcl_EvalFile( interp, MATVIEW_RC_FILE );

		Tk_MainLoop();

		Tcl_Eval( interp, "exit" );

		exit( 1 );
	}

	/* Return to tkMain to Start Tk_MainLoop() */

	return( TCL_OK );
}


int
program_init()
{
	char *getenv();

	char *matview_dir;
	char *defdir;

	char tmp[1024];

	/* Seed Random Generator */

	srandom( (int) time( (time_t *) NULL ) );

	/* Get MatView Main Directory */

	matview_dir = getenv( "MATVIEW_ROOT" );

	if ( matview_dir == NULL )
	{
		defdir = matview_default_dir();

		sprintf( tmp, "%s", defdir );

		FREE( defdir );
	}

	else
		sprintf( tmp, "%s", matview_dir );

	MATVIEW_DIR = copy_str( tmp );

	/* Initialize Tcl / Tk Globals */

	SET_TCL_GLOBAL( interp, "network_display", "FALSE" );

	MENUVAR_AUTO_COLOR = (char *) NULL;
	MENUVAR_GREYSCALE = (char *) NULL;
	MENUVAR_FILTER_VALUES = (char *) NULL;
	MENUVAR_MIN_MAX = (char *) NULL;
	MENUVAR_CUMULATIVE = (char *) NULL;
	MENUVAR_ABSVALS = (char *) NULL;
	MENUVAR_PATTERN = (char *) NULL;
	MENUVAR_TRANSPOSE = (char *) NULL;
	MENUVAR_NOLOADVALUES = (char *) NULL;

	PTR_VALID = FALSE;

	PTR_X = -1;
	PTR_Y = -1;

	/* Initialize MatView Interface Globals */

	PRINT_VIEW = create_view();

	THUMB_VIEW = create_view();

	MAIN_VIEW = create_view();

	SEP_VIEW = create_view();

	MATVIEW_MATFILE = (char *) NULL;

	init_region( MAX_REGION );

	PLANE_CHOICE = PLANE_XY;

	REGIONS_SAVE = (REGION_STACK) NULL;
	REGIONS = (REGION_STACK) NULL;

	AUTO_COLOR_MAP = TRUE;  /* XXX eem2314 */

	GREYSCALE = FALSE;

	FILTER_VALUES = FALSE;

	VIEW_ABS_VALS = FALSE;

	VIEW_PATTERN = FALSE;

	TRANSPOSE_MATRIX = FALSE;

	SHOW_MIN_MAX = FALSE;

	SHOW_CUMULATIVE = FALSE;

	CUMULATIVE_MIN = 0.0;
	CUMULATIVE_MAX = 0.0;
	CUMULATIVE_SET = 0;

	NOLOAD_VALUES = FALSE;

	REDUCTION_FUNCTION = RF_AVG;

	DATA_MODE = DATA_MODE_MATRIX;

	VIS_MODE = VIS_MODE_NORMAL;

	COMPLEX_VIEW = COMPLEX_RE;

	PROGRESS_TICKS = -1;

	COLOR_MIN_VALUE = 0.0;
	COLOR_MAX_VALUE = 0.0;

	ZERO_COLOR = (char *) NULL;

	COLOR_DIST = -1.0;

	COLORIDS = (int *) NULL;
	NCOLORIDS = 0;

	GRID_STATUS = FALSE;
	GRID_ROWS = GRID_COLS = -1;
	GRID_ROFF = GRID_COFF = 0;

	THUMB_BORDER = -1;

	USE_SEP_VIEW = -1;

	REGION_ID = -1;

	/* Initialize MatView Matrix Globals */

	MAT = create_matrix();

	MAT->symm = MS_UNSYMMETRIC;

	MAT->loaded = ML_UNLOADED;

	SAVE_FIELD = (double *) NULL;
	SAVE_FIELD_NELEMS = -1;

#ifdef FILE_BUFFERING
	FILE_BUF_SIZE = 1024; /* Kbytes */
#endif

	NO_SHELL = TRUE;  /* XXX eem2314 */

	vflag = 0;

	return( TCL_OK );
}


int
read_args()
{
	char *argv_str;

	int i, j, k;
	int do_usage;
	int len;

	argv_str = GET_TCL_GLOBAL( interp, "argv" );

	/* XXX gcc says arg 4 is an incompatible pointer */
	if ( Tcl_SplitList( interp, argv_str, &Argc, &Argv ) == TCL_ERROR )
		return( TCL_ERROR );

	do_usage = 0;

	for ( i=0 ; i < Argc ; i++ )
	{
		/* Short Circuit to Usage */
		if ( !strcmp( Argv[i], "-help" )
				|| !strcmp( Argv[i], "-Help" )
				|| !strcmp( Argv[i], "-h" )
				|| !strcmp( Argv[i], "-H" ) )
			usage();

		/* Short Circuit to "-network" Option */
		if ( !strcmp( Argv[i], "-network" ) )
		{
			SET_TCL_GLOBAL( interp, "network_display", "TRUE" );

			PROGRESS_TICKS = 10;
		}

		/* Short Circuit to "-filter" Option */
		else if ( !strcmp( Argv[i], "-filter" ) )
			FILTER_VALUES = TRUE;

		/* Short Circuit to "-noshell" Option */
		else if ( !strcmp( Argv[i], "-noshell" ) )
			NO_SHELL = TRUE;

		/* Short Circuit to "-abs" Option */
		else if ( !strcmp( Argv[i], "-abs" ) )
			VIEW_ABS_VALS = TRUE;

		/* Standard Options */
		else if ( Argv[i][0] == '-' )
		{
			k = i + 1;

			len = strlen( Argv[i] );

			for ( j=1 ; j < len ; j++ )
			{
				switch ( Argv[i][j] )
				{
					case 'F':
					{
						if ( MATVIEW_MATFILE != NULL )
							FREE( MATVIEW_MATFILE );

						if ( k < Argc )
							MATVIEW_MATFILE = copy_str( Argv[ k++ ] );

						else
							do_usage++;

						break;
					}

					case 'R':
					{
						if ( k < Argc )
						{
							if ( !strcmp( Argv[ k ], "avg" ) )
								REDUCTION_FUNCTION = RF_AVG;

							else if ( !strcmp( Argv[ k ], "min" ) )
								REDUCTION_FUNCTION = RF_MIN;

							else if ( !strcmp( Argv[ k ], "max" ) )
								REDUCTION_FUNCTION = RF_MAX;

							else if ( !strcmp( Argv[ k ], "sum" ) )
								REDUCTION_FUNCTION = RF_SUM;

							else if ( !strcmp( Argv[ k ], "midpoint" ) )
								REDUCTION_FUNCTION = RF_MIDPOINT;

							else if ( !strcmp( Argv[ k ], "count" ) )
								REDUCTION_FUNCTION = RF_COUNT;

							else if ( !strcmp( Argv[ k ], "std" ) )
								REDUCTION_FUNCTION = RF_STD;

							else if ( !strcmp( Argv[ k ], "var" ) )
								REDUCTION_FUNCTION = RF_VAR;

							else
							{
								printf(
								"Unknown Reduction Function \"%s\".\n",
									Argv[ k ] );
							
								do_usage++;
							}

							k++;
						}

						else
							do_usage++;

						break;
					}

					case 'D':
					{
						if ( k < Argc )
						{
							if ( !strcmp( Argv[ k ], "matrix" ) )
								DATA_MODE = DATA_MODE_MATRIX;

							else if ( !strcmp( Argv[ k ], "field" ) )
								DATA_MODE = DATA_MODE_FIELD;

							else
							{
								printf( "Unknown Data Mode \"%s\".\n",
									Argv[ k ] );
							
								do_usage++;
							}

							k++;
						}

						else
							do_usage++;

						break;
					}

					case 'V':
					{
						if ( k < Argc )
						{
							if ( !strcmp( Argv[ k ], "normal" ) )
								VIS_MODE = VIS_MODE_NORMAL;

							else if ( !strcmp( Argv[ k ], "diff" ) )
								VIS_MODE = VIS_MODE_DIFF;

							else
							{
								printf( "Unknown Vis Mode \"%s\".\n",
									Argv[ k ] );
							
								do_usage++;
							}

							k++;
						}

						else
							do_usage++;

						break;
					}

					case 'A':
					{
						AUTO_COLOR_MAP = TRUE;

						break;
					}

					case 'G':
					{
						GREYSCALE = TRUE;

						break;
					}

					case 'P':
					{
						VIEW_PATTERN = TRUE;

						break;
					}

					case 'X':
					{
						if ( k < Argc )
						{
							if ( !strcmp( Argv[ k ], "re" ) )
								COMPLEX_VIEW = COMPLEX_RE;

							else if ( !strcmp( Argv[ k ], "im" ) )
								COMPLEX_VIEW = COMPLEX_IM;

							else if ( !strcmp( Argv[ k ], "magn" ) )
								COMPLEX_VIEW = COMPLEX_MAGN;

							else if ( !strcmp( Argv[ k ], "phase" ) )
								COMPLEX_VIEW = COMPLEX_PHASE;

							else
							{
								printf(
								"Unknown Complex View Type \"%s\".\n",
									Argv[ k ] );
							
								do_usage++;
							}

							complex_was_set++;

							k++;
						}

						else
							do_usage++;

						break;
					}

					case 'Z':
					{
						if ( k < Argc )
						{
							if ( ZERO_COLOR != NULL )
								FREE( ZERO_COLOR );

							ZERO_COLOR = copy_str( Argv[ k++ ] );
						}

						else
							do_usage++;

						break;
					}

					case 'T':
					{
						TRANSPOSE_MATRIX = TRUE;

						break;
					}

					case 'M':
					{
						SHOW_MIN_MAX = TRUE;

						break;
					}

					case 'C':
					{
						SHOW_CUMULATIVE = TRUE;

						break;
					}

					case 'N':
					{
						NOLOAD_VALUES = TRUE;

						break;
					}

#ifdef FILE_BUFFERING
					case 'B':
					{
						if ( k < Argc )
							FILE_BUF_SIZE = atoi( Argv[k++] );
						
						else
							do_usage++;

						break;
					}
#endif

					case 'I':
					{
						if ( k + 1 < Argc )
						{
							if ( !strcmp( Argv[k],
									"progress_ticks" ) )
								PROGRESS_TICKS = atoi( Argv[k+1] );
							
							else
							{
								printf(
								"Unknown Interface Parameter \"%s\".\n",
									Argv[k] );
								
								do_usage++;
							}

							k += 2;
						}

						else
							do_usage++;

						break;
					}

					case 'v': vflag++; break;

					default:
					{
						printf( "Unknown Option -%c\n", Argv[i][j] );

						do_usage++;

						break;
					}
				}
			}

			i = k - 1;
		} else {
			// Mod by eem2314: allow spec filename without -F:
			if ( MATVIEW_MATFILE != NULL )
				FREE( MATVIEW_MATFILE );

			MATVIEW_MATFILE = copy_str( Argv[ i ] );
		}
	}

	if ( do_usage )
		usage();

	return( TCL_OK );
}


void
usage()
{
	printf( "\nusage:  matview [-F <mat file>] " );
	printf( "[-R <rf>] [-D <mode>] [-V <mode>]\n" );
	printf( "\t\t[-AGPTMC] [-X <view>] [-abs] [-Z <color>] [-N]\n" );
	printf( "\t\t[-I <p> <val>] [-network] [-noshell] [-v]\n\n" );

	printf( "where:\n" );
	printf( "------\n" );

	printf( "-F <mat file> = Set Matrix File.\n" );

	printf( "-R <rf>       = Set Reduction Function:\n" );
	printf( "                <rf> =\n" );
	printf( "                   \"avg\"      - Average Values.\n" );
	printf( "                   \"min\"      - Minimum of Values.\n" );
	printf( "                   \"max\"      - Maximum of Values.\n" );
	printf( "                   \"sum\"      - Sum Values.\n" );
	printf( "                   \"midpoint\" - Midpoint of Values.\n" );
	printf( "                   \"count\"    - " );
	printf( "Count Number of Non-Zero Values.\n" );
	printf( "                   \"std\"      - Standard Deviation.\n" );
	printf( "                   \"var\"      - Variance.\n" );

	printf( "-D <mode>     = Set Data Mode:\n" );
	printf( "                <mode> =\n" );
	printf( "                   \"matrix\"   - " );
	printf( "Matrix Semantics Mode.\n" );
	printf( "                   \"field\"    - " );
	printf( "Field Semantics Mode.\n" );

	printf( "-V <mode>     = Set Visualization Mode:\n" );
	printf( "                <mode> =\n" );
	printf( "                   \"normal\"   - Normal Mode.\n" );
	printf( "                   \"diff\"     - Differential Mode.\n" );

	printf( "-A            = Turn On Auto Color Map Ranging.\n" );
	printf( "-G            = Turn On Greyscale View Drawing.\n" );
	printf( "-P            = Display Matrix Using Pattern Only " );
	printf( "(No Values).\n" );

	printf( "-X <view>     = Set Complex Value View Type:\n" );
	printf( "                <view> =\n" );
	printf( "                   \"re\"    - View Real Component.\n" );
	printf( "                   \"im\"    - View Imag Component.\n" );
	printf( "                   \"magn\"  - View Vector Magnitude.\n" );
	printf( "                   \"phase\" - View Vector Phase.\n" );

	printf( "-T            = Display Matrix in Transposed Form.\n" );
	printf( "-M            = Show Min & Max Data Values.\n" );
	printf( "-C            = Show Cumulative Min & Max Values.\n" );

	printf( "-abs          = Turn On Absolute Value Viewing.\n" );

	printf( "-Z <color>    = Set Default Zero Value Color.\n" );

	printf( "-N            = Don't Load Matrix Data Values, " );
	printf( "Just Pattern.\n" );

#ifdef FILE_BUFFERING
	printf( "-B <size>     = Internal File Buffer Size " );
	printf( "(in Kbytes).\n" );
#endif

	printf( "-I <p> <val>  = Interface Paramater Tuning.\n" );
	printf( "                <p> <val> =\n" );
	printf( "                   \"progress_ticks\" <nticks>\n" );

	printf( "-network      = Optimize for display over network.\n" );

	printf( "-noshell      = Do Not Start Tk Application Shell.\n" );

	printf( "-v            = Verbose Operation.\n" );

	printf( "\n" );

	exit( 0 );
}


int
tcl_init()
{
	/* Initialize TCL / TK */

	if ( Tcl_Init( interp ) == TCL_ERROR )
		return( TCL_ERROR );

	if ( Tk_Init( interp ) == TCL_ERROR )
		return( TCL_ERROR );

	/* Create TCL Commands for Action Routines */

	Tcl_CreateCommand( interp, "canvas_handle",
		canvas_handle_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "thumb_nav",
		thumb_nav_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_plane_choice",
		set_plane_choice_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_color_function",
		set_color_function_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "get_complex_view",
		get_complex_view_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_complex_view",
		set_complex_view_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_color_range",
		set_color_range_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "draw_color_map",
		draw_color_map_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_vis_region",
		set_vis_region_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_zero_color",
		set_zero_color_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "get_data_mode",
		get_data_mode_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_data_mode",
		set_data_mode_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "get_reduction",
		get_reduction_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_reduction",
		set_reduction_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "get_mat_name",
		get_mat_name_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_mat_name",
		set_mat_name_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_sep_view",
		set_sep_view_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "get_vis_mode",
		get_vis_mode_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_vis_mode",
		set_vis_mode_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "do_print",
		do_print_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_grid",
		set_grid_proc, (ClientData) NULL, (vfp) NULL );

	/* Create TCL Utility Commands */

	Tcl_CreateCommand( interp, "get_matview_default_dir",
		get_matview_default_dir_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "interfaceParamHandle",
		interface_param_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "register_menuvar",
		register_menuvar_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "load_bitmap_file",
		load_bitmap_file_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "fill_about_text",
		fill_about_text_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "fix_help_line",
		fix_help_line_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "double_click",
		double_click_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "strip_label",
		strip_label_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "reload_mat",
		reload_mat_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "resize_mat",
		resize_mat_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "save_field",
		save_field_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "title_info",
		title_info_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "set_opt",
		set_opt_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "quit",
		(Tcl_CmdProc *) quit_proc, (ClientData) NULL, (vfp) NULL );

	Tcl_CreateCommand( interp, "pad",
		pad_proc, (ClientData) NULL, (vfp) NULL );

	return( TCL_OK );
}


int
window_init()
{
	static char fname[1024];

	FILE *fptest;

	char *tmp;

	/* Create Interface - Run TCL Script */

	fptest = fopen( "matview.tcl", "r" );

	if ( fptest != NULL )
	{
		strcpy( fname, "matview.tcl" );

		fclose( fptest );
	}

	else
		sprintf( fname, "%s/matview.tcl", MATVIEW_DIR );

	if ( Tcl_EvalFile( interp, fname ) == TCL_ERROR )
		return( TCL_ERROR );

	/* Read in TCL Global Constants */

	MV = GET_TCL_GLOBAL( interp, "MV" );
	MV_M = GET_TCL_GLOBAL( interp, "MV_M" );
	MV_G = GET_TCL_GLOBAL( interp, "MV_G" );

	PROGRESS_HEIGHT = GET_TCL_GLOBAL( interp, "progress_height" );
	PROGRESS_WIDTH = GET_TCL_GLOBAL( interp, "progress_width" );

	/* Get TCL Spacing Constants */

	tmp = GET_TCL_GLOBAL( interp, "navigate_size" );
	NAVIGATE_SIZE = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "border_space" );
	BORDER_SPACE = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "frame_border" );
	FRAME_BORDER = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "row_height" );
	ROW_HEIGHT = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "col_width" );
	COL_WIDTH = atoi( tmp );

	/* Initialize Main View Window Display Object */

	MAIN_VIEW->win = create_win();

	MAIN_VIEW->win->toplevel = GET_TCL_GLOBAL( interp, "MV" );

	MAIN_VIEW->win->canvas = GET_TCL_GLOBAL( interp, "MV_C" );

	MAIN_VIEW->image = GET_TCL_GLOBAL( interp, "MV_I" );

	tmp = GET_TCL_GLOBAL( interp, "image_height" );
	MAIN_VIEW->ht = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "image_width" );
	MAIN_VIEW->wt = atoi( tmp );

	MAIN_VIEW->progress_canvas = GET_TCL_GLOBAL( interp, "MV_P" );

	tmp = GET_TCL_GLOBAL( interp, "progress_bar" );
	MAIN_VIEW->progress_bar = atoi( tmp );

	MAIN_VIEW->win->htvar = "matview_height";
	MAIN_VIEW->win->wtvar = "matview_width";

	MAIN_VIEW->win->xoff = MAIN_VIEW->win->yoff =
		BORDER_SPACE + FRAME_BORDER;

	tmp = GET_TCL_GLOBAL( interp, "image_id" );
	MAIN_VIEW->image_id = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "rect_id" );
	MAIN_VIEW->rect_id = atoi( tmp );

	MAIN_VIEW->do_auto_color = TRUE;

	MAIN_VIEW->do_diff_mode = TRUE;

	MAIN_VIEW->do_grid = TRUE;

	MAIN_VIEW->draw_type = DRAW_NORMAL | DRAW_PORTRAIT;

	/* Initialize Thumb View Window Display Object */

	THUMB_VIEW->win = create_win();

	THUMB_VIEW->win->toplevel = GET_TCL_GLOBAL( interp, "MV" );

	THUMB_VIEW->win->canvas = GET_TCL_GLOBAL( interp, "MV_TC" );

	THUMB_VIEW->image = GET_TCL_GLOBAL( interp, "MV_TI" );

	tmp = GET_TCL_GLOBAL( interp, "thumb_height" );
	THUMB_VIEW->ht = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "thumb_width" );
	THUMB_VIEW->wt = atoi( tmp );

	THUMB_VIEW->progress_canvas = GET_TCL_GLOBAL( interp, "MV_P" );

	tmp = GET_TCL_GLOBAL( interp, "progress_bar" );
	THUMB_VIEW->progress_bar = atoi( tmp );

	THUMB_VIEW->win->htvar = "thumb_height";
	THUMB_VIEW->win->wtvar = "thumb_width";

	THUMB_VIEW->win->xoff = THUMB_VIEW->win->yoff = 2 * FRAME_BORDER;

	tmp = GET_TCL_GLOBAL( interp, "thumb_image_id" );
	THUMB_VIEW->image_id = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "thumb_rect_id" );
	THUMB_VIEW->rect_id = atoi( tmp );

	THUMB_VIEW->do_auto_color = FALSE;

	THUMB_VIEW->do_diff_mode = FALSE;

	THUMB_VIEW->do_grid = FALSE;

	THUMB_VIEW->draw_type = DRAW_NORMAL | DRAW_PORTRAIT;

	/* Initialize Separate View Window Display Object */

	SEP_VIEW->win = create_win();

	SEP_VIEW->win->toplevel = GET_TCL_GLOBAL( interp, "VF" );

	SEP_VIEW->win->canvas = GET_TCL_GLOBAL( interp, "VF_C" );

	SEP_VIEW->image = GET_TCL_GLOBAL( interp, "MV_I" );

	tmp = GET_TCL_GLOBAL( interp, "vf_height" );
	SEP_VIEW->ht = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "vf_width" );
	SEP_VIEW->wt = atoi( tmp );

	SEP_VIEW->progress_canvas = GET_TCL_GLOBAL( interp, "VF_P" );

	tmp = GET_TCL_GLOBAL( interp, "sep_progress_bar" );
	SEP_VIEW->progress_bar = atoi( tmp );

	SEP_VIEW->win->htvar = "view_field_height";
	SEP_VIEW->win->wtvar = "view_field_width";

	SEP_VIEW->win->xoff = SEP_VIEW->win->yoff =
		BORDER_SPACE + FRAME_BORDER;

	tmp = GET_TCL_GLOBAL( interp, "vf_image_id" );
	SEP_VIEW->image_id = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "vf_rect_id" );
	SEP_VIEW->rect_id = atoi( tmp );

	SEP_VIEW->do_auto_color = TRUE;

	SEP_VIEW->do_diff_mode = TRUE;

	SEP_VIEW->do_grid = TRUE;

	SEP_VIEW->draw_type = DRAW_NORMAL | DRAW_PORTRAIT;

	/* Initialize Print View Window Display Object */

	PRINT_VIEW->win = create_win();

	PRINT_VIEW->win->toplevel = GET_TCL_GLOBAL( interp, "PV" );

	PRINT_VIEW->win->canvas = GET_TCL_GLOBAL( interp, "PV_C" );

	PRINT_VIEW->image = GET_TCL_GLOBAL( interp, "PV_I" );

	tmp = GET_TCL_GLOBAL( interp, "pv_height" );
	PRINT_VIEW->ht = atoi( tmp );

	tmp = GET_TCL_GLOBAL( interp, "pv_width" );
	PRINT_VIEW->wt = atoi( tmp );

	PRINT_VIEW->progress_canvas = GET_TCL_GLOBAL( interp, "PD_P" );

	tmp = GET_TCL_GLOBAL( interp, "pd_progress_bar" );
	PRINT_VIEW->progress_bar = atoi( tmp );

	PRINT_VIEW->win->htvar = "pview_field_height";
	PRINT_VIEW->win->wtvar = "pview_field_width";

	PRINT_VIEW->win->xoff = PRINT_VIEW->win->yoff =
		BORDER_SPACE + FRAME_BORDER;

	tmp = GET_TCL_GLOBAL( interp, "pv_image_id" );
	PRINT_VIEW->image_id = atoi( tmp );

	PRINT_VIEW->do_auto_color = FALSE;

	PRINT_VIEW->do_diff_mode = TRUE;

	PRINT_VIEW->do_grid = TRUE;

	PRINT_VIEW->draw_type = DRAW_CANVAS | DRAW_PORTRAIT;

	/* Set Option Default Menuvars */

	if ( MENUVAR_AUTO_COLOR )
		SET_TCL_GLOBAL( interp, MENUVAR_AUTO_COLOR,
			( AUTO_COLOR_MAP ) ? "ON" : "OFF" );

	if ( MENUVAR_GREYSCALE )
		SET_TCL_GLOBAL( interp, MENUVAR_GREYSCALE,
			( GREYSCALE ) ? "ON" : "OFF" );

	if ( MENUVAR_FILTER_VALUES )
		SET_TCL_GLOBAL( interp, MENUVAR_FILTER_VALUES,
			( FILTER_VALUES ) ? "ON" : "OFF" );

	if ( MENUVAR_MIN_MAX )
		SET_TCL_GLOBAL( interp, MENUVAR_MIN_MAX,
			( SHOW_MIN_MAX ) ? "ON" : "OFF" );

	if ( MENUVAR_CUMULATIVE )
		SET_TCL_GLOBAL( interp, MENUVAR_CUMULATIVE,
			( SHOW_CUMULATIVE ) ? "ON" : "OFF" );

	if ( MENUVAR_ABSVALS )
		SET_TCL_GLOBAL( interp, MENUVAR_ABSVALS,
			( VIEW_ABS_VALS ) ? "ON" : "OFF" );

	if ( MENUVAR_PATTERN )
		SET_TCL_GLOBAL( interp, MENUVAR_PATTERN,
			( VIEW_PATTERN ) ? "ON" : "OFF" );

	if ( MENUVAR_TRANSPOSE )
		SET_TCL_GLOBAL( interp, MENUVAR_TRANSPOSE,
			( TRANSPOSE_MATRIX ) ? "ON" : "OFF" );

	if ( MENUVAR_NOLOADVALUES )
		SET_TCL_GLOBAL( interp, MENUVAR_NOLOADVALUES,
			( NOLOAD_VALUES ) ? "ON" : "OFF" );

	/* Set Complex View No-Default Flag */

	if ( complex_was_set )
		SET_TCL_GLOBAL( interp, "nodefault_complex_view", "TRUE" );

	/* Set Other Option Defaults */

	if ( ZERO_COLOR )
	{
		SET_TCL_GLOBAL( interp, "nodefault_zero_color", "TRUE" );

		UPDATE_ENTRY( "zero_entry", ZERO_COLOR );
		SET_TCL_GLOBAL( interp, "fixed_zero_color", "ON" );
	}

	/* Set Interface Paramater Option Defaults */

	if ( PROGRESS_TICKS > 0 )
		SET_TCL_GLOBAL( interp, "nodefault_progress_ticks", "TRUE" );

	else
		PROGRESS_TICKS = 100;

	/* Re-Layout Main Panel */

	if ( Tcl_Eval( interp, "layout_main_panel" ) == TCL_ERROR )
		return( TCL_ERROR );

	return( TCL_OK );
}


char *
matview_default_dir()
{
	struct passwd *pw;

	char *getenv();

	char tmp[1024];

	char *home;

	char* froot = getenv("FROOT");					/* XXX eem2314 */
	if (froot) {   										/* XXX eem2314 */
		sprintf(tmp, "%s/share/matview", froot);   /* XXX eem2314 */
		return( copy_str( tmp ) );   					/* XXX eem2314 */
	}   														/* XXX eem2314 */

#ifdef INSTALL_ROOT
	sprintf( tmp, "%s", INSTALL_ROOT );
#else
	sprintf( tmp, "%s", MATVIEW_DEFAULT_DIR );
#endif

	if ( tmp[0] != '/' )
	{
		home = getenv( "HOME" );

		if ( home == NULL )
		{
			if ( (pw = getpwuid( getuid() )) != NULL )
				home = pw->pw_dir;
		
			else
				home = "/";
		}

#ifdef INSTALL_ROOT
		sprintf( tmp, "%s/%s", home, INSTALL_ROOT );
#else
		sprintf( tmp, "%s/%s", home, MATVIEW_DEFAULT_DIR );
#endif
	}

	return( copy_str( tmp ) );
}


/* ARGSUSED */
int
get_matview_default_dir_proc( clientData, itp, argc, argv )
ClientData clientData;
Tcl_Interp *itp;
int argc;
char **argv;
{
	Tcl_SetResult( itp, matview_default_dir(), TCL_DYNAMIC );

	return( TCL_OK );
}

