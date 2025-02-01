/* $Id: matview.h,v 1.58 1999/03/25 19:50:39 kohl Exp $ */


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


/* System Header Files */

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/time.h>
#include <pwd.h>
#include <string.h>
#include <math.h>


#ifndef PI
#define PI (M_PI)
#endif


#ifndef NO_BZERO
void bzero();
#endif

#ifndef NO_BCOPY
void bcopy();
#endif


/* TCL / TK Header Files */

#include <tk.h>


/* MatView Constants & Macros */

#include "const.h"


/* Function Types */

typedef void (*vfp)();

typedef int (*ifp)();


/* Cell Declarations */

typedef int CELL[ MAX_DIM ];

#define CELL_SIZE( _dim )			( _dim )

#define CELL_OF_DIM( _cell, _dim )	( (_cell)[_dim] )

/* Region Declarations */

typedef int REGION[2 * MAX_DIM];

#define REGION_SIZE( _dim )			( 2 * ( _dim ) )

#define REGION_MIN( _region, _dim )	( (_region)[ 2 * (_dim) ] )
#define REGION_MAX( _region, _dim )	( (_region)[ ( 2 * (_dim) ) + 1 ] )

struct region_stack_struct
{
	REGION region;
	/* REGION_STACK */ struct region_stack_struct *next;
};

typedef struct region_stack_struct *REGION_STACK;


/* Matrix Structures */

struct complex_struct
{
	int x, y;
	double val;
	double jval;
};

typedef struct complex_struct *COMPLEX;


struct data_struct
{
	int x, y;
	double val;
};

typedef struct data_struct *DATA;


struct template_struct
{
	int x, y;
};

typedef struct template_struct *TEMPLATE;


struct matrix_struct
{
	TEMPLATE template;
	COMPLEX complex;
	DATA data;
	int alloc_size;
	int nelems;
	int dim;
	int glb[MAX_DIM];
	int gub[MAX_DIM];
	int symm;
	FILE *fpsave;
	int is_compressed;
	ifp read_mat_func;
	int load_inst;
	int first_loaded;
	int last_loaded;
	int loaded;
};

typedef struct matrix_struct *MATRIX;


/* Windowing Structures */

struct rect_struct
{
	int x1, y1;
	int x2, y2;
	int drawn;
	double value;
};

typedef struct rect_struct *RECT;

struct grid_struct
{
	int id;
	int x1, y1;
	int x2, y2;
};

typedef struct grid_struct *GRID;

struct win_struct
{
	const char *toplevel;
	const char *canvas;
	char *htvar, *wtvar;
	double scale;
	int csx, csy;
	REGION region;
	int xoff, yoff;
	int border_id;
};

typedef struct win_struct *WIN;


/* View Structure :-Q */

struct view_struct
{
	WIN win;

	REGION visregion;
	CELL viscell;

	RECT rectids;
	int rectidsize;
	int nrectids;
	int nx, ny;

	GRID grids;
	int gridsize;
	int ngrids;

	int do_auto_color;
	int do_diff_mode;
	int do_grid;

	double *field;
	double *ftmp1;
	double *ftmp2;
	int field_nelems;

	const char *image;
	int wt, ht;

	double scale;
	int iscale;

	const char *progress_canvas;
	int progress_bar;

	int need_reduction;
	int need_redraw;
	int hide_image;
	int draw_type;
	int instance;
	int sequence;
	int fileinst;
	int last_rf;

	int image_id;
	int rect_id;
};

typedef struct view_struct *VIEW;


/* Routines */

MATRIX	create_matrix();

VIEW	create_view();

RECT	find_rect();

WIN		create_win();

FILE	*fopen();

double	distance();
double	angle_of();
double	atof();

void	**make_args();

char	*matview_default_dir();
char	*read_line();
char	*upper_str();
char	*trunc_str();
char	*suffix_of();
char	*copy_str();
char	*pad_num();

void	check_need_redraw();
void	clear_view();
void	color_value();
void	compute_view_min_max();
void	copy_cell();
void	copy_region();
void	do_load_mat();
int	load_mat(int);
void	do_mat_cleanup();
void	draw_color_graph();
void	draw_thumb_border();
void	draw_thumb_border_region();
void	free_mat();
void	free_matrix();
void	free_region_stack();
void	free_view();
void	free_win();
void	handle_grid();
void	handle_min_max();
void	handle_reduce_status();
void	hide_image();
void	init_cell();
void	init_region();
void	load_current_field();
void	memcheck();
void	nudge_view();
void	prepare_print();
void 	push_region();
void	recolor_matrix();
void	set_current_field();
void	free_current_field();
void	save_current_field();
void	save_image_to_file();
void	show_image();
void	landscape_coords();
void	swap_views();
void	update_color_range();
void	update_progress();
void	update_vis_region();

int		get_matview_default_dir_proc();

int		read_harwell_boeing();
int		read_mat_file();
int		read_matlab_mat();
int		read_mm_coordinate();
int		read_mm_coordinate_real();
int		read_mm_coordinate_complex();
int		read_mm_coordinate_complex();

int		alloc_mat();
int		append_cmd_str();
int		check_in_bounds();
int		check_region();
int		compare();
int		double_click();
int		filecheck();
int		generate_view_field();
int		get_plane_indices();
int		handle_frame();
int		incr_field_index();
int		is_bogus_value();
int		load_rest_of_mat();
int		match_cell();
int		match_region();
int		pop_region();
int		program_init();
int		repush_region();
int		split_line();
int		transpose_matrix();

