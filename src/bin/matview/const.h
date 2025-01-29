/* $Id: const.h,v 1.61 2008/12/03 20:31:24 kohl Exp $ */


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


/* MatView Version */

#define MATVIEW_VERSION "1.0b6"


/* Define Default Location of MatView Main Directory */

#define MATVIEW_DEFAULT_DIR "MatView"


/* Define Default Location of MatView User-Specific Startup File */

#define MATVIEW_RC_FILE "~/.matviewrc"


/* Define Max Matrix Size to Read from File (in MBs) */

#ifndef MAX_MAT_SIZE
#define MAX_MAT_SIZE		32
#endif


/* Matrix Values Constants */

#define MATRIX_PATTERN		0
#define MATRIX_REAL			1
#define MATRIX_COMPLEX		2


/* Complex Value Viewing Options */

#define COMPLEX_RE			0
#define COMPLEX_IM			1
#define COMPLEX_MAGN		2
#define COMPLEX_PHASE		3


/* Matrix Symmetry Constants */

#define MS_UNSYMMETRIC		0
#define MS_SYMMETRIC		1
#define MS_SKEW_SYMMETRIC	2
#define MS_HERMITIAN		3


/* Matrix Loaded Constants */

#define ML_UNLOADED			0
#define ML_LOADED			1
#define ML_PARTIAL			2


/* Reduction Function Constants */

#define RF_AVG		0
#define RF_MIN		1
#define RF_MAX		2
#define RF_SUM		3
#define RF_MIDPOINT	4
#define RF_COUNT	5
#define RF_STD		6
#define RF_VAR		7


/* Data Mode Constants */

#define DATA_MODE_MATRIX	1
#define DATA_MODE_FIELD		2


/* Vis Mode Constants */

#define VIS_MODE_NORMAL		1
#define VIS_MODE_DIFF		2


/* Plane Choice Constants */

#define PLANE_XY	1
#define PLANE_YZ	2
#define PLANE_XZ	3

#define MAX_DIM		3


/* Storage Order */

#define COLUMN_MAJOR	0
#define ROW_MAJOR		1


/* Draw Types */

#define DRAW_NORMAL		1
#define DRAW_IMAGE		2
#define DRAW_CANVAS		4
#define DRAW_PORTRAIT	8
#define DRAW_LANDSCAPE	16


/* True / False */

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif


/* View Instance Counter Checks */
/* NOTE: Beware of the {}, if invoked with a ; at the end... */

#define INST_CHECK_ALT_RETURN( _inst, _Inst, _cmd, _ret ) \
if ( (_inst) != (_Inst) ) \
{ \
	(_cmd); \
	return( _ret ); \
}

#define INST_CHECK_RETURN( _v, _inst, _Inst, _unhide, _ret ) \
if ( (_inst) != (_Inst) ) \
{ \
	if ( _unhide ) show_image( _v ); \
	return( _ret ); \
}

#define INST_CHECK( _v, _inst, _Inst, _unhide ) \
if ( (_inst) != (_Inst) ) \
{ \
	if ( _unhide ) show_image( _v ); \
	return; \
}


/* Divisibility Macro */

#define DIVISIBLE( _dvd, _dvr ) \
	( ( ( (_dvd) / (_dvr) ) * (_dvr) ) == (_dvd) )


/* Floating Point Equivalence Macro */

#define EQUIV( _d1, _d2, _delta ) \
	( ( (_d1) != 0.0 ) ? \
		( ( fabs( (_d1) - (_d2) ) / fabs( _d1 ) ) < (_delta) ) \
		: ( fabs( (_d1) - (_d2) ) < (_delta) ) )


/* Screen to Region Coordinate Macros */

#define SCREEN_TO_REGION( _xy, _off, _min, _max, _n, _iscale ) \
	( ( (((_xy) - (_off)) / (_iscale)) \
			* ( (_max) - (_min) + 1 ) / (_n) ) \
		+ (_min) )

#define REGION_TO_SCREEN( _xy, _off, _min, _max, _n, _iscale ) \
	( ( (_iscale) * ( (_n) * ( (_xy) - (_min) ) \
		/ ( (_max) - (_min) + 1 ) ) ) + (_off) )


/* TCL Globals Macros */

#define GET_TCL_GLOBAL( _itp, _name ) \
	Tcl_GetVar( _itp, _name, TCL_GLOBAL_ONLY )

#define SET_TCL_GLOBAL( _itp, _name, _str ) \
	Tcl_SetVar( _itp, _name, _str, TCL_GLOBAL_ONLY )


/* Entry Update Macro */

#define UPDATE_ENTRY( _entry, _value ) \
{ \
	sprintf( TMP_CMD, "%s.%s delete 0 end", MV, _entry ); \
	Tcl_Eval( interp, TMP_CMD ); \
\
	sprintf( TMP_CMD, "%s.%s insert 0 %s", MV, _entry, _value ); \
	Tcl_Eval( interp, TMP_CMD ); \
} \
\


#define UPDATE_IENTRY( _entry, _ival ) \
{ \
	sprintf( TMP_CMD, "%s.%s delete 0 end", MV, _entry ); \
	Tcl_Eval( interp, TMP_CMD ); \
\
	sprintf( TMP_CMD, "%s.%s insert 0 %d", MV, _entry, _ival ); \
	Tcl_Eval( interp, TMP_CMD ); \
} \
\


#define UPDATE_FENTRY( _entry, _fval ) \
{ \
	sprintf( TMP_CMD, "%s.%s delete 0 end", MV, _entry ); \
	Tcl_Eval( interp, TMP_CMD ); \
\
	sprintf( TMP_CMD, "%s.%s insert 0 %lg", MV, _entry, _fval ); \
	Tcl_Eval( interp, TMP_CMD ); \
} \
\


/* Drawing Macros */


/* Images */

#define DRAW_IMAGE_POINT( _itp, _img, _x, _y, _color ) \
(	sprintf( TMP_CMD, "%s put %s -to %d %d", _img, _color, _x, _y ), \
	Tcl_Eval( _itp, TMP_CMD ) \
)

#define DRAW_IMAGE_RECT( _itp, _img, _x1, _y1, _x2, _y2, _color ) \
(	sprintf( TMP_CMD, "%s put %s -to %d %d %d %d", \
		_img, _color, _x1, _y1, _x2, _y2 ), \
	Tcl_Eval( _itp, TMP_CMD ) \
)

#define DRAW_IMAGE_OUTLINE( _itp, _img, _x1, _y1, _x2, _y2, _color ) \
(	DRAW_IMAGE_RECT( _itp, _img, _x1, _y1, _x2, (_y1) + 1, _color ), \
	DRAW_IMAGE_RECT( _itp, _img, _x2, _y1, (_x2) + 1, _y2, _color ), \
	DRAW_IMAGE_RECT( _itp, _img, _x2, _y2, _x1, (_y2) + 1, _color ), \
	DRAW_IMAGE_RECT( _itp, _img, _x1, _y2, (_x1) + 1, _y1, _color ) \
)


/* Rectangles */

#define DRAW_RECT( _itp, _c, _x1, _y1, _x2, _y2, _fill, _cl, _wt, _t ) \
(	sprintf( TMP_COORD1, "-fill \"%s\"", _fill ), \
	sprintf( TMP_COORD2, "-outline \"%s\"", _cl ), \
	sprintf( TMP_COORD3, "-width %d", _wt ), \
	sprintf( TMP_COORD4, "-tag %s", _t ), \
	sprintf( TMP_CMD, "%s create rectangle %d %d %d %d %s %s %s %s", \
		_c, _x1, _y1, _x2, _y2, \
		TMP_COORD1, TMP_COORD2, TMP_COORD3, TMP_COORD4 ), \
	Tcl_Eval( _itp, TMP_CMD ), \
	atoi( Tcl_GetStringResult(_itp) ) \
)

#define RECT_COORDS( _itp, _c, _id, _x1, _y1, _x2, _y2 ) \
(	sprintf( TMP_CMD, "%s coords %d %d %d %d %d", \
		_c, _id, _x1, _y1, _x2, _y2 ), \
	Tcl_Eval( _itp, TMP_CMD ) \
)

#define RECT_COLOR( _itp, _c, _id, _fill, _color ) \
(	sprintf( TMP_CMD, \
		"%s itemconfigure %d -fill \"%s\" -outline \"%s\"", \
		_c, _id, _fill, _color ), \
	Tcl_Eval( _itp, TMP_CMD ) \
)


/* Grid Lines */

#define DRAW_LINE( _itp, _c, _x1, _y1, _x2, _y2, _fill, _wt ) \
(	sprintf( TMP_COORD1, "-fill \"%s\"", _fill ), \
	sprintf( TMP_COORD2, "-width %d", _wt ), \
	sprintf( TMP_CMD, "%s create line %d %d %d %d %s %s", \
		_c, _x1, _y1, _x2, _y2, TMP_COORD1, TMP_COORD2 ), \
	Tcl_Eval( _itp, TMP_CMD ), \
	atoi( Tcl_GetStringResult(_itp)) \
)

#define LINE_COORDS( _itp, _c, _id, _x1, _y1, _x2, _y2 ) \
(	sprintf( TMP_CMD, "%s coords %d %d %d %d %d", \
		_c, _id, _x1, _y1, _x2, _y2 ), \
	Tcl_Eval( _itp, TMP_CMD ) \
)

#define LINE_COLOR( _itp, _c, _id, _fill ) \
(	sprintf( TMP_CMD, \
		"%s itemconfigure %d -fill \"%s\"", _c, _id, _fill ), \
	Tcl_Eval( _itp, TMP_CMD ) \
)


/* Text */

#define DRAW_TEXT( _itp, _c, _x, _y, _str, _color, _font ) \
(	sprintf( TMP_CMD, \
		"%s create text %d %d -text {%s} -fill {%s} -font {%s} %s", \
		_c, _x, _y, _str, _color, _font, "-anchor nw" ), \
	Tcl_Eval( _itp, TMP_CMD ), \
	atoi( Tcl_GetStringResult(_itp)) \
)


/* Coordinate / Angle / Scale Macros */

#define X_OF_ANGLE( _x, _y, _angle, _scale ) \
	( (int) ( (_scale) * ( ( (double) (_y) * sin( _angle ) ) \
		- ( (double) (_x) * cos( _angle ) ) ) ) )

#define Y_OF_ANGLE( _x, _y, _angle, _scale ) \
	( (int) ( (_scale) * ( ( (double) (_y) * cos( _angle ) ) \
		+ ( (double) (_x) * sin( _angle ) ) ) ) )

#define ANGLE_COORDS( _x, _y, _x0, _y0, _angle, _scale ) \
	(_x0) + X_OF_ANGLE( (_x) - (_x0), (_y) - (_y0), _angle, _scale ), \
	(_y0) + Y_OF_ANGLE( (_x) - (_x0), (_y) - (_y0), _angle, _scale )


#define RAND_RAD( _rad ) \
	( ( 3 * (_rad) / 4 ) \
		+ ( (int) random() \
			% ( ( (_rad) / 2 ) ? ( (_rad) / 2 ) : 1 ) ) )

#define RADX( _rad, _angle )    ( (_rad) * cos( _angle ) )

#define RADY( _rad, _angle )    ( (_rad) * sin( _angle ) )


/* Raise Macro */

#define RAISE_ITEM( _itp, _c, _id ) \
(	sprintf( TMP_CMD, "%s raise %d", _c, _id ), \
	Tcl_Eval( _itp, TMP_CMD ) \
)


/* Un-Drawing Macro */

#define DELETE_ITEM( _itp, _c, _id ) \
(	sprintf( TMP_CMD, "%s delete %d", _c, _id ), \
	Tcl_Eval( _itp, TMP_CMD ) \
)


#ifdef MEM_DEBUG

/******* Debug Versions ************/

#define MALLOC( _d ) \
	( printf( "malloc size = %lu\n", (unsigned long) _d ), \
		TMP_PTR = (void *) Tcl_Alloc( _d ), \
		printf( "\tmem returned = 0x%lx\n", (long) TMP_PTR ), \
		fflush( stdout ), \
		TMP_PTR )

#define FREE( _c ) \
	( printf( "free ptr = 0x%lx\n", (long) _c ), Tcl_Free( (char*)_c ) )

#else

/****** Non-Debug Versions ********/

#define MALLOC( _d )    Tcl_Alloc( _d )

#define FREE( _c )    Tcl_Free( (char*)_c )

#endif

