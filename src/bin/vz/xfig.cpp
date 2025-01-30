//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
/*
 * xfig Output
 *
 * Beorn Johnson
 * Alan Kramer
 * David Harrison
 */


#include "config.h"
#include "xgout.h"

#define HEIGHT	792
#define FIX(X)	X = HEIGHT - X;

typedef struct {
	char	*title_font;
	char	*axis_font;
	int	title_size;
	int	axis_size;
	FILE	*strm;
} Info;

const char	*xfig_prologue[ ] = {
	"#FIG 3.2 Produced by Flaps vz",
	"Landscape",
	"Center",
	"Inches",
	"Letter",
	"100.0",
	"Single",
	"-2",
	"1200 2",       // fig_units/inch coord_system
	0
};

/*
 * Hardcopy Interface for Xgraph
 *
 * Major differences from first version:
 *   Four new parameters are passed to the device initialization routine:
 *   title_family, title_size, axis_family, and axis_size.  See the 
 *   description of xg_init() for details.
 *   
 *   Clipping is done automatically by xgraph.  The xg_clip() routine
 *   is obsolete.
 *
 *   The xg_line() routine has become the xg_seg() routine.  It now
 *   draws segments rather than a series of lines.  
 * 
 *   A new field (max_segs) in the device structure now specifies 
 *   the maximum number of segments the device can handle in a group.
 */


/*
 * Adding an output device to xgraph
 *
 * Step 1
 *   Write versions of the following routines for your device:
 *   xg_init(), xg_text(), xg_seg(), xg_dot(), and xg_end().
 *   The interface and function of these routines are described
 *   in detail below.  These routines should be named according
 *   to your device.  For example,  the initialization routine
 *   for the Postscript output device is psInit().  Also,  name
 *   your source file after your device (e.g. the postscript
 *   routines are in the file ps.c).  Instructions continue
 *   after the description of the interface routines.
 */


void xfigText(char*,int,int,char const*,int,int);
void xfigDot(char*,int,int,int,int,int);
void xfigSeg(char*,int,XSegment*,int,int,int,int);
void xfigEnd(char*);

int
xfigInit(
FILE *strm,			/* Output stream              */
int width, int height,		/* Size of space (microns)    */
char *title_family,		/* Name of title font family  */
double title_size,		/* Title font height (points) */
char *axis_family,		/* Name of axis font family   */
double axis_size,		/* Axis font height (points)  */
int flags,		        /* Flags                      */
xgOut *out_info,		/* Device info (RETURN)       */
char errmsg[ERRBUFSIZE])	/* Error message area         */
/*
 * This routine is called by xgraph just before drawing is to
 * begin.  The desired size of the plot is given by `width'
 * and `height'.  The parameters `title_family', `title_size',
 * `axis_family', and `axis_size' specify the names of the
 * title and axis fonts and their vertical sizes (in points).
 * These parameters can be ignored if your device does not
 * support multiple fonts.  Binary flags are specified in
 * the `flags' field.  These include:
 *  D_DOCU:
 *      If this flag is set,  it indicates the user has specified that
 *	the output will be included in some larger document.  Devices
 *	may choose to use this information to produce output that
 *	can be integrated into documents with less effort.  For example,
 *	the Postscript output routines produce bounding box information
 *	when this flag is set.
 * The routine should fill in all of the fields of `out_info' with 
 * appropriate values.  The values are described below:
 *  area_w, area_h:  
 * 	Size of the drawing space in device coordinates.
 *	This should take in account the requested area
 *	given by `width', and `height'.
 *  bdr_pad:  
 * 	Xgraph will leave this number of device coordinates around
 *	all of the outer edges of the graph.
 *  axis_pad: 
 *	Additional space around axis labels (in devcoords)
 *	so that the labels do not appear crowded.
 *  legend_pad:
 *	Space (in devcoords) from the top of legend text to
 *	the representative line drawn above the legend text.
 *  tick_len:    
 *	Size of a tick mark placed on axis (in devcoords)
 *  axis_width:  
 *	An estimate of the width of a large character in
 *      the axis font (in devcoords).  This can be an overestimate.  An
 *      underestimate may produce bad results.
 *  axis_height: 
 *	An estimate of the height of a large character in
 *      the axis labeling font (in devcoords).
 *  title_width, title_height:  
 *	Same as above except for the title font.
 *  max_segs:
 *	Due to buffering constraints,  some devices may not be able to
 *	handle massive segment lists.  This parameter tells xgraph not
 *	to send more than `max_segs' segments in one request.
 * Output to the device should be written to the stream `strm'.
 * The functions are described individually below.  After filling
 * in the parameters and setting the function pointers,  the routine
 * should initialize its drawing state and store any extra needed
 * information in `user_state'.  This value will be passed to all
 * other routines during the drawing sequence.  If the device
 * cannot initialize,  it should return a zero status and fill
 * `errmsg' with an informative error message.
 */
{
	Info	*xfig_info;
	// const char	**l;
	double	scx, scy;

	xfig_info = (Info *) calloc(1, sizeof(*xfig_info));

	for (const char** l = xfig_prologue; *l; l++)
		fprintf(strm, "%s\n", *l);

	out_info -> dev_flags = 0;
	//scx = width / 612;
	//scy = height / 792.0;
	scx = width / 200;
	scy = height / 300.0;
	if (scx > scy) {
		scy /= scx;
		scx = 1;
	}
	else {
		scx /= scy;
		scy = 1;
	}
	out_info -> bdr_pad = (int)title_size/4;
	out_info -> axis_pad = (int)(2.0 * axis_size);
	out_info -> legend_pad = 0;

	//out_info -> area_w = (int)(width * 0.00283); /* pts per micron */
	//out_info -> area_h = (int)(height * 0.00283);
	// make figures very large in xfig to get smooth curves:
	out_info -> area_w = (int)(width * 0.0283); /* pts per micron */
	out_info -> area_h = (int)(height * 0.0283);

	out_info -> tick_len = (int)axis_size;
	out_info -> axis_height = (int)axis_size;
	out_info -> title_height = (int)title_size;
	out_info -> axis_width = (int)((axis_size*5.0) / 12.0);
	out_info -> title_width = (int)((title_size*5.0) / 12.0);
	out_info -> max_segs = 1000000;
	out_info -> xg_text = xfigText;
	out_info -> xg_seg = xfigSeg;
	out_info -> xg_dot = xfigDot;
	out_info -> xg_end = xfigEnd;

	out_info -> user_state = (char *) xfig_info;

	xfig_info -> title_font = title_family;
	xfig_info -> axis_font = axis_family;
	xfig_info -> title_size = (int)title_size;
	xfig_info -> axis_size = (int)axis_size;
	xfig_info -> strm = strm;
	return 1;
}

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


void
xfigText(
char *user_state,		/* Value set in xg_init   */
int x, int y,			/* Text position (pixels) */
char const* text,			/* Null terminated text   */
int just,			/* Justification (above)  */
int style)			/* Text style (above)     */
/*
 * This routine should draw text at the indicated position using
 * the indicated justification and style.  The justification refers
 * to the location of the point in reference to the text.  For example,
 * if just is T_LOWERLEFT,  (x,y) should be located at the lower left
 * edge of the text string.
 */
{
}

/* Line Styles */
#define L_AXIS		0
#define L_ZERO		1
#define L_VAR		2

void xfigSeg(
char *user_state,		/* Value set in xg_init */
int ns,				/* Number of segments   */
XSegment *seglist,		/* X array of segments  */
int width,			/* Width of lines       */
int style,			/* See above            */
int lappr,			/* Line appearence      */
int color)			/* Line color (if any)  */
/*
 * This routine draws a number of line segments at the points
 * given in `seglist'.  Note that contiguous segments need not share
 * endpoints but often do.  All segments should be `width' devcoords wide
 * and drawn in style `style'.  If `style' is L_VAR,  the parameters
 * `color' and `lappr' should be used to draw the line.  Both
 * parameters vary from 0 to 7.  If the device is capable of
 * color,  `color' varies faster than `style'.  If the device 
 * has no color,  `style' will vary faster than `color' and
 * `color' can be safely ignored.  However,  if the
 * the device has more than 8 line appearences,  the two can
 * be combined to specify 64 line style variations.
 * Xgraph promises not to send more than the `max_segs' in the
 * xgOut structure passed back from xg_init().
 */
{
	T_(Trace trc(1,"xfigSeg");)
	Info	*xfig_info = (Info *) user_state;
	int	i, j, k;
	FILE* stream = xfig_info->strm;

	T_(trc.dprint(ns," segments, style ",style,", color ",color,", width ",width);)

	// I don't know why but this fcn gets called with width=0
	if (width == 0) width = 1;

	// polyline header:
	// fprintf(stream, "2 1 0 %d 0 7 50 -1 -1 0.0 0 0 -1 0 0 %d\n", width, ns+1);
	fprintf (stream, "2 1 %d %d %d 7 50 -1 -1 0.0 0 0 -1 0 0 %d\n",
			style, width, color, ns+1);

	for (i = 0; i < ns; i = j) {
		for (j = i + 1; j < ns
		    && seglist[j - 1].x2 == seglist[j].x1
		    && seglist[j - 1].y2 == seglist[j].y1;
		    j++) ;

		for (k = i; k < j; k++)
			fprintf(stream, "%d %d\n", seglist[k].x1, seglist[k].y1);

		fprintf(stream, "%d %d\n", seglist[k-1].x2, seglist[k-1].y2);
	}
}

/* Marker styles */
#define P_PIXEL		0
#define P_DOT		1
#define P_MARK		2

void xfigDot(
char *user_state,		/* Value set in xg_init    */
int x, int y,			/* Location in pixel units */
int style,			/* Dot style               */
int type,			/* Type of marker          */
int color)			/* Marker color (if any)   */
/*
 * This routine should draw a marker at location `x,y'.  If the
 * style is P_PIXEL,  the dot should be a single pixel.  If
 * the style is P_DOT,  the dot should be a reasonably large
 * dot.  If the style is P_MARK,  it should be a distinguished
 * mark which is specified by `type' (0-7).  If the output
 * device is capable of color,  the marker should be drawn in
 * `color' (0-7) which corresponds with the color for xg_line.
 */
{
/*
 * Description of xfig data for circles:
 *
 *  type	name			(brief description)
 *  ----	----			-------------------
 *  int	object_code		(always 1)
 *  int	sub_type		(1: ellipse defined by radii
 *                    2: ellipse defined by diameters
 *   					    3: circle defined by radius
 *   					    4: circle defined by diameter)
 *  int	line_style		(enumeration type)
 *  int	thickness		(1/80 inch)
 *  int	pen_color		(enumeration type, pen color)
 *  int	fill_color		(enumeration type, fill color)
 *  int	depth			(enumeration type)
 *  int	pen_style		(pen style, not used)
 *  int	area_fill		(enumeration type, -1 = no fill)
 *  float	style_val		(1/80 inch)
 *  int	direction		(always 1)
 *  float	angle			(radians, the angle of the x-axis)
 *  int	center_x, center_y	(Fig units)
 *  int	radius_x, radius_y	(Fig units)
 *  int	start_x, start_y	(Fig units; the 1st point entered)
 *  int	end_x, end_y		(Fig units; the last point entered)
 */

#ifdef NEVER // disable for now
	Info	*xfig_info = (Info *) user_state;
	FILE* stream = xfig_info->strm;
	int sub_type = 3;
	int line_style = 0;
	int thickness = 1;
	int pen_color = 0;
	int fill_color = 0;
	int depth = 50;
	int pen_style = -1;
	int area_fill = 20;
	float style_val = 0.0;
	int direction = 1;
	float angle = 0.0;
	int radius_x = 15;
	int radius_y = 15;

	//fprintf(stream, "1 3 0 1 0 0 50 -1 20 0.0 1 0.0 %d %d 10 10 %d %d %d %d\n",
	//		x, y, x, y, x, y, x, y);
	fprintf(stream, "1 %d %d %d %d %d %d %d %d %f %d %f %d %d %d %d %d %d %d %d\n",
			sub_type,
			line_style,
			thickness,
			pen_color,
			fill_color,
			depth,
			pen_style,
			area_fill,
			style_val,
			direction,
			angle,
			x, y,
			radius_x,
			radius_y,
			x, y,
			x, y);
#endif // NEVER // disable for now

}

void xfigEnd(char* user_state)
/*
 * This routine is called after a drawing sequence is complete.
 * It can be used to clean up the user state and set the device
 * state appropriately.  This routine is optional in the structure.
 */
{
	// Info	*xfig_info = (Info *) user_state;

	// fprintf(idraw -> strm, "End %%I eop\n");
	// fclose(xfig_info -> strm);
}

/*
 * Adding an output device to xgraph
 *
 * Step 2
 *   Edit the file hard_devices.c.  Declare your initialization
 *   function and add your device to the list of devices,
 *   hard_devices[].  The structure hard_dev is described below:
 */


/*
 * dev_spec:
 *    The dev_spec field should be a command that directly outputs to
 *    your device.  The command should contain one %s directive that
 *    will be filled in with the name of the device from the hardcopy
 *    dialog.
 * dev_file:
 *    The default file to write output to if the user selects `To File'.
 * dev_printer:
 *    The default printer to write output to if the user selects
 *    `To Device'.
 * dev_max_dim:
 *    The default maximum dimension for the device in centimeters.
 * dev_title_font, dev_title_size:
 *    The default title font and size.  Sizes are specified in
 *    points (1/72 inch).
 * dev_axis_font, dev_axis_size:
 *    The default axis font and size.
 */

/*
 * Adding an output device to xgraph
 *
 * Step 3
 *   Edit the file Makefile.  Add your source file to the SRC variable
 *   and the corresponding object file to the OBJ variable.  Finally,
 *   remake xgraph.  Your device should now be available in the
 *   hardcopy dialog.
 */

