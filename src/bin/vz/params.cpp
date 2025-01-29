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
 * Xgraph Parameters
 *
 * This file contains routines for setting and retrieving
 * various X style display parameters for xgraph.
 */

#include <cstdlib>		/* added 1/31/92 EEM  for atof() */

#include "config.h"
#include "params.h"

using namespace std;

//unused:     static int strihash(char* string, int modulus);

/* For use by convenience macros */
params param_temp, *param_temp_ptr;
XColor param_null_color = { 
	0, 0, 0, 0, 0, 0 };
param_style param_null_style = { 
	STYLE, 0, (char *) 0 };

class param_full {
public:
	param_types type;
	string text_form;
	params* real_form;

	param_full(param_types t, string const& tf, params* rf) :
		type(t), text_form(tf), real_form(rf) {}
}; // param_full

// static st_table *param_table = (st_table *) 0;
map<string,param_full*> param_table;
typedef map<string,param_full*>::iterator PTIter;



static Display *param_disp = 0;
static Colormap param_cmap = 0;
static int param_scrn = 0;

static void free_resource (params *val);
static params *resolve_entry( string const& name, param_types type, string const& form);

static int do_color (char const* name, XColor *color);
static int do_font (char const* name, XFontStruct **font_info);
static int do_style (char const* list, param_style *val);
static int do_bool (char const* name, int *val);

#define DEF_INT		"0"
#define DEF_STR		""
#define DEF_FONT	"fixed"
#define DEF_PIXEL	"black"
#define DEF_STYLE	"1"
#define DEF_BOOL	"false"
#define DEF_DBL		"0.0"

#define DEF_MAX_FONT	1024
#define DEF_MAX_NAMES	10

#define DUP(str) strdup(str)


void
param_init (Display *disp, Colormap cmap) {
/*
 * Initializes parameter package.  The display and colormap arguments
 * are used to resolve font and pixel values.
 */
	// XXX param_table = st_init_table(stricmp, strihash);
	param_disp = disp;
	param_cmap = cmap;
	/* This could also be a parameter for greater generality */
	param_scrn = DefaultScreen(disp);
}



void
param_set(
const char *nm,			/* Name of parameter   */
param_types type,		/* Type                */
const char *val)			/* Text form for value */
/*
 * Sets the parameter with the given name to have the type
 * `type' and the text value `value'.  This will be evaluated
 * to its full form the first time it is referenced using
 * param_get().  If it is already filled, the old value
 * will be reclaimed.
 */
{
	param_full *entry;
	string name(nm);

//XXX	if (!param_table) {
//XXX		(void) fprintf(stderr, "Parameter table not initialized\n");
//XXX		return;
//XXX	}
	PTIter pos = param_table.find(name);
	if (pos == param_table.end()) {
		entry = new param_full(type, val, 0);
		param_table.insert(make_pair(name,entry));
	} else {
		entry = pos->second;
		if (entry->real_form) free_resource(entry->real_form);
		entry->real_form = (params *) 0;
		entry->type = type;
		entry->text_form = val;
	}
}


void param_reset(
const char *name,			/* Name of parameter   */
const char *val)			/* Text form for value */
/*
 * This routine sets the value of an existing parameter to a new
 * value.  The type of the parameter remains the same.  Changes
 * in type should be done by using param_set() directly.
 */
{
	param_full *entry;

//XXX	if (!param_table) {
//XXX		(void) fprintf(stderr, "Parameter table not initialized\n");
//XXX		return;
//XXX	}
	PTIter pos = param_table.find(name);
	if (pos == param_table.end()) {
		(void) fprintf(stderr, "Cannot reset unknown parameter `%s'\n", name);
	} else {
		entry = pos->second;
		param_set(name, entry->type, val);
	}
}




params*
param_get(
char const* nm,			/* Name of parameter */
params *val)			/* Result value      */
/*
 * Retrieves a value from the parameter table.  The value
 * is placed in `val'.  If successful, the routine will
 * return `val'.  Otherwise, it will return zero.
 */
{
	param_full *entry;
	string name(nm);

//XXX	if (!param_table) {
//XXX		(void) fprintf(stderr, "Parameter table not initialized\n");
//XXX		return (params *) 0;
//XXX	}
	PTIter pos = param_table.find(name);
	if (pos == param_table.end()) {
		return 0;
	} else {
		entry = pos->second;
		if (!entry->real_form) {
			entry->real_form = resolve_entry(name, entry->type, entry->text_form);
		}
		*val = *(entry->real_form);
		return val;
	}
}


static void
free_resource (params *val)
/*
 * Reclaims a resource based on its type.
 */
{
	switch (val->type) {
	case INT:
	case STR:
	case BOOL:
	case DBL:
		/* No reclaimation necessary */
		break;
	case PIXEL:
		if ((val->pixv.value.pixel != WhitePixel(param_disp, param_scrn)) &&
		    (val->pixv.value.pixel != BlackPixel(param_disp, param_scrn))) {
			XFreeColors(param_disp, param_cmap, &(val->pixv.value.pixel), 1, 0);
		}
		break;
	case FONT:
		XFreeFont(param_disp, val->fontv.value);
		break;
	case STYLE:
		(void) free(val->stylev.dash_list);
		break;
	}
	(void) free((char *) val);
}



static params *
resolve_entry( string const& name, param_types type, string const& form) {
/*
 * Allocates and returns an appropriate parameter structure
 * by translating `form' into its native type as given by `type'.
 * If it can't translate the given form, it will fall back onto
 * the default.
 */
	params *result = (params *) calloc(1, sizeof(params));

	result->type = type;
	switch (type) {
	case INT:
		if (sscanf(form.c_str(), "%d", &result->intv.value) != 1) {
			(void) fprintf(stderr,
			    "Parameter %s: can't translate `%s' into an integer (defaulting to `%s')\n",
			    name.c_str(), form.c_str(), DEF_INT);
			result->intv.value = atoi(DEF_INT);
		}
		break;
	case STR:
		result->strv.value = strdup(form.c_str());
		break;
	case PIXEL:
		if (!do_color(form.c_str(), &result->pixv.value)) {
			(void) fprintf(stderr, "Parameter %s: can't translate `%s' into a color (defaulting to `%s')\n",
			    name.c_str(), form.c_str(), DEF_PIXEL);
			(void) do_color(DEF_PIXEL, &result->pixv.value);
		}
		break;
	case FONT:
		if (!do_font(form.c_str(), &result->fontv.value)) {
			(void) fprintf(stderr, "Parameter %s: can't translate `%s' into a font (defaulting to `%s')\n",
			    name.c_str(), form.c_str(), DEF_FONT);
			(void) do_font(DEF_FONT, &result->fontv.value);
		}
		break;
	case STYLE:
		if (!do_style(form.c_str(), &result->stylev)) {
			(void) fprintf(stderr, "Parameter %s: can't translate `%s' into a line style (defaulting to `%s')\n",
			    name.c_str(), form.c_str(), DEF_STYLE);
			(void) do_style(DEF_STYLE, &result->stylev);
		}
		break;
	case BOOL:
		if (!do_bool(form.c_str(), &result->boolv.value)) {
			(void) fprintf(stderr, "Parameter %s: can't translate `%s' into a boolean flag (defaulting to `%s')\n",
			    name.c_str(), form.c_str(), DEF_BOOL);
			(void) do_bool(DEF_BOOL, &result->boolv.value);
		}
		break;
	case DBL:
		if (sscanf(form.c_str(), "%lf", &result->dblv.value) != 1) {
			(void) fprintf(stderr,
			    "Parameter %s: can't translate `%s' into a double (defaulting to `%s')\n",
			    name.c_str(), form.c_str(), DEF_DBL);
			result->dblv.value = atof(DEF_DBL);
		}
		break;
	}
	return result;
}



static int
do_color(
char const* name,			/* Name for color */
XColor *color)			/* Returned color */
/*
 * Translates `name' into a color and attempts to get the pixel
 * for the color using XAllocColor().
 */
{
	int result = 1;

	if (XParseColor(param_disp, param_cmap, name, color)) {
		if (stricmp(name, "black") == 0) {
			color->pixel = BlackPixel(param_disp, param_scrn);
			XQueryColor(param_disp, param_cmap, color);
		} else if (stricmp(name, "white") == 0) {
			color->pixel = WhitePixel(param_disp, param_scrn);
			XQueryColor(param_disp, param_cmap, color);
		} else {
			result = XAllocColor(param_disp, param_cmap, color);
		}
	} else {
		result = 0;
	}
	return result;
}



static int
do_font(
char const* name,			/* Name of desired font      */
XFontStruct **font_info)	/* Returned font information */
/*
 * This routine translates a font name into a font structure.  The
 * font name can be in two forms.  The first form is <family>-<size>.
 * The family is a family name (like helvetica) and the size is
 * in points (like 12).  If the font is not in this form, it
 * is assumed to be a regular X font name specification and
 * is looked up using the standard means.
 */
{
	char name_copy[DEF_MAX_FONT], query_spec[DEF_MAX_FONT];
	char *font_family, *font_size, **font_list;
	int font_size_value, font_count = 0, i;

	/* First attempt to interpret as font family/size */
	(void) strcpy(name_copy, name);
	if ((font_size = strchr(name_copy, '-'))) {
		*font_size = '\0';
		font_family = name_copy;
		font_size++;
		font_size_value = atoi(font_size);
		if (font_size_value > 0) {
			/* Still a little iffy -- what about weight and roman vs. other */
			(void) sprintf(query_spec,
			    "*-*-%s-medium-r-normal-*-*-%d-*-*-*-*-iso8859-*",
			    font_family, font_size_value * 10);
			font_list = XListFonts(param_disp, query_spec,
			    DEF_MAX_NAMES, &font_count);

			/* Load first one that you can */
			for (i = 0;  i < font_count;  i++) {
				if ((*font_info = XLoadQueryFont(param_disp, font_list[i]))) {
					break;
				}
			}
			if (*font_info) return 1;
		}
	}
	/* Assume normal font name */
	*font_info = XLoadQueryFont(param_disp, name);
	if (*font_info) return 1;
	return 0;
}


static int do_style(
char const* list,			/* List of ones and zeros */
param_style *val)		/* Line style returned    */
/*
 * Translates a string representation of a dash specification into
 * a form suitable for use in XSetDashes().  Assumes `list'
 * is a null terminated string of ones and zeros.
 */
{
	char const* i;
	char* spot;
	char last_char;
	int count;

	for (i = list;  *i;  i++) {
		if ((*i != '0') && (*i != '1')) break;
	}
	if (!*i) {
		val->len = 0;
		last_char = '\0';
		for (i = list;  *i;  i++) {
			if (*i != last_char) {
				val->len += 1;
				last_char = *i;
			}
		}
		val->dash_list = (char *) calloc(1, (unsigned) (sizeof(char)*val->len + 1));
		last_char = *list;
		spot = val->dash_list;
		count = 0;
		for (i = list;  *i;  i++) {
			if (*i != last_char) {
				*spot++ = (char) count;
				last_char = *i;
				count = 1;
			} else {
				count++;
			}
		}
		*spot = (char) count;
		return 1;
	} else {
		return 0;
	}
}


static const char *positive[] = { 
	"on", "yes", "true", "1", "affirmative", (char *) 0 };
static const char *negative[] = { 
	"off", "no", "false", "0", "negative", (char *) 0 };

static int
do_bool(
char const* name,			/* String representation */
int *val)			/* Returned value        */
/*
 * Translates a string representation into a suitable binary value.
 * Can parse all kinds of interesting boolean type words.
 */
{
	const char **term;

	for (term = positive;  *term;  term++) {
		if (stricmp(name, *term) == 0) break;
	}
	if (*term) {
		*val = 1;
		return 1;
	}
	for (term = negative;  *term;  term++) {
		if (stricmp(name, *term) == 0) break;
	}
	if (*term) {
		*val = 0;
		return 1;
	}
	return 0;
}


void param_dump()
/*
 * Dumps all of the parameter values to standard output.
 */
{
	//XXX st_foreach(param_table, dump_it, (char *) 0);
}



int
stricmp(char const* a, char const* b)
/*
 * This routine compares two strings disregarding case.
 */
{
	int value;

	if ((a == (char *) 0) || (b == (char *) 0)) {
		return a - b;
	}

	for ( /* nothing */;
	    ((*a | *b) &&
	    !(value = ((isupper(*a) ? *a - 'A' + 'a' : *a) -
	    (isupper(*b) ? *b - 'A' + 'a' : *b))));
	    a++, b++)
		/* Empty Body */;

	return value;
}

char*
pm_str(char const* name) {
	param_temp_ptr = param_get(name, &param_temp);
	if (!param_temp_ptr) {
		ostringstream os;
		os << name << " has not been defined";
		throw runtime_error(os.str());
	}
	return param_temp_ptr->strv.value;
}

Pixel
PM_PIXEL(char const* name) {
	param_temp_ptr = param_get(name, &param_temp);
	if (!param_temp_ptr) {
		ostringstream os;
		os << name << " has not been defined";
		throw runtime_error(os.str());
	}
	return param_temp_ptr->pixv.value.pixel;
}

XFontStruct*
PM_FONT(char const* name) {
	param_temp_ptr = param_get(name, &param_temp);
	if (!param_temp_ptr) {
		ostringstream os;
		os << name << " has not been defined";
		throw runtime_error(os.str());
	}
	return param_temp_ptr->fontv.value;
}

double
PM_DBL(char const* name) {
	param_temp_ptr = param_get(name, &param_temp);
	if (!param_temp_ptr) {
		ostringstream os;
		os << name << " has not been defined";
		throw runtime_error(os.str());
	}
	return param_temp_ptr->dblv.value;
}
