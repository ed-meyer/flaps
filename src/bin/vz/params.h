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
 * Xgraph parameters
 */

#ifndef _PARAMS_H_
#define _PARAMS_H_

#include "plotcurve.h"

typedef enum param_types_defn {
    INT, STR, PIXEL, FONT, STYLE, BOOL, DBL 
} param_types;

typedef struct params_int_defn {
    param_types type;		/* INT */
    int value;
} param_int;

typedef struct params_str_defn {
    param_types type;		/* STR */
    char *value;
} param_str;

typedef struct params_pix_defn {
    param_types type;		/* PIXEL */
    XColor value;
} param_pix;

typedef struct params_font_defn {
    param_types type;		/* FONT */
    XFontStruct *value;
} param_font;

typedef struct params_style_defn {
    param_types type;		/* STYLE */
    int len;
    char *dash_list;
} param_style;

typedef struct params_bool_defn {
    param_types type;		/* BOOL */
    int value;
} param_bool;

typedef struct params_dbl_defn {
    param_types type;		/* DBL */
    double value;
} param_dbl;

typedef union params_defn {
    param_types type;
    param_int intv;		/* INT */
    param_str strv;		/* STR */
    param_pix pixv;		/* PIXEL */
    param_font fontv;		/* FONT */
    param_style stylev;		/* STYLE */
    param_bool boolv;		/* BOOL */
    param_dbl dblv;		/* DBL */
} params;

void param_init(Display *disp, Colormap cmap);
void param_set(const char *name, param_types type, const char *val);
void param_reset(const char *name, const char *val);
params* param_get(char const* name, params *val);
void param_dump();

int stricmp (char const* a, char const* b);

/* Some convenience macros */

extern params param_temp, *param_temp_ptr;
extern XColor param_null_color;
extern param_style param_null_style;

#define PM_INT(name)	\
((param_temp_ptr = param_get(name, &param_temp)) ? \
 param_temp_ptr->intv.value : \
 (abort(), (int) 0))

char* pm_str (char const* name);

#define PM_COLOR(name)	\
((param_temp_ptr = param_get(name, &param_temp)) ? \
 param_temp_ptr->pixv.value : \
 (abort(), param_null_color))

Pixel PM_PIXEL (char const* name);
XFontStruct* PM_FONT (char const* name);
double PM_DBL (char const* name);

#define PM_STYLE(name)	\
((param_temp_ptr = param_get(name, &param_temp)) ? \
 param_temp_ptr->stylev : \
 (abort(), param_null_style))

#define PM_BOOL(name)	\
((param_temp_ptr = param_get(name, &param_temp)) ? \
 param_temp_ptr->boolv.value : \
 (abort(), 0))


#endif /* _PARAMS_H_ */
