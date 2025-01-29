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
 * Hardcopy Devices
 *
 * This file contains the basic output device table.  The hardcopy
 * dialog is automatically constructed from this table.
 */

#include "config.h"
#include "hard_devices.h"
#include "xfig.h"
#include "params.h"
#include "ps.h"


struct hard_dev hard_devices[] = {
	{ "xfig", xfigInit, 0,
	"vz.fig", "xfig", 19.0, "Times-Bold", 18.0,
	"Times-Roman", 12.0, NONE },
	{ "Postscript", psInit, "%s", "vz.ps", "lpr",
	19.0, "Times-Bold", 18.0, "Times-Roman", 12.0, NO }
};


int hard_count = sizeof(hard_devices)/sizeof(struct hard_dev);

#define CHANGE_D(name, field) \
if (param_get(name, &val)) { \
    if (val.type == DBL) { \
       hard_devices[idx].field = val.dblv.value; \
    } \
}

#define CHANGE_S(name, field) \
if (param_get(name, &val)) { \
    if (val.type == STR) { \
       (void) strcpy(hard_devices[idx].field, val.strv.value); \
    } \
}


void hard_init()
/*
 * Changes values in hard_devices structures in accordance with
 * parameters set using the parameters module.
 */
{
	char newname[1024];
	int idx;
	params val;

	for (idx = 0;  idx < hard_count;  idx++) {
		(void) sprintf(newname, "%s.Dimension", hard_devices[idx].dev_name);
		CHANGE_D(newname, dev_max_dim);
		(void) sprintf(newname, "%s.OutputTitleFont", hard_devices[idx].dev_name);
		CHANGE_S(newname, dev_title_font);
		(void) sprintf(newname, "%s.OutputTitleSize", hard_devices[idx].dev_name);
		CHANGE_D(newname, dev_title_size);
		(void) sprintf(newname, "%s.OutputAxisFont", hard_devices[idx].dev_name);
		CHANGE_S(newname, dev_axis_font);
		(void) sprintf(newname, "%s.OutputAxisSize", hard_devices[idx].dev_name);
		CHANGE_D(newname, dev_axis_size);
	}
}
