//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
#ifndef HARD_DEVICES_H
#define HARD_DEVICES_H
/*
 * Hardcopy Device Header
 *
 * This file declares the types required for the hardcopy table
 * found in hard_devices.c.
 */

#include "xgout.h"

#define MFNAME	25

typedef enum hard_dev_docu_defn { NONE, NO, YES } hard_dev_docu;

struct hard_dev {
    const char *dev_name;		/* Device name                */
    int (*dev_init)(FILE*,int,int,char*,double,char*,double,
	 	int,xgOut*,char*);		/* Initialization function    */
    const char *dev_spec;		/* Default pipe program       */
    char dev_file[MFNAME];	/* Default file name          */
    char dev_printer[MFNAME];	/* Default printer name       */
    double dev_max_dim;		/* Default maximum dimension (cm)    */
    char dev_title_font[MFNAME];/* Default name of title font        */
    double dev_title_size;	/* Default size of title font (pnts) */
    char dev_axis_font[MFNAME];	/* Default name of axis font         */
    double dev_axis_size;	/* Default size of axis font (pnts)  */
    hard_dev_docu dev_docu;	/* Document predicate                */
};

extern int hard_count;
extern struct hard_dev hard_devices[];

extern void hard_init();


#endif /* HARD_DEVICES_H */
