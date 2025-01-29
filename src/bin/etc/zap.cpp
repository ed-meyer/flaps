//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#include <cctype>
#include <cerrno>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>

char *Progname;
char const* Pscmd = "ps";

static int ttyin(void);
static int strindex(char *s, char *t);
static FILE *efopen(char const* file, char const* mode);

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
/*------------------------------------------------------------------
 * zap:   kill processes by name
 *        uses SIGTERM unless there is a -n on the cmdline
 *------------------------------------------------------------------*/

int
main (int argc, char **argv) {
   FILE *fin;
   char buf[BUFSIZ];
	char cmd[128];
   int pid;
	int start = 1;
	char *cp;
	char* group = 0;
	int sig = SIGTERM;
	int force = FALSE;

   Progname = argv[0];

	if (argc == 1) {
		fprintf (stderr,"Usage: %s [-n][-f] name [name ...]\n", Progname);
		fprintf (stderr,"   -n:  signal number (see /usr/include/sys/csignal\n");
		fprintf (stderr,"   -f:  do not prompt for permission to kill\n");
		return(0);
	}

	cp = argv[start];
	while (*cp == '-') {
		start++;
		cp++;
		if (isdigit(*cp)) {
			sig = atoi(cp);
		} else if (*cp == 'f') {
			force = TRUE;
		} else if (*cp == 'g') {
			group = argv[start++];
			sprintf (cmd, "%s -g %s", Pscmd, group);
			fprintf (stderr,"killing process group %s (%s)\n",group,cmd);
		} else {
			fprintf (stderr,"%s: illegal option (%s)\n", argv[0], argv[1]);
			return(1);
		}
		cp = argv[start];
	}

	if (group) {
#ifdef LINUX
		sprintf (cmd, "%s g %s", Pscmd, group);
#else
		sprintf (cmd, "%s -g %s", Pscmd, group);
#endif
	} else {
#ifdef LINUX
		sprintf (cmd, "%s e", Pscmd);
#else
		sprintf (cmd, "%s -e", Pscmd);
#endif
	}

	fprintf (stderr,"killing: %s\n",cmd);
		
	/*
	 * Open a pipe to the ps command and read its output
	 */
   if ((fin = popen(cmd, "r")) == NULL) {
      fprintf (stderr,"%s: can't run %s\n", Progname, cmd);
      return(1);
   }

   fgets(buf, sizeof buf, fin);
   fprintf (stderr,"%s", buf);
   while (fgets(buf, sizeof buf, fin) != NULL) {
      if (strindex(buf, argv[start]) >= 0) {
			buf[strlen(buf)-1] = '\0';     /* suppress \n */
			if (force) {
				sscanf(buf, "%d", &pid);
				kill(pid, sig);
				fprintf (stderr,"%s\n",buf);
			} else {
				fprintf (stderr,"%s? ", buf);
				if (ttyin() == 'y') {
					sscanf(buf, "%d", &pid);
					kill(pid, sig);
				}
			}
      }
	}
   return(0);
}

static int
ttyin() {
   char buf[BUFSIZ];
   static FILE *tty = NULL;

   if (tty == NULL)
      tty = efopen("/dev/tty", "r");
   for (;;) {
      if (fgets(buf, BUFSIZ, tty) == NULL || buf[0] == 'q')
         return(0);
      else if (buf[0] == '!') {
         system(buf+1);
         printf ("!\n");
      }
      else
         return buf[0];
   }
}

static int
strindex(char *s, char *t) {
   int i, n;

   n = strlen(t);
   for (i = 0; s[i] != '\0'; i++)
      if (strncmp(s+i, t, n) == 0)
         return i;
   return -1;
}

static FILE *
efopen(char const* file, char const* mode) {
   FILE *fp;

   fp = fopen(file, mode);
   if (fp == NULL) {
		fprintf (stderr,"%s: can't open file %s mode %s: %s\n",
					Progname, file, mode, strerror(errno));
		exit(1);
	}
	return fp;
}
