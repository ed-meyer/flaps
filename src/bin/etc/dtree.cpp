//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
// Copyright (C) 2024 Edward E. Meyer
// This file is part of Flaps; Flaps is free software: you can redistribute
// it and/or modify it under the terms of the GNU General Public License.
// See the file COPYING in the root directory.
// Flaps is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

#define _ALL_SOURCE 1			/* for getopt() */
/*-----------------------  dtree.c  ------------------------
 *
 *  Directory tree drawing program
 *  The original program appeared in Dr. Dobbs' Journal
 *  (Oct 86?) C Chest.
 *  Modified by Ed Meyer to work under Unix
 *  8/90 - modified to use readdir etc.
 *---------------------------------------------------------*/

#include "config.h"

#include <cassert>
#include <sys/stat.h>
#ifdef HAVE_VNODE_H
#endif
#ifdef HAVE_DIRENT_H
#include <dirent.h>
#endif
#include <cerrno>
#ifdef HAVE_LIMITS_H
#endif
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <unistd.h>

#define ERR strerror(errno)

#define D1(x)
#define D2(x)
#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define MAXFILES  1000      /* max files per subdirectory */

#ifndef PATH_MAX
#define PATH_MAX 1023
#endif

typedef struct {
	char *name;
	size_t size;
} Dir;

Dir* Largest = 0;
int NLargest = 0;
int KLargest = 0;
int MaxLargest = 32768;

char const* Branch[] = { "|", "+-----", "+-----" };
#define     VERT        Branch[0]
#define     ELL         Branch[1]
#define     T_RIGHT     Branch[2]


/*------------------------------------------------------
 *  Globals
 *-----------------------------------------------------*/
static char *Startdir;			/* starting directory */
char *Root = (char*)NULL;		/* root directory */
static char Map[64/8];          /* bit map */
size_t    Size;               /* directory size */
static int Print_size = FALSE;     /* default: no sizes */
static int Print_tree = TRUE;     /* default: print tree */
static int  Short_pname = 1;    /* default: short pathname  */
char *Progname;                   /* name of this program (argv[0]) */

/*------------------------------------------------------------------
 * local (static) routines
 *------------------------------------------------------------------*/
static void usage(void);
static void setBit(int c, int val);
static void pline (int depth, int terminate);
static void pname (char *dname);
static void prnt (char *dname, int others);
static char **sub_dir (char *dname, int *count);
static char *strsav(char *string);
static char *dodot (char *str);
static void cmdline(int argc, char **argv);

#define testbit(x)  (Map[x >> 3] & (1 << (x & 0x07)))

static void
setBit (int c, int val) {
    if (val)
        Map[c >> 3] |= 1 << (c & 0x07);
    else
        Map[c >> 3] &= ~(1 << (c & 0x07));
}

static void
pline (int depth, int terminate) {
    int i;
    
    for (i=0; i<depth-1; i++)
        printf (testbit(i) ? "%s     " : "      ",VERT);

    if (terminate) printf ("\n");
}

static void
pname (char *dname) {
    char    *name;

    if (!Short_pname || (dname[0] == '/' && !dname[1]))
        name = dname;
    else if ((name = strrchr(dname,'/'))) name++;
    else
        name = dname;

    printf ("%s",name);
    if (Print_size) printf ("(%dk)",(int)(Size/(long)1000));
    printf ("\n");
}

static void
prnt (char *dname, int others) {
   char            **vects;
   int             count;
   static int      depth = -1;
         
   D1(E("prnt(dname<%s>,others<%d>)\n", dname, others);)

   if (++depth && Print_tree) {
      pline (depth,0);
      printf ("%s",others ? T_RIGHT : ELL);
   }

   vects = sub_dir (dname, &count);

   D1(T("%d directories on %s\n",count,dname);)

	if (Print_tree)
		pname (dname);

   while (--count >= 0) {
      setBit (depth,count);
      D1(T("call prnt with <%s> <0x%x>",*vects,*vects);)
      prnt (*vects++,count);
   }
    
   if (!others && Print_tree)
		pline (depth,1);

   --depth;
   D1(X("prnt\n");)
}

/*-------------------------------------------------
 *  sub_dir:  set up a list of subdirectories
 *---------------------------------------------------*/
static char **
sub_dir (char *dname, int *count) {
   struct dirent   *entry;
   struct stat     stbuf;
   DIR             *dp;
   char            **dir, buf[PATH_MAX];

   D1(E("sub_dir to search %s\n",dname);)

   *count = 0;

   dir = (char **)malloc(sizeof(char *)*MAXFILES);
   assert(dir != (char**)NULL);

   if ((dp=opendir(dname)) == NULL) {
       fprintf (stderr,"can't open directory %s: %s",dname, ERR);
       return 0;
   }

    Size = 0;
			/*
			 * For each directory entry, call lstat(). lstat() differs
			 * from stat() if the entry is a symbolic link - it gives
			 * info on the link instead of the file the link references
			 */
    while ((entry = readdir(dp)) != NULL) {
      D2(T("dir entry: <%s>\n",entry->d_name);)
      if (entry->d_name[0] != '.' && strcmp(entry->d_name,"lost+found") != 0) {
         strcpy (buf,dname);
         strcat (buf,"/");
         strcat (buf,entry->d_name);
         if (lstat (buf,&stbuf) == -1) {
            D1(T("can't access %s\n",buf);)
         } else {
				D2(T("%s: mode %#o\n",buf,stbuf.st_mode);)
				switch(stbuf.st_mode & S_IFMT) {
					case S_IFDIR:
						D2(T("is link? %d\n",S_ISLNK(stbuf.st_mode));)
						dir[*count] = strsav(buf);
						D2(T("found directory %d <%s>\n",*count,dir[*count]);)
						(*count)++;
					case S_IFREG:
						Size += (long)stbuf.st_size;
						D2(T("Size = %ld\n",Size);)
				}
         }
      }
   }
   closedir(dp);
	if (Largest && KLargest < MaxLargest) {
		Largest[KLargest].name = strsav(dname);
		Largest[KLargest].size = Size;
		KLargest++;
	}
   D1(X("sub_dir: count<%d>\n",*count);)
   return dir;
}

static char *
strsav(char *string) {
    char *cp;

    cp = (char *)malloc(strlen(string)+1);
    assert(cp != (char *)NULL);
    strcpy (cp,string);
    return cp;
}
                     
/*-----------------------------------------------------
 *  dodot - expand dots in a pathname
 *-----------------------------------------------------*/
static char *
dodot (char *str) {
   char *rval, root[PATH_MAX];

   if (!strchr(str,'.')) return str;

   if (chdir(str) || !getcwd(root,PATH_MAX)) {
      fprintf (stderr,"%s: Can't find %s\n",Progname,str);
      exit(1);
   }
    
   chdir (Startdir);
	rval = strsav (root);
   return rval;
}
                                                        
/*---------------------------------------------------------
 *  optstring is a string of characters which determine
 *  what command line options ('-' followed by a letter)
 *  are accepted. The sys5 routine getopt is used to
 *  get the letters, so all options must appear first in
 *  the command line (before the plotfile name)
 *  A letter followed by a colon in optstring means that
 *  the next word on the command line will be taken as
 *  a value associated with that option.
 *--------------------------------------------------------*/
char const* optstring = "hlr:sS:";

/*-------------------------------------------------------------
 *  cmdline - parse the command line
 *------------------------------------------------------------*/
static void
cmdline(int argc, char *argv[]) {
    int             c;
    int             errflg=0;
    extern char     *optarg;
    extern int      optind;
    extern int      opterr;
    extern int      optopt;

    opterr = 0;                     /* set to 1 to quiet getopt */

   while ((c = getopt(argc,argv,optstring)) != EOF) {
		D1(T("got option \"%c\" (%c) arg<%s>\n",optopt,c,optarg?optarg:"none");)
      switch (c) {
			case 'h':
				usage();
				exit(0);
         case 'l':
            Short_pname = 0;
            break;
         case 'r':
            Root = strsav(optarg);
            break;
			case 'S':
				NLargest = atoi(optarg);
				Print_tree = FALSE;
         case 's':
            Print_size = 1;
				if (NLargest == 0)
					NLargest = 20;
				Largest = (Dir*)calloc(MaxLargest, sizeof(Dir));
				D1(T("print %d largest directories\n",NLargest);)
            break;
         case '?':
            fprintf (stderr,"%s: unrecognized cmd line arg (%c)\n",
					argv[0],optopt);
            errflg++;
      }
	}

	D1(T("cmdline: argc %d, optind %d\n",argc,optind);)
	if (optind < argc) {
		Root = strsav(argv[optind]);
		D1(T("got Root <%s>\n",Root);)
		if (*Root == '?' || *Root == '0' || strcmp(Root, "help") == 0) {
			usage();
			exit(0);
		}
	}

	if (errflg) {
		usage();
		exit(1);
	}
}

static void
usage() {
    fprintf (stderr,"Usage: %s [-l][-s][-S n][path]\n",Progname);
    fprintf (stderr,"  -l        - full pathnames\n");
    fprintf (stderr,"  -s        - give directory size (kbytes)\n");
    fprintf (stderr,"  -S n      - list n largest directories only\n");
    fprintf (stderr,"  path      - root directory (default: cwd)\n");
}

int
cmpDir (const void* a, const void* b) {
	Dir *ap = (Dir*)a;
	Dir *bp = (Dir*)b;
	return bp->size - ap->size;
}

int
main(int argc, char **argv) {
	char *cp;
	char buf[PATH_MAX];
              
   Progname = strsav(argv[0]);
	cp = strrchr(argv[0], '/');
	if (cp)
		Progname = strsav(cp+1);
	else
		Progname = strsav(argv[0]);

	D1(E("%s",Progname);)

   cmdline (argc,argv);

	if ((Startdir = getcwd(buf, sizeof(buf))) == NULL) {
		fprintf (stderr,"cannot get cwd: %s\n", ERR);
		exit(1);
    }
    if (!Root) Root = Startdir;
    Root = dodot(Root);
	 if (Print_tree)
		 printf ("\n\t\tDirectory Tree with root at %s\n",Root);
    prnt (Root,0);
		/*
		 * Print a summary of the NLargest directorys
		 */
	if (Largest) {
		int i;
		size_t tot = 0;
		D1(for (i=0;i<KLargest;i++)\
			T(" %10d %s\n", Largest[i].size, Largest[i].name);)
		if (KLargest > 1)
			qsort(Largest, KLargest, sizeof(Dir), cmpDir);
		NLargest = MIN(NLargest, KLargest);
		if (NLargest > 0) {
			printf ("  %d largest directories from %s:\n", NLargest, Root);
			for (i=0; i<NLargest; i++) {
				printf (" %10ld %s\n", Largest[i].size, Largest[i].name);
			}
			for (i=0; i<KLargest; i++)
				tot += Largest[i].size;
			if (tot > 500000000) {
				printf (" ----- Total: %ld bytes (%g Gb)\n",
					tot, (double)tot/1.0e+9);
			} else {
				printf (" ----- Total: %ld bytes (%g Mb)\n",
					tot, (double)tot/1.0e+6);
			}
		}
	}
	return 0;
}
