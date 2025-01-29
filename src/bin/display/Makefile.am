
FORCE:

display_SOURCES = \
main.cpp \
specs.cpp

noinst_HEADERS = \
specs.h

bin_PROGRAMS = display

if PROFILING
AM_CFLAGS = -pg
endif

LDADD = -L$(FROOT)/lib -lflaps -lrt

main.o:		version.c

version.c:		FORCE
		sed "s,REPLACE WITH DATE,`date`," ../../lib/version.tmpl >version.c

CLEANFILES = version.c
