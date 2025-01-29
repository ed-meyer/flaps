
FORCE:

parameters_SOURCES = \
main.cpp

bin_PROGRAMS = parameters

if PROFILING
AM_CFLAGS = -pg
endif

LDADD=-L$(FROOT)/lib -lflaps -lrt

main.o:		version.cpp

version.cpp:		FORCE
		sed "s,REPLACE WITH DATE,`date`," ../../lib/version.tmpl >version.cpp

CLEANFILES = \
version.cpp
