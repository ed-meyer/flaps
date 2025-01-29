
# Makefile for fem, whirl-flutter beam modelling

fem_SOURCES = \
dlm.cpp \
main.cpp \
modal.cpp \
ss.cpp \
vzdata.cpp

include_HEADERS = fem.h

bin_PROGRAMS = fem

LDADD = -L$(FROOT)/lib -lflaps -lrt

FORCE:
fem_DEPENDENCIES = FORCE

main.o:		version.cpp

version.cpp:		FORCE
		sed "s,REPLACE WITH DATE,`date`," ../../lib/version.tmpl >version.cpp

CLEANFILES = \
version.cpp
