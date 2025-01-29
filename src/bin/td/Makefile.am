
FORCE:

CC = $(CXX)
CFLAGS = $(CXXFLAGS)

if PROFILING
AM_CFLAGS = -pg
endif

AM_CPPFLAGS = -I$(FROOT)/include

LDADD = $(Libraries) -lpthread

td_SOURCES = \
main.cpp \
hopf.cpp \
region.cpp \
specs.cpp

noinst_HEADERS = \
hopf.h \
region.h \
specs.h

bin_PROGRAMS = td

main.o:		version.cpp

version.cpp:		FORCE
		sed "s,REPLACE WITH DATE,`date`," ../../lib/version.tmpl >version.cpp


# test programs built with -DMAIN

CLEANFILES = \
version.cpp
