# src/bin/flut Makefile

FORCE:

# asan stuff:
# CXXFLAGS += -fsanitize=address
# LIBS += -lasan

LIBS += -L$(FROOT)/lib -lflaps -lrt
LIBS += -lpthread


flut_SOURCES = \
bifurcation.cpp \
concurrent.cpp \
flutcurve.cpp \
startpts.cpp \
main.cpp \
process.cpp \
region.cpp \
specs.cpp \
target.cpp

noinst_HEADERS = \
concurrent.h \
Exception.h \
flutcurve.h \
startpts.h \
process.h \
region.h \
specs.h \
target.h

bin_PROGRAMS = flut

if PROFILING
AM_CFLAGS = -pg
endif

main.o:		version.cpp

version.cpp:		FORCE
		sed "s,REPLACE WITH DATE,`date`," ../../lib/version.tmpl >version.cpp

CLEANFILES = \
version.cpp
