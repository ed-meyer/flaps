
# this package is from csm.ornl.gov/~kohl/MatView
# src/bin/matview Makefile

FORCE:

matview_SOURCES = \
commands.c \
file.c \
matview.c \
matview_glob.c \
util.c \
view.c

# these are needed at runtime: install in froot/share/matview
pkgdata_DATA = \
compass.xbm \
matview.tcl \
procs.tcl \
util.tcl

noinst_HEADERS = \
const.h \
matview_glob.h \
matview.h

CPPFLAGS += -I/usr/include/tcl

Libtcl = -ltcl
Libtk = -ltk

X11_LIBS = \
$(X_PRE_LIBS) \
$(X_LIBS) \
$(X_EXTRA_LIBS) \
-lX11

LDADD = -L$(FROOT)/lib $(Libtcl) $(Libtk) -lm

bin_PROGRAMS = matview

matview_DEPENDENCIES = FORCE

EXTRA_DIST = $(noinst_HEADERS) $(pkgdata_DATA)
