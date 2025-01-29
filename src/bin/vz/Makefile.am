
# src/bin/vz Makefile

FORCE:

vz_SOURCES = \
dialog.cpp \
hard_devices.cpp \
main.cpp \
params.cpp \
plotcurve.cpp \
ps.cpp \
specs.cpp \
xfig.cpp \
xgX.cpp \
xtb.cpp \
dot.11 \
mark1.11 \
mark2.11 \
mark3.11 \
mark4.11 \
mark5.11 \
mark6.11 \
mark7.11 \
mark8.11 \
tagged.11


Amvz = \
amdata.cpp \
amvz.cpp \
swf.cpp \
trackball.cpp \
prefdialog.cpp

amvz_SOURCES = $(Amvz) plotcurve.cpp specs.cpp params.cpp

amvz_CXXFLAGS = $(CXXFLAGS) -DSTANDALONE

vz_SOURCES += $(Amvz)

# gl2ps.h

noinst_HEADERS = \
amdata.h \
amvz.h \
amvz.xpm \
copyright.h \
idraw.h \
xfig.h \
swf.h \
hard_devices.h \
params.h \
pkparam.h \
plotcurve.h \
plotter.h \
prefdialog.h \
ps.h \
specs.h \
trackball.h \
xgout.h \
xtb.h

bin_PROGRAMS = vz

CPPFLAGS += `wx-config --cxxflags`

wxlibs = `wx-config --libs gl,core,base`
# LDADD = $(Libraries) $(FLIBS) $(wxlibs) $(X_PRE_LIBS) $(X_LIBS) $(X_EXTRA_LIBS) -lX11 -lGL -lGLU -lgl2ps
# AM_LDFLAGS = $(FLIBS) $(wxlibs) $(X_PRE_LIBS) $(X_LIBS) $(X_EXTRA_LIBS) -lX11 -lGL -lGLU -lgl2ps
# AM_LDFLAGS = `wx-config --libs gl,core,base` -lGL -lGLU -lgl2ps

# LIBS = $(Libraries) $(FLIBS) $(wxlibs) $(X_PRE_LIBS) $(X_LIBS) $(X_EXTRA_LIBS) -lX11 -lGL -lGLU -lgl2ps

LIBS += -L$(FROOT)/lib -lflaps -lrt
LIBS += `wx-config --libs gl,core,base`
LIBS += `pkg-config --libs opengl`
LIBS += `pkg-config --libs gl glu`
LIBS += -lX11

bin_PROGRAMS += amvz

vz_DEPENDENCIES = FORCE
