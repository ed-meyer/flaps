
# src/bin/viz Makefile

FORCE:

# pre-processor flags
# AM_CPPFLAGS = -I$(FROOT)/include
CPPFLAGS += `wx-config --cxxflags`
# uncomment this for threads:
# CPPFLAGS += -DTHREADS=1

LIBS += -L$(FROOT)/lib -lflaps -lrt
# LIBS = $(Libraries)
LIBS += `wx-config --libs gl,core,base`
LIBS += `pkg-config --libs opengl`
LIBS += `pkg-config --libs gl glu`
LIBS += -lX11
LIBS += -lrt
LIBS += -lpng
# LIBS += `pkg-config --libs ftgl`

viz_SOURCES = \
viz.cpp \
amviz.cpp \
mathplot.cpp \
plotcurve.cpp \
specs.cpp \
trackball.cpp \
util.cpp

amviz_SOURCES = \
amain.cpp \
amviz.cpp \
plotcurve.cpp \
specs.cpp \
trackball.cpp \
util.cpp

colors_SOURCES = util.cpp


noinst_HEADERS = \
amviz.h \
amviz.xpm \
mathplot.h \
plotcurve.h \
specs.h \
trackball.h \
viz.h \
util.h

bin_PROGRAMS = viz amviz colors

viz_DEPENDENCIES = FORCE

colors:		util.cpp util.h
		$(CXX) $(CPPFLAGS) $(CFLAGS) -DMAIN -o colors -g util.cpp $(LIBS)
