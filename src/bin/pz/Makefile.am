
FORCE:

pz_SOURCES = \
main.cpp \
specs.cpp

noinst_HEADERS = \
specs.h

bin_PROGRAMS = pz

if PROFILING
AM_CFLAGS = -pg
endif

LDADD=-L$(FROOT)/lib -lflaps -lrt

main.o:		version.cpp

version.cpp:		FORCE
		sed "s,REPLACE WITH DATE,`date`," ../../lib/version.tmpl >version.cpp

# build a userdf
userdf:	userdf.cpp FORCE
	$(CXX) $(CFLAGS) -I. -I$(FROOT)/share/pz $(AM_CPPFLAGS) $(AM_CFLAGS) -o userdf -g userdf.cpp $(LDADD) $(FLIBS) $(LIBS)

CLEANFILES = version.cpp

