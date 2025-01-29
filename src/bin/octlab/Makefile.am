
FORCE:

octlab_SOURCES = \
octlab.cpp \
specs.cpp

noinst_HEADERS = \
specs.h

bin_PROGRAMS = octlab

LDADD=-L$(FROOT)/lib -lflaps

main.o:		version.cpp

version.cpp:		FORCE
		sed "s,REPLACE WITH DATE,`date`," ../../lib/version.tmpl >version.cpp

CLEANFILES = \
version.cpp
