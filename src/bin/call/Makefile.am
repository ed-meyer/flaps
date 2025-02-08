
FORCE:

call_SOURCES = \
main.cpp

bin_PROGRAMS = call

LDADD=-L$(FROOT)/lib -lflaps -lrt -llapack -lblas

main.o:		version.cpp

version.cpp:		FORCE
		sed "s,REPLACE WITH DATE,`date`," ../../lib/version.tmpl >version.cpp

CLEANFILES = \
version.cpp
