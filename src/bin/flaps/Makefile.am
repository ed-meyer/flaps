
# Makefile for flaps, the Flaps executive

flaps_SOURCES = \
catalog.cpp \
cmds.cpp \
export.cpp \
functions.cpp \
gyro.cpp \
import.cpp \
lexcmds.cpp \
load.cpp \
lti.cpp \
main.cpp \
motd.cpp \
save_restore.cpp \
specs.cpp \
symmetrize.cpp

noinst_HEADERS = \
cmds.h \
functions.h \
load.h \
specs.h

LDADD = -L$(FROOT)/lib -lflaps -lrt

bin_PROGRAMS = flaps

FORCE:
flaps_DEPENDENCIES = FORCE

load$(EXEEXT):	FORCE load.cpp
		$(COMPILE) -DMAIN -o load load.cpp $(LDADD)

save_restore$(EXEEXT):	FORCE save_restore.cpp
		$(COMPILE) -DMAIN -o save_restore save_restore.cpp $(Libax) $(Liblapack) $(Libblas) $(LIBS)
