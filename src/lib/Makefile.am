# src/lib Makefile

SUFFIXES = .cpp .h

# CC = $(CXX)
# CFLAGS = $(CXXFLAGS) -fPIC

ACLOCAL_AMFLAGS = -I m4

if PROFILING
AM_CFLAGS = -pg
endif

if ANALYZE
CXXFLAGS += -fanalyzer
CXXFLAGS += -Wno-analyzer-malloc-leak
endif

libflaps_la_SOURCES = \
Ad.cpp \
atmos.cpp \
adeqn.cpp \
blas.cpp \
csignals.cpp \
Curve.cpp \
custom.cpp \
df.cpp \
einspline.cpp \
exim.cpp \
fma.cpp \
fio.cpp \
fptype.cpp \
interp.cpp \
interval.cpp \
lapack.cpp \
lexer.cpp \
lti.cpp \
matrix.cpp \
message.cpp \
nrc.cpp \
pac.cpp \
Pz.cpp \
Par.cpp \
pset.cpp \
qr.cpp \
Regex.cpp \
run_program.cpp \
settings.cpp \
Stack.cpp \
svd.cpp \
text.cpp \
trace.cpp \
tyme.cpp \
vastr.cpp

include_HEADERS = \
Ad.h \
atmos.h \
adeqn.h \
blas.h \
conv.h \
csignals.h \
Curve.h \
custom.h \
df.h \
einspline.h \
exim.h \
extract.h \
fma.h \
fio.h \
fptype.h \
interp.h \
interval.h \
lapack.h \
lexer.h \
lti.h \
main.h \
matrix.h \
message.h \
nrc.h \
ObjectCounter.h \
pac.h \
Pz.h \
Par.h \
pset.h \
qr.h \
Regex.h \
run_program.h \
settings.h \
Stack.h \
svd.h \
text.h \
trace.h \
tyme.h \
transpose.h \
vastr.h \
version.tmpl

lib_LTLIBRARIES = libflaps.la

pkgdata_DATA = custom.cpp config.h

LDADD = -L$(prefix)/lib -lflaps -lrt

# Test programs (compiled with -DMAIN)

Ad$(EXEEXT):	Ad.cpp Ad.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o Ad -g Ad.cpp $(LDADD) $(FLIBS) $(LIBS)

blas$(EXEEXT):	blas.cpp blas.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o blas -g blas.cpp $(LDADD) $(FLIBS) $(LIBS)

adeqn$(EXEEXT):	adeqn.cpp libflaps.la
		$(CXX) -I. $(AM_CPPFLAGS) $(CFLAGS) -DMAIN -o adeqn -g adeqn.cpp $(LDADD) $(FLIBS) $(LIBS)

Curve$(EXEEXT):	Curve.cpp Curve.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o Curve -g Curve.cpp $(LDADD) $(FLIBS) $(LIBS)

custom$(EXEEXT):	custom.cpp custom.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o custom -g custom.cpp $(LDADD) $(FLIBS) $(LIBS)

einspline$(EXEEXT):	einspline.cpp einspline.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o einspline -g einspline.cpp $(LDADD) $(FLIBS) $(LIBS)

interp$(EXEEXT):	interp.cpp interp.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o interp -g interp.cpp $(LDADD) $(FLIBS) $(LIBS)


fma$(EXEEXT):	fma.cpp fma.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o fma -g fma.cpp -I. $(LDADD) $(FLIBS) $(LIBS)

fptype$(EXEEXT):	fptype.cpp fptype.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o fptype -g fptype.cpp -I. $(LDADD) $(FLIBS) $(LIBS)

debug$(EXEEXT):	FORCE debug.cpp debug.h
		$(CXX) $(CFLAGS) -DMAIN -o debug -g debug.cpp -I. $(LDADD) -L$(FROOT)/lib $(FLIBS) $(LIBS)

denom$(EXEEXT):	denom.cpp denom.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o denom -g denom.cpp -I. $(LDADD) $(FLIBS) $(LIBS)

exim$(EXEEXT):	exim.cpp exim.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o exim -g exim.cpp -I. $(LDADD) $(FLIBS) $(LIBS)

fdir$(EXEEXT):	fdir.cpp fdir.h
		$(CXX) $(CFLAGS) -DMAIN -o fdir -g fdir.cpp -I. $(LDADD) $(FLIBS) $(LIBS)

fio$(EXEEXT):	fio.cpp fio.h
		$(CXX) $(CFLAGS) -DMAIN -o fio -g fio.cpp -I. $(LDADD) $(FLIBS) $(LIBS)

lapack$(EXEEXT):	lapack.cpp libflaps.la
		$(CXX) $(CFLAGS) -I. $(AM_CPPFLAGS) -DMAIN -o lapack -g lapack.cpp $(LDADD) $(FLIBS) $(LIBS)

lexer$(EXEEXT):	lexer.cpp lexer.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -I. -o lexer -g lexer.cpp $(LDADD) $(FLIBS) $(LIBS)

message$(EXEEXT):	message.cpp message.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o message -g message.cpp $(LDADD) $(FLIBS) $(LIBS)

matrix$(EXEEXT):	matrix.cpp libflaps.la
		$(CXX) $(CFLAGS) -I. $(AM_CPPFLAGS) -DMAIN -o matrix -g matrix.cpp $(LDADD) $(FLIBS) $(LIBS)

Par$(EXEEXT):	Par.cpp Par.h libflaps.la
		$(CXX) $(CFLAGS) -I. $(AM_CPPFLAGS) -DMAIN -o Par -g Par.cpp $(LDADD) $(FLIBS) $(LIBS)

pset$(EXEEXT):	pset.cpp pset.h
		$(CXX) $(CFLAGS) -I. $(AM_CPPFLAGS) -DMAIN -o pset -g pset.cpp $(LDADD) $(FLIBS)


nrc$(EXEEXT):	nrc.cpp nrc.h
		$(CXX) $(CFLAGS) -DMAIN -o nrc -g nrc.cpp $(LDADD) $(FLIBS) $(LIBS)

atmos$(EXEEXT):	atmos.cpp atmos.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o atmos -g atmos.cpp $(LDADD) $(FLIBS) $(LIBS)

text$(EXEEXT):	text.cpp text.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o text -g text.cpp $(LDADD) $(FLIBS)

Regex$(EXEEXT):	Regex.cpp Regex.h libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o Regex -g Regex.cpp $(LDADD) $(LIBS)

settings$(EXEEXT):	settings.cpp libflaps.la
		$(CXX) $(CFLAGS) -DMAIN -o settings -g settings.cpp $(LDADD) $(FLIBS) $(LIBS)

df$(EXEEXT):	FORCE df.cpp libflaps.la
		$(CXX) $(CFLAGS) -I. $(AM_CPPFLAGS) $(AM_CFLAGS) -DMAIN -o df -g df.cpp $(LDADD) -lfftw3 $(FLIBS) $(LIBS)

csignals$(EXEEXT):	FORCE csignals.cpp libflaps.la
		$(CXX) $(CFLAGS) -I. $(AM_CPPFLAGS) $(AM_CFLAGS) -DMAIN -o csignals -g csignals.cpp $(LDADD) $(FLIBS) $(LIBS)

tyme$(EXEEXT):	FORCE tyme.cpp libflaps.la
		$(CXX) $(CFLAGS) -I. $(AM_CPPFLAGS) $(AM_CFLAGS) -DMAIN -o tyme -g tyme.cpp $(LDADD) $(FLIBS) $(LIBS)


pac$(EXEEXT):		FORCE pac.cpp pac.h
		$(CXXCOMPILE) -DMAIN -o pac pac.cpp -L$(FROOT)/lib -lflaps $(LIBS) -lrt -lpthread
#		$(CXXCOMPILE) -fsanitize=address -DMAIN -o pac pac.cpp -lasan -L$(FROOT)/lib -lflaps $(LIBS)

qr$(EXEEXT):	FORCE qr.cpp
		$(COMPILE) -DMAIN -o qr qr.cpp lapack.o $(Libraries)

vastr$(EXEEXT):	FORCE vastr.cpp
		$(CXXCOMPILE) -DMAIN -o vastr vastr.cpp lapack.o $(Libraries)
 
FORCE:

CLEANFILES = \
apf$(EXEEXT) \
Ad$(EXEEXT) \
adeqn$(EXEEXT) \
denom$(EXEEXT) \
df$(EXEEXT) \
fdir$(EXEEXT) \
einspline$(EXEEXT) \
lapack$(EXEEXT) \
lexer$(EXEEXT) \
message$(EXEEXT) \
matrix$(EXEEXT) \
exim$(EXEEXT) \
mid$(EXEEXT) \
fio$(EXEEXT) \
fma$(EXEEXT) \
Par$(EXEEXT) \
pset$(EXEEXT) \
fptype$(EXEEXT) \
nrc$(EXEEXT) \
Regex$(EXEEXT) \
settings$(EXEEXT) \
stdatm$(EXEEXT) \
text$(EXEEXT) \
tyme$(EXEEXT) \
csignals$(EXEEXT) \
pac$(EXEEXT) \
qr$(EXEEXT) \
vastr$(EXEEXT)
