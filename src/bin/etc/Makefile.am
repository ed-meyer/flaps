
# src/bin/etc

SUFFIXES = .cpp .h

# XXX convert linpackd to c++
# linpackd
# iwyu

bin_PROGRAMS = \
catalog \
cores \
logger \
dtree \
exitstatus \
resource_limits \
sizeof \
wildcard \
zap

cores_SOURCES = cores.cpp
logger_SOURCES = logger.cpp
dtree_SOURCES = dtree.cpp
exitstatus_SOURCES = exitstatus.cpp
catalog_SOURCES = catalog.cpp
iwyu_SOURCES = iwyu.cpp
resource_limits_SOURCES = resource_limits.cpp
sizeof_SOURCES = sizeof.cpp
wildcard_SOURCES = wildcard.cpp
zap_SOURCES = zap.cpp
# linpackd_SOURCES = linpackd.f

LDADD=-L$(FROOT)/lib -lflaps -lrt

FORCE:

# linpackd$(EXEEXT):		linpackd.f
# 		gfortran -o linpackd -O2 linpackd.f

