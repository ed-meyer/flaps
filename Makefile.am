
SUBDIRS = @subdirs@

FORCE:

EXTRA_DIST = \
COPYING \
doc/manual.pdf \
header \
Profile \
PACKAGES \
README \
Version

CLEANFILES = \
include/* \
bin/* \
lib/*

backup:	FORCE
	make dist
	if [ -d $(HOME)/Drive/backup ]; then \
		mv $(distdir).tar.gz $(HOME)/Drive/backup/$(distdir).`date +%s`.tar.gz; \
	else \
		mv $(distdir).tar.gz $(distdir).`date +%s`.tar.gz; \
	fi
	
