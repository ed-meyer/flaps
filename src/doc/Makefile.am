
# src/doc Makefile
# NOTES:
#   to extract pages from manual.pdf and create a new pdf file:
#      texexec --pdfselect --selection=78,79,80 manual.pdf

SUFFIXES = .fig .tex .eps .ps .pdf

# subdirs is defined in configure.ac
SUBDIRS = @subdirs@

docdir=${prefix}/doc

install-pdf:	manual.pdf
	mkdir -p ${prefix}/doc
	cp -u manual.pdf ${prefix}/doc

install:
	make install-pdf

# We must ignore errors from pdflatex because it always
# returns with an exit status of 1 - I believe this is
# due to the fact that the manual is using some extensions
# evidenced by the fact that pdflatex prints "This is pdfTeXk"
# which according to the man page is because we are using extensions.

SRC = \
flaps.bib \
manual.tex

manual.pdf:	manual.tex
#	cd figures && make
	pdflatex manual.tex
	bibtex manual
	makeindex manual.idx
	pdflatex manual.tex
	pdflatex manual.tex
	detex -l manual.tex | wc
	detex -l manual.tex | diction -s | fmt >manual.diction
	hunspell -t -p manual.words -l manual.tex | sort | uniq -i
	@ echo "words listed by hunspell that are legit should be added to manual.words"

FORCE:

noinst_DATA = \
flaps.bib \
manual.pdf

refs:		FORCE
		sed -n 's/^.*label{eqn:\([^}]*\)}.*/\1/p' manual.tex >eqn.ref
		sed -n 's/^.*label{sect:\([^}]*\)}.*/\1/p' manual.tex >sect.ref
		sed -n 's/^.*label{chap:\([^}]*\)}.*/\1/p' manual.tex >chap.ref
		sed -n 's/^.*label{table:\([^}]*\)}.*/\1/p' manual.tex >table.ref
		sed -n 's/^.*label{ex:\([^}]*\)}.*/\1/p' manual.tex >ex.ref
		sed -n 's/^.*label{opt:\([^}]*\)}.*/\1/p' manual.tex >opt.ref

# Update the manual in MyDrive
drive:		manual.pdf
	cp -u manual.pdf $(HOME)/MyDrive

# tarfile:		FORCE
# 	tar zcf manual.tgz manual.tex fig *.eps Makefile

EXTRA_DIST = $(SRC) manual.pdf

CLEANFILES = \
manual.out \
manual.lof \
manual.lot \
manual.diction \
manual.ps \
manual.pdf \
manual.aux \
manual.bbl \
manual.blg \
manual.dvi \
manual.idx \
manual.ilg \
manual.ind \
manual.log \
manual.ps \
manual.bib \
manual.toc
