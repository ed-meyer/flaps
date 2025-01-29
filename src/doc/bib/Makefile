
articles_to_fetch.pdf:	FORCE articles_to_fetch.tex
	make clean
	pdflatex articles_to_fetch.tex
	bibtex articles_to_fetch
	pdflatex articles_to_fetch.tex
	pdflatex articles_to_fetch.tex

articles_to_fetch.tex:		articles_to_fetch.bib
		echo '\documentclass{article}' >articles_to_fetch.tex
		echo '\usepackage{url}' >> articles_to_fetch.tex
		echo '\\begin{document}' >> articles_to_fetch.tex
		sed -n 's/^@.*{\([^,]*\),/\\cite{\1}/p' articles_to_fetch.bib >>articles_to_fetch.tex
		echo '\\bibliography{flaps}' >> articles_to_fetch.tex
		echo '\\bibliographystyle{plain}' >> articles_to_fetch.tex
		echo '\\end{document}' >> articles_to_fetch.tex

bib.pdf:	FORCE bib.tex
	make clean
	pdflatex bib.tex
	bibtex bib
	pdflatex bib.tex
	pdflatex bib.tex

bib.tex:		flaps.bib
		echo '\documentclass{article}' >bib.tex
		echo '\usepackage{url}' >> bib.tex
		echo '\\begin{document}' >> bib.tex
		sed -n 's/^@.*{\([^,]*\),/\\cite{\1}/p' flaps.bib >>bib.tex
		echo '\\bibliography{flaps}' >> bib.tex
		echo '\\bibliographystyle{plain}' >> bib.tex
		echo '\\end{document}' >> bib.tex

# Make check.aiaa.bbl and check.mpdi.bbl for inclusion
# in publications

aiaa.bbl:	FORCE flaps.bib
		../bin/checkbib ./flaps.bib aiaa
		cp checkdir/check.aiaa.bbl aiaa.bbl

mdpi.bbl:	FORCE flaps.bib
		../bin/checkbib ./flaps.bib mdpi
		cp checkdir/check.mdpi.bbl mdpi.bbl

FORCE:

clean:
	rm -f \
	bib.aux \
	bib.bbl \
	bib.blg \
	bib.log \
	bib.out \
	bib.pdf
