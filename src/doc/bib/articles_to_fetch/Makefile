
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


FORCE:

clean:
	rm -f \
	*.aux \
	*.bbl \
	*.blg \
	*.log \
	*.out \
	*.pdf \
	*.tex
