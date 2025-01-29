
slides.pdf:	slides.tex
	pdflatex slides.tex
	hunspell -t -p slides.words -l slides.tex | sort | uniq -i
	@ echo "words listed by hunspell that are legit should be added to slides.words"

