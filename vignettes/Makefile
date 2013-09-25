all:	pdf

pdf:	SimBindProfiles.tex
	pdflatex SimBindProfiles
	bibtex SimBindProfiles
	pdflatex SimBindProfiles
	pdflatex SimBindProfiles
	make clean

clean:
	rm -f *.out *.bbl *.aux *.log *.blg *.brf *.toc
	rm -f *.nfs*
	rm -f *.bed
	rm -f SimBindProfiles-*.pdf
	rm -f *.png

