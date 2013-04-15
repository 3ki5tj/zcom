deps = Makefile

all: dih.pdf

dih.pdf: dih.tex dihfig.pdf $(deps)
	pdflatex $<

dihfig.ps: dihfig.mp $(deps)
	mpost $<
	mv dihfig.1 $@

dihfig.pdf: dihfig.ps $(deps)
	epstopdf dihfig.ps

.PHONY: clean

clean:
	rm -f *.log *~ *.aux *.mpx *.dvi *.ps *.pdf
