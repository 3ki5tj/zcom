deps = Makefile

all: cvg.pdf

cvg.pdf: cvg.tex $(deps)
	pdflatex $<

.PHONY: clean

clean:
	rm -f *.log *~ *.aux *.mpx *.dvi *.ps *.pdf
