.PHONY: clean

dih.pdf: dih.tex dihfig.ps
	pdflatex $<
	pdflatex $<

dihfig.ps: dihfig.mp
	mpost $<
	mv dihfig.1 $@

dih.pdf: dih.ps
	ps2pdf dih.ps

dih.ps: dih.dvi
	dvips dih.dvi

clean:
	rm -f *.log *~ *.aux *.mpx *.dvi *.ps
