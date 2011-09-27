.PHONY: clean

dih.dvi: dih.tex dihfig.ps
	latex $<
	latex $<

dihfig.ps: dihfig.mp
	mpost $<
	mv dihfig.1 $@

clean:
	rm -f *.log *~ *.aux *.mpx *.dvi *.ps
