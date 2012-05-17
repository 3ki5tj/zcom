prj=distr

$(prj).pdf: $(prj).tex  $(prj)text.pdf Makefile \
  phi.pdf ene.pdf vol.pdf rdf.pdf bb.png ref.pdf phiaj.pdf
	pdflatex $(prj).tex

# the integral sign in phi.pdf needs fonts, the option -dPDFSETTINGS=/print is necessary

phi.pdf: phi.ps
	ps2pdf -dPDFSETTINGS=/printer -dEmbedAllFonts=true $<

ene.pdf: ene.ps
	ps2pdf $<

vol.pdf: vol.ps
	ps2pdf $<

rdf.pdf: rdf.ps
	ps2pdf $<

phiaj.pdf: phiaj.ps
	ps2pdf $<

ref.pdf: ref.ps
	ps2pdf $<

ref.ps: ref.gp
	gpload $<

.PHONY: clean

clean:
	rm -rf $(prj).log *.aux *~

