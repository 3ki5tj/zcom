prj=distr

$(prj).pdf: $(prj).tex  $(prj)text.pdf \
  phi.pdf ene.pdf vol.pdf rdf.pdf bb.png ref.pdf
	pdflatex $(prj).tex

phi.pdf: phi.ps
	ps2pdf $<

ene.pdf: ene.ps
	ps2pdf $<

vol.pdf: vol.ps
	ps2pdf $<

rdf.pdf: rdf.ps
	ps2pdf $<

ref.pdf: ref.ps
	ps2pdf $<

ref.ps: ref.gp
	gpload $<

.PHONY: clean

clean:
	rm -rf $(prj).log *.aux *~

