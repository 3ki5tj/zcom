CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) GNUmakefile
bin = ene vol rdf

all: $(bin)

ene: ene.c $(dep)
	$(CX) -o $@ $< 

vol: vol.c $(dep)
	$(CX) -o $@ $< 

rdf: rdf.c $(dep)
	$(CX) -o $@ $< 

zip: zcom.h
	zip dsprog ene.c ene.cfg vol.c vol.cfg rdf.c rdf.cfg zcom.h Makefile README

zcom.h: zcom.full.h GNUmakefile
	./zcompick.py -i zcom.full.h -o zcom.h -k "def util ss rng opt cfg distr av rc rv2 rv3 md specfunc lj" 

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos dsprog.zip zcom.h
	rm -rf bin obj *.depend *.layout ds.dat dsb.dat fe.dat

.PHONY: clean all

