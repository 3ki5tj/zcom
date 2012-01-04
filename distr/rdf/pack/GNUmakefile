CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = rdf

all: $(bin)

rdf: rdf.c $(dep)
	$(CX) -o $@ $< 

run: rdf
	python run.py

zcom.h::
	./zcompick.py -i ../zcom.h -o zcom.h -k "def util ss rng opt cfg distr av rc rv2 rv3 md specfunc lj" 

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos
	rm -rf bin obj *.depend *.layout ds.dat dsb.dat

.PHONY: clean all

