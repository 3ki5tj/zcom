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
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos
	rm -rf bin obj *.depend *.layout
	rm -rf data/dsaj*.dat data/dsbaj*.dat data/ds.dat data/dsb.dat \
               data/rdf.ps data/rdfref.ps

.PHONY: clean all

