CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) GNUmakefile
bin = ljmc ljmd sqmc is isent

all: $(bin)

ljmc: ljmc.c $(dep)
	$(CX) -o $@ $< 

ljmd: ljmd.c $(dep)
	$(CX) -o $@ $< 

sqmc: sqmc.c $(dep)
	$(CX) -o $@ $< 

is: is.c $(dep)
	$(CX) -o $@ $< 

isent: isent.c $(dep)
	$(CX) -o $@ $< 

zip: zcom.h
	zip rptprog *.c *.ma zcom.h islogdos*.txt Makefile README

zcom.h: zcom.full.h GNUmakefile
	./zcompick.py -i zcom.full.h -o zcom.h -k "def util ss rng opt argopt rpt av hist rv2 rv3 md lj ising2" 

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos rptprog.zip zcom.h
	rm -rf bin obj *.depend *.layout *.dat

.PHONY: clean all

