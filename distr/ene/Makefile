CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
#CX = icc -Wall -O3 -g -vec_report0
CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops
LM=-lm 
src = zcom.h
dep = $(src) Makefile
bin = ene err merge

all: $(bin)

ene: ene.c $(dep)
	$(CX) -o $@ $< $(LM)

err: err.c $(dep)
	$(CX) -o $@ $< $(LM)

merge: merge.c $(dep)
	$(CX) -o $@ $< $(LM)

run: ene err merge
	python run.py

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos
	rm -rf bin obj *.depend *.layout fe*.dat ds.dat dsb.dat err.dat
	rm -f data/ene.ps data/eneref.ps data/ds.dat data/dsb.dat \
	  data/*_aj.dat data/*_mfav.dat data/*_mfii.dat data/*_c_ix.dat \
	  data/dsb_c.dat data/ds_c.dat data/dsb_m.dat data/ds_m.dat \
	  data/*merge*.dat data/err.dat data[lh]/dsb.dat data[lh]/ds.dat \
	  data[lh]/dsb_c_ii.dat data[lh]/ds_c_ii.dat

.PHONY: clean all

