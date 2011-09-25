CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = cago.c
dep = $(src) cago.h def.h rv3.h md.c md.h pdb.h pdb.c dihc.h dihc.c util.c util.h rotfit.c rotfit.h Makefile
bin = testmd testmd1 testcvg

all: $(bin)

testmd: testmd.c $(dep)
	$(CX) -o $@ $< 

testmd1: testmd1.c $(dep)
	$(CX) -o $@ $< 

testcvg: testcvg.c $(dep)
	$(CX) -o $@ $<

clean:
	rm -f *~ *.out MTSEED *.dat *.pos $(bin) final.pdb ref.pdb

.PHONY: clean all

