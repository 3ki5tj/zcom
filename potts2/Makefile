CC = gcc -Wall -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = potts2.c
dep = $(src) potts2.h rng.c rng.h util.h Makefile
bin = testem testmc

all: $(bin)

testem: testem.c $(dep)
	$(CC) -o $@ $< 

testmc: testmc.c $(dep)
	$(CX) -o $@ $<

clean:
	rm -f *~ *.out MTSEED *.dat $(bin)

.PHONY: clean all

