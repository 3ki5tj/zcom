CC = gcc -Wall -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -lm
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = is2.c
dep = $(src) is2.h rng.c rng.h Makefile
bin = testmc testem testexact

all: $(bin)

testmc: testmc.c $(dep)
	$(CX) -o $@ $< $(src) 

testem: testem.c $(dep)
	$(CC) -o $@ $< $(src) 

testexact: testexact.c $(dep)
	$(CC) -o $@ $< $(src) 

clean:
	rm -f *~ *.out MTSEED *.dat $(bin)

.PHONY: clean all

