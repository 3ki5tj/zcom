CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = abpro.c
dep = $(src) abpro.h rv3.h rv2.h util.c util.h Makefile
bin = testef testcn

all: $(bin)

testef: testef.c $(dep)
	$(CC) -o $@ $< 

testcn: testcn.c $(dep)
	$(CC) -o $@ $< 

clean:
	rm -f *~ *.out MTSEED *.dat $(bin)

.PHONY: clean all

