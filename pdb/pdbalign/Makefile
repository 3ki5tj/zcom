CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = pdbalign

all: $(bin)

pdbalign: pdbalign.c $(dep)
	$(CX) -o $@ $< 

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

install: $(bin)
	cp pdbalign ~/app/bin/

clean:
	rm -f *~ *.out MTSEED $(bin) fit*.pdb

.PHONY: clean all

