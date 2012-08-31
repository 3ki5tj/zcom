CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = ljmd ljmc sqmc sqmc1 is isent

all: $(bin)

ljmc: ljmc.c $(dep)
	$(CX) -o $@ $< 

ljmd: ljmd.c $(dep)
	$(CX) -o $@ $< 

sqmc: sqmc.c $(dep)
	$(CX) -o $@ $< 

sqmc1: sqmc1.c $(dep)
	$(CX) -o $@ $< 

is: is.c $(dep)
	$(CX) -o $@ $< 

isent: isent.c $(dep)
	$(CX) -o $@ $< 

run: $(bin)
	python sizemcl.py
	python rcmcl.py
	python others.py

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos
	rm -rf bin obj *.depend *.layout *.dat 

.PHONY: clean all

