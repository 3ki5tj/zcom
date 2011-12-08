CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = ljmatch

all: $(bin)

ljmatch: ljmatch.c $(dep)
	$(CX) -o $@ $< 

zcom.h::
	cd ../../lj && python ./mk2d.py && cd ljmatch
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) shakefail* rattlefail* tmh*.e tmh*.t tmh*.ehis tmh*.tr

.PHONY: clean all

