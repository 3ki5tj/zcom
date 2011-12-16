CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = vol

all: $(bin)

vol: vol.c $(dep)
	$(CX) -o $@ $< 

run: vol
	python run.py

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos
	rm -rf bin obj *.depend *.layout fe.dat

.PHONY: clean all

