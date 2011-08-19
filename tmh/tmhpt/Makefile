CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = tmhpt

all: $(bin)

tmhpt: tmhpt.c $(dep)
	$(CX) -o $@ $< 

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) shakefail* rattlefail* tmh*.e tmh*.t tmh*.ehis TRACE

.PHONY: clean all

