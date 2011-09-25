CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = tmhab

all: $(bin)

tmhab: tmhab.c $(dep)
	$(CC) -o $@ $< 

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) shakefail* rattlefail* tmh*.e tmh*.t tmh*.ehis TRACE ab.pos

.PHONY: clean all

