CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = tmhab tmhab_d gtmhab

all: $(bin)

tmhab: tmhab.c $(dep)
	$(CX) -o $@ $< 

tmhab_d: tmhab.c $(dep)
	$(CC) -o $@ $< 

gtmhab: gtmhab.c $(dep) gtmh.h
	$(CX) -o $@ $< -lglut

zcom.h::
	cd ../../abpro && python ./mk2d.py && cd tmhab
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) shakefail* rattlefail* tmh*.e tmh*.t tmh*.ehis tmh*.tr ab*.pos
	rm -rf bin obj
	cd test && ./clean

.PHONY: clean all

