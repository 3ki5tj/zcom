CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm
src = zcom.h
dep = $(src) Makefile
bin = tmhcago

all: $(bin)

tmhcago: tmhcago.c $(dep)
	$(CX) -o $@ $<

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) tmh*.e tmh*.t tmh*.ehis TRACE *.pos guess.range

.PHONY: clean all

