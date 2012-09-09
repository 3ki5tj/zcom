# CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CC = icc -Wall -O3 -g -vec_report0
src = zcom.h
dep = $(src) GNUmakefile
bin = $(patsubst %.c,%, $(wildcard *.c))

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CC) -o $@ $<

zip: zcom.h
	zip dsprog *.c *.cfg zcom.h Makefile README

zcom.h: zcom.full.h GNUmakefile
	./zcompick.py -i zcom.full.h -o zcom.h -k "def util ss rng opt cfg distr av rc rv2 rv3 md specfunc hist lj" 

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos dsprog.zip zcom.h
	rm -rf bin obj *.depend *.layout ds.dat dsb.dat fe.dat

.PHONY: clean all

