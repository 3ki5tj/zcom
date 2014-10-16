CC = gcc -Wall -Wextra -O0 -g3 -pedantic
LM = -lm
# optimized compiling
CX = icc -Wall -Wremarks -Wcheck -O3 -g -vec_report0 -fno-alias -fargument-noalias
LMX =
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm -Wstrict-aliasing -fstrict-aliasing
#LMX = -lm

deps = lj.c lj.h ljeos.h ljmc.h ljmd.h ljrdf.h \
       rv2.h rv3.h util.h md.h rng.h Makefile
bins = $(patsubst %.c,%, $(wildcard test*.c))
bins_d = $(patsubst %.c,%_d, $(wildcard test*.c))
autosrc = lj.h ljmc.h

debug: $(bins_d)

build: $(bins)

all: $(bins) $(bins_d) ljeosKN

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

ljeosKN: ljeosKN.f KN.f
	gfortran $^ -o $@

genhdr : $(autosrc)

$(autosrc) : %.h : %.0.h add2d.py
	python add2d.py

clean:
	rm -f *~ *.out MTSEED *.dat \
	  $(bins) $(bins_d) ljeosKN \
	  lj.pos *.his
	rm -rf bin obj *.layout *.depend

.PHONY: clean all debug build

