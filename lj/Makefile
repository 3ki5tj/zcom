CC = gcc -Wall -Wextra -O0 -g3 -pedantic
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0 -fno-alias -fargument-noalias
LMX =
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm -Wstrict-aliasing -fstrict-aliasing
#LMX = -lm

deps = lj.c lj.h ljmc.h ljmd.h ljrdf.h \
       def.h rv2.h rv3.h util.h md.h rng.h Makefile
bins = $(patsubst %.c,%, $(wildcard test*.c))
bins_d = $(patsubst %.c,%_d, $(wildcard test*.c))

debug: $(bins_d)

build: $(bins)

all: $(bins) $(bins_d)

$(bins) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

lj.h: lj.0.h mklj2d.py
	python mklj2d.py

ljmc.h: ljmc.0.h mklj2d.py
	python mklj2d.py

clean:
	rm -f *~ *.out MTSEED *.dat $(bins) $(bins_d) shakefail* rattlefail* lj.pos *.his
	rm -rf bin obj *.layout *.depend

.PHONY: clean all debug build

