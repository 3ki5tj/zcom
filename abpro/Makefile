CC = gcc -Wall -Wextra -Wno-unknown-pragmas -O0 -g3 -pedantic
LM = -lm
# optimized compiling
# -wd161 to disable unrecognized #pragma
CX = icc -Wall -O3 -march=native -vec_report0
LMX =
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops
#LMX =-lm
CM = $(CX) -openmp -D_OPENMP -openmp-report0
LMM =
#CM = $(CX) -openmp-profile -D_OPENMP -openmp-report0
#CM = $(CC) -fopenmp

deps = abpro.h rv2.h rv3.h rm3.h util.h md.h Makefile
bins = $(patsubst %.c,%, $(wildcard test*.c))
bins_d = $(patsubst %.c,%_d, $(wildcard test*.c))
bins_m = $(patsubst %.c,%_m, $(wildcard test*.c))

debug: $(bins_d)

build: $(bins)

omp: $(bins_m)

all: $(bins) $(bins_d) $(bins_m)

$(bins_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bins) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bins_m) : %_m : %.c $(deps)
	$(CM) -o $@ $< $(LMM)

abpro.h: abpro.0.h mkab2d.py
	python mkab2d.py $@

clean:
	rm -f *~ *.out MTSEED *.dat $(bins) $(bins_d) $(bins_m) \
	  shakefail* rattlefail* ab.pos guide.gvs
	$(MAKE) clean -C abemin

.PHONY: clean all debug build omp

