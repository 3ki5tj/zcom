CC = gcc -Wall -Wextra -Wno-unknown-pragmas -O0 -g3 -pedantic
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -march=native -vec_report0 -wd161
LMX =
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops
#LMX =-lm
CM = $(CX) -openmp -D_OPENMP -openmp-report0
LMM =
#CM = $(CX) -openmp-profile -D_OPENMP -openmp-report0
#CM = $(CC) -fopenmp

deps = abpro.c abpro.h rv3.h rv2.h util.h md.c md.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c))
bin_d = $(patsubst %.c,%_d, $(wildcard test*.c))
bin_m = $(patsubst %.c,%_m, $(wildcard test*.c))

all: $(bin) $(bin_d) $(bin_m)

$(bin_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bin) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bin_m) : %_m : %.c $(deps)
	$(CM) -o $@ $< $(LMM)

abpro.c: abpro.0.c mk2d.py
	python mk2d.py $@

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) $(bin_d) $(bin_m) \
	  shakefail* rattlefail* ab.pos guide.gvs
	$(MAKE) clean -C tmhab
	$(MAKE) clean -C abemin

.PHONY: clean all

