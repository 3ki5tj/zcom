CC = gcc -Wall -Wextra -Wno-unknown-pragmas -O0 -g3 -pedantic -lm
# optimized compiling
CX = icc -Wall -O3 -vec_report0 -wd161
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
CM = $(CX) -openmp -D_OPENMP -openmp-report0 
#CM = $(CX) -openmp-profile -D_OPENMP -openmp-report0
#CM = $(CC) -fopenmp

deps = abpro.c abpro.h rv3.h rv2.h util.h md.c md.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c)) \
      $(patsubst %.c,%_d, $(wildcard test*.c)) \
      $(patsubst %.c,%_m, $(wildcard test*.c))

all: $(bin)

$(bin): $(deps)

%_d: %.c
	$(CC) -o $@ $< -lm 

%: %.c
	$(CX) -o $@ $< 

%_m: %.c
	$(CM) -o $@ $<

abpro.c: abpro.0.c mk2d.py
	python mk2d.py $@

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) shakefail* rattlefail* ab.pos guide.gvs
	$(MAKE) clean -C tmhab
	$(MAKE) clean -C abemin

.PHONY: clean all

