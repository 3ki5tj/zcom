# generic make file for new modules

CC = gcc -O0 -g3 -Wall -Wextra -pedantic
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -march=native -vec_report0
LMX =
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops
#LMX =-lm

deps = cfg.h ss.h opt.h def.h util.h Makefile
bins = $(patsubst %.c,%, $(wildcard test*.c))
bins_d = $(patsubst %.c,%_d, $(wildcard test*.c))

debug: $(bins_d)

build: $(bins)

all: $(bins) $(bins_d)

$(bins_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bins) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

subdirs =

clean:
	$(RM) -f *~ *.out MTSEED *.dat $(bins) $(bins_d)
	-for d in $(subdirs); do (cd $$d; $(MAKE) clean ); done
	rstrip.py -lRv

.PHONY: clean all debug build

