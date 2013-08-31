# generic make file for new modules

CC = gcc -O0 -g3 -Wall -Wextra -pedantic
LM = -lm
# optimized compiling
# -wd161: disable unknown #pragma
CX = icc -Wall -Wremarks -Wcheck -O3 -march=native -vec_report0 -wd161
LMX =
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops
#LMX =-lm

deps = util.h svd.h Makefile
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
	rstrip.py -lRv

.PHONY: clean all debug build

