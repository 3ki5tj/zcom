src=rng.c
hdr=rng.h
dep=$(src) $(hdr) Makefile
CC = gcc -Wall -Wextra -Wno-unknown-pragmas -O0 -g3 -pedantic -lm
# optimized compiling
CX = icc -Wall -Wremarks -Wcheck -wd161 -O3 -g -vec_report0 -fno-alias -fargument-noalias
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm -Wstrict-aliasing -fstrict-aliasing

deps = Makefile rng.h util.h
bin = $(patsubst %.c,%, $(wildcard test*.c)) \
      $(patsubst %.c,%_d, $(wildcard test*.c)) \
      mptest

all: $(bin)

$(bin): $(deps)

%_d: %.c
	$(CC) -o $@ $< -lm

%: %.c
	$(CX) -o $@ $<

mptest: mptest.c $(deps)
	$(CC) -fopenmp -o $@ $< -lm

clean:
	rm -f $(bin) out* MTSEED mt.dat *~

.PHONY: all clean

