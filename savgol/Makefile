bin=test
hdr=savgol.h lu.h
src=savgol.c lu.c
dep=$(src) $(hdr) Makefile

all: $(bin)

CC = gcc -O0 -g3 -Wall -pedantic -Wextra

test: test.c $(dep)
	$(CC) -o $@ $< $(src) -lm

clean:
	rm -f $(bin) *~ savgol[12]d.dat

.PHONY: all clean

