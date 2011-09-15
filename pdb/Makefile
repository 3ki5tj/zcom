bin=test1
hdr=pdb.h
src=pdb.c
dep=$(src) $(hdr) Makefile

all: $(bin)

CC = gcc -O0 -g3 -Wall -pedantic -Wextra

test1: test1.c $(dep)
	$(CC) -o $@ $< $(src) -lm

clean:
	rm -f $(bin) *~ out.pdb a.out

.PHONY: all clean

