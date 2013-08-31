bin=test1
hdr=pdb.h
src=pdb.c
dep=$(src) $(hdr) Makefile

all: $(bin)

CC = gcc -O0 -g3 -Wall -pedantic -Wextra
#CC = icc -O3 -Wall -Wremarks -Wcheck

test1: test1.c $(dep)
	$(CC) -o $@ $< -lm

subdirs = pdbalign

clean:
	rm -f $(bin) *~ out.pdb a.out
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done

.PHONY: all clean

