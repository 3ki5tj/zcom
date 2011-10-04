CC = gcc -O0 -g3 -Wall -pedantic -Wextra -lm
hdr=argopt.h util.h def.h
src=argopt.c util.c 
dep=$(src) $(hdr) Makefile

bin=test1
all: $(bin)


test1: test1.c $(dep)
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ a.out

.PHONY: all clean

