CC = gcc -O0 -g3 -Wall -pedantic -Wextra -lm
dep=def.h util.h argopt.h argopt.c Makefile

bin=test1
all: $(bin)


test1: test1.c $(dep)
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ a.out

.PHONY: all clean

