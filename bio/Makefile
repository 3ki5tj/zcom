all: test

CC = gcc -O0 -g3 -Wall -Wextra -pedantic

test: test.c bio.h endn.h
	$(CC) -DENDN_DBG -o $@ $<

clean:
	rm -f io.bin *~ test *.orig *.spr

.PHONY: clean all
