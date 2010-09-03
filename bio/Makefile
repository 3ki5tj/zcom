all: test

CC = gcc -O0 -g3 -Wall

test: test.c bio.h endn.c
	$(CC) -DENDN_DBG -o $@ $< endn.c

clean:
	rm -f io.bin *~ test *.orig *.spr

.PHONY: clean all
