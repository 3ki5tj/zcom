bin=test1
hdr=cfg.h ss.h
src=cfg.c ss.c
dep=$(src) $(hdr) Makefile

all: $(bin)

CC = gcc -O0 -g3 -Wall -pedantic -Wextra

%: %.c $(dep)
	$(CC) -o $@ $< $(src) -lm

clean:
	rm -f $(bin) *~

.PHONY: all clean

