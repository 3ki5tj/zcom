CC=gcc -std=c99 -O0 -g3 -Wall -Wextra -pedantic
LM=-lm

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): %: %.c util.h util.c Makefile
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bin) *~ *.o

.PHONY: all clean

