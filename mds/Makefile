CC=gcc -O0 -g3 -Wall -Wextra -pedantic
LM=-lm
CC=icc -O2 -Wall -Wremarks -Wcheck
deps=util.h mds.h

bins = $(patsubst %.c,%, $(wildcard test*.c))

all: $(bins)

test1: test1.c $(deps)
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bin) *~ *.o

.PHONY: all clean

