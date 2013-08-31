#CC=gcc -O0 -g3 -Wall -Wextra -pedantic
LM=-lm
CC=icc -Wall -Wremarks -Wcheck

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): specfunc.c specfunc.h Makefile

%: %.c
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bin) *~ *.o

.PHONY: all clean

