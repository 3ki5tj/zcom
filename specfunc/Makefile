CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): specfunc.c specfunc.h Makefile

%: %.c
	$(CC) -o $@ $< -lm

clean:
	rm -f $(bin) *~ *.o

.PHONY: all clean

