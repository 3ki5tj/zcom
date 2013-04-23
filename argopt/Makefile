CC=gcc -O0 -g3 -Wall -Wextra -pedantic
LM=-lm

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin) : % : %.c def.h argopt.h argopt.c opt.h Makefile
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bin) *~ *.o

.PHONY: all clean

