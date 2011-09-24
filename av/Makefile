CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm
deps=av.h av.c

bin=test1
all: $(bin)

test1: test1.c $(deps)
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

