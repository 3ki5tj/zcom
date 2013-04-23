CC=gcc -O0 -g3 -Wall -Wextra -pedantic

all: test

test: test.c lu.c
	$(CC) -Wall -o $@ $^

clean:
	rm -f test *~ *.o

.PHONY: all clean

