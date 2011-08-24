CC=gcc -O0 -g3 -Wall -Wextra
deps=trace.c ss.c

all: test test0

test0: test0.c $(deps)
	$(CC) -Wall -o $@ $^

test: test.c $(deps)
	$(CC) -Wall -o $@ $^ 

clean:
	rm -f test *~ *.o TRACE

.PHONY: all clean

