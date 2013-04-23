CC=gcc -O0 -g3 -Wall -Wextra -pedantic
LM=-lm
deps=def.h util.h

bin=test1 testsqr testsvdsolve
all: $(bin)

test1: test1.c $(deps)
	$(CC) -o $@ $< $(LM)

testsqr: testsqr.c $(deps)
	$(CC) -o $@ $< $(LM)

testsvdsolve: testsvdsolve.c $(deps)
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bin) *~ *.o

.PHONY: all clean

