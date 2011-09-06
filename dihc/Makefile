CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm
deps=def.h rv3.h dihc.c dihc.h

bin=test1 test2
all: $(bin)

test1: test1.c $(deps)
	$(CC) -o $@ $<

test2: test2.c $(deps)
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

