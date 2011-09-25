CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm
deps=def.h rv3.h

bin=testang testprinc testvdist testdet testeig testsvd
all: $(bin)

.c: $(deps)
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

