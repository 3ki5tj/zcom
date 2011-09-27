CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm
deps=def.h rv3.h

bin=testang testprinc testvdist testdet testeig testsvd testrot testdih1 testdih2
all: $(bin)

%: %.c $(deps)
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

