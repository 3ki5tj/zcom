CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm
deps=rv3.c rv3.h def.h

bin=testang testprinc testvdist
all: $(bin)

testang: testang.c $(deps)
	$(CC) -o $@ $<

testprinc: testprinc.c $(deps)
	$(CC) -o $@ $< 

testvdist: testvdist.c $(deps)
	$(CC) -o $@ $< 

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

