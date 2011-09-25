CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm
deps=def.h rv3.h

bin=testang testprinc testvdist testdet
all: $(bin)

testang: testang.c $(deps)
	$(CC) -o $@ $<

testprinc: testprinc.c $(deps)
	$(CC) -o $@ $< 

testvdist: testvdist.c $(deps)
	$(CC) -o $@ $< 

testdet: testdet.c $(deps)
	$(CC) -o $@ $< 

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

