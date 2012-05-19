CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm
deps=def.h util.h rng.h rng.c mds.c mds.h eig.c eig.h clus.c clus.h

bin=test1 test2 testload
all: $(bin)

test1: test1.c $(deps)
	$(CC) -o $@ $<

test2: test2.c $(deps)
	$(CC) -o $@ $<

testload: testload.c $(deps)
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o x.txt clus.txt 

.PHONY: all clean

