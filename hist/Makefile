all: test1 testrdf test2

src=hist.c hist.h

CC=gcc -Wall -lm

test1: test1.c $(src)
	$(CC) -o $@ $< 

testrdf: testrdf.c $(src)
	$(CC) -o $@ $< 

test2: test2.c $(src)
	$(CC) -o $@ $< 

clean:
	rm -f test1 testrdf test2 a.out *~ *.o *.log HIST* RDF* rdf*.dat

.PHONY: all clean

