all: test testrdf test2

src=hist.c hist.h

test: test.c $(src)
	$(CC) -o $@ $< 

testrdf: testrdf.c $(src)
	$(CC) -o $@ $< 

test2: test2.c $(src)
	$(CC) -o $@ $< 

clean:
	rm -f test testrdf test2 a.out *~ *.o *.log HIST* RDF* rdf*.dat

.PHONY: all clean

