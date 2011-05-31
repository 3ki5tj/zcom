all: test testrdf

test: test.c
	$(CC) -o $@ $^ 

testrdf: testrdf.c
	$(CC) -o $@ $^ 

clean:
	rm -f test testrdf a.out *~ *.o *.log HIST* RDF* rdf*.dat

.PHONY: all clean

