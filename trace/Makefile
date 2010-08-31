all: test

test: test.c trace.c ss.c
	$(CC) -o $@ $^ 

clean:
	rm -f test *~ *.o TRACE

.PHONY: all clean

