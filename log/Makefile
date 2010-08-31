all: test

test: test.c log.c ss.c
	$(CC) -o $@ $^ 

clean:
	rm -f test *~ *.o *.log

.PHONY: all clean

