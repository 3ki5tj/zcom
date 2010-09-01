all: test

test: test.c cfg.c ss.c
	$(CC) -o $@ $^ -lm

clean:
	rm -f test *~

.PHONY: all clean

