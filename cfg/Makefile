binf=test test2
all: $(binf)

CC = gcc -O0 -g3 -Wall

test test2: cfg.h

test: test.c cfg.c ss.c
	$(CC) -o $@ $^ -lm

test2: test2.c cfg.c ss.c
	$(CC) -o $@ $^ -lm

clean:
	rm -f $(binf) *~

.PHONY: all clean

