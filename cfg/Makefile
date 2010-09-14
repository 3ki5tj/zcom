binf=test 
all: $(binf)

CC = gcc -O0 -g3 -Wall

hdrs=cfg.h ss.h die.h
src=cfg.c ss.c die.c

test: test.c $(src) $(hdrs)
	$(CC) -o $@ $< $(src) -lm

clean:
	rm -f $(binf) *~

.PHONY: all clean

