bin=test1 test2
hdr=cfg.h ss.h opt.h def.h util.h
src=cfg.c ss.c
dep=$(src) $(hdr) Makefile

all: $(bin)

CC = gcc -O0 -g3 -Wall -pedantic -Wextra

$(bin): $(dep)

%:: %.c
	$(CC) -o $@ $< -lm

clean:
	rm -f $(bin) *~ a.out *.o

.PHONY: all clean

