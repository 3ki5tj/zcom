INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm $(INCLUDE)

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): util.h lp.h lp.c Makefile

%: %.c
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o a.out

.PHONY: all clean

