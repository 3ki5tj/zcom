INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm $(INCLUDE)

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): log.h log.c Makefile

%: %.c
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o a.out my.log LOG TRACE

.PHONY: all clean

