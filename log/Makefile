INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic $(INCLUDE)

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin) : %: %.c log.h log.c Makefile
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o a.out my.log LOG TRACE

.PHONY: all clean

