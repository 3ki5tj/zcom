INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic $(INCLUDE)
LM=-lm

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): def.h rc.h Makefile

%: %.c
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bin) *~ *.o a.out

.PHONY: all clean

