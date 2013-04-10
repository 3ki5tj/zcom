INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm $(INCLUDE)

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): def.h rv3.h rotfit.c rotfit.h Makefile

%: %.c
	$(CC) -o $@ $< -lm

clean:
	rm -f $(bin) *~ *.o a.outi *.pos

.PHONY: all clean

