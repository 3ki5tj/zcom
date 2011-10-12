CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm -lglut

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): def.h rv3.h glez.h glez.c Makefile

%: %.c
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

