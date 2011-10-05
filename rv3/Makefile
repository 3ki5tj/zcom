CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): def.h rv3.h Makefile

%: %.c
	$(CC) -o $@ $<

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

