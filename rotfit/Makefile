INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic $(INCLUDE)
LM=-lm

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin) : % : %.c def.h rv3.h rotfit.c rotfit.h Makefile
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bin) *~ *.o a.outi *.pos

.PHONY: all clean

