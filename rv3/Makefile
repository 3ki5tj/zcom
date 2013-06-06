INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic $(INCLUDE)
LM=-lm

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin) : % : %.c def.h rv3.h Makefile
	$(CC) -o $@ $< $(LM)

rv3.h: rv3.0.h addfd.py
	python addfd.py

clean:
	rm -f $(bin) *~ *.o a.out
	$(MAKE) -C doc clean

.PHONY: all clean

