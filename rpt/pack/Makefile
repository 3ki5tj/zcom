# CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CC = icc -Wall -O3 -g -vec_report0
src = zcom.h
dep = $(src) Makefile
bin = $(patsubst %.c,%, $(wildcard *.c))

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CC) -o $@ $<

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos 
	rm -rf bin obj *.depend *.layout *.dat

.PHONY: clean all

