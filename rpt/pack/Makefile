# CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CC = gcc -Wall -Wextra -O3 -g3 
src = zcom.h
dep = $(src) Makefile
bin = $(patsubst %.c,%, $(wildcard *.c))

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CC) $< -o $@ -lm

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos 
	rm -rf bin obj *.depend *.layout *.dat

.PHONY: clean all

