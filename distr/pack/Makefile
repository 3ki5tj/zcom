# debug compiling
#CC = gcc -Wall -Wextra -O0 -g3 -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
#LM = -lm
# optimized compiling
CC = gcc -Wall -Wextra -O3 -march=native -g -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
LM = -lm
#CC = icc -Wall -O3 -g -vec_report0
#LM = 
src = zcom.h
dep = $(src) Makefile
bin = $(patsubst %.c,%, $(wildcard *.c))

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CC) -o $@ $< $(LM)

clean:
	rm -f *~ *.out MTSEED $(bin) lj.pos
	rm -rf bin obj *.depend *.layout ds.dat dsb.dat

.PHONY: clean all

