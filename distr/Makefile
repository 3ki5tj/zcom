INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic $(INCLUDE)
LM=-lm
CX=icc -O3 -march=native -g  -Wall -vec_report0 $(INCLUDE)
LMX=

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin) : % : %.c def.h distr.h Makefile
	$(CX) -o $@ $< $(LMX)

test2: distr2d.h

testft: distr2d.h

subdirs = ene vol rdf bb

clean:
	rm -f $(bin) *~ *.o a.out
	-for d in $(subdirs); do ( $(MAKE) -C $$d clean ); done

.PHONY: all clean

