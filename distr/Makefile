INCLUDE=-I./include
CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm $(INCLUDE)
CX=icc -O3 -g  -Wall -vec_report0 $(INCLUDE)

bin=$(patsubst %.c,%, $(wildcard test*.c))

all: $(bin)

$(bin): def.h distr.h Makefile

%: %.c
	$(CX) -o $@ $<

test2: distr2d.h

testft: distr2d.h

subdirs = ene vol rdf bb

clean:
	rm -f $(bin) *~ *.o a.out
	-for d in $(subdirs); do (cd $$d; $(MAKE) clean ); done

.PHONY: all clean

