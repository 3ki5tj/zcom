INCLUDE = -I./include
CC = gcc -O0 -g3 -Wall -Wextra -pedantic -DRV3_DEBUG $(INCLUDE)
LM = -lm
CX = icc -O3 -Wall $(INCLUDE)
LMX =

deps = def.h rv3.h rm3.h rv3calcdih.h rv3gi.h Makefile
bin_d = $(patsubst %.c,%_d, $(wildcard test*.c))
bin = $(patsubst %.c,%, $(wildcard test*.c))

all: $(bin_d) $(bin)

$(bin_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bin) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

rv3.h: rv3.0.h addfd.py
	python addfd.py

rm3.h: rm3.0.h addfd.py
	python addfd.py

clean:
	rm -f $(bin) $(bin_d) *~ *.o a.out
	$(MAKE) -C doc clean

.PHONY: all clean

