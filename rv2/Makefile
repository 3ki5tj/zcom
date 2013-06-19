INCLUDE =
CC = gcc -O0 -g3 -Wall -Wextra -pedantic -DRV2_DEBUG $(INCLUDE)
LM = -lm
CX = icc -O3 -Wall -Wextra $(INCLUDE)
LMX =

deps = def.h rv2.h Makefile
bin_d = $(patsubst %.c,%_d, $(wildcard test*.c))
bin = $(patsubst %.c,%, $(wildcard test*.c))

$(bin_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bin) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

rv2.h: rv2.0.h addfd.py
	python addfd.py

clean:
	rm -f $(bin) $(bin_d) *~ *.o a.out
	$(MAKE) -C doc clean

.PHONY: all clean

