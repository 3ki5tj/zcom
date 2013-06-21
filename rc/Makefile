INCLUDE =
CC = gcc -O0 -g3 -Wall -Wextra -pedantic -DRV2_DEBUG $(INCLUDE)
LM = -lm
CX = icc -O3 -Wall -Wextra $(INCLUDE)
LMX =

autosrc = rc.h
deps = def.h $(autosrc) Makefile
bin_d = $(patsubst %.c,%_d, $(wildcard test*.c))
bin = $(patsubst %.c,%, $(wildcard test*.c))

$(bin_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bin) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

genhdr: $(autosrc)

$(autosrc) : %.h : %.0.h addfd.py
	python addfd.py

clean:
	rm -f $(bin) $(bin_d) *~ *.o a.out

.PHONY: all clean

