CX  = icc -O3 -Wall -Wremarks -Wcheck
LMX =

CC = gcc -O0 -g3 -Wall -Wextra -pedantic -Wno-long-long
LM = -lm -lquadmath

bin = $(patsubst %.c, %, $(wildcard test*.c))
bin_d = $(patsubst %, %_d, $(bin))

src = argopt.h argopt.h opt.h util.h

all: $(bin) $(bin_d)

$(bin) : % : %.c $(src) Makefile
	$(CX) -o $@ $< $(LMX)

$(bin_d) : %_d : %.c $(src) Makefile
	$(CC) -o $@ $< $(LM)

clean:
	rm -f $(bin) $(bin_d) *~ *.o

.PHONY: all clean

