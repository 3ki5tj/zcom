CC = gcc -Wall -Wextra -pedantic -O0 -g3
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -march=native -g -vec_report0
LMX =
#CX = gcc -Wall -O3 -g -march=native -mcpu=native -mfpmath=sse -funroll-all-loops
#LMX = -lm
deps = ss.c ss.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c))
bin_d = $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin) $(bin_d)

$(bin): $(deps)

$(bin_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bin) : %: %.c $(deps)
	$(CX) -o $@ $< $(LMX)

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) $(bin_d)

.PHONY: clean all

