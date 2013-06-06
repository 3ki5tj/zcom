CC = gcc -Wall -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -march=native -vec_report0
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops -lm
deps = potts2.c potts2.h rng.h util.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c))
bin_d = $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin) $(bin_d)

$(bin_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< -lm

$(bin) : % : %.c $(deps)
	$(CX) -o $@ $<

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) $(bin_d)

.PHONY: clean all

