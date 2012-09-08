CC = gcc -Wall -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
dep = potts2.c potts2.h rng.c rng.h util.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c)) \
      $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin)

$(bin): $(deps)

%_d: %.c
	$(CC) -o $@ $< 

%: %.c
	$(CX) -o $@ $< 

clean:
	rm -f *~ *.out MTSEED *.dat $(bin)

.PHONY: clean all

