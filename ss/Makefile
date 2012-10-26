CC = gcc -Wall -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
deps = ss.c ss.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c)) \
      $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin)

$(bin): $(deps)

%_d: %.c
	$(CC) -o $@ $< -lm

%: %.c
	$(CX) -o $@ $< 

clean:
	rm -f *~ *.out MTSEED *.dat $(bin)

.PHONY: clean all

