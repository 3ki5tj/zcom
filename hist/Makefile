CC = gcc -Wall -O0 -g3 -lm -pedantic
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
LMX =
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
#LMX = -lmx
dep = util.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c))
bin_d = $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin) $(bin_d)

$(bin_d) : %_d: %.c $(deps)
	$(CC) -o $@ $< $(LM)

$(bin) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) $(bin_d) HIST*

.PHONY: all clean

