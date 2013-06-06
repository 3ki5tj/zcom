CC = gcc -Wall -Wextra -O0 -g3 -pedantic
LM = -lm
# optimized compiling
CX = icc -Wall -O3
LMX =
#CX = gcc -Wall -O3 -march=native -mfpmath=sse -funroll-all-loops
#LMX = -lm
deps = hmc.h md.c md.h rv3.h rv2.h rng.h def.h util.h Makefile
bin = $(patsubst %.c,%,$(wildcard test*.c))
bin_d = $(patsubst %.c,%_d,$(wildcard test*.c))

all: $(bin) $(bin_d)

$(bin) : % : %.c $(deps)
	$(CX) -o $@ $< $(LM)

$(bin_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

clean:
	rm -f *~ *.out MTSEED *.dat *.his $(bin) $(bin_d)

.PHONY: clean all

