OPT = -Iinclude
CC = gcc -Wall -Wextra -O0 -g3 -pedantic $(OPT)
LM = -lm
# optimized compiling
CX = icc -Wall -Wremarks -Wcheck -O3 -march=native $(OPT)
LMX =
#CX = gcc -Wall -O3 -march=native -mfpmath=sse -funroll-all-loops
#LMX = -lm
prj = md
src = $(prj).h
deps = $(src) $(prj).h rv3.h rv2.h rng.h util.h Makefile
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

