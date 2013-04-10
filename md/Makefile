CC = gcc -Wall -Wextra -O0 -g3 -pedantic -lm
# optimized compiling
CX = icc -Wall -O3
#CX = gcc -Wall -O3 -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
prj = md
src = $(prj).c
dep = $(src) $(prj).h rv3.h rv2.h def.h Makefile
bin = $(patsubst %.c,%,$(wildcard test*.c)) \
      $(patsubst %.c,%_d,$(wildcard test*.c))

all: $(bin)

%:: %.c $(dep)
	$(CX) -o $@ $< 

%_d:: %.c $(dep)
	$(CC) -o $@ $< -lm

clean:
	rm -f *~ *.out MTSEED *.dat *.his $(bin)

.PHONY: clean all

