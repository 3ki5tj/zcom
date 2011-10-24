CC = gcc -Wall -Wextra -O0 -g3 -pedantic -lm
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 

deps = abpro.c abpro.h rv3.h rv2.h util.h md.c md.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c)) \
      $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin)

$(bin): $(deps)

%_d: %.c
	$(CC) -o $@ $< 

%: %.c
	$(CX) -o $@ $< 

abpro.c: abpro.0.c mk2d.py
	./mk2d.py $@

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) shakefail* rattlefail* ab.pos
	$(MAKE) clean -C tmhab
	$(MAKE) clean -C abemin

.PHONY: clean all

