CC = gcc -Wall -Wextra -O0 -g3 -pedantic -lm
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0 -fno-alias -fargument-noalias
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm -Wstrict-aliasing -fstrict-aliasing

deps = lj.c lj.h rv3.h rv2.h util.h md.c md.h Makefile
bin = $(patsubst %.c,%, $(wildcard test*.c)) \
      $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin)

$(bin): $(deps)

%_d: %.c
	$(CC) -o $@ $< 

%: %.c
	$(CX) -o $@ $< 

lj.c: lj.0.c mk2d.py
	./mk2d.py $@

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) shakefail* rattlefail* ab.pos
	rm -rf bin obj *.layout *.depend

.PHONY: clean all

