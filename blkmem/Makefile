src=blkmem.c
hdr=blkmem.h
dep=$(src) $(hdr) Makefile
CC = gcc -Wall -Wextra -Wno-unknown-pragmas -O0 -g3 -pedantic -lm
# optimized compiling
CX = icc -Wall -Wremarks -Wcheck -wd161 -O3 -g -vec_report0 -fno-alias -fargument-noalias
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm -Wstrict-aliasing -fstrict-aliasing

deps = Makefile blkmem.h
bin = $(patsubst %.c,%, $(wildcard test*.c)) \
      $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin)

$(bin): $(deps)

%_d: %.c
	$(CC) -o $@ $< -lm

%: %.c
	$(CX) -o $@ $<

clean:
	rm -f $(bin) out* *~

.PHONY: all clean

