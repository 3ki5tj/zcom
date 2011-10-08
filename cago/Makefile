incl = -I./include
CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic $(incl)
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0 -wn1 $(incl)
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 

deps = cago.c cago.h def.h rv3.h md.c md.h pdb.h pdb.c util.c util.h rotfit.c rotfit.h \
       include/argopt.c include/argopt.h Makefile
bin := $(patsubst %.c,%, $(wildcard test*.c)) \
       $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin)

$(bin): $(deps)

% :: %.c
	$(CX) -o $@ $<

%_d :: %.c
	$(CC) -o $@ $<

clean:
	rm -f *~ *.out MTSEED *.dat *.pos $(bin) final.pdb ref.pdb

.PHONY: clean all

