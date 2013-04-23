incl = -I./include
CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic $(incl)
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -g -march=native -vec_report0 -wn1 $(incl)
LMX =
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops -lm

deps = cago.c cago.h def.h rv3.h md.c md.h pdb.h pdb.c util.h rotfit.c rotfit.h \
       include/argopt.c include/argopt.h Makefile
bin := $(patsubst %.c,%, $(wildcard test*.c))
bin_d := $(patsubst %.c,%_d, $(wildcard test*.c))

all: $(bin) $(bin_d)

$(bin) : % : %.c $(deps)
	$(CX) -o $@ $< $(LMX)

$(bin_d) : %_d : %.c $(deps)
	$(CC) -o $@ $< $(LM)

subdirs=gl example

clean:
	rm -f *~ *.out MTSEED *.dat *.pos $(bin) $(bin_d) final.pdb ref.pdb
	-for d in $(subdirs); do ($(MAKE) -C $$d clean); done

.PHONY: clean all

