incl = -I./include
CC = gcc -Wall -Wextra -pedantic -O0 -g3 $(incl)
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -g -march=native -vec_report0 -wn1 $(incl)
LMX =
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops -lm

deps = cago.h rv3.h rm3.h md.h pdb.h util.h \
       include/argopt.h Makefile
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

