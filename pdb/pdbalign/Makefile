CC = gcc -Wall -Wextra -O0 -g3 -pedantic
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -g -march=native -vec_report0
LMX =
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops
#LMX = -lm
src = zcom.h
dep = $(src) Makefile
bin = pdbalign

all: $(bin)

pdbalign: pdbalign.c $(dep)
	$(CX) -o $@ $< $(LMX)

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

install: $(bin)
	cp pdbalign ~/app/bin/

clean:
	rm -f *~ *.out MTSEED $(bin) fit*.pdb

.PHONY: clean all

