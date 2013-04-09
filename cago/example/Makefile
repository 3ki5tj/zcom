CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99
# optimized compiling
#CX = icc -Wall -O3 -g -vec_report0
#LM =
CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops -lm 
LM = -lm
src = zcom.h
dep = $(src) Makefile
bin = go

all: $(bin)

go: go.c $(dep)
	$(CX) -o $@ $< $(LM) 

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) *.his *.pos

.PHONY: clean all

