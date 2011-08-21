CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = abemin

all: $(bin)

abemin: abemin.c $(dep)
	$(CX) -o $@ $< 

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

install::
	cp abemin ~/app/bin/

clean:
	rm -f *~ *.out MTSEED $(bin) shakefail* rattlefail* tmh*.e tmh*.t tmh*.ehis TRACE ab.pos

.PHONY: clean all

