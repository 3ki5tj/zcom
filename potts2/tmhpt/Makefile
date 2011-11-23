CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic
# optimized compiling
CX = icc -Wall -O3 -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = tmhpt tmhpt0 tmhpt_e

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CX) -o $@ $< 

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -rf *~ *.out MTSEED *.dat $(bin) shakefail* rattlefail* tmh*.e tmh*.t tmh*.ehis TRACE pt.pos tmh*.tr bin obj pt2.depend

.PHONY: clean all

