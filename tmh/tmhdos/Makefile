CC = gcc -Wall -Wextra -O0 -g3 -pedantic
LM = -lm
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
LMX =
#CX = gcc -Wall -O3 -g -march=native -mfpmath=sse -funroll-all-loops -lm
#LMX = -lm
src = zcom.h
dep = $(src) Makefile
bin = tmhdos tmhdos_d

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CX) -o $@ $< $(LMX)

%_d: %.c
	$(CC) -o $@ $< $(LM)

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

install:: $(bin)
	cp tmhdos ~/app/bin/

clean:
	rm -f *~ *.out MTSEED *.dat $(bin) shakefail* rattlefail* tmh*.e tmh*.t tmh*.ehis TRACE

.PHONY: clean all

