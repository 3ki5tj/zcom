CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -lglut -lGLU
# optimized compiling
CX = icc -Wall -O3 -g -lm -vec_report0 -lglut -lGLU
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
dep = zcom.h Makefile
bin = glgo

all: $(bin)

$(bin): $(dep)

%: %.c
	$(CX) -o $@ $< 

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) shakefail* rattlefail* tmh*.e tmh*.t tmh*.ehis TRACE ab.pos

.PHONY: clean all

