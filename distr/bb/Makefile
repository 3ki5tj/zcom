CC = gcc -Wall -Wextra -O0 -g3 -lm -pedantic -std=c99 -Wno-unknown-pragmas -Wno-variadic-macros
# optimized compiling
CX = icc -Wall -O3 -g -vec_report0
#CX = gcc -Wall -O3 -g -march=pentium4 -mfpmath=sse -funroll-all-loops -lm 
src = zcom.h
dep = $(src) Makefile
bin = bb

all: $(bin)

bb: bb.c $(dep)
	$(CX) -o $@ $< 

# run bb from bbcompat.txt
run: bb.png
	ls bb.png

bb.png: bbii.txt bb.gp
	gpload bb.gp

bbii.txt: $(bin) bbcompat.txt
	./bb

# forced run
frun::
	./bb
	gpload bb.gp

zcom.h::
	cd ../.. && python ./assemble.py -a -v -1

clean:
	rm -f *~ *.out MTSEED $(bin) 
	rm -rf bin obj *.depend *.layout bbii.txt bb.png 

.PHONY: clean all

