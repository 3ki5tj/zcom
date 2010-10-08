src=rng.c
hdr=rng.h
dep=$(src) $(hdr) Makefile
CC = icc -Wall -O0 -g3
CX = icc -Wall -O3 
#CX = gcc -Wall -O3 -g -march=pentium4 -funroll-all-loops -mfpmath=sse

all: cmp1 cmp2

# compares the result from test and testref
cmp1: test testref
	rm -f out1 out2
	./test > out1
	./testref > out2
	diff out1 out2

# test if two short runs produce the same result
# as a long run
cmp2: testsv test
	rm -f MTSEED out1 out3
	./test > out1
	./testsv > out3
	./testsv >> out3
	diff out1 out3

test: test.c $(dep) 
	$(CC) -o $@ $< $(src) -lm

testref: test.c zcomref.h
	$(CC) -DREF -o $@ $< -lm

testsv: testsv.c $(dep)
	$(CC) -o $@ $< $(src) -lm

testsp: testsv.c $(dep) 
	$(CX) -DTESTSP -DN=2000000000 -o $@ $< $(src) -lm

clean:
	rm -f test testsv testsp testref out* MTSEED mt.dat *~

.PHONY: all clean

