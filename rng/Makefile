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

test: test.c rng.c ss.c
	$(CC) -o $@ $^ -lm

testsv: test.c rng.c ss.c
	$(CC) -DSAVE -o $@ $^ -lm

testref: test.c zcomref.h
	$(CC) -DREF -o $@ $< -lm

clean:
	rm -f test testref testsv out* MTSEED *~

.PHONY: all clean

