CC=gcc -O0 -g3 -Wall -Wextra -pedantic -lm
deps=rv2.c rv2.h def.h

bin=testang 
all: $(bin)

testang: testang.c $(deps)
	$(CC) -o $@ $< -lm

clean:
	rm -f $(bin) *~ *.o 

.PHONY: all clean

