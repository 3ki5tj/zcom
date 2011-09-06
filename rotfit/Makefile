CC = gcc -Wall -Wextra -pedantic -O0 -g3 -lm
#CC = icc -Wall -lm -O0 -g3

test1: test1.c rotfit.c rotfit.h Makefile 
	$(CC) $< -o $@

clean: 
	rm -rf *~ *.o *.obj *.exe *.ncb *.plg *.dat *.opt test1 Debug
