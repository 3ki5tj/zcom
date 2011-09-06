CC = gcc -O0 -g3 -Wall -Wextra -pedantic
CX = icc -O3 -Wall

zcom.h::
	python assemble.py -a -v1

zcom.zip::
	git archive --format=zip HEAD > zcom.zip

# test object file
zcom.o: zcom.c zcom.h Makefile
	$(CC) -c -DZCOM_XFUNCS $< -o $@
	wc $@

clean:
	rm -f *~ zcom.o */*~ */a.out
