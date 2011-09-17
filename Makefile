CC = gcc -O0 -g3 -Wall -Wextra -pedantic
CX = icc -O3 -Wall
usbdir = /media/RII
prj = zcom

$(prj).h::
	cd abpro && ./mk2d.py && cd ..
	python assemble.py -a -v1

$(prj).zip::
	git archive --format=zip HEAD > $@

# test object file
$(prj).o: $(prj).c $(prj).h Makefile
	$(CC) -c -DZCOM_XFUNCS $< -o $@
	wc $@

clean:
	rm -f *~ $(prj).o $(prj).zip */*~ */a.out *.tmp

pack: $(prj).zip
	wc $<
	gnome-open $<

sync: $(prj).h $(prj).zip
	mv $(prj).zip $(usbdir)/
	cp $(prj).h $(usbdir)/

