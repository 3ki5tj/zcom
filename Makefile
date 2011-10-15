CC = gcc -O0 -g3 -Wall -Wextra -pedantic
CX = icc -O3 -Wall
usbdir = /media/RII
prj = zcom

subdirs = def util ss endn bio rng rv2 rv3 eig lu svd rotfit savgol specfunc \
	  argopt cfg log av hist mds pdb clus ising2 potts2 md lj abpro cago tmh

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
	$(RM) -f *~ $(prj).o $(prj).zip */*~ */*/*~ */a.out *.tmp
	-for d in $(subdirs); do (cd $$d; $(MAKE) clean ); done

pack: $(prj).zip
	wc $<
	gnome-open $<

usb: $(prj).h $(prj).zip
	mv $(prj).zip $(usbdir)/
	cp $(prj).h $(usbdir)/

dodep:
	git add [a-z0-9]*/*.h
	git add [a-z0-9]*/*/zcom.h
