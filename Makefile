CC = gcc -O0 -g3 -Wall -Wextra -pedantic -Wno-variadic-macros
LM = -lm
usbdir = /media/C3
prj = zcom

$(prj).h::
	$(MAKE) -C lj lj.c
	$(MAKE) -C abpro abpro.c
	python assemble.py -a -v 1

zip: $(prj).zip

$(prj).zip::
	git archive --format=zip -9 HEAD > $@

# test object file
$(prj).o: $(prj).c $(prj).h Makefile
	$(CC) -c -DZCOM_XFUNCS $< -o $@ $(LM)
	wc $@

subdirs = def util ss endn bio rng rc rv2 rv3 eig lu svd rotfit savgol specfunc \
	  argopt cfg log av hist mds pdb clus ising2 potts2 md lj \
	  abpro cago objgen

clean:
	$(RM) -f *~ $(prj).o $(prj).zip */*~ */*/*~ */a.out *.tmp
	-for d in $(subdirs); do (cd $$d; $(MAKE) clean ); done
	-rstrip.py -R *.[ch] *.py *.ma *.txt *.cfg README* *akefile

pack: $(prj).zip
	wc $<
	gnome-open $<

usb: $(prj).h $(prj).zip
	mv $(prj).zip $(usbdir)/
	cp $(prj).h $(usbdir)/

usball::
	$(MAKE) clean
	$(MAKE) usb
	zip -r --symlinks --filesync -9 $(usbdir)/zcomall.zip *

dodep:
	git add [a-z0-9]*/*.h
	git add [a-z0-9]*/*/zcom.h
