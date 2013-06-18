CC = gcc -O0 -g3 -Wall -Wextra -pedantic -Wno-variadic-macros
LM = -lm
usbdir = /media/C3/code
prj = zcom

$(prj).h::
	$(MAKE) -C rv3 rv3.h rm3.h
	$(MAKE) -C lj lj.h ljmc.h
	$(MAKE) -C abpro abpro.h
	python assemble.py -a -v 1

zip: $(prj).zip

$(prj).zip::
	git archive --format=zip -9 HEAD > $@

# test object file
$(prj).o: $(prj).c $(prj).h Makefile
	$(CC) -c -DZCOM_XFUNCS $< -o $@ $(LM)
	wc $@

subdirs = def util ss endn bio rng rc rv2 rv3 eig lu svd savgol specfunc \
	  argopt cfg log av hist mds pdb md ising2 potts2 lj abpro cago

clean:
	$(RM) -f *~ $(prj).o $(prj).zip */*~ */*/*~ */a.out *.tmp
	-for d in $(subdirs); do (cd $$d; $(MAKE) clean ); done
	-rstrip.py -Rv

pack: $(prj).zip
	wc $<
	gnome-open $<

usb: $(prj).h $(prj).zip
	mv $(prj).zip $(usbdir)/
	cp $(prj).h $(usbdir)/

usball::
	$(MAKE) clean
	$(MAKE) usb
	zip -r --symlinks --filesync -9 $(usbdir)/zcomall.zip * \
	  --exclude "*.swp" "*~"

# add symbolic links of header files that are referenced elsewhere
dodep::
	git add [a-z0-9]*/_*.h
	git add [a-z0-9]*/*/zcom.h

# run the code beautifier
cspacer::
	python cspacer/cspacer.py -RLwv

.PHONY: zip clean pack usb usball dodep cspacer

