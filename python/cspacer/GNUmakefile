subdirs =

clean:
	$(RM) -f *~ *.tmp *.m.py
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

# create self-contained python scripts
mpy:
	python pyembed.py -v -a -o cspacer.m.py cspacer.py
	python pyembed.py -v -a -o cindent.m.py cindent.py
	python pyembed.py -v -a -o rstrip.m.py rstrip.py

.PHONY: clean mpy

