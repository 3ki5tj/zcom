subdirs =

clean:
	$(RM) -f *~ *.tmp *.sc.py
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

# create self-contained python scripts
scpy:
	python pyembed.py -v -a -o cspacer.sc.py cspacer.py
	python pyembed.py -v -a -o cindent.sc.py cindent.py
	python pyembed.py -v -a -o rstrip.sc.py rstrip.py

.PHONY: clean pack

