subdirs = cspacer pyembed objgen

clean:
	$(RM) -f *~ */*~ */*/*~ *.tmp *.m.py
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

.PHONY: clean

