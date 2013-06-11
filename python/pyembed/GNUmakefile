subdirs =

clean:
	$(RM) -f *~ */*~ */*/*~ */a.out *.tmp *.m.py
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

mpy:
	python pyembed.py pyembed.py --mod="pyrmfunc" -o pyembed.m.py

.PHONY: clean mpy

