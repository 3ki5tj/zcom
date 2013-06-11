subdirs =

clean:
	$(RM) -f *~ */*~ */*/*~ */a.out *.tmp *.sc.py
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

.PHONY: clean pack

