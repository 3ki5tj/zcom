subdirs =

clean:
	$(RM) -f *~ *.tmp *.bak* tmp* *.o *.pyc at.c mb.c spb.c bb.c test.c
	-for d in $(subdirs); do ($(MAKE) -C $$d clean ); done
	-rstrip.py -Rv

.PHONY: clean pack

