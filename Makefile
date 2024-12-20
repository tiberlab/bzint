SHELL = /bin/sh

TOPDIR = $(CURDIR)

SRCDIR = src
TESTDIR = test
SUBDIRS = $(SRCDIR) $(TESTDIR)


.PHONY :  clean $(SUBDIRS)


all : $(SRCDIR)

$(SUBDIRS) :
	$(MAKE) -C $@


$(TESTDIR) : $(SRCDIR)


clean :
	for dir in $(SUBDIRS); do \
		$(MAKE) -C $$dir clean; \
	done
