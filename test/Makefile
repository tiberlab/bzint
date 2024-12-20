SHELL = /bin/sh


SRCDIR = ../src
FC = ifort
#FC = g95

FFLAGS = -I$(SRCDIR)
LDFLAGS = -L$(SRCDIR) -lBZint
OBJECTS = $(SRCDIR)/simplex2D.o $(SRCDIR)/geometry.o \
	  $(SRCDIR)/realroots.o $(SRCDIR)/qsurface.o \
	  $(SRCDIR)/lintrafo.o

TEST = test_theta test_theta2

TARGETS = test_ellipse test_line test_hyperbola \
	 test_parabola test_deg test_qsband test_qsband_theta

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@

.PHONY :  clean


all : $(TARGETS)


$(TARGETS) : % : %.o $(OBJECTS)
	$(FC) $(OBJECTS) $@.o -o $@

$(TEST) : % : %.o $(OBJECTS)
	$(FC) $(OBJECTS) $@.o -o $@

#ifndef DYNAMIC
#	$(FC) $(OBJECTS) $@.o -o $@
#else
#	$(FC) -shared $(LDFLAGS) $@.o -o $@
#endif

clean :
	rm -f *.o *.mod *~
	rm -f $(TARGETS) $(TEST)
