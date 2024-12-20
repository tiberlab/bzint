SHELL = /bin/sh

FC = ifort
#FC = g95
FFLAGS =

LIBNAME = BZint

SRC = simplex2D.f90 qsurface.f90 lintrafo.f90 \
      realroots.f90 geometry.f90
OBJECTS = $(SRC:.f90=.o)

%.o : %.f90
	$(FC) -c $(FFLAGS) $< -o $@



all : $(OBJECTS)

shared : all
	ifort -shared -o lib$(LIBNAME).so $(OBJECTS)

qsurface.o : realroots.o geometry.o lintrafo.o
simplex2D.o : qsurface.o


.PHONY :  clean
clean :
	rm -f *.o *.so *.mod *~
