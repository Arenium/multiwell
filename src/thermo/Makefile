# absoft flags
#FC = f77
#FFLAGS =  -O -s -lm -lfio  -ffortran-bounds-check

# GNU FORTRAN compiler
#FC = gfortran
##FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics
#For de-bugging:
#FFLAGS = -O -fno-automatic -fbounds-check -std=f95 -Wextra -Wall -pedantic -fall-intrinsics -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans
#FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics -g -fbacktrace -ffpe-trap=zero,overflow,underflow -std=f95

FC = gfortran
FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics

OTHERM = datetime.o rotunits.o qhinder.o classrot.o alngam.o \
         lenstr.o quantrot.o qmorse.o thermo.o qboxmod.o qelect.o \
         qdrso.o MaxVeff.o read_dat.o qsymtop.o thermochem.o \
         hindrance_fit.o write_out.o element.o ghrlev.o \
         shrlev.o uhrlev.o eis.o chemformula.o tuneck.o qvbint.o \
         ucase.o qtrans.o


thermo:    $(OTHERM)
	$(FC)   $(FFLAGS) -o thermo $(OTHERM);
	cp thermo ../../bin/

clean:
	rm -f *.o thermo ../../bin/thermo
