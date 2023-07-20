# absoft flags
#FC = f77
#FFLAGS =  -O -s -lm -lfio  -ffortran-bounds-check

# GNU FORTRAN compiler
#FC = gfortran
#For de-bugging:
#FFLAGS = -O -fno-automatic -fbounds-check -std=f95 -Wextra -Wall -pedantic
#Standard:
#FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics

FC = gfortran
FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics

OTHERM = CMdistance.o MomInert.o rotate.o eis.o elemass.o


mominert:    $(OTHERM)
	$(FC)   $(FFLAGS) -o mominert $(OTHERM) ;
	cp mominert ../../bin

clean:
	rm -f *.o mominert ../../bin/mominert
