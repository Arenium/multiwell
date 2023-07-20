# GNU FORTRAN-77 compiler 
#FC = g77
#FFLAGS =  -O -ffortran-bounds-check

# GNU FORTRAN compiler
#FC = gfortran
#For de-bugging:
#FFLAGS = -O -fno-automatic -fbounds-check -std=f95 -Wextra -Wall -pedantic
#Standard:
#FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics

FC = gfortran
FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics

Osctst = sctst.o datetime.o energyxyz.o energy.o lenstr.o calderixyz.o \
           nvmax.o nvmaxxyz.o ran1.o shrlev.o uhrlev.o ghrlev.o eis.o \
           ckderiv.o ckderivxyz.o convib.o READ_WXYZ.o rotlev.o krotlev.o  \
           trialmax.o boxlev.o morlev.o rotunits.o crotlev.o gam.o \
           alngam.o  omega.o calpn.o stoplev.o ucase.o

sctst:    $(Osctst)
	$(FC) $(FFLAGS) -o sctst $(Osctst) ;
	cp sctst ../../bin

clean:
	rm -f *.o sctst ../../bin/sctst
