# GNU FORTRAN-77 compiler 
#FC = g77
#FFLAGS =  -O -ffortran-bounds-check

# GNU FORTRAN compiler
#FC = gfortran
#For de-bugging:
#FFLAGS = -O -fno-automatic -fbounds-check -std=f95 -Wextra -Wall -pedantic
#Standard:
#FFLAGS = -O -fno-automatic -fbounds-check

FC = gfortran
FFLAGS = -O -fno-automatic -fbounds-check

Obdens = vibcalcs.o HRMOD.o SEPSUBS.o bdens.o datetime.o energyxyz.o energy.o lenstr.o calderixyz.o \
           nvmax.o nvmaxxyz.o ran1.o ckderiv.o ckderivxyz.o eis.o \
           trialmax.o \
           output.o calein.o energydoloops.o\
           energyxyzdoloops.o nvmaxdoloops.o nvmaxxyzdoloops.o\
           ucase.o gxsolver.o dgxsolver.o

bdens:    $(Obdens)
	$(FC) $(FFLAGS) -o bdens $(Obdens) ;
	cp bdens ../../bin
clean:
	rm -f *.o *.mod bdens main ../../bin/bdens
