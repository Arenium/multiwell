# absoft flags
#FC = f77
#FFLAGS =  -O -s -lm -lfio  -ffortran-bounds-check

# GNU FORTRAN compiler
#FC = gfortran
#For de-bugging:
#FFLAGS = -O -fno-automatic -fbounds-check -std=f95 -Wextra -Wall -pedantic
#Standard:
#FFLAGS = -O -fno-automatic -fbounds-check

FC = gfortran
FFLAGS = -O -fno-automatic -fbounds-check

ODEN = HRMOD.o SEPSUBS.o densum.o sterab.o lenstr.o datetime.o ucase.o eis.o

densum:    $(ODEN)
	$(FC) $(FFLAGS) -o densum $(ODEN);
	cp densum ../../bin 

clean: 
	rm -f *.o *.mod densum ../../bin/densum

