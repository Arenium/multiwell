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

RIB = lamm.o read_in.o cmxyz.o congi.o invg.o write_out.o

lamm:    $(RIB)
	$(FC) $(FFLAGS) -o lamm $(RIB);
	cp lamm ../../bin

clean:
	rm -f *.o lamm ../../bin/lamm

