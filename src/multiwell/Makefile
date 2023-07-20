# GNU FORTRAN compiler
#FC = gfortran
#For de-bugging:
#FFLAGS = -O -fno-automatic -fbounds-check -std=f95 -Wextra -Wall -pedantic
#FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics -g -O0 -fbacktrace
#FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics 
#FFLAGS = -O -fno-automatic -g3 -fcheck=all -Wall -fbacktrace
#Standard:
#FFLAGS = -O -fno-automatic -fall-intrinsics 


FC = gfortran
FFLAGS = -O -fno-automatic -g3 -fcheck=all -fall-intrinsics

OMULTI =  declare_mod.o utili_mod.o subs_mod.o I_O_mod.o \
          bookstep_mod.o main.o
    
multiwell:    $(OMULTI)
	$(FC)   $(FFLAGS) -o multiwell $(OMULTI); 
	cp -f multiwell ../../bin 

clean:   
	rm -f *.o multiwell *.mod ../../bin/multiwell

