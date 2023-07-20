FC = mpif90
FFLAGS = -O3
FLAGZERO = -finit-local-zero -fallow-argument-mismatch
OBJECTS = decl_alloc.o paradensum.o alngam.o boxlev.o calderixyz.o \
          ckderiv.o ckderivxyz.o convib.o crotlev.o datetime.o eis.o \
          new_energy.o energyxyz.o gam.o ghrlev.o krotlev.o lenstr.o \
          morlev.o ran1.o READ_WXYZ.o rotlev.o rotunits.o shrlev.o \
          uhrlev.o ucase.o stoplev.o wind_setup.o

paradensum:	$(OBJECTS)
	        $(FC) $(FFLAGS) $(OBJECTS) -o paradensum;
	        cp paradensum ../../bin/

%.o:	%.F90
	$(FC) $(FFLAGS) $(FLAGZERO) -c $< -o $@

%.o:	%.f
	$(FC) $(FFLAGS) $(FLAGZERO) -c $< -o $@

%.o:	%.f90
	$(FC) $(FFLAGS) $(FLAGZERO) -c $< -o $@

clean:
	rm -f *.o *.mod paradensum ../../bin/paradensum


                                                                   		
