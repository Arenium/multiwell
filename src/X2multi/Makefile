# absoft flags
#FC = f77
#FFLAGS =  -O -s -lm -lfio  -ffortran-bounds-check

# GNU FORTRAN-77 compiler
#FC = g77
#FFLAGS =  -O -ffortran-bounds-check

# GNU FORTRAN compiler
#FC = gfortran
#FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics

FC = gfortran
FFLAGS = -O -fno-automatic -fbounds-check -fall-intrinsics

 SUB1 = extractdatamominert.o \
        calcrotor.o writefreq.o\
        writemominert.o writedensum.o writethermo.o \
	X2multiwell.o writemultiwelldat.o writeCFG.o \
	OpSystem.o

GAUSS= ./gaussian/extractdataOUT.o \
       ./gaussian/readXYZforMominert.o

        all: gauss2multi

        gauss2multi:    $(SUB1) $(GAUSS)
	$(FC)   $(FFLAGS) -o gaussian2multi $(SUB1) $(GAUSS);
	cp gaussian2multi ../../bin;
	$(FC)   $(FFLAGS) -o install_gauss2multi ./gaussian/install_gauss2multi.f;
	./install_gauss2multi;
	cp gaussian2multi ../../bin;
	cp gauss2multi ../../bin;
	
        clean: 
	rm -f *.o ./*/*.o gauss2multi install_gauss2multi gaussian2multi ../../bin/gauss2multi ../../bin/gaussian2multi

