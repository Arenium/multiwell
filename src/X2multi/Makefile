# Compiler and flags for Linux/MacOS
FC = gfortran
FFLAGS = -O3 -fno-automatic -fbounds-check -fall-intrinsics
SUB1 = extractdatamominert.o \
       calcrotor.o writefreq.o\
       writemominert.o writedensum.o writethermo.o \
	X2multiwell.o writemultiwelldat.o writeCFG.o \
	OpSystem.o

GAUSS = ./gaussian/extractdataOUT.o \
        ./gaussian/readXYZforMominert.o

EXECUTABLE = gauss2multi

# Compiler and flags for Windows cross-compilation
FC_WINDOWS = x86_64-w64-mingw32-gfortran
FFLAGS_WINDOWS = -O3 -fno-automatic -fbounds-check -fall-intrinsics -static
SUB1_W = $(SUB1)
GAUSS_W = $(GAUSS)
EXECUTABLE_WINDOWS = $(EXECUTABLE).exe

# Default target
all: unix windows

# Target for Unix
unix: $(EXECUTABLE)

$(EXECUTABLE): $(SUB1) $(GAUSS)
	$(FC) $(FFLAGS) -o $(EXECUTABLE) $(SUB1) $(GAUSS)
	cp $(EXECUTABLE) ../../bin/
	$(FC) $(FFLAGS) -o install_gauss2multi ./gaussian/install_gauss2multi.f
	./install_gauss2multi
	cp gaussian2multi ../../bin

# Target for Windows cross-compilation
windows: FC = $(FC_WINDOWS)
windows: FFLAGS = $(FFLAGS_WINDOWS)
windows: LIBS = -static-libgcc -static-libgfortran
windows: $(EXECUTABLE_WINDOWS)

$(EXECUTABLE_WINDOWS): $(SUB1_W) $(GAUSS_W)
	$(FC) $(FFLAGS) -o $(EXECUTABLE_WINDOWS) $(SUB1_W) $(GAUSS_W) $(LIBS)
	cp $(EXECUTABLE_WINDOWS) ../../bin/windows/
	$(FC) $(FFLAGS) -o install_gauss2multi ./gaussian/install_gauss2multi.f
	./install_gauss2multi
	cp gaussian2multi ../../bin/windows/

%.f: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.wo: %.f
	$(FC_WINDOWS) $(FFLAGS_WINDOWS) -c $< -o $@

clean:
	rm -f *.f *.wo $(EXECUTABLE) $(EXECUTABLE_WINDOWS)
	rm -f *.mod $(EXECUTABLE) $(EXECUTABLE_WINDOWS)
	rm -f ../../bin/$(EXECUTABLE) ../../bin/windows/$(EXECUTABLE_WINDOWS)

