# Compiler and flags for Linux/MacOS
FC = gfortran
FFLAGS = -O3 -fno-automatic -fbounds-check -fall-intrinsics 
SOURCES = sctst.f datetime.f energyxyz.f energy.f lenstr.f calderixyz.f \
          nvmax.f nvmaxxyz.f ran1.f shrlev.f uhrlev.f ghrlev.f eis.f \
          ckderiv.f ckderivxyz.f convib.f READ_WXYZ.f rotlev.f krotlev.f  \
          trialmax.f boxlev.f morlev.f rotunits.f crotlev.f gam.f \
          alngam.f  omega.f calpn.f stoplev.f ucase.f
OBJECTS = $(SOURCES:.f=.o)
EXECUTABLE = sctst

# Compiler and flags for Windows cross-compilation
FC_WINDOWS = x86_64-w64-mingw32-gfortran
FFLAGS_WINDOWS = -O3 -fno-automatic -fbounds-check -fall-intrinsics -static
SOURCES_WINDOWS = $(SOURCES)
OBJECTS_WINDOWS = $(SOURCES_WINDOWS:.f=.wo)
EXECUTABLE_WINDOWS = $(EXECUTABLE).exe

# Default target
all: unix windows

# Target for Unix
unix: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FFLAGS) $(OBJECTS) -o $@
	cp $(EXECUTABLE) ../../bin/

# Target for Windows cross-compilation
windows: FC = $(FC_WINDOWS)
windows: FFLAGS = $(FFLAGS_WINDOWS)
windows: LIBS = -static-libgcc -static-libgfortran
windows: $(EXECUTABLE_WINDOWS)

$(EXECUTABLE_WINDOWS): $(OBJECTS_WINDOWS)
	$(FC) $(FFLAGS) $(OBJECTS_WINDOWS) -o $@ $(LIBS)
	cp $(EXECUTABLE_WINDOWS) ../../bin/windows

%.f: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.wo: %.f
	$(FC_WINDOWS) $(FFLAGS_WINDOWS) -c $< -o $@

clean:
	rm -f *.f *.wo $(EXECUTABLE) $(EXECUTABLE_WINDOWS)
	rm -f *.mod $(EXECUTABLE) $(EXECUTABLE_WINDOWS)
	rm -f ../../bin/$(EXECUTABLE) ../../bin/windows/$(EXECUTABLE_WINDOWS)

