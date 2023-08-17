# Compiler and flags for Linux/MacOS
FC = gfortran
FFLAGS = -O3 -fno-automatic -fbounds-check -fall-intrinsics
SOURCES = datetime.f rotunits.f qhinder.f classrot.f alngam.f \
          lenstr.f quantrot.f qmorse.f thermo.f qboxmod.f qelect.f \
          qdrso.f MaxVeff.f read_dat.f qsymtop.f thermochem.f \
          hindrance_fit.f write_out.f element.f ghrlev.f \
          shrlev.f uhrlev.f eis.f chemformula.f tuneck.f qvbint.f \
          ucase.f qtrans.f
OBJECTS = $(SOURCES:.f=.o)
EXECUTABLE = thermo

# Compiler and flags for Windows cross-compilation
FC_WINDOWS = x86_64-w64-mingw32-gfortran
FFLAGS_WINDOWS = -O3 -fno-automatic -fbounds-check -fall-intrinsics -static
SOURCES_WINDOWS = $(SOURCES)
OBJECTS_WINDOWS = $(SOURCES_WINDOWS:.f=.wo)
EXECUTABLE_WINDOWS = thermo.exe

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
	cp $(EXECUTABLE_WINDOWS) ../../bin/windows/

%.o: %.f
	$(FC) $(FFLAGS) -c $< -o $@

%.wo: %.f
	$(FC_WINDOWS) $(FFLAGS_WINDOWS) -c $< -o $@

clean:
	rm -f *.o *.wo $(EXECUTABLE) $(EXECUTABLE_WINDOWS)
	rm -f *.mod $(EXECUTABLE) $(EXECUTABLE_WINDOWS)
	rm -f ../../bin/$(EXECUTABLE) ../../bin/windows/$(EXECUTABLE_WINDOWS)

