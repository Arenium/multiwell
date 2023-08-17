# Compiler and flags for Linux/MacOS
FC = gfortran
FFLAGS = -O3 -fno-automatic -fbounds-check -mcmodel=medium
SOURCES = DECLARE_MOD.f MICRO_MOD.f bsort.f canon_min.f canon_rates.f canon_writeout.f check_input.f chemformula.f chkdens.f covsrt.f dnt.f dsort.f eis.f element.f finalprint.f find_jmax.f find_tmins.f funcs.f gaussj.f get_viblo.f hsort.f indexx.f lfit.f lowerc.f main.f multiwell_write.f nkinfint.f nrgconvert.f parabfit.f japrewrite.f qele.f qhinf.f qroq.f qvibf.f rank.f read_input.f sestamp.f sort_input.f stamp.f super_mol.f tghrlev.f tqstop.f tshrlev.f veff.f write_input.f
OBJECTS = $(SOURCES:.f=.o)
EXECUTABLE = ktools

# Compiler and flags for Windows cross-compilation
FC_WINDOWS = x86_64-w64-mingw32-gfortran
FFLAGS_WINDOWS = -O3 -fno-automatic -fbounds-check -static
SOURCES_WINDOWS = $(SOURCES)
OBJECTS_WINDOWS = $(SOURCES_WINDOWS:.f=.wo)
EXECUTABLE_WINDOWS = ktools.exe

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

