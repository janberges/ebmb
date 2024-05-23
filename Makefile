FC = gfortran

flags_gfortran = -std=f2003 -pedantic -Wall -Wno-maybe-uninitialized
flags_ifort = -O0 -stand f03 -warn all

FFLAGS = ${flags_$(FC)}

bin/critical: LDLIBS = -llapack -lblas

needless = .DS_Store ebmb.pyc manual/ebmb.aux manual/.ebmb.lb manual/ebmb.log manual/ebmb.out manual/ebmb.synctex.gz ~temporary.dat

# generated by makemake90 bin=bin mod=modules:

modules_gfortran = -Jmodules
modules_ifort = -module modules

override FFLAGS += ${modules_$(FC)}

needless += src/critical.o src/dos.o src/ebmb.o src/eigenvalues.o src/eliashberg_eigenvalue.o src/eliashberg_eigenvalue_cdos.o src/eliashberg_self_energy.o src/eliashberg_self_energy_cdos.o src/eliashberg_spectral_function.o src/formatting.o src/global.o src/io_load.o src/io_store.o src/io_tell.o src/pade.o src/real_axis.o src/tc.o src/tools.o modules/*.mod

programs = bin/critical bin/ebmb bin/tc

.PHONY: all clean cleaner

all: $(programs)

clean:
	rm -f $(needless)

cleaner: clean
	rm -f $(programs)

$(programs):
	$(FC) $(FFLAGS) -o $@ $^ $(LDLIBS)

%.o: %.f90
	$(FC) $(FFLAGS) -c $< -o $@

bin/critical: src/critical.o src/eigenvalues.o src/eliashberg_eigenvalue.o src/eliashberg_eigenvalue_cdos.o src/eliashberg_self_energy.o src/eliashberg_spectral_function.o src/global.o src/io_load.o src/tools.o
bin/ebmb: src/dos.o src/ebmb.o src/eliashberg_self_energy.o src/eliashberg_self_energy_cdos.o src/eliashberg_spectral_function.o src/formatting.o src/global.o src/io_load.o src/io_store.o src/io_tell.o src/pade.o src/real_axis.o src/tools.o
bin/tc: src/eliashberg_self_energy.o src/eliashberg_self_energy_cdos.o src/eliashberg_spectral_function.o src/formatting.o src/global.o src/io_load.o src/tc.o src/tools.o

src/critical.o: src/eliashberg_eigenvalue.o src/eliashberg_eigenvalue_cdos.o src/global.o src/io_load.o
src/dos.o: src/eliashberg_self_energy.o src/global.o
src/ebmb.o: src/dos.o src/eliashberg_self_energy.o src/eliashberg_self_energy_cdos.o src/global.o src/io_load.o src/io_store.o src/io_tell.o src/real_axis.o
src/eigenvalues.o: src/global.o src/tools.o
src/eliashberg_eigenvalue.o: src/eigenvalues.o src/eliashberg_self_energy.o src/global.o
src/eliashberg_eigenvalue_cdos.o: src/eigenvalues.o src/eliashberg_spectral_function.o src/global.o
src/eliashberg_self_energy.o: src/eliashberg_spectral_function.o src/global.o src/tools.o
src/eliashberg_self_energy_cdos.o: src/eliashberg_spectral_function.o src/global.o
src/eliashberg_spectral_function.o: src/global.o src/tools.o
src/formatting.o: src/global.o
src/io_load.o: src/eliashberg_spectral_function.o src/global.o src/tools.o
src/io_store.o: src/global.o
src/io_tell.o: src/formatting.o src/global.o
src/pade.o: src/global.o
src/real_axis.o: src/global.o src/pade.o src/tools.o
src/tc.o: src/eliashberg_self_energy.o src/eliashberg_self_energy_cdos.o src/formatting.o src/global.o src/io_load.o
src/tools.o: src/global.o
