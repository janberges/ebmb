compiler = gfortran
mode = validate

ifeq ($(compiler), gfortran)
  ifeq ($(mode), validate)
    options = -std=f2003 -Wall -Wno-maybe-uninitialized -pedantic
  else ifeq ($(mode), optimize)
    options = -O3
  endif
  override options += -J mod
else ifeq ($(compiler), g95)
  ifeq ($(mode), validate)
    options = -std=f2003 -Wall -pedantic
  else ifeq ($(mode), optimize)
    options = -O3
  endif
  override options += -fmod=mod
else ifeq ($(compiler), ifort)
  ifeq ($(mode), validate)
    options = -O0 -warn all
  else ifeq ($(mode), optimize)
    options = -O3
  endif
  override options += -module mod
else ifeq ($(compiler), f95)
  ifeq ($(mode), validate)
    options = -w4
  else ifeq ($(mode), optimize)
    options = -O3
  endif
  override options += -moddir=mod
endif

needless = *.aux *.dat .*.lb *.log mod/*.mod *.out *.pyc *.synctex *.synctex.gz ~temporary.* .DS_Store

# generated by makemake.py:

programs = critical ebmb tc

.PHONY: all clean cleaner

all: $(programs)

clean:
	@rm -f $(needless) eliashberg/eigenvalue.o eliashberg/self_energy.o eliashberg/self_energy_cdos.o global.o io/formatting.o io/load.o io/store.o io/tell.o programs/critical.o programs/ebmb.o programs/tc.o real_axis/pade.o real_axis/real_axis.o

cleaner: clean
	@rm -f $(programs)

$(programs):
	@echo link $@
	@$(compiler) -o $@ $^ $(external)

%.o: %.f90
	@echo compile $*
	@$(compiler) $(options) -c $< -o $@

critical: eliashberg/eigenvalue.o global.o io/load.o programs/critical.o
ebmb: eliashberg/self_energy.o eliashberg/self_energy_cdos.o global.o io/formatting.o io/load.o io/store.o io/tell.o programs/ebmb.o real_axis/pade.o real_axis/real_axis.o
tc: eliashberg/self_energy.o eliashberg/self_energy_cdos.o global.o io/formatting.o io/load.o programs/tc.o

eliashberg/eigenvalue.o: global.o
eliashberg/self_energy.o: global.o
eliashberg/self_energy_cdos.o: global.o
io/formatting.o: global.o
io/load.o: global.o
io/store.o: global.o
io/tell.o: global.o io/formatting.o
programs/critical.o: eliashberg/eigenvalue.o global.o io/load.o
programs/ebmb.o: eliashberg/self_energy.o eliashberg/self_energy_cdos.o global.o io/load.o io/store.o io/tell.o real_axis/real_axis.o
programs/tc.o: eliashberg/self_energy.o eliashberg/self_energy_cdos.o global.o io/formatting.o io/load.o
real_axis/pade.o: global.o
real_axis/real_axis.o: global.o real_axis/pade.o
