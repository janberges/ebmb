compiler = gfortran

ifeq ($(compiler), gfortran)
  optional = -std=f2003 -Wall -pedantic
endif

optimize = true

ifeq ($(optimize), true)
  optional += -O3
endif

needless = *.aux *.dat .*.lb *.log *.out *.pyc *.synctex *.synctex.gz ~temporary.* .DS_Store

# generated by makemake.py:

programs = dos eb_local

.PHONY: all clean cleaner

all: $(programs)

clean:
	@rm -f $(needless) arguments.mod arguments.o dos.o eliashberg/constant_dos.o eliashberg/variable_dos.o eliashberg_constant_dos.mod eliashberg_variable_dos.mod filenames.mod filenames.o global.mod global.o integration.mod integration.o intervals.mod intervals.o io.mod io.o main.o pade.mod pade.o realaxis.mod realaxis.o tc.mod tc.o

cleaner: clean
	@rm -f $(programs)

$(programs):
	@echo link $@
	@$(compiler) -o $@ $^ $(external)

%.o: %.f90
	@echo compile $*
	@$(compiler) $(optional) -c $< -o $@

dos: arguments.o dos.o global.o intervals.o
eb_local: arguments.o eliashberg/constant_dos.o eliashberg/variable_dos.o filenames.o global.o integration.o intervals.o io.o main.o pade.o realaxis.o tc.o

dos.o: arguments.o global.o intervals.o
eliashberg/constant_dos.o: global.o
eliashberg/variable_dos.o: global.o
integration.o: global.o
intervals.o: global.o
io.o: filenames.o global.o integration.o
main.o: arguments.o eliashberg/constant_dos.o eliashberg/variable_dos.o global.o io.o realaxis.o tc.o
pade.o: global.o
realaxis.o: global.o intervals.o pade.o
tc.o: eliashberg/constant_dos.o eliashberg/variable_dos.o global.o
