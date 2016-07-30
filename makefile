compiler = gfortran
mode = validate

ifeq ($(compiler), gfortran)
  ifeq ($(mode), validate)
    options = -std=f2003 -Wall -pedantic
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

programs = eb_local

.PHONY: all clean cleaner

all: $(programs)

clean:
	@rm -f $(needless) eliashberg/constant_dos.o eliashberg/variable_dos.o global.o io/arguments.o io/load.o io/store_data.o io/store_text.o main.o real_axis/intervals.o real_axis/pade.o real_axis/realaxis.o tc.o

cleaner: clean
	@rm -f $(programs)

$(programs):
	@echo link $@
	@$(compiler) -o $@ $^ $(external)

%.o: %.f90
	@echo compile $*
	@$(compiler) $(options) -c $< -o $@

eb_local: eliashberg/constant_dos.o eliashberg/variable_dos.o global.o io/arguments.o io/load.o io/store_data.o io/store_text.o main.o real_axis/intervals.o real_axis/pade.o real_axis/realaxis.o tc.o

eliashberg/constant_dos.o: global.o
eliashberg/variable_dos.o: global.o
io/load.o: global.o io/arguments.o
io/store_data.o: global.o
io/store_text.o: global.o
main.o: eliashberg/constant_dos.o eliashberg/variable_dos.o global.o io/arguments.o io/load.o io/store_data.o io/store_text.o real_axis/realaxis.o tc.o
real_axis/intervals.o: global.o
real_axis/pade.o: global.o
real_axis/realaxis.o: global.o real_axis/intervals.o real_axis/pade.o
tc.o: eliashberg/constant_dos.o eliashberg/variable_dos.o global.o
