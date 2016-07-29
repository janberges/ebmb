E B _ l o c a l

This program solves the multiband Eliashberg equations for local self-energies
and Einstein phonon spectra.

For further information typeset and read 'manual.tex'.


I n s t a l l a t i o n

By default, the makefile builds the program with the GNU Fortran compiler and in
validation mode, i.e. with all warnings for the Fortran 2003 standard turned on.
On the contrary, to obtain an optimized executable run

  $ make mode=optimize

Other compilers and options may also be chosen, e.g.

  $ make compiler=ifort options=-O3


A c k n o w l e d g m e n t

Parts of the program are inspired by the EPW code and work of Malte Rösner.
