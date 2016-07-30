E B _ l o c a l

This program solves the multiband Eliashberg equations for local self-energies
and Einstein phonon spectra.

For further information typeset and read 'manual.tex'.


I n s t a l l a t i o n

By default, the makefile builds the program with the GNU Fortran compiler and in
validation mode, i.e. with all warnings for the Fortran 2003 standard turned on.
There are two command-line arguments to change this behavior:

  $ make compiler=gfortran|g95|ifort|f95 mode=validate|optimize


A c k n o w l e d g m e n t

Parts of the program are inspired by the EPW code and work of Malte Rösner.
