E B _ i s o

This  programs  solves  the  isotropic Eliashberg  equations  for  the  Hubbard-
Holstein model. The names  of the parameter files to be  processed are passed as
command-line arguments. The required file format is exemplified by 'example.in'.


I n s t a l l a t i o n

Simply  run  'make', provided  recent  versions  of  Fortran compiler  and  Make
utility are installed. Compiler preferences may be adjusted in the Makefile.


M a n u a l

Format of the bare (.dat) output files:

     double  mu*               |  Only if c is T:
    integer  status(Z, Delta)  |
    integer  n                 |     double  Delta0 (eV)
                               |    integer  status(Delta0)
  n doubles  omega (eV)        |    integer  N
  n doubles  Z                 |
  n doubles  Delta (eV)        |  N doubles  omega (eV)
                               |  N doubles  Re[Delta] (eV)
  character  c (T or F)        |  N doubles  Im[Delta] (eV)


A c k n o w l e d g m e n t

Parts of the program are inspired by the EPW code and work of Malte Rösner.
