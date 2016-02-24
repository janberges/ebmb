E B _ i s o

This  programs  solves  the  isotropic Eliashberg  equations  for  the  Hubbard-
Holstein model. The names  of the parameter files to be  processed are passed as
command-line arguments. The required file format is exemplified by 'example.in'.


I n s t a l l a t i o n

Simply  run  'make', provided  recent  versions  of  Fortran compiler  and  Make
utility are installed. Compiler preferences may be adjusted in the Makefile.


M a n u a l

     double  mu*               |     double  Delta0 (eV)
    integer  status(Z, Delta)  |    integer  status(Delta0)
    integer  n                 |    integer  N
                               |
  n doubles  omega (eV)        |  N doubles  omega (eV)
  n doubles  Z                 |  N doubles  Re[Delta] (eV)
  n doubles  Delta (eV)        |  N doubles  Im[Delta] (eV)

  Format of the bare (.dat) output files. The right column is optional.


A c k n o w l e d g m e n t

Parts of the program are inspired by the EPW code and work of Malte Rösner.
