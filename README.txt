E B _ i s o

This  programs  solves  the  isotropic Eliashberg  equations  for  the  Hubbard-
Holstein model. The names  of the parameter files to be  processed are passed as
command-line arguments. The required file format is exemplified by 'example.in'.


I n s t a l l a t i o n

Simply  run  'make', provided  recent  versions  of  Fortran compiler  and  Make
utility are installed. Compiler preferences may be adjusted in the Makefile.


M a n u a l

Format of the bare (.dat) output files:

    integer  status(Z, Delta)  |  character  'T' (continue) or 'F' (EOF)
    integer  n                 |
                               |    integer  N
  n doubles  omega (eV)        |
  n doubles  Z                 |  N doubles  omega (eV)
  n doubles  Delta (eV)        |  N doubles  Re[Delta] (eV)
                               |  N doubles  Im[Delta] (eV)
     double  phiC (eV)         |
                               |     double  Delta0 (eV)
     double  Tc (K)            |    integer  status(Delta0)


A c k n o w l e d g m e n t

Parts of the program are inspired by the EPW code and work of Malte Rösner.
