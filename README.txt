E B _ i s o

This  programs  solves  the  isotropic Eliashberg  equations  for  the  Hubbard-
Holstein model. The names  of the parameter files to be  processed are passed as
command-line arguments. The required file format is exemplified by 'example.in'.


I n s t a l l a t i o n

Simply  run  'make', provided  recent  versions  of  Fortran compiler  and  Make
utility are installed. Compiler preferences may be adjusted in the Makefile.


M a n u a l

The Einstein frequency  'omegaE', the Electron-phonon coupling  'lambda' and the
Coulomb pseudo-potential 'mu*',  which must be given in the  parameter file, are
understood as variables of the critical temperature

       omegaE           1.04 (1 + lambda)
  Tc = ------ exp ------------------------------
       1.2 kB     mu* + 0.62 lambda mu* - lambda

according to McMillan (1968) and Dynes (1972).

Format of the bare (.dat) output files:

    integer  status(Z, Delta)  |  Only if c is 'T':
    integer  n                 |
                               |    integer  N
  n doubles  omega (eV)        |
  n doubles  Z                 |  N doubles  omega (eV)
  n doubles  Delta (eV)        |  N doubles  Re[Z]
                               |  N doubles  Im[Z]
     double  phiC (eV)         |  N doubles  Re[Delta] (eV)
     double  mu*EB             |  N doubles  Im[Delta] (eV)
     double  Tc (K)            |
                               |     double  Delta0 (eV)
  character  c ('T' or 'F')    |    integer  status(Delta0)


A c k n o w l e d g m e n t

Parts of the program are inspired by the EPW code and work of Malte Rösner.
