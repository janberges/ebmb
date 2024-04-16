! Copyright (C) 2016-2024 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15) ! double precision (8 B)
   integer, parameter :: qp = selected_real_kind(30) ! quad precision (16 B)
   integer, parameter :: i4 = selected_int_kind(9)   ! signed integer (4 B)

   real(dp), parameter :: pi = 4 * atan(1.0_dp) ! 3.14159...
   real(dp), parameter :: kB = 8.61733e-5_dp    ! Boltzmann constant (meV/K)

   integer, parameter :: unit = 11 ! file unit number

   type parameters
      character(1024) :: file = 'none'   ! name of output file
      character(50)   :: form = 'F16.12' ! number format

      logical :: tell = .true. ! use standard output?

      real(dp) :: T = 10.0_dp ! temperature (K)

      real(dp) :: omegaE  = 0.02_dp ! Einstein frequency (eV)
      real(dp) :: cutoff  = 15.0_dp ! overall cutoff frequency (omegaE)
      real(dp) :: cutoffC = -1.0_dp ! Coulomb cutoff frequency (omegaE)

      integer(i4) :: bands = 1 ! number of electronic bands

      real(dp), allocatable :: lambda(:, :) ! electron-phonon coupling
      real(dp), allocatable :: muStar(:, :) ! Coulomb pseudo-potential

      real(dp), allocatable :: energy(:) ! free-electron energy (eV)
      real(dp), allocatable :: dos(:, :) ! density of Bloch states (a.u.)

      real(dp), allocatable :: omega(:) ! phonon energy (frequency argument)
      real(dp), allocatable :: a2F(:, :, :) ! Eliashberg spectral function

      logical :: ldos = .false. ! density of states given?
      logical :: la2F = .false. ! Eliashberg spectral function given?

      real(dp) ::  n = 0.0_dp ! initial occupancy number
      real(dp) :: mu = 0.0_dp ! initial chemical potential (eV)

      logical :: conserve = .true. ! conserve particle number?
      logical :: chi      = .true. ! consider energy shift?

      integer(i4) :: limit = 250000 ! maximum number of iterations

      real(dp) :: error = 1e-05_dp ! bisection error (a.u.)
      real(dp) ::  zero = 1e-10_dp ! negligible gap at critical temperature (eV)
      real(dp) ::  rate = 1e-01_dp ! growth rate for bound search

      real(dp) :: lower =  0.0_dp ! minimum real-axis frequency (eV)
      real(dp) :: upper = -1.0_dp ! maximum real-axis frequency (eV)

      real(dp) :: eta = 0.0_dp ! infinitesimal energy (eV)

      integer(i4) :: resolution = 0       ! real-axis resolution
      logical     :: measurable = .false. ! find measurable gap?

      logical :: unscale = .true.  ! unscale Coulomb pseudo-potential?
      logical :: rescale = .true.  ! rescale Coulomb pseudo-potential?
      logical :: imitate = .false. ! cut off renormalization function?

      logical :: normal = .false. ! enforce normal state?

      logical :: power = .true. ! use power method for single band?
   end type parameters

   type matsubara
      real(dp), allocatable :: omega(:)    ! frequency (eV)
      real(dp), allocatable :: Z    (:, :) ! renormalization
      real(dp), allocatable :: chi  (:, :) ! energy shift (eV)
      real(dp), allocatable :: Delta(:, :) ! gap (eV)
      real(dp), allocatable :: phi  (:, :) ! order parameter (eV)
      real(dp), allocatable :: phiC (:)    ! constant Coulomb contribution (eV)

      integer(i4) :: status ! convergence status
   end type matsubara

   type continued
      real   (dp), allocatable :: omega (:)    ! frequency (eV)
      complex(dp), allocatable :: Z     (:, :) ! renormalization
      complex(dp), allocatable :: chi   (:, :) ! energy shift (eV)
      complex(dp), allocatable :: Delta (:, :) ! gap (eV)
      real   (dp), allocatable :: Delta0(:)    ! measurable gap (eV)
      real   (dp), allocatable :: dos   (:, :) ! quasiparticle density (a.u.)

      integer(i4), allocatable :: status(:) ! convergence status
   end type continued

   type occupancy
      real(dp) :: n0, n   ! initial and final occupancy number
      real(dp) :: mu0, mu ! initial and final chemical potential (eV)
   end type occupancy

   real(dp) :: epsilon = 1e-15_dp ! negligible float difference (a.u.)

   interface operator(.ap.)
      module procedure ap
   end interface

   interface operator(.na.)
      module procedure na
   end interface

contains

   elemental function ap(lhs, rhs)
      logical :: ap
      real(dp), intent(in) :: lhs, rhs

      ap = abs(lhs - rhs) .le. epsilon
   end function ap

   elemental function na(lhs, rhs)
      logical :: na
      real(dp), intent(in) :: lhs, rhs

      na = abs(lhs - rhs) .gt. epsilon
   end function na

end module global
