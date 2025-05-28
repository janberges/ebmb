! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module globals
   implicit none
   public

   integer, parameter :: dp = selected_real_kind(15) ! double precision (8 B)
   integer, parameter :: qp = selected_real_kind(30) ! quad precision (16 B)
   integer, parameter :: i4 = selected_int_kind(9)   ! signed integer (4 B)

   real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp) ! 3.14159...

   real(dp), parameter :: eVSI = 1.602176634e-19_dp ! electronvolt (J)
   real(dp), parameter :: kBSI = 1.380649e-23_dp    ! Boltzmann constant (J/K)
   real(dp), parameter :: kB = kBSI / eVSI          ! Boltzmann constant (eV/K)

   integer, parameter :: fun = 11 ! file unit number

   type parameters
      character(1024) :: output = 'none'   ! name of output file
      character(50)   :: flomat = 'F16.12' ! number format

      logical :: tell = .true. ! use standard output?

      real(dp) :: T = 10.0_dp ! temperature (K)

      real(dp) :: omegaE  = 0.02_dp ! Einstein frequency (eV)
      real(dp) :: cutoff  = 15.0_dp ! overall cutoff frequency (omegaE)
      real(dp) :: cutoffC = -1.0_dp ! Coulomb cutoff frequency (omegaE)

      real(dp) :: omegaLog = 0.0_dp ! logarithmic avarage frequency (eV)
      real(dp) :: omega2nd = 0.0_dp ! second-moment avarage frequency (eV)

      integer(i4) :: bands = 1 ! number of electronic bands

      real(dp), allocatable :: lambda(:, :) ! electron-phonon coupling
      real(dp), allocatable :: muStar(:, :) ! Coulomb pseudo-potential

      real(dp), allocatable :: energy(:) ! free-electron energy (eV)
      real(dp), allocatable :: dos(:, :) ! density of Bloch states (1/eV)

      real(dp), allocatable :: omega(:) ! phonon energy (frequency argument)
      real(dp), allocatable :: a2F(:, :, :) ! Eliashberg spectral function

      logical :: ldos = .false. ! density of states given?
      logical :: la2F = .false. ! Eliashberg spectral function given?

      real(dp) ::  n = -1.0_dp ! initial occupancy number
      real(dp) :: mu =  0.0_dp ! initial chemical potential (eV)

      logical :: conserve = .true.  ! conserve particle number?
      logical :: chi      = .true.  ! consider energy shift?
      logical :: chiC     = .false. ! consider Coulomb part of energy shift?
      logical :: Sigma    = .false. ! calculate normal self-energy?

      integer(i4) :: steps = 250000 ! maximum number of iterations

      real(dp) :: error = 1e-05_dp ! bisection error (a.u.)
      real(dp) :: toln  = 1e-10_dp ! tolerance for occupancy number
      real(dp) ::  zero = 1e-10_dp ! negligible gap at critical temperature (eV)
      real(dp) ::  rate = 1e-01_dp ! growth rate for bound search

      real(dp) :: lower =  0.0_dp ! minimum real-axis frequency (eV)
      real(dp) :: upper = -1.0_dp ! maximum real-axis frequency (eV)

      integer(i4) :: points = 0 ! number of real-axis frequencies

      real(dp) :: logscale = 1.0_dp ! scaling of logarithmic sampling

      real(dp) :: eta = 1e-3_dp ! broadening of retarded objects (eV)

      logical :: measurable = .false. ! find measurable gap?

      logical :: unscale = .true.  ! unscale Coulomb pseudo-potential?
      logical :: rescale = .true.  ! rescale Coulomb pseudo-potential?
      logical :: imitate = .false. ! cut off renormalization function?

      logical :: divdos = .true.  ! divide by DOS at Fermi level?
      logical :: stable = .false. ! calculate quasiparticle DOS differently?
      logical :: normal = .false. ! enforce normal state?
      logical :: realgw = .false. ! do real-axis GW0 calculation?
      logical :: eta0Im = .true.  ! send broadening to zero in imaginary part?
      logical :: noZchi = .false. ! skip normal self-energy components?

      logical :: power = .true. ! use power method for single band?
   end type parameters

   type matsubara
      real(dp), allocatable :: omega(:)    ! frequency (eV)
      real(dp), allocatable :: Z    (:, :) ! renormalization
      real(dp), allocatable :: chi  (:, :) ! energy shift (eV)
      real(dp), allocatable :: Delta(:, :) ! gap (eV)
      real(dp), allocatable :: phi  (:, :) ! order parameter (eV)
      real(dp), allocatable :: chiC (:)    ! Coulomb part of energy shift (eV)
      real(dp), allocatable :: phiC (:)    ! Coulomb part of order parameter (eV)

      complex(dp), allocatable :: Sigma(:, :) ! normal self-energy (eV)

      integer(i4) :: steps ! steps until convergence
   end type matsubara

   type continued
      real   (dp), allocatable :: omega (:)    ! frequency (eV)
      complex(dp), allocatable :: Z     (:, :) ! renormalization
      complex(dp), allocatable :: chi   (:, :) ! energy shift (eV)
      complex(dp), allocatable :: Sigma (:, :) ! normal self-energy (eV)
      complex(dp), allocatable :: Delta (:, :) ! gap (eV)
      real   (dp), allocatable :: Delta0(:)    ! measurable gap (eV)
      real   (dp), allocatable :: dos   (:, :) ! quasiparticle density (1/eV)

      integer(i4), allocatable :: steps(:) ! steps until convergence
   end type continued

   type occupancy
      real(dp) :: states  ! integral of density of states
      real(dp) :: inspect ! integral of spectral function
      real(dp) :: n0, n   ! initial and final occupancy number
      real(dp) :: mu0, mu ! initial and final chemical potential (eV)
   end type occupancy

   real(dp) :: eps = 1e-13_dp ! negligible float difference (a.u.)

   interface operator(.ap.)
      module procedure ap, apc
   end interface

   interface operator(.na.)
      module procedure na, nac
   end interface

contains

   elemental function ap(lhs, rhs)
      logical :: ap
      real(dp), intent(in) :: lhs, rhs

      ap = abs(lhs - rhs) .le. eps
   end function ap

   elemental function na(lhs, rhs)
      logical :: na
      real(dp), intent(in) :: lhs, rhs

      na = abs(lhs - rhs) .gt. eps
   end function na

   elemental function apc(lhs, rhs)
      logical :: apc
      complex(dp), intent(in) :: lhs, rhs

      apc = abs(lhs - rhs) .le. eps
   end function apc

   elemental function nac(lhs, rhs)
      logical :: nac
      complex(dp), intent(in) :: lhs, rhs

      nac = abs(lhs - rhs) .gt. eps
   end function nac
end module globals
