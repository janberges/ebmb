module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: qp = selected_real_kind(30)

   real(dp), parameter :: pi = 4 * atan(1.0_dp)
   real(dp), parameter :: kB = 8.61733e-05_dp ! Boltzmann constant (eV/K)

   type universal
      character(:), allocatable :: name

      real(dp) :: T = 10.0_dp ! temperature (K)

      real(dp) :: TcMD ! McMillan's critical temperature (eV)

      real(dp), allocatable :: TcEB(:) ! Eliashberg's critical temperature (eV)

      logical :: critical = .false. ! find critical temperature?

      real(dp) :: error = 1e-03_dp ! valid error of critical temperature (K)
      real(dp) :: bound = 1e+00_dp ! lower bound of critical temperature (K)
      real(dp) :: small = 1e-10_dp ! maximum gap at critical temperature (eV)

      integer :: bands = 1 ! number of electronic bands

      real(dp) :: omegaE = 0.02_dp ! Einstein frequency (eV)

      real(dp), allocatable :: lambda(:, :) ! electron-phonon coupling
      real(dp), allocatable :: muStar(:, :) ! Coulomb pseudo-potential

      logical :: DOS = .false. ! consider full density of states?

      real(dp) :: upper = 10.0_dp ! overall cutoff (omegaE)
      real(dp) :: lower = -1.0_dp ! Coulomb cutoff (omegaE)

      integer :: limit = 100000 ! maximum number of fixed-point steps

      logical :: measurable = .false. ! find measurable gap?
      integer :: resolution = 0       ! real axis resolution

      character(4) :: form = 'both'     ! output format
      character(9) :: edit = 'ES15.6E3' ! number format

      logical :: standalone = .false. ! include parameters in output file?
      logical ::    rescale = .true.  ! rescale Coulomb pseudo-potential?

      real(dp), allocatable :: energy(:) ! free-electron energy (eV)
      real(dp), allocatable :: density(:, :) ! density of Bloch states (a.u.)
      real(dp), allocatable :: weight(:, :) ! integration weight (eV)
   end type universal

   type matsubara
      integer :: status ! convergence status

      real(dp), allocatable :: omega(:) ! frequency (eV)
      real(dp), allocatable :: Delta(:, :) ! gap (eV)
      real(dp), allocatable :: phi(:, :) ! order parameter (eV)
      real(dp), allocatable :: chi(:, :) ! energy shift (eV)
      real(dp), allocatable :: Z(:, :) ! renormalization

      real(dp), allocatable :: phiC(:) ! constant Coulomb contribution (eV)
   end type matsubara

   type continued
      real(dp), allocatable :: omega(:) ! frequency (eV)
      complex(dp), allocatable :: Delta(:, :) ! gap (eV)
      complex(dp), allocatable :: chi(:, :) ! energy shift (eV)
      complex(dp), allocatable :: Z(:, :) ! renormalization

      real(dp), allocatable :: Delta0(:) ! measurable gap (eV)
      integer, allocatable :: status(:) ! convergence status
   end type continued

   real(dp) :: negligible_difference = 1e-15_dp

   interface operator(.ap.)
      module procedure ap
   end interface

contains

   elemental function ap(lhs, rhs)
      logical :: ap
      real(dp), intent(in) :: lhs, rhs

      ap = abs(lhs - rhs) .le. negligible_difference
   end function ap
end module global
