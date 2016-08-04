module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15) ! double precision
   integer, parameter :: qp = selected_real_kind(30) ! quad precision

   real(dp), parameter :: pi = 4 * atan(1.0_dp) ! 3.14159...
   real(dp), parameter :: kB = 8.61733e-05_dp   ! Boltzmann constant (meV/K)

   integer, parameter :: unit = 11 ! file unit number

   type universal
      character(99) :: file = 'none'   ! name of output file
      character(50) :: form = 'F16.12' ! number format

      logical :: tell = .true. ! use standard output?

      real(dp), pointer :: variable => null() ! parameter to be optimized

      real(dp) :: T = 10.0_dp ! temperature (K)

      real(dp) :: omegaE = 0.02_dp ! Einstein frequency (eV)
      real(dp) :: cutoff = 15.0_dp  ! overall cutoff frequency (omegaE)
      real(dp) :: cutout = -1.0_dp  ! Coulomb cutoff frequency (omegaE)

      integer :: bands = 1 ! number of electronic bands

      real(dp), allocatable :: lambda(:, :) ! electron-phonon coupling
      real(dp), allocatable :: muStar(:, :) ! Coulomb pseudo-potential

      real(dp), allocatable :: energy(:) ! free-electron energy (eV)
      real(dp), allocatable :: dos(:, :) ! density of Bloch states (a.u.)

      logical :: chi = .false. ! find energy shift?

      integer :: limit = 250000 ! maximum number of iterations

      real(dp) :: error = 1e-05_dp ! bisection error (a.u.)
      real(dp) ::  zero = 1e-10_dp ! negligible gap at critical temperature (eV)
      real(dp) ::  rate = 1e-01_dp ! growth rate for bound search
      real(dp) :: shift = 1e+02_dp ! eigenvalue shift for power method

      integer :: resolution = 0       ! real-axis resolution
      logical :: measurable = .false. ! find measurable gap?

      logical :: rescale = .true.  ! rescale Coulomb pseudo-potential?
      logical :: imitate = .false. ! cut off renormalization function?
   end type universal

   type matsubara
      real(dp), allocatable :: omega(:)    ! frequency (eV)
      real(dp), allocatable :: Z    (:, :) ! renormalization
      real(dp), allocatable :: chi  (:, :) ! energy shift (eV)
      real(dp), allocatable :: Delta(:, :) ! gap (eV)
      real(dp), allocatable :: phi  (:, :) ! order parameter (eV)
      real(dp), allocatable :: phiC (:)    ! constant Coulomb contribution (eV)

      integer :: status ! convergence status
   end type matsubara

   type continued
      real   (dp), allocatable :: omega (:)    ! frequency (eV)
      complex(dp), allocatable :: Z     (:, :) ! renormalization
      complex(dp), allocatable :: chi   (:, :) ! energy shift (eV)
      complex(dp), allocatable :: Delta (:, :) ! gap (eV)
      real   (dp), allocatable :: Delta0(:)    ! measurable gap (eV)

      integer, allocatable :: status(:) ! convergence status
   end type continued

   real(dp) :: epsilon = 1e-15_dp ! negligible float difference (a.u.)

   interface operator(.ap.)
      module procedure ap
   end interface

contains

   elemental function ap(lhs, rhs)
      logical :: ap
      real(dp), intent(in) :: lhs, rhs

      ap = abs(lhs - rhs) .le. epsilon
   end function ap
end module global
