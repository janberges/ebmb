module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: qp = selected_real_kind(30)

   real(dp), parameter :: pi = 4 * atan(1.0_dp) ! 3.14159...
   real(dp), parameter :: kB = 8.61733e-05_dp   ! Boltzmann constant (meV/K)

   type universal
      character(:), allocatable :: name

      character(11) :: mode = 'self-energy'

      real(dp), pointer :: variable => null() ! parameter to be optimized

      real(dp) :: T = 0 ! temperature (K)

      real(dp) :: TcMD ! McMillan's critical temperature (eV)

      real(dp), allocatable :: TcEB(:) ! Eliashberg's critical temperature (eV)

      real(dp) :: error = 1e-05_dp ! valid error of critical temperature (K)
      real(dp) :: bound = 1e+00_dp ! lower bound of critical temperature (K)
      real(dp) :: small = 1e-10_dp ! maximum gap at critical temperature (eV)

      integer :: bands = 1 ! number of electronic bands

      real(dp) :: omegaE = 0.02_dp ! Einstein frequency (eV)

      real(dp), allocatable :: lambda(:, :) ! electron-phonon coupling
      real(dp), allocatable :: muStar(:, :) ! Coulomb pseudo-potential

      logical :: chi = .false. ! consider full DOS and calculate chi?

      real(dp) :: upper = 15.0_dp ! overall cutoff (omegaE)
      real(dp) :: lower = -1.0_dp ! Coulomb cutoff (omegaE)

      integer :: limit = 250000 ! maximum number of fixed-point steps

      logical :: measurable = .false. ! find measurable gap?
      integer :: resolution = 0       ! real axis resolution

      character(4) :: form = 'both'     ! output format
      character(9) :: edit = 'ES15.6E3' ! number format

      logical :: standalone = .false. ! include parameters in output file?
      logical ::    rescale = .true.  ! rescale Coulomb pseudo-potential?
      logical ::    cutoffZ = .false. ! cut off renormalization function?

      real(dp), allocatable :: energy(:) ! free-electron energy (eV)
      real(dp), allocatable :: dos(:, :) ! density of Bloch states (a.u.)

      real(dp) ::  rate = 1e-01_dp ! growth rate for bound search
      real(dp) :: shift = 1e+02_dp ! eigenvalue shift for power method
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

   real(dp) :: epsilon = 1e-15_dp

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
