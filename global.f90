module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: qp = selected_real_kind(30)

   real(dp), parameter :: pi = 4 * atan(1.0_dp)
   real(dp), parameter :: kB = 8.61733e-05_dp ! Boltzmann constant (eV/K)

   type universal
      character(:), allocatable :: name

      real(dp) :: T ! temperature (eV)

      real(dp) :: TcMD ! McMillan's critical temperature (eV)

      real(dp), allocatable :: TcEB(:) ! Eliashberg's critical temperature (eV)

      logical :: critical ! find critical temperature?

      real(dp) :: error ! valid error of critical temperature (eV)
      real(dp) :: bound ! lower bound of critical temperature (eV)
      real(dp) :: small ! maximum gap at critical temperature (eV)

      integer :: bands ! number of electronic bands

      real(dp) :: omegaE ! Einstein frequency (eV)

      real(dp), allocatable :: lambda(:, :) ! electron-phonon coupling
      real(dp), allocatable :: muStar(:, :) ! Coulomb pseudo-potential

      logical :: DOS ! consider full density of states?

      real(dp) :: upper ! overall cutoff (eV)
      real(dp) :: lower ! Coulomb cutoff (eV)

      integer :: limit ! maximum number of fixed-point steps

      logical :: measurable ! find measurable gap?
      integer :: resolution ! real axis resolution

      character(4) :: form ! output format
      character(9) :: edit ! number format

      logical :: standalone ! include parameters in output file?
      logical :: rescale ! rescale Coulomb pseudo-potential?

      real(dp), allocatable :: energy(:) ! free-electron energy (eV)
      real(dp), allocatable :: density(:, :) ! density of Bloch states (a.u.)
      real(dp), allocatable :: weight(:, :) ! integration weight (eV)
   end type universal

   type matsubara
      real(dp), allocatable :: muStar(:, :) ! rescaled Coulomb pseudo-potential

      integer :: status ! convergence status

      integer :: u ! index of overall cutoff frequency
      integer :: l ! index of Coulomb cutoff frequency

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
