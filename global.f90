module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: qp = selected_real_kind(30)

   real(dp), parameter :: pi = 4 * atan(1.0_dp)

   type universal
      character(:), allocatable :: name

      real(dp) :: T ! temperature (eV)
      real(dp) :: Tc ! McMillan's critical temperature (eV)

      logical :: critical ! find critical temperature?

      real(dp) :: error ! valid error of critical temperature (eV)
      real(dp) :: bound ! lower bound of critical temperature (eV)
      real(dp) :: small ! negligible gap (eV)

      real(dp) :: omegaE ! Einstein frequency (eV)
      real(dp) :: lambda ! electron-phonon coupling
      real(dp) :: muStar ! Coulomb pseudo-potential

      logical :: DOS ! consider full density of states?

      real(dp) :: upper ! overall cutoff frequency (eV)
      real(dp) :: lower ! Coulomb cutoff frequency (eV)

      integer :: limit ! maximum number of fixed-point steps

      logical :: measurable ! find measurable gap?
      integer :: resolution ! real axis resolution

      character(4) :: form ! output format
      logical :: standalone ! include parameters in output file?

      real(dp), allocatable :: energy(:) ! free-electron energy (eV)
      real(dp), allocatable :: density(:) ! density of Bloch states (a.u.)
      real(dp), allocatable :: weight(:) ! integration weight (eV)
   end type universal

   type matsubara
      real(dp) :: muStar ! rescaled Coulomb pseudo-potential

      real(dp), allocatable :: lambda(:) ! electron-phonon coupling
      real(dp), allocatable :: mu(:) ! Coulomb pseudo-potential

      integer :: status ! convergence status

      integer :: u ! index of overall cutoff frequency
      integer :: l ! index of Coulomb cutoff frequency

      real(dp), allocatable :: omega(:) ! frequency (eV)
      real(dp), allocatable :: Delta(:) ! gap (eV)
      real(dp), allocatable :: phi(:) ! order parameter (eV)
      real(dp), allocatable :: chi(:) ! energy shift (eV)
      real(dp), allocatable :: Z(:) ! renormalization

      real(dp) :: phiC ! constant Coulomb contribution (eV)
   end type matsubara

   type continued
      real(dp), allocatable :: omega(:) ! frequency (eV)
      complex(dp), allocatable :: Delta(:) ! gap (eV)
      complex(dp), allocatable :: chi(:) ! energy shift (eV)
      complex(dp), allocatable :: Z(:) ! renormalization

      real(dp) :: Delta0 ! measurable gap (eV)
      integer :: status ! convergence status
   end type continued

   real(dp) :: negligible_difference = 1e-15_dp

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

      ap = abs(lhs - rhs) .le. negligible_difference
   end function ap

   elemental function na(lhs, rhs)
      logical :: na
      real(dp), intent(in) :: lhs, rhs

      na = abs(lhs - rhs) .gt. negligible_difference
   end function na
end module global
