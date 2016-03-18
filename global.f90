module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: qp = selected_real_kind(30)

   real(dp), parameter :: pi = 4 * atan(1.0_dp)
   real(dp), parameter :: qe = 1.60217662e-19_dp ! elementary charge (C)
   real(dp), parameter :: kB = 1.38064852e-23_dp ! Boltzmann constant (J/K)

   type universal
      character(:), allocatable :: name

      real(dp) :: kT ! temperature (eV)
      real(dp) :: Tc ! McMillan's critical temperature (K)

      logical :: critical ! find critical temperature?

      real(dp) :: small ! negligible gap (eV)
      real(dp) :: error ! error of critical temperature (eV)

      real(dp) :: omegaE ! Einstein frequency (eV)
      real(dp) :: lambda ! electron-phonon coupling
      real(dp) :: muStar ! Coulomb pseudo-potential

      real(dp) :: upper ! overall cutoff frequency (eV)
      real(dp) :: lower ! Coulomb cutoff frequency (eV)

      integer :: limit ! maximum number of fixed-point steps

      logical :: measurable ! find measurable gap?
      integer :: resolution ! real axis resolution

      character(4) :: form ! output format
   end type universal

   type matsubara
      real(dp) :: muStar ! rescaled Coulomb pseudo-potential

      integer :: status ! convergence status

      real(dp), allocatable :: omega(:) ! frequency (eV)
      real(dp), allocatable :: Delta(:) ! gap (eV)
      real(dp), allocatable :: Z(:) ! renormalization

      real(dp) :: phiC ! constant Coulomb contribution (eV)
   end type matsubara

   type continued
      real(dp), allocatable :: omega(:) ! frequency (eV)
      complex(dp), allocatable :: Delta(:) ! gap (eV)
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
