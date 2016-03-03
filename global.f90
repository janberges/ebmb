module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: qp = selected_real_kind(30)

   real(dp), parameter :: pi = 4 * atan(1.0_dp)
   real(dp), parameter :: qe = 1.60217662e-19_dp ! elementary charge (C)
   real(dp), parameter :: kB = 1.38064852e-23_dp ! Boltzmann constant (J/K)

   type info
      character(:), allocatable :: name

      real(dp) :: kT ! temperature (eV)

      real(dp) :: omegaE ! Einstein frequency (eV)
      real(dp) :: lambda ! electron-phonon coupling
      real(dp) :: muStar ! Coulomb pseudo-potential

      real(dp) :: upper ! overall cutoff frequency (eV)
      real(dp) :: lower ! Coulomb cutoff frequency (eV)

      logical :: continue ! continue to real axis?
      integer :: resolution ! real axis resolution

      integer :: limit ! maximum number of steps
      real(dp) :: tiny ! negligible difference (a.u.)

      character(4) :: form ! output format

      real(dp), allocatable :: omega(:) ! Matsubara frequency (eV)
      real(dp), allocatable :: Delta(:) ! Matsubara gap (eV)
      real(dp), allocatable :: Z(:) ! Matsubara renormalization

      real(dp) :: phiC ! constant Coulomb contribution (eV)
      real(dp) :: muStarEB ! rescaled Coulomb pseudo-potential
      real(dp) :: Tc ! McMillan's critical temperature (K)

      real(dp), allocatable :: omega_(:) ! frequency (eV)
      complex(dp), allocatable :: Delta_(:) ! gap (eV)
      complex(dp), allocatable :: Z_(:) ! renormalization

      real(dp) :: Delta0 ! leading gap (eV)

      integer :: status, statusDelta0 ! convergence status
   end type info
end module global
