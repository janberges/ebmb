module global
   implicit none

   integer, parameter :: dp = selected_real_kind(15)
   integer, parameter :: qp = selected_real_kind(30)

   real(dp), parameter :: pi = 4 * atan(1.0_dp)
   real(dp), parameter :: qe = 1.60217662e-19 ! elementary charge (C)
   real(dp), parameter :: kB = 1.38064852e-23 ! Boltzmann constant (J/K)

   type info
      character(:), allocatable :: name

      real(dp) :: kT ! temperature (eV)

      real(dp) :: omega_E ! Einstein frequency (eV)
      real(dp) :: g ! electron-phonon coupling (eV)
      real(dp) :: U ! on-site Coulomb repulsion (eV)

      real(dp) :: DOS ! density of states per spin and unit cell (1/eV)

      real(dp) :: upper ! overall cutoff frequency (eV)
      real(dp) :: lower ! Coulomb cutoff frequency (eV)

      real(dp) :: mu ! Coulomb pseudo-potential

      logical :: continue ! continue to real axis?
      integer :: resolution ! real axis resolution

      integer :: limit ! maximum number of steps
      real(dp) :: tiny ! negligible difference (a.u.)

      character(4) :: form ! output format

      real(dp), allocatable :: omega(:) ! Matsubara frequencies (eV)
      real(dp), allocatable :: Delta(:) ! imaginary-axis gap (eV)

      real(dp), allocatable :: Z(:) ! renormalization

      real(dp), allocatable :: energy(:) ! real-axis energies (eV)
      complex(dp), allocatable :: gap(:) ! real axis gap (eV)

      real(dp) :: Delta0 ! leading gap (eV)

      integer :: status, statusDelta0 ! convergence status
   end type info
end module global
