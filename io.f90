module io
   use filenames
   use global
   implicit none

   private
   public :: load, save_text, save_data

   integer, parameter :: unit = 11

contains

   subroutine load(file, i)
      character(*), intent(in) :: file
      type(info), intent(out) :: i

      i%name = stem(file)

      open (unit, file=file, action='read', status='old')

      read (unit, *) i%kT ! temperature (K)

      i%kT = i%kT * kB / qe ! (eV)

      read (unit, *) i%omegaE ! Einstein frequency (eV)
      read (unit, *) i%lambda ! Electron-phonon coupling
      read (unit, *) i%muStar ! Coulomb pseudo-potential

      read (unit, *) i%upper ! general cutoff frequency (omegaE)
      read (unit, *) i%lower ! Coulomb cutoff frequency (omegaE)

      if (i%lower .lt. 0) i%lower = i%upper

      i%upper = i%upper * i%omegaE ! (eV)
      i%lower = i%lower * i%omegaE ! (eV)

      read (unit, *) i%limit ! maximum number of fixed-point steps

      read (unit, *) i%continue ! continue to real axis?
      read (unit, *) i%resolution ! real axis resolution

      read (unit, *) i%form ! output format

      read (unit, *) negligible_difference ! negligible float difference

      close (unit)
   end subroutine load

   subroutine save_text(i)
      type(info), intent(in) :: i

      integer :: n

      open (unit, file=i%name // '.out', action='write', status='replace')

      write (unit, "('Rescaled Coulomb pseudo-potential:')")
      write (unit, '(/, ES23.14E3)') i%muStarEB

      write (unit, "(/, 'McMillan and Dynes'' critical temperature:')")
      write (unit, "(/, ES23.14E3, ' K')") i%Tc

      write (unit, "(/, 'Imaginary-axis solution (', I0, '):')") i%status
      write (unit, '(/, 3A23)') 'omega/eV', 'Z', 'Delta/eV'

      do n = lbound(i%omega, 1), ubound(i%omega, 1)
         write (unit, '(3ES23.14E3)') i%omega(n), i%Z(n), i%Delta(n)
      end do

      write (unit, '(/, A23)') 'phiC/eV'
      write (unit, '(ES23.14E3)') i%phiC

      if (i%continue) then
         write (unit, "(/, 'Measurable gap (', I0, '):')") i%statusDelta0
         write (unit, "(/, ES15.6E3, ' eV')") i%Delta0

         write (unit, "(/, 'Real-axis solution:')")

         write (unit, '(/, 5A15)') &
            'omega/eV', 'Re[Z]', 'Im[Z]', 'Re[Delta]/eV', 'Im[Delta]/eV'

         do n = 1, i%resolution
            write (unit, '(5ES15.6E3)') i%omega_(n), i%Z_(n), i%Delta_(n)
         end do
      end if

      close (unit)
   end subroutine save_text

   subroutine save_data(i)
      type(info), intent(in) :: i

      open (unit, file=i%name // '.dat', action='write', status='replace', &
         form='unformatted', access='stream')

      write (unit) 'MUEB', i%muStarEB
      write (unit) 'TCMD', i%Tc

      write (unit) 'IMAG'
      write (unit) i%status, size(i%omega), i%omega, i%Z, i%Delta, i%phiC

      if (i%continue) then
         write (unit) 'EDGE', i%statusDelta0, i%Delta0

         write (unit) 'REAL'
         write (unit) i%resolution, i%omega_
         write (unit) real(i%Z_), aimag(i%Z_)
         write (unit) real(i%Delta_), aimag(i%Delta_)
      end if

      close (unit)
   end subroutine save_data
end module io
