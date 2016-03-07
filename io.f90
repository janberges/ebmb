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

      i%critical = i%kT .lt. 0 ! find critical temperature?

      i%kT = i%kT * kB / qe ! (eV)

      read (unit, *) i%small ! negligible gap (eV)
      read (unit, *) i%error ! error of critical temperature (K)

      i%error = i%error * kB / qe ! (eV)

      read (unit, *) i%omegaE ! Einstein frequency (eV)
      read (unit, *) i%lambda ! Electron-phonon coupling
      read (unit, *) i%muStar ! Coulomb pseudo-potential

      read (unit, *) i%upper ! general cutoff frequency (omegaE)
      read (unit, *) i%lower ! Coulomb cutoff frequency (omegaE)

      if (i%lower .lt. 0) i%lower = i%upper

      i%upper = i%upper * i%omegaE ! (eV)
      i%lower = i%lower * i%omegaE ! (eV)

      read (unit, *) i%limit ! maximum number of fixed-point steps

      read (unit, *) i%measurable ! find measurable gap?
      read (unit, *) i%resolution ! real axis resolution

      read (unit, *) i%form ! output format

      read (unit, *) negligible_difference ! negligible float difference

      close (unit)
   end subroutine load

   subroutine save_text(i, im, re)
      type(info), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      integer :: n

      open (unit, file=i%name // '.out', action='write', status='replace')

      write (unit, "('Rescaled Coulomb pseudo-potential:')")
      write (unit, '(/, ES23.14E3)') im%muStar

      write (unit, "(/, 'McMillan and Dynes'' critical temperature:')")
      write (unit, "(/, ES23.14E3, ' K')") i%Tc

      if (i%critical) then
         write (unit, "(/, 'Eliashberg''s critical temperature:')")
         write (unit, "(/, ES23.14E3, ' K')") i%kT * qe / kB
         close (unit)
         return
      end if

      write (unit, "(/, 'Imaginary-axis solution (', I0, '):')") im%status
      write (unit, '(/, 3A23)') 'omega/eV', 'Z', 'Delta/eV'

      do n = lbound(im%omega, 1), ubound(im%omega, 1)
         write (unit, '(3ES23.14E3)') im%omega(n), im%Z(n), im%Delta(n)
      end do

      write (unit, '(/, A23)') 'phiC/eV'
      write (unit, '(ES23.14E3)') im%phiC

      if (i%measurable) then
         write (unit, "(/, 'Measurable gap (', I0, '):')") re%status
         write (unit, "(/, ES15.6E3, ' eV')") re%Delta0
      end if

      if (i%resolution .gt. 0) then
         write (unit, "(/, 'Real-axis solution:')")

         write (unit, '(/, 5A15)') &
            'omega/eV', 'Re[Z]', 'Im[Z]', 'Re[Delta]/eV', 'Im[Delta]/eV'

         do n = 1, i%resolution
            write (unit, '(5ES15.6E3)') re%omega(n), re%Z(n), re%Delta(n)
         end do
      end if

      close (unit)
   end subroutine save_text

   subroutine save_data(i, im, re)
      type(info), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      open (unit, file=i%name // '.dat', action='write', status='replace', &
         form='unformatted', access='stream')

      write (unit) 'MUEB', im%muStar
      write (unit) 'TCMD', i%Tc

      if (i%critical) then
         write (unit) 'TCEB', i%kT * qe / kB
         close (unit)
         return
      end if

      write (unit) 'IMAG'
      write (unit) im%status, size(im%omega), im%omega, im%Z, im%Delta, im%phiC

      if (i%measurable) then
         write (unit) 'EDGE', re%status, re%Delta0
      end if

      if (i%resolution .gt. 0) then
         write (unit) 'REAL'
         write (unit) i%resolution, re%omega
         write (unit) real(re%Z), aimag(re%Z)
         write (unit) real(re%Delta), aimag(re%Delta)
      end if

      close (unit)
   end subroutine save_data
end module io
