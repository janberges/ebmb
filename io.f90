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

      read (unit, *) i%upper ! general cutoff frequency (eV)
      read (unit, *) i%lower ! Coulomb cutoff frequency (eV)

      read (unit, *) i%continue ! continue to real axis?
      read (unit, *) i%resolution ! real axis resolution

      read (unit, *) i%limit ! maximum number of steps
      read (unit, *) i%tiny ! negligible difference (a.u.)

      read (unit, *) i%form ! output format

      close (unit)
   end subroutine load

   subroutine save_text(i)
      type(info), intent(in) :: i

      integer :: n

      open (unit, file=i%name // '.out', action='write', status='replace')

      write (unit, '(3A23)') 'phiC/eV', 'mu*EB', 'Tc/K'
      write (unit, '(3ES23.14E3)') i%phiC, i%muStarEB, i%Tc

      write (unit, "(/, 'Imaginary-axis solution (', I0, '):')") i%status

      write (unit, '(/, 3A23)') 'omega/eV', 'Z', 'Delta/eV'

      do n = lbound(i%omega, 1), ubound(i%omega, 1)
         write (unit, '(3ES23.14E3)') i%omega(n), i%Z(n), i%Delta(n)
      end do

      if (i%continue) then
         write (unit, "(/, 'Real-axis solution:')")

         write (unit, '(/, A23)') 'Delta0/eV'
         write (unit, "(ES23.14E3, ' (', I0, ')')") i%Delta0, i%statusDelta0

         write (unit, '(/, 3A23)') 'omega/eV', 'Re[Delta]/eV', 'Im[Delta]/eV'

         do n = 1, i%resolution
            write (unit, '(3ES23.14E3)') i%energy(n), i%gap(n)
         end do
      end if

      close (unit)
   end subroutine save_text

   subroutine save_data(i)
      type(info), intent(in) :: i

      open (unit, file=i%name // '.dat', action='write', status='replace', &
         form='unformatted', access='stream')

      write (unit) i%status, size(i%omega)
      write (unit) i%omega, i%Z, i%Delta
      write (unit) i%phiC, i%muStarEB, i%Tc

      write (unit) merge('T', 'F', i%continue)

      if (i%continue) then
         write (unit) i%resolution
         write (unit) i%energy, real(i%gap), aimag(i%gap)
         write (unit) i%Delta0, i%statusDelta0
      end if

      close (unit)
   end subroutine save_data
end module io
