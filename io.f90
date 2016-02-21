module io
   use filenames
   use global
   implicit none

   private
   public :: load, store

   integer, parameter :: unit = 11

contains

   subroutine load(file, i)
      character(*), intent(in) :: file
      type(info), intent(out) :: i

      i%name = stem(file)

      open (unit, file=file, action='read', status='old')

      read (unit, *) i%kT ! temperature (K)

      i%kT = i%kT * kB / qe ! (eV)

      read (unit, *) i%omega_E ! Einstein frequency (eV)
      read (unit, *) i%g ! electron-phonon coupling (eV)
      read (unit, *) i%U ! on-site Coulomb repulsion (eV)

      read (unit, *) i%DOS ! density of states per spin and unit cell (1/eV)

      read (unit, *) i%upper ! general cutoff frequency (eV)
      read (unit, *) i%lower ! Coulomb cutoff frequency (eV)

      read (unit, *) i%mu ! Coulomb pseudo-potential

      read (unit, *) i%continue ! continue to real axis?
      read (unit, *) i%resolution ! real axis resolution

      read (unit, *) i%limit ! maximum number of steps
      read (unit, *) i%tiny ! negligible difference (a.u.)

      close (unit)
   end subroutine load

   subroutine store(i)
      type(info), intent(in) :: i

      integer :: n

      open (unit, file=i%name // '.out', action='write', status='replace')

      write (unit, '(2A23)') 'mu*'
      write (unit, '(2ES23.14E3)') i%mu

      write (unit, "(/, 'Imaginary-axis solution (', I0, '):')") i%status

      write (unit, '(/, 4A23)') 'omega/eV', 'Z', 'Delta/eV'

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
   end subroutine store
end module io
