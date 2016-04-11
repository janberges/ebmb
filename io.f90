module io
   use filenames
   use global
   use integration
   implicit none

   private
   public :: load, store

   integer :: n, m
   integer, parameter :: unit = 11

contains

   subroutine load(file, i)
      character(*), intent(in) :: file
      type(universal), intent(out) :: i

      character(50) :: DOSfile
      real(dp), allocatable :: DOS(:)

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

      read (unit, *) DOSfile ! file with density of states

      i%DOS = DOSfile .ne. 'none' ! consider full density of states?

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

      if (i%DOS) then
         open (unit, file=DOSfile, action='read', status='old')

         read (unit, *) n ! density-of-states resolution

         allocate(DOS(n)) ! density of states (a.u.)

         allocate(i%energy(n)) ! free-electron energy (eV)
         allocate(i%weight(n)) ! integration weight (eV)

         do m = 1, n
            read(unit, *) i%energy(m), DOS(m)
         end do

         close (unit)

         call differential(i%energy, i%weight)

         n = minloc(abs(i%energy), 1) ! index of Fermi level

         i%weight = i%weight * DOS / DOS(n)
      end if
   end subroutine load

   subroutine store(i, im, re)
      type(universal), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      if (i%form .eq. 'text' .or. i%form .eq. 'both') call store_text(i, im, re)
      if (i%form .eq. 'data' .or. i%form .eq. 'both') call store_data(i, im, re)
   end subroutine store

   subroutine store_text(i, im, re)
      type(universal), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      open (unit, file=i%name // '.out', action='write', status='replace')

      write (unit, "('Rescaled Coulomb pseudo-potential:')")
      write (unit, '(/, ES23.14E3)') im%muStar

      write (unit, "(/, 'McMillan and Dynes'' critical temperature:')")
      write (unit, "(/, ES23.14E3, ' K')") i%Tc

      if (i%critical) then
         write (unit, "(/, 'Eliashberg''s critical temperature:')")
         write (unit, "(/, ES23.14E3, ' K')") i%kT * qe / kB
      end if

      write (unit, "(/, 'Imaginary-axis solution (', I0, '):')") im%status

      if (i%DOS) then
         write (unit, '(/, 4A23)') 'omega/eV', 'Z', 'chi/eV', 'Delta/eV'

         do n = 0, im%n - 1
            write (unit, '(4ES23.14E3)') &
               im%omega(n), im%Z(n), im%chi(n), im%Delta(n)
         end do
      else
         write (unit, '(/, 3A23)') 'omega/eV', 'Z', 'Delta/eV'

         do n = 0, im%n - 1
            write (unit, '(3ES23.14E3)') im%omega(n), im%Z(n), im%Delta(n)
         end do
      end if

      write (unit, '(/, A23)') 'phiC/eV'
      write (unit, '(ES23.14E3)') im%phiC

      if (i%measurable) then
         write (unit, "(/, 'Measurable gap (', I0, '):')") re%status
         write (unit, "(/, ES15.6E3, ' eV')") re%Delta0
      end if

      if (i%resolution .gt. 0) then
         write (unit, "(/, 'Real-axis solution:')")

         if (i%DOS) then
            write (unit, '(/, 7A15)') 'omega/eV', 'Re[Z]', 'Im[Z]', &
               'Re[chi]', 'Im[chi]', 'Re[Delta]/eV', 'Im[Delta]/eV'

            do n = 1, i%resolution
               write (unit, '(7ES15.6E3)') &
                  re%omega(n), re%Z(n), re%chi(n), re%Delta(n)
            end do
         else
            write (unit, '(/, 5A15)') &
               'omega/eV', 'Re[Z]', 'Im[Z]', 'Re[Delta]/eV', 'Im[Delta]/eV'

            do n = 1, i%resolution
               write (unit, '(5ES15.6E3)') re%omega(n), re%Z(n), re%Delta(n)
            end do
         end if
      end if

      close (unit)
   end subroutine store_text

   subroutine store_data(i, im, re)
      type(universal), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      open (unit, file=i%name // '.dat', action='write', status='replace', &
         form='unformatted', access='stream')

      write (unit) 'REAL:DIM:', 0

      write (unit) 'mu*EB:', im%muStar
      write (unit) 'TcMD:', i%Tc

      if (i%critical) write (unit) 'TcEB:', i%kT * qe / kB

      write (unit) 'INT:status:', im%status

      write (unit) 'REAL:DIM:', 1, im%n

      write (unit) 'iomega:', im%omega
      write (unit) 'Z:', im%Z

      if (i%DOS) write (unit) 'chi:', im%chi

      write (unit) 'Delta:', im%Delta

      write (unit) 'DIM:', 0

      write (unit) 'phiC:', im%phiC

      if (i%measurable) then
         write (unit) 'INT:status[Delta0]:', re%status
         write (unit) 'REAL:Delta0:', re%Delta0
      end if

      if (i%resolution .gt. 0) then
         write (unit) 'DIM:', 1, i%resolution

         write (unit) 'omega:', re%omega

         write (unit) 'Re[Z]:', real(re%Z)
         write (unit) 'Im[Z]:', aimag(re%Z)

         if (i%DOS) then
            write (unit) 'Re[chi]:', real(re%chi)
            write (unit) 'Im[chi]:', aimag(re%chi)
         end if

         write (unit) 'Re[Delta]:', real(re%Delta)
         write (unit) 'Im[Delta]:', aimag(re%Delta)
      end if

      close (unit)
   end subroutine store_data
end module io
