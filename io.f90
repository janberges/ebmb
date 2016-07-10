module io
   use filenames
   use global
   use integration
   implicit none

   private
   public :: load, store

   real(dp), parameter :: k = 8.61733e-05_dp ! Boltzmann constant (eV/K)

   integer :: n, m
   integer, parameter :: unit = 11

contains

   subroutine load(file, i)
      character(*), intent(in) :: file
      type(universal), intent(out) :: i

      character(50) :: DOSfile

      i%name = stem(file)

      open (unit, file=file, action='read', status='old')

      read (unit, *) i%T ! temperature (K)

      i%critical = i%T .lt. 0 ! find critical temperature?

      i%T = k * i%T ! (eV)

      read (unit, *) i%error ! valid error of critical temperature (K)
      read (unit, *) i%bound ! lower bound of critical temperature (K)
      read (unit, *) i%small ! negligible gap (eV)

      i%error = k * i%error ! (eV)
      i%bound = k * i%bound ! (eV)

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
      read (unit, *) i%standalone ! include parameters in output file?

      read (unit, *) i%rescale ! rescale Coulomb pseudo-potential?

      read (unit, *) negligible_difference ! negligible float difference

      close (unit)

      if (i%DOS) then
         open (unit, file=DOSfile, action='read', status='old')

         read (unit, *) n ! density-of-states resolution

         allocate(i%energy(n)) ! free-electron energy (eV)
         allocate(i%density(n)) ! density of states (a.u.)
         allocate(i%weight(n)) ! integration weight (eV)

         do m = 1, n
            read(unit, *) i%energy(m), i%density(m)
         end do

         close (unit)

         call differential(i%energy, i%weight)

         n = minloc(abs(i%energy), 1) ! index of Fermi level

         i%weight = i%weight * i%density / i%density(n)
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
      write (unit, "(/, ES23.14E3, ' K')") i%Tc / k

      if (i%critical) then
         write (unit, "(/, 'Eliashberg''s critical temperature:')")
         write (unit, "(/, ES23.14E3, ' K')") i%T / k
      end if

      write (unit, "(/, 'Imaginary-axis solution (', I0, '):')") im%status

      if (i%DOS) then
         write (unit, '(/, 4A23)') 'omega/eV', 'Z', 'chi/eV', 'Delta/eV'

         do n = 0, im%u - 1
            write (unit, '(4ES23.14E3)') &
               im%omega(n), im%Z(n), im%chi(n), im%Delta(n)
         end do
      else
         write (unit, '(/, 3A23)') 'omega/eV', 'Z', 'Delta/eV'

         do n = 0, im%u - 1
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

      if (i%standalone) then
         write (unit, "(/, 'Parameters:', /)")

         write (unit, "(F9.3, 2X, A, /)") i%T / k, 'temperature (K)'

         if (i%critical) then
            write (unit, '(ES9.1E3, 2X, A)') i%error / k, &
               'valid error of critical temperature (K)'

            write (unit, '(ES9.1E3, 2X, A)') i%bound / k, &
               'lower bound of critical temperature (K)'

            write (unit, '(ES9.1E3, 2X, A, /)') i%small, 'negligible gap (eV)'
         end if

         write (unit, '(F9.3, 2X, A)') i%omegaE, 'Einstein frequency (eV)'
         write (unit, '(F9.3, 2X, A)') i%lambda, 'electron-phonon coupling'
         write (unit, '(F9.3, 2X, A, /)') i%muStar, 'Coulomb pseudo-potential'

         if (im%l .lt. im%u) write (unit, '(I9, 2X, A)') im%l, &
            'index of Coulomb cutoff frequency'

         write (unit, '(I9, 2X, A, /)') i%limit, &
            'maximum number of fixed-point steps'

         write (unit, '(ES9.1E3, 2X, A)') negligible_difference, &
            'negligible float difference (a.u.)'

         if (i%DOS) then
            write (unit, "(/, 'Density of Bloch states:')")
            write (unit, '(/, 2A23)') 'E/eV', 'DOS/a.u.'

            do n = 1, size(i%energy)
               write (unit, '(2ES23.14E3)') i%energy(n), i%density(n)
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
      write (unit) 'TcMD:', i%Tc / k

      if (i%critical) write (unit) 'TcEB:', i%T / k

      write (unit) 'INT:status:', im%status

      write (unit) 'REAL:DIM:', 1, im%u

      write (unit) 'iomega:', im%omega
      write (unit) 'Z:', im%Z

      if (i%DOS) write (unit) 'chi:', im%chi

      write (unit) 'Delta:', im%Delta

      write (unit) 'DIM:', 0

      write (unit) 'phiC:', im%phiC

      if (i%measurable) then
         write (unit) 'INT:status0:', re%status
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

      if (i%standalone) then
         write (unit) 'DIM:', 0

         write (unit) 'T:', i%T / k

         if (i%critical) then
            write (unit) 'error:', i%error / k
            write (unit) 'bound:', i%bound / k
            write (unit) 'small:', i%small
         end if

         write (unit) 'omegaE:', i%omegaE
         write (unit) 'lambda:', i%lambda
         write (unit) 'mu*MD:', i%muStar

         write (unit) 'INT:'

         if (im%l .lt. im%u) write (unit) 'cutoff:', im%l

         write (unit) 'limit:', i%limit

         write (unit) 'REAL:epsilon:', negligible_difference

         if (i%DOS) then
            write (unit) 'DIM:', 1, size(i%energy)

            write (unit) 'energy:', i%energy
            write (unit) 'density:', i%density
         end if
      end if

      close (unit)
   end subroutine store_data
end module io
