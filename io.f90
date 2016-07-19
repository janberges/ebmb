module io
   use filenames
   use global
   use integration
   implicit none

   private
   public :: load, store

   integer :: p, n, m
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

      i%T = kB * i%T ! (eV)

      read (unit, *) i%error ! valid error of critical temperature (K)
      read (unit, *) i%bound ! lower bound of critical temperature (K)
      read (unit, *) i%small ! maximum gap at critical temperature (eV)

      i%error = kB * i%error ! (eV)
      i%bound = kB * i%bound ! (eV)

      read (unit, *) i%bands ! number of electronic bands

      read (unit, *) i%omegaE ! Einstein frequency (eV)

      allocate(i%lambda(i%bands, i%bands))
      allocate(i%muStar(i%bands, i%bands))

      do p = 1, i%bands
         read (unit, *) i%lambda(:, p) ! Electron-phonon coupling
      end do

      do p = 1, i%bands
         read (unit, *) i%muStar(:, p) ! Coulomb pseudo-potential
      end do

      read (unit, *) DOSfile ! file with density of states

      i%DOS = DOSfile .ne. 'none' ! consider full density of states?

      read (unit, *) i%upper ! general cutoff (Einstein frequency)
      read (unit, *) i%lower ! Coulomb cutoff (Einstein frequency)

      if (i%lower .lt. 0) i%lower = i%upper

      i%upper = i%upper * i%omegaE ! (eV)
      i%lower = i%lower * i%omegaE ! (eV)

      read (unit, *) i%limit ! maximum number of fixed-point steps

      read (unit, *) i%measurable ! find measurable gap?
      read (unit, *) i%resolution ! real axis resolution

      if (i%critical) then
         i%measurable = .false.
         i%resolution = 0
      end if

      read (unit, *) i%form ! output format
      read (unit, *) i%edit ! number format

      read (unit, *) i%standalone ! include parameters in output file?
      read (unit, *) i%rescale ! rescale Coulomb pseudo-potential?

      read (unit, *) negligible_difference ! negligible float difference

      close (unit)

      if (i%DOS) then
         open (unit, file=DOSfile, action='read', status='old')

         read (unit, *) n ! density-of-states resolution

         allocate(i%energy(n)) ! free-electron energy (eV)
         allocate(i%density(n, i%bands)) ! density of states (a.u.)
         allocate(i%weight(n, i%bands)) ! integration weight (eV)

         do m = 1, n
            read(unit, *) i%energy(m), i%density(m, :)
         end do

         close (unit)

         call differential(i%energy, i%weight(:, 1))

         do p = 2, i%bands
            i%weight(:, p) = i%weight(:, 1)
         end do

         n = minloc(abs(i%energy), 1) ! index of Fermi level

         do p = 1, i%bands
            i%weight(:, p) = i%weight(:, p) * i%density(:, p) / i%density(n, p)
         end do
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

      character(20) :: float, count, matrix, head, body

      character(100) :: test
      integer :: width

      write (test, "(" // i%edit // ", '|')") pi
      width = index(test, '|') - 1

      write (float, "('(', A, ')')") i%edit
      write (count, "('(I', I0, ')')") width
      write (matrix, "('(', I0, A, ')')") i%bands, i%edit
      write (head, "('(7A', I0, ')')") width
      write (body, "('(7', A, ')')") i%edit

      open (unit, file=i%name // '.out', action='write', status='replace')

      write (unit, "('McMillan and Dynes'' critical temperature (K):', /)")
      write (unit, float) i%TcMD / kB

      if (i%critical) then
         write (unit, "(/, 'Eliashberg''s critical temperature (K):', /)")
         write (unit, float) i%TcEB / kB
      else
         write (unit, "(/, 'imaginary-axis solution [', I0, ']:', /)") im%status

         if (i%DOS) then
            write (unit, head) 'omega/eV', 'Z', 'Delta/eV', 'chi/eV'

            do p = 1, i%bands
               call rule(4)

               do n = 0, im%u - 1
                  write (unit, body) &
                     im%omega(n), im%Z(n, p), im%Delta(n, p), im%chi(n, p)
               end do
            end do
         else
            write (unit, head) 'omega/eV', 'Z', 'Delta/eV'

            do p = 1, i%bands
               call rule(3)

               do n = 0, im%u - 1
                  write (unit, body) &
                     im%omega(n), im%Z(n, p), im%Delta(n, p)
               end do
            end do
         end if

         write (unit, "(/, 'constant Coulomb contribution (eV):', /)")
         write (unit, float) im%phiC

         if (i%measurable) then
            write (unit, "(/, 'measurable gap (eV):', /)")

            do p = 1, i%bands
               write (unit, float, advance='no') re%Delta0(p)
               write (unit, "(' [', I0, ']')") re%status(p)
            end do
         end if

         if (i%resolution .gt. 0) then
            write (unit, "(/, 'real-axis solution:', /)")

            if (i%DOS) then
               write (unit, head) 'omega/eV', 'Re[Z]', 'Im[Z]', &
                  'Re[Delta]/eV', 'Im[Delta]/eV', 'Re[chi]', 'Im[chi]'

               do p = 1, i%bands
                  call rule(7)

                  do n = 1, i%resolution
                     write (unit, body) &
                        re%omega(n), re%Z(n, p), re%Delta(n, p), re%chi(n, p)
                  end do
               end do
            else
               write (unit, head) &
                  'omega/eV', 'Re[Z]', 'Im[Z]', 'Re[Delta]/eV', 'Im[Delta]/eV'

               do p = 1, i%bands
                  call rule(5)

                  do n = 1, i%resolution
                     write (unit, body) &
                        re%omega(n), re%Z(n, p), re%Delta(n, p)
                  end do
               end do
            end if
         end if
      end if

      if (i%standalone) then
         if (i%critical) then
            write (unit, "(/, 'valid error of critical temperature (K)', /)")
            write (unit, float) i%error / kB

            write (unit, "(/, 'lower bound of critical temperature (K)', /)")
            write (unit, float) i%bound / kB

            write (unit, "(/, 'maximum gap at critical temperature (eV)', /)")
            write (unit, float) i%small
         else
            write (unit, "(/, 'temperature (K):', /)")
            write (unit, float) i%T / kB
         end if

         write (unit, "(/, 'number of electronic bands:', /)")
         write (unit, count) i%bands

         write (unit, "(/, 'Einstein frequency (eV):', /)")
         write (unit, float) i%omegaE

         write (unit, "(/, 'electron-phonon coupling:', /)")
         write (unit, matrix) i%lambda

         write (unit, "(/, 'Coulomb pseudo-potential:', /)")
         write (unit, matrix) i%muStar
      end if

      if (.not. i%critical) then
         write (unit, "(/, 'rescaled Coulomb pseudo-potential:', /)")
         write (unit, matrix) im%muStar
      end if

      if (i%standalone) then
         write (unit, "(/, 'overall cutoff (eV):', /)")
         write (unit, float) i%upper

         if (i%lower .lt. i%upper) then
            write (unit, "(/, 'Coulomb cutoff (eV):', /)")
            write (unit, float) i%lower
         end if

         write (unit, "(/, 'maximum number of fixed-point steps:', /)")
         write (unit, count) i%limit

         write (unit, "(/, 'negligible float difference (a.u.):', /)")
         write (unit, float) negligible_difference

         if (i%DOS) then
            write (unit, "(/, 'density of Bloch states:', /)")
            write (unit, head) 'E/eV', 'DOS/a.u.'

            call rule(i%bands + 1)

            do n = 1, size(i%energy)
               write (unit, body) i%energy(n), i%density(n, :)
            end do
         end if
      end if

      close (unit)

   contains

      subroutine rule(n)
         integer, intent(in) :: n

         write (unit, '(A)') repeat('~', n * width)
      end subroutine rule

   end subroutine store_text

   subroutine store_data(i, im, re)
      type(universal), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      open (unit, file=i%name // '.dat', action='write', status='replace', &
         form='unformatted', access='stream')

      write (unit) 'REAL:DIM:', 0
      write (unit) 'TcMD:', i%TcMD / kB

      if (i%critical) then
         if (i%bands .gt. 1) write (unit) 'DIM:', 1, i%bands

         write (unit) 'TcEB:', i%TcEB / kB
      else
         write (unit) 'INT:DIM:', 0
         write (unit) 'status:', im%status

         write (unit) 'REAL:DIM:', 1, im%u
         write (unit) 'iomega:', im%omega

         if (i%bands .gt. 1) write (unit) 'DIM:', 2, i%bands, im%u

         write (unit) 'Z:', im%Z
         write (unit) 'Delta:', im%Delta

         if (i%DOS) write (unit) 'chi:', im%chi

         write (unit) 'DIM:'

         if (i%bands .gt. 1) then
            write (unit) 1, i%bands
         else
            write (unit) 0
         end if

         write (unit) 'phiC:', im%phiC

         if (i%measurable) then
            write (unit) 'INT:status0:', re%status
            write (unit) 'REAL:Delta0:', re%Delta0
         end if

         if (i%resolution .gt. 0) then
            write (unit) 'DIM:', 1, i%resolution

            write (unit) 'omega:', re%omega

            if (i%bands .gt. 1) write (unit) 'DIM:', 2, i%bands, i%resolution

            write (unit) 'Re[Z]:', real(re%Z)
            write (unit) 'Im[Z]:', aimag(re%Z)

            write (unit) 'Re[Delta]:', real(re%Delta)
            write (unit) 'Im[Delta]:', aimag(re%Delta)

            if (i%DOS) then
               write (unit) 'Re[chi]:', real(re%chi)
               write (unit) 'Im[chi]:', aimag(re%chi)
            end if
         end if
      end if

      if (i%standalone) then
         write (unit) 'DIM:', 0

         if (i%critical) then
            write (unit) 'error:', i%error / kB
            write (unit) 'bound:', i%bound / kB
            write (unit) 'small:', i%small
         else
            write (unit) 'T:', i%T / kB
         end if

         write (unit) 'INT:bands:', i%bands
         write (unit) 'REAL:omegaE:', i%omegaE

         if (i%bands .gt. 1) write (unit) 'DIM:', 2, i%bands, i%bands

         write (unit) 'lambda:', i%lambda
         write (unit) 'mu*MD:', i%muStar
      end if

      if (.not. i%critical) write (unit) 'mu*EB:', im%muStar

      if (i%standalone) then
         write (unit) 'DIM:', 0
         write (unit) 'upper:', i%upper

         if (i%lower .lt. i%upper) write (unit) 'lower:', i%lower

         write (unit) 'INT:limit:', i%limit
         write (unit) 'REAL:epsilon:', negligible_difference

         if (i%DOS) then
            write (unit) 'DIM:', 1, size(i%energy)
            write (unit) 'energy:', i%energy

            if (i%bands .gt. 1) write (unit) 'DIM:', 2, i%bands, size(i%energy)

            write (unit) 'density:', i%density
         end if
      end if

      close (unit)
   end subroutine store_data
end module io
