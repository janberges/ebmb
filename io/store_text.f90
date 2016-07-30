module io_store_text
   use global
   implicit none

   private
   public :: store_text

contains

   subroutine store_text(i, im, re)
      type(universal), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      integer, parameter :: unit = 11
      integer :: p, n, width

      character( 20) :: float, count, matrix, head, body
      character(100) :: test

      write (test, "(" // i%edit // ", '|')") pi
      width = index(test, '|') - 1

      write (float, "('(', A, ')')") i%edit
      write (count, "('(I', I0, ')')") width
      write (matrix, "('(', I0, A, ')')") i%bands, i%edit
      write (head, "('(7A', I0, ')')") width
      write (body, "('(7', A, ')')") i%edit

      open (unit, file=i%name // '.out', action='write', status='replace')

      write (unit, "('McMillan and Dynes'' critical temperature (K):', /)")
      write (unit, float) i%TcMD

      if (i%critical) then
         write (unit, "(/, 'Eliashberg''s critical temperature (K):', /)")
         write (unit, float) i%TcEB
      else
         write (unit, "(/, 'imaginary-axis solution [', I0, ']:', /)") im%status

         if (i%chi) then
            write (unit, head) 'omega/eV', 'Z', 'Delta/eV', 'chi/eV'

            do p = 1, i%bands
               call rule(4)

               do n = 0, size(im%omega) - 1
                  write (unit, body) &
                     im%omega(n), im%Z(n, p), im%Delta(n, p), im%chi(n, p)
               end do
            end do
         else
            write (unit, head) 'omega/eV', 'Z', 'Delta/eV'

            do p = 1, i%bands
               call rule(3)

               do n = 0, size(im%omega) - 1
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

            if (i%chi) then
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
            write (unit, float) i%error

            write (unit, "(/, 'lower bound of critical temperature (K)', /)")
            write (unit, float) i%bound

            write (unit, "(/, 'maximum gap at critical temperature (eV)', /)")
            write (unit, float) i%small
         else
            write (unit, "(/, 'temperature (K):', /)")
            write (unit, float) i%T
         end if

         write (unit, "(/, 'number of electronic bands:', /)")
         write (unit, count) i%bands

         write (unit, "(/, 'Einstein frequency (eV):', /)")
         write (unit, float) i%omegaE

         write (unit, "(/, 'electron-phonon coupling:', /)")
         write (unit, matrix) i%lambda

         write (unit, "(/, 'Coulomb pseudo-potential:', /)")
         write (unit, matrix) i%muStar

         write (unit, "(/, 'overall cutoff (Einstein frequency):', /)")
         write (unit, float) i%upper

         if (i%lower .lt. i%upper) then
            write (unit, "(/, 'Coulomb cutoff (Einstein frequency):', /)")
            write (unit, float) i%lower
         end if

         write (unit, "(/, 'maximum number of fixed-point steps:', /)")
         write (unit, count) i%limit

         write (unit, "(/, 'negligible float difference (a.u.):', /)")
         write (unit, float) negligible_difference

         if (i%chi) then
            write (unit, "(/, 'density of Bloch states:', /)")
            write (unit, head) 'E/eV', 'DOS/a.u.'

            call rule(i%bands + 1)

            do n = 1, size(i%energy)
               write (unit, body) i%energy(n), i%dos(n, :)
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
end module io_store_text
