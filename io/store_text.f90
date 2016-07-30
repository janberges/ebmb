module io_store_text
   use global
   implicit none

   private
   public :: store_text

contains

   subroutine store_text(x, im, re)
      type(universal), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      integer, parameter :: unit = 11
      integer :: i, n, width

      character( 20) :: float, count, matrix, head, body
      character(100) :: test

      write (test, "(" // x%edit // ", '|')") pi
      width = index(test, '|') - 1

      write (float, "('(', A, ')')") x%edit
      write (count, "('(I', I0, ')')") width
      write (matrix, "('(', I0, A, ')')") x%bands, x%edit
      write (head, "('(7A', I0, ')')") width
      write (body, "('(7', A, ')')") x%edit

      open (unit, file=x%name // '.out', action='write', status='replace')

      write (unit, "('McMillan and Dynes'' critical temperature (K):', /)")
      write (unit, float) x%TcMD

      if (x%critical) then
         write (unit, "(/, 'Eliashberg''s critical temperature (K):', /)")
         write (unit, float) x%TcEB
      else
         write (unit, "(/, 'imaginary-axis solution [', I0, ']:', /)") im%status

         if (x%chi) then
            write (unit, head) 'omega/eV', 'Z', 'Delta/eV', 'chi/eV'

            do i = 1, x%bands
               call rule(4)

               do n = 0, size(im%omega) - 1
                  write (unit, body) &
                     im%omega(n), im%Z(n, i), im%Delta(n, i), im%chi(n, i)
               end do
            end do
         else
            write (unit, head) 'omega/eV', 'Z', 'Delta/eV'

            do i = 1, x%bands
               call rule(3)

               do n = 0, size(im%omega) - 1
                  write (unit, body) &
                     im%omega(n), im%Z(n, i), im%Delta(n, i)
               end do
            end do
         end if

         write (unit, "(/, 'constant Coulomb contribution (eV):', /)")
         write (unit, float) im%phiC

         if (x%measurable) then
            write (unit, "(/, 'measurable gap (eV):', /)")

            do i = 1, x%bands
               write (unit, float, advance='no') re%Delta0(i)
               write (unit, "(' [', I0, ']')") re%status(i)
            end do
         end if

         if (x%resolution .gt. 0) then
            write (unit, "(/, 'real-axis solution:', /)")

            if (x%chi) then
               write (unit, head) 'omega/eV', 'Re[Z]', 'Im[Z]', &
                  'Re[Delta]/eV', 'Im[Delta]/eV', 'Re[chi]', 'Im[chi]'

               do i = 1, x%bands
                  call rule(7)

                  do n = 1, x%resolution
                     write (unit, body) &
                        re%omega(n), re%Z(n, i), re%Delta(n, i), re%chi(n, i)
                  end do
               end do
            else
               write (unit, head) &
                  'omega/eV', 'Re[Z]', 'Im[Z]', 'Re[Delta]/eV', 'Im[Delta]/eV'

               do i = 1, x%bands
                  call rule(5)

                  do n = 1, x%resolution
                     write (unit, body) &
                        re%omega(n), re%Z(n, i), re%Delta(n, i)
                  end do
               end do
            end if
         end if
      end if

      if (x%standalone) then
         if (x%critical) then
            write (unit, "(/, 'valid error of critical temperature (K)', /)")
            write (unit, float) x%error

            write (unit, "(/, 'lower bound of critical temperature (K)', /)")
            write (unit, float) x%bound

            write (unit, "(/, 'maximum gap at critical temperature (eV)', /)")
            write (unit, float) x%small
         else
            write (unit, "(/, 'temperature (K):', /)")
            write (unit, float) x%T
         end if

         write (unit, "(/, 'number of electronic bands:', /)")
         write (unit, count) x%bands

         write (unit, "(/, 'Einstein frequency (eV):', /)")
         write (unit, float) x%omegaE

         write (unit, "(/, 'electron-phonon coupling:', /)")
         write (unit, matrix) x%lambda

         write (unit, "(/, 'Coulomb pseudo-potential:', /)")
         write (unit, matrix) x%muStar

         write (unit, "(/, 'overall cutoff (Einstein frequency):', /)")
         write (unit, float) x%upper

         if (x%lower .lt. x%upper) then
            write (unit, "(/, 'Coulomb cutoff (Einstein frequency):', /)")
            write (unit, float) x%lower
         end if

         write (unit, "(/, 'maximum number of fixed-point steps:', /)")
         write (unit, count) x%limit

         write (unit, "(/, 'negligible float difference (a.u.):', /)")
         write (unit, float) negligible_difference

         if (x%chi) then
            write (unit, "(/, 'density of Bloch states:', /)")
            write (unit, head) 'E/eV', 'DOS/a.u.'

            call rule(x%bands + 1)

            do n = 1, size(x%energy)
               write (unit, body) x%energy(n), x%dos(n, :)
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
