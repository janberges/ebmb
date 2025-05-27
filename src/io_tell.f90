! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module io_tell
   use formatting
   use globals
   implicit none

   private
   public :: tell

contains

   subroutine tell(x, im, re, oc)
      type(parameters), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re
      type(occupancy), intent(in) :: oc

      integer :: i, n ! band and Matsubara indices

      character(:), allocatable :: head, body, more ! edit descriptors

      call measure(x%flomat)

      head = edit('(99Aw)') ! 99 replaces Fortran-2008 unlimited format control.
      body = edit('(99x)')

      print "('imaginary-axis solution [', I0, ']:', /)", im%steps

      if (x%Sigma) then
         print head, 'omega/eV', 'Re[Sigma]/eV', 'Im[Sigma]/eV', 'Delta/eV'

         do i = 1, x%bands
            print rule(4)

            do n = 0, size(im%omega) - 1
               print body, im%omega(n), im%Sigma(n, i), im%Delta(n, i)
            end do
         end do
      else if (x%ldos) then
         print head, 'omega/eV', 'Z', 'chi/eV', 'Delta/eV'

         do i = 1, x%bands
            print rule(4)

            do n = 0, size(im%omega) - 1
               print body, im%omega(n), im%Z(n, i), im%chi(n, i), im%Delta(n, i)
            end do
         end do
      else
         print head, 'omega/eV', 'Z', 'Delta/eV'

         do i = 1, x%bands
            print rule(3)

            do n = 0, size(im%omega) - 1
               print body, im%omega(n), im%Z(n, i), im%Delta(n, i)
            end do
         end do
      end if

      more = edit('(x)')

      if (x%ldos) then
         print "(/, 'integral of density of states:', /)"
         print more, oc%states

         if (x%points .gt. 0) then
            print "(/, 'integral of spectral function:', /)"
            print more, oc%inspect
         end if

         print "(/, 'initial and final occupancy number:', /)"
         print more, oc%n0, oc%n

         print "(/, 'initial and final chemical potential (eV):', /)"
         print more, oc%mu0, oc%mu
      end if

      if (x%la2F) then
         print "(/, 'effective electron-phonon coupling:', /)"
         do i = 1, x%bands
            print body, x%lambda(:, i)
         end do

         print "(/, 'effective phonon frequency (eV):', /)"
         print more, x%omegaE

         print "(/, 'logarithmic average phonon frequency (eV):', /)"
         print more, x%omegaLog

         print "(/, 'second-moment average phonon frequency (eV):', /)"
         print more, x%omega2nd
      end if

      if (x%ldos) then
         print "(/, 'Coulomb part of energy shift (eV):', /)"
         print more, im%chiC
      end if

      print "(/, 'Coulomb part of order parameter (eV):', /)"
      print more, im%phiC

      if (x%measurable) then
         print "(/, 'measurable gap (eV):', /)"

         more = edit("(x, ' [', I0, ']')")

         do i = 1, x%bands
            print more, re%Delta0(i), re%steps(i)
         end do
      end if

      if (x%points .gt. 0) then
         print "(/, 'real-axis solution:', /)"

         if (x%Sigma) then
            if (x%ldos) then
               print head, 'omega/eV', 'Re[Sigma]', 'Im[Sigma]', &
                  'Re[Delta]/eV', 'Im[Delta]/eV', 'DOS/(1/eV)'

               do i = 1, x%bands
                  print rule(6)

                  do n = 1, x%points
                     print body, re%omega(n), re%Sigma(n, i), re%Delta(n, i), &
                        re%dos(n, i)
                  end do
               end do
            else
               print head, 'omega/eV', 'Re[Sigma]', 'Im[Sigma]', &
                  'Re[Delta]/eV', 'Im[Delta]/eV'

               do i = 1, x%bands
                  print rule(5)

                  do n = 1, x%points
                     print body, re%omega(n), re%Sigma(n, i), re%Delta(n, i)
                  end do
               end do
            end if
         else
            if (x%ldos) then
               print head, 'omega/eV', 'Re[Z]', 'Im[Z]', &
                  'Re[chi]', 'Im[chi]', 'Re[Delta]/eV', 'Im[Delta]/eV', &
                  'DOS/(1/eV)'

               do i = 1, x%bands
                  print rule(8)

                  do n = 1, x%points
                     print body, re%omega(n), re%Z(n, i), re%chi(n, i), &
                        re%Delta(n, i), re%dos(n, i)
                  end do
               end do
            else
               print head, &
                  'omega/eV', 'Re[Z]', 'Im[Z]', 'Re[Delta]/eV', 'Im[Delta]/eV'

               do i = 1, x%bands
                  print rule(5)

                  do n = 1, x%points
                     print body, re%omega(n), re%Z(n, i), re%Delta(n, i)
                  end do
               end do
            end if
         end if
      end if
   end subroutine tell
end module io_tell
