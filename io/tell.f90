module io_tell
   use formatting
   use global
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

      character(:), allocatable :: head, body, form ! edit descriptors

      call measure(x%form)

      head = edit('(7Aw)')
      body = edit('(7x)')

      print "('imaginary-axis solution [', I0, ']:', /)", im%status

      if (x%chi) then
         print head, 'omega/eV', 'Z', 'Delta/eV', 'chi/eV'

         do i = 1, x%bands
            print rule(4)

            do n = 0, size(im%omega) - 1
               print body, im%omega(n), im%Z(n, i), im%Delta(n, i), im%chi(n, i)
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

      form = edit('(x)')

      if (x%chi) then
         print "(/, 'initial and final occupancy number:', /)"
         print edit(form), oc%n0, oc%n

         print "(/, 'initial and final chemical potential (eV):', /)"
         print edit(form), oc%mu0, oc%mu
      end if

      print "(/, 'constant Coulomb contribution (eV):', /)"
      print edit(form), im%phiC

      if (x%measurable) then
         print "(/, 'measurable gap (eV):', /)"

         form = edit("(x, ' [', I0, ']')")

         do i = 1, x%bands
            print form, re%Delta0(i), re%status(i)
         end do
      end if

      if (x%resolution .gt. 0) then
         print "(/, 'real-axis solution:', /)"

         if (x%chi) then
            print head, 'omega/eV', 'Re[Z]', 'Im[Z]', &
               'Re[Delta]/eV', 'Im[Delta]/eV', 'Re[chi]', 'Im[chi]'

            do i = 1, x%bands
               print rule(7)

               do n = 1, x%resolution
                  print body, re%omega(n), re%Z(n, i), &
                     re%Delta(n, i), re%chi(n, i)
               end do
            end do
         else
            print head, &
               'omega/eV', 'Re[Z]', 'Im[Z]', 'Re[Delta]/eV', 'Im[Delta]/eV'

            do i = 1, x%bands
               print rule(5)

               do n = 1, x%resolution
                  print body, re%omega(n), re%Z(n, i), re%Delta(n, i)
               end do
            end do
         end if
      end if
   end subroutine tell
end module io_tell
