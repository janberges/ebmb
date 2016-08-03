module io_store
   use global
   implicit none

   private
   public :: store

contains

   subroutine store(x, im, re)
      type(universal), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      open (unit, &
         file=x%file, action='write', status='replace', access='stream')

      write (unit) 'INT:DIM:', 0
      write (unit) 'status:', im%status

      write (unit) 'REAL:DIM:', 1, size(im%omega)
      write (unit) 'iomega:', im%omega

      if (x%bands .gt. 1) write (unit) 'DIM:', 2, x%bands, size(im%omega)

      write (unit) 'Z:', im%Z
      write (unit) 'Delta:', im%Delta

      if (x%chi) write (unit) 'chi:', im%chi

      write (unit) 'DIM:'

      if (x%bands .gt. 1) then
         write (unit) 1, x%bands
      else
         write (unit) 0
      end if

      write (unit) 'phiC:', im%phiC

      if (x%measurable) then
         write (unit) 'INT:status0:', re%status
         write (unit) 'REAL:Delta0:', re%Delta0
      end if

      if (x%resolution .gt. 0) then
         write (unit) 'DIM:', 1, x%resolution

         write (unit) 'omega:', re%omega

         if (x%bands .gt. 1) write (unit) 'DIM:', 2, x%bands, x%resolution

         write (unit) 'Re[Z]:', real(re%Z)
         write (unit) 'Im[Z]:', aimag(re%Z)

         write (unit) 'Re[Delta]:', real(re%Delta)
         write (unit) 'Im[Delta]:', aimag(re%Delta)

         if (x%chi) then
            write (unit) 'Re[chi]:', real(re%chi)
            write (unit) 'Im[chi]:', aimag(re%chi)
         end if
      end if

      close (unit)
   end subroutine store
end module io_store
