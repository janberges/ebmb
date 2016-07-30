module io_store_data
   use global
   implicit none

   private
   public :: store_data

contains

   subroutine store_data(x, im, re)
      type(universal), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      integer, parameter :: unit = 11

      open (unit, file=x%name // '.dat', action='write', status='replace', &
         form='unformatted', access='stream')

      write (unit) 'REAL:DIM:', 0
      write (unit) 'TcMD:', x%TcMD

      if (x%critical) then
         if (x%bands .gt. 1) write (unit) 'DIM:', 1, x%bands

         write (unit) 'TcEB:', x%TcEB
      else
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
      end if

      if (x%standalone) then
         write (unit) 'DIM:', 0

         if (x%critical) then
            write (unit) 'error:', x%error
            write (unit) 'bound:', x%bound
            write (unit) 'small:', x%small
         else
            write (unit) 'T:', x%T
         end if

         write (unit) 'INT:bands:', x%bands
         write (unit) 'REAL:omegaE:', x%omegaE

         if (x%bands .gt. 1) write (unit) 'DIM:', 2, x%bands, x%bands

         write (unit) 'lambda:', x%lambda
         write (unit) 'muStar:', x%muStar

         write (unit) 'DIM:', 0
         write (unit) 'upper:', x%upper

         if (x%lower .lt. x%upper) write (unit) 'lower:', x%lower

         write (unit) 'INT:limit:', x%limit
         write (unit) 'REAL:epsilon:', negligible_difference

         if (x%chi) then
            write (unit) 'DIM:', 1, size(x%energy)
            write (unit) 'energy:', x%energy

            if (x%bands .gt. 1) write (unit) 'DIM:', 2, x%bands, size(x%energy)

            write (unit) 'density:', x%dos
         end if
      end if

      close (unit)
   end subroutine store_data
end module io_store_data
