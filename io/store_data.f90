module io_store_data
   use global
   implicit none

   private
   public :: store_data

contains

   subroutine store_data(i, im, re)
      type(universal), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re

      integer, parameter :: unit = 11

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

         write (unit) 'REAL:DIM:', 1, size(im%omega)
         write (unit) 'iomega:', im%omega

         if (i%bands .gt. 1) write (unit) 'DIM:', 2, i%bands, size(im%omega)

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
         write (unit) 'muStar:', i%muStar

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
end module io_store_data
