module realaxis
   use global
   use intervals
   use pade
   implicit none

contains

   subroutine realize(i, im, re)
      type(universal), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(out) :: re

      integer :: n
      real(dp) :: Delta0

      if (i%measurable .or. i%resolution .gt. 0) then
         call coefficients(im%omega, im%Delta)
      end if

      if (i%measurable) then
         re%Delta0 = 1
         re%status = -1

         do n = 1, i%limit
            Delta0 = real(continuation(re%Delta0), dp)

            if (re%Delta0 .ap. Delta0) re%status = n

            re%Delta0 = Delta0

            if (n .eq. re%status) exit
         end do
      end if

      if (i%resolution .gt. 0) then
         allocate(re%omega(i%resolution))

         call interval(re%omega, 0.0_dp, i%upper, lower=.true., upper=.false.)

         allocate(re%Delta(i%resolution))

         do n = 1, i%resolution
            re%Delta(n) = continuation(re%omega(n))
         end do

         call coefficients(im%omega, im%Z)

         allocate(re%Z(i%resolution))

         do n = 1, i%resolution
            re%Z(n) = continuation(re%omega(n))
         end do

         if (i%DOS) then
            call coefficients(im%omega, im%chi)

            allocate(re%chi(i%resolution))

            do n = 1, i%resolution
               re%chi(n) = continuation(re%omega(n))
            end do
         end if
      end if
   end subroutine realize
end module realaxis
