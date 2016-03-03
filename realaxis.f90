module realaxis
   use global
   use intervals
   use pade
   implicit none

contains

   subroutine realize(i)
      type(info), intent(inout) :: i

      integer :: n
      real(dp) :: Delta0

      call coefficients(i%omega, i%Delta)

      allocate(i%omega_(i%resolution))
      allocate(i%Delta_(i%resolution))

      call interval(i%omega_, 0.0_dp, i%upper, lower=.true., upper=.true.)

      do n = 1, i%resolution
         i%Delta_(n) = continuation(i%omega_(n))
      end do

      i%Delta0 = 1
      i%statusDelta0 = -1

      do n = 1, i%limit
         Delta0 = real(continuation(i%Delta0), dp)

         if (abs(i%Delta0 - Delta0) .le. i%tiny) i%statusDelta0 = n

         i%Delta0 = Delta0

         if (n .eq. i%statusDelta0) exit
      end do

      call coefficients(i%omega, i%Z)

      allocate(i%Z_(i%resolution))

      do n = 1, i%resolution
         i%Z_(n) = continuation(i%omega_(n))
      end do
   end subroutine realize
end module realaxis
