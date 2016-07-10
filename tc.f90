module tc
   use eliashberg
   use global
   implicit none

contains

   subroutine estimate(i)
      type(universal), intent(inout) :: i

      i%Tc = i%omegaE / 1.2_dp * exp(-1.04_dp * (1 + i%lambda) &
         / (i%lambda - 0.62_dp * i%lambda * i%muStar - i%muStar))
   end subroutine estimate

   subroutine bisection(i, im)
      type(universal), intent(inout) :: i
      type(matsubara), intent(out) :: im

      real(dp), parameter :: ratio = 0.1_dp

      real(dp) :: lower, upper

      i%T = i%Tc

      call solve(i, im)

      if (im%Delta(0) .gt. i%small) then
         do while (im%Delta(0) .gt. i%small)
            lower = i%T

            i%T = i%T * (1 + ratio)

            call solve(i, im)
         end do

         upper = i%T
      else
         do while (im%Delta(0) .le. i%small)
            if (i%T .lt. i%bound) return

            upper = i%T

            i%T = i%T * (1 - ratio)

            call solve(i, im)
         end do

         lower = i%T
      end if

      do
         i%T = (lower + upper) / 2

         call solve(i, im)

         if (im%Delta(0) .gt. i%small) then
            lower = i%T
         else if ((upper - lower) / 2 .gt. i%error) then
            upper = i%T
         else
            exit
         end if
      end do
   end subroutine bisection
end module tc
