module tc
   use eliashberg
   use global
   implicit none

contains

   subroutine estimate(i)
      type(info), intent(inout) :: i

      i%Tc = qe / kB * i%omegaE / 1.2_dp * exp(-1.04_dp * (1 + i%lambda) &
         / (i%lambda - 0.62_dp * i%lambda * i%muStar - i%muStar))
   end subroutine estimate

   subroutine bisection(i)
      type(info), intent(inout) :: i

      real(dp), parameter :: ratio = 0.1_dp

      real(dp) :: lower, upper

      i%kT = i%Tc * kB / qe

      call solve(i)

      if (i%Delta(0) .gt. i%small) then
         lower = i%kT

         do while (i%Delta(0) .gt. i%small)
            i%kT = i%kT * (1 + ratio)
            call solve(i)
         end do

         upper = i%kT
      else
         upper = i%kT

         do while (i%Delta(0) .le. i%small)
            i%kT = i%kT * (1 - ratio)
            call solve(i)
         end do

         lower = i%kT
      end if

      do
         i%kT = (lower + upper) / 2

         if ((upper - lower) / 2 .le. i%error) exit

         call solve(i)

         if (i%Delta(0) .gt. i%small) then
            lower = i%kT
         else
            upper = i%kT
         end if
      end do
   end subroutine bisection
end module tc
