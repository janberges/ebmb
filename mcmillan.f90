module mcmillan
   use global
   implicit none

contains

   subroutine estimate(i)
      type(info), intent(inout) :: i

      i%Tc = qe / kB * i%omegaE / 1.2_dp * exp(-1.04_dp * (1 + i%lambda) &
         / (i%lambda - 0.62_dp * i%lambda * i%muStar - i%muStar))
   end subroutine estimate
end module mcmillan
