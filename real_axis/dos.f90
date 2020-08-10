module dos
   use global
   use eliashberg_self_energy, only: weight
   implicit none

   private
   public :: density_of_states

contains

   subroutine density_of_states(x, re, oc)
      type(parameters), intent(in) :: x
      type(continued), intent(inout) :: re
      type(occupancy), intent(in) :: oc

      integer :: i, n

      complex(dp) :: omg, phi
      complex(dp) :: eps(size(x%energy))

      allocate(re%dos(x%resolution, x%bands))

      do i = 1, x%bands
         do n = 1, x%resolution
            omg = re%Z(n, i) * cmplx(re%omega(n), x%eta, dp)
            phi = re%Z(n, i) * re%Delta(n, i)

            eps(:) = x%energy - oc%mu + re%chi(n, i)

            re%dos(n, i) = -sum(weight(:, i) / pi &
               * aimag((omg + eps) / (omg ** 2 - eps ** 2 - phi ** 2)))
         end do
      end do

      re%dos(:, :) = abs(re%dos)
   end subroutine density_of_states
end module dos
