! Copyright (C) 2016-2024 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module dos
   use global
   use eliashberg_self_energy, only: weight
   use pade
   implicit none

   private
   public :: density_of_states, density_of_states_stable

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

   subroutine density_of_states_stable(x, im, re, oc)
      type(parameters), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(inout) :: re
      type(occupancy), intent(in) :: oc

      integer :: i, n

      real(dp) :: omg, eps(size(x%energy))
      complex(dp) :: green(size(im%omega))

      allocate(re%dos(x%resolution, x%bands))

      do i = 1, x%bands
         do n = 1, size(im%omega)
            omg = im%Z(n, i) * im%omega(n)

            eps(:) = x%energy - oc%mu + im%chi(n, i)

            green(n) = -sum(weight(:, i) * cmplx(eps, omg, dp) &
               / (omg ** 2 + eps ** 2 + im%phi(n, i) ** 2))
         end do

         call coefficients(im%omega, green)

         re%dos(:, i) = -aimag(continuation(cmplx(re%omega, x%eta, dp))) / pi
      end do
   end subroutine density_of_states_stable
end module dos
