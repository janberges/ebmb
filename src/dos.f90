! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module dos
   use globals
   use eliashberg_self_energy, only: weight_dos => weight
   use pade
   use tools, only: differential
   implicit none

   private
   public :: density_of_states, density_of_states_stable

contains

   subroutine density_of_states(x, re, oc)
      type(parameters), intent(in) :: x
      type(continued), intent(inout) :: re
      type(occupancy), intent(inout) :: oc

      integer :: i, n

      real(dp) :: weight(x%points)
      complex(dp) :: omg, phi
      complex(dp) :: eps(size(x%energy))

      allocate(re%dos(x%points, x%bands))

      do i = 1, x%bands
         !$omp parallel do private(omg, phi, eps)
         do n = 1, x%points
            omg = re%Z(n, i) * cmplx(re%omega(n), x%eta, dp)
            phi = re%Z(n, i) * re%Delta(n, i)

            eps(:) = x%energy - oc%mu + re%chi(n, i)

            re%dos(n, i) = -sum(weight_dos(:, i) / pi &
               * aimag((omg + eps) / (omg ** 2 - eps ** 2 - phi ** 2)))
         end do
         !$omp end parallel do
      end do

      re%dos(:, :) = abs(re%dos)

      call differential(re%omega, weight)

      oc%inspect = sum(weight * sum(re%dos, 2))
   end subroutine density_of_states

   subroutine density_of_states_stable(x, im, re, oc)
      type(parameters), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(inout) :: re
      type(occupancy), intent(inout) :: oc

      integer :: i, n

      real(dp) :: omg, eps(size(x%energy)), weight(x%points)
      complex(dp) :: green(0:size(im%omega) - 1)

      allocate(re%dos(x%points, x%bands))

      do i = 1, x%bands
         !$omp parallel do private(omg, eps)
         do n = 0, size(im%omega) - 1
            omg = im%Z(n, i) * im%omega(n)

            eps(:) = x%energy - oc%mu + im%chi(n, i)

            green(n) = -sum(weight_dos(:, i) * cmplx(eps, omg, dp) &
               / (omg ** 2 + eps ** 2 + im%phi(n, i) ** 2))
         end do
         !$omp end parallel do

         call coefficients(im%omega, green)

         !$omp parallel do
         do n = 1, x%points
            re%dos(n, i) = -aimag(continuation(cmplx(re%omega(n), x%eta, dp))) &
               / pi
         end do
         !$omp end parallel do
      end do

      call differential(re%omega, weight)

      oc%inspect = sum(weight * sum(re%dos, 2))
   end subroutine density_of_states_stable
end module dos
