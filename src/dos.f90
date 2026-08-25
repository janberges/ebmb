! Copyright (C) 2016-2026 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module dos
   use globals
   use eliashberg_self_energy, only: weight_dos => weight
   use pade
   use tools, only: differential
   implicit none

   private
   public :: density_of_states, density_of_states_pade

contains

   subroutine density_of_states(x, re, oc)
      type(parameters), intent(in) :: x
      type(continued), intent(inout) :: re
      type(occupancy), intent(inout) :: oc

      integer :: i, n

      real(dp) :: weight(x%points)
      complex(dp) :: omg
      complex(dp) :: eps(size(x%energy))

      allocate(re%dos(x%points, x%bands))

      do i = 1, x%bands
         !$omp parallel do private(omg, eps)
         do n = 1, x%points
            omg = cmplx(re%omega(n), x%eta, dp) - re%Sigma(n, i) + re%chi(n, i)

            eps(:) = x%energy - oc%mu + re%chi(n, i)

            re%dos(n, i) = -sum(weight_dos(:, i) / pi &
               * aimag((omg + eps) / (omg ** 2 - eps ** 2 - re%phi(n, i) ** 2)))
         end do
         !$omp end parallel do
      end do

      re%dos(:, :) = abs(re%dos)

      call differential(re%omega, weight)

      oc%inspect = sum(weight * sum(re%dos, 2))
   end subroutine density_of_states

   subroutine density_of_states_pade(x, im, re, oc)
      type(parameters), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(inout) :: re
      type(occupancy), intent(inout) :: oc

      integer :: i, n, nP

      real(dp) :: nE, omg, eps(size(x%energy)), weight(x%points)
      complex(dp), allocatable :: green(:)

      nE = x%omegaE / (2.0_dp * pi * kB * x%T)

      nP = min(max(1, ceiling(x%cutoffP * nE - 0.5_dp)), size(im%omega))

      allocate(green(0:nP - 1))
      allocate(re%dos(x%points, x%bands))

      do i = 1, x%bands
         !$omp parallel do private(omg, eps)
         do n = 0, nP - 1
            omg = im%Z(n, i) * im%omega(n)

            eps(:) = x%energy - oc%mu + im%chi(n, i)

            green(n) = -sum(weight_dos(:, i) * cmplx(eps, omg, dp) &
               / (omg ** 2 + eps ** 2 + im%phi(n, i) ** 2))
         end do
         !$omp end parallel do

         call coefficients(im%omega(:nP - 1), green)

         !$omp parallel do
         do n = 1, x%points
            re%dos(n, i) = -aimag(continuation(cmplx(re%omega(n), x%eta, dp))) &
               / pi
         end do
         !$omp end parallel do
      end do

      call differential(re%omega, weight)

      oc%inspect = sum(weight * sum(re%dos, 2))
   end subroutine density_of_states_pade
end module dos
