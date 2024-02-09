! Copyright (C) 2016-2024 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module eliashberg_spectral_function
   use global
   use tools, only: differential
   implicit none

   private
   public :: integrate_a2F, lambda_from_a2F

   logical :: initial = .true.
   integer :: i, j

   real(dp), allocatable :: weight(:, :, :)

contains

   subroutine integrate_a2F(x)
      type(parameters), intent(inout) :: x

      call lambda_from_a2F(x, x%lambda, 0)

      ! matrix-element average best way to obtain scalar frequency?

      x%omegaE = 0.0_dp

      do i = 1, x%bands
         do j = 1, x%bands
            x%omegaE = x%omegaE + exp(2 / x%lambda(j, i) &
               * sum(weight(:, j, i) * log(x%omega) / x%omega))
         end do
      end do

      x%omegaE = x%omegaE / x%bands ** 2
   end subroutine integrate_a2F

   subroutine lambda_from_a2F(x, lambda, n)
      type(parameters), intent(in) :: x
      real(dp), intent(out) :: lambda(:, :)
      integer, intent(in) :: n

      if (initial) call initialize(x)

      do i = 1, x%bands
         do j = 1, x%bands
            lambda(j, i) = 2 * sum(weight(:, j, i) * x%omega &
               / (x%omega ** 2 + (2 * pi * n * kB * x%T) ** 2))
         end do
      end do
   end subroutine lambda_from_a2F

   subroutine initialize(x)
      type(parameters), intent(in) :: x

      initial = .false.

      if (allocated(weight)) deallocate(weight)
      allocate(weight(size(x%omega), x%bands, x%bands))

      call differential(x%omega, weight(:, 1, 1))

      do i = 1, x%bands
         do j = 1, x%bands
            if (i .ne. 1 .or. j .ne. 1) weight(:, j, i) = weight(:, 1, 1)
         end do
      end do

      weight(:, :, :) = weight * x%a2F
   end subroutine initialize
end module eliashberg_spectral_function
