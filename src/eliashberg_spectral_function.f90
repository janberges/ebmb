! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module eliashberg_spectral_function
   use global
   use tools, only: differential
   implicit none

   private
   public :: integrate_a2F, lambda_from_a2F, initialize, weight

   logical :: initial = .true.
   integer :: i, j

   real(dp), allocatable :: weight(:, :, :)

contains

   subroutine integrate_a2F(x)
      type(parameters), intent(inout) :: x

      real(dp) :: omegaMax

      call lambda_from_a2F(x, x%lambda, 0)

      x%omegaLog = exp(2.0_dp / sum(x%lambda) &
         * sum(sum(sum(weight, 2), 2) * log(x%omega) / x%omega))

      x%omega2nd = sqrt(2.0_dp / sum(x%lambda) &
         * sum(sum(sum(weight, 2), 2) * x%omega))

      do i = size(x%omega), 1, -1
         if (any(x%a2F(i, :, :) .gt. 0.0_dp)) then
            omegaMax = x%omega(i)
            exit
         end if
      end do

      !x%omegaE = x%omegaLog
      x%omegaE = x%omega2nd ! choice by Allen and Dynes
      !x%omegaE = x%omegaMax
   end subroutine integrate_a2F

   subroutine lambda_from_a2F(x, lambda, n)
      type(parameters), intent(in) :: x
      real(dp), intent(out) :: lambda(:, :)
      integer, intent(in) :: n

      if (initial) call initialize(x)

      do i = 1, x%bands
         do j = 1, x%bands
            lambda(j, i) = 2.0_dp * sum(weight(:, j, i) * x%omega &
               / (x%omega ** 2 + (2 * n * pi * kB * x%T) ** 2))
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
