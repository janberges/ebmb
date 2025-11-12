! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module real_axis
   use globals
   use pade
   use tools, only: interval
   implicit none

   private
   public :: realize

contains

   subroutine realize(x, im, re)
      type(parameters), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(out) :: re

      integer :: i, n
      real(dp) :: Delta0
      complex(dp), allocatable :: omega(:)

      if (x%measurable) then
         allocate(re%Delta0(x%bands))
         allocate(re%steps(x%bands))
      end if

      if (x%points .gt. 0) then
         allocate(re%omega(x%points))
         allocate(omega(x%points))
         allocate(re%Delta(x%points, x%bands))
         allocate(re%Z(x%points, x%bands))

         if (x%ldos) allocate(re%chi(x%points, x%bands))
         if (x%Sigma) allocate(re%Sigma(x%points, x%bands))

         call interval(re%omega, x%lower, x%upper, lower=.true., upper=.true., &
            logscale=x%logscale)

         omega(:) = cmplx(re%omega, x%eta, dp)
      end if

      if (x%measurable .or. x%points .gt. 0) then
         do i = 1, x%bands
            call coefficients(im%omega, cmplx(im%Delta(:, i), kind=dp))

            if (x%measurable) then
               re%Delta0(i) = 1.0_dp
               re%steps(i) = -1

               do n = 1, x%steps
                  Delta0 = real(continuation(cmplx(re%Delta0(i), kind=dp)))

                  if (re%Delta0(i) .ap. Delta0) re%steps(i) = n

                  re%Delta0(i) = Delta0

                  if (n .eq. re%steps(i)) exit
               end do
            end if

            if (x%points .gt. 0) then
               !$omp parallel do
               do n = 1, x%points
                  re%Delta(n, i) = continuation(omega(n))
               end do
               !$omp end parallel do

               call coefficients(im%omega, cmplx(im%Z(:, i), kind=dp))

               !$omp parallel do
               do n = 1, x%points
                  re%Z(n, i) = continuation(omega(n))
               end do
               !$omp end parallel do

               if (x%ldos) then
                  call coefficients(im%omega, cmplx(im%chi(:, i), kind=dp))

                  !$omp parallel do
                  do n = 1, x%points
                     re%chi(n, i) = continuation(omega(n))
                  end do
                  !$omp end parallel do
               end if

               if (x%Sigma) then
                  call coefficients(im%omega, im%Sigma(:, i))

                  !$omp parallel do
                  do n = 1, x%points
                     re%Sigma(n, i) = continuation(omega(n))
                  end do
                  !$omp end parallel do
               end if
            end if
         end do
      end if
   end subroutine realize
end module real_axis
