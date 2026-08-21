! Copyright (C) 2016-2026 Jan Berges
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

      integer :: i, n, nP
      real(dp) :: Delta0
      complex(dp), allocatable :: omega(:)

      nP = ceiling(x%cutoffP * x%omegaE / (2.0_dp * pi * kB * x%T) - 0.5_dp)

      if (nP .lt. 1) nP = 1
      if (nP .gt. size(im%omega)) nP = size(im%omega)

      if (x%gap) then
         allocate(re%Delta0(x%bands))
         allocate(re%steps(x%bands))
      end if

      if (x%points .gt. 0) then
         allocate(re%omega(x%points))
         allocate(omega(x%points))
         allocate(re%Z(x%points, x%bands))
         allocate(re%Delta(x%points, x%bands))
         allocate(re%phi(x%points, x%bands))

         if (x%ldos) allocate(re%chi(x%points, x%bands))

         allocate(re%Sigma(x%points, x%bands))

         call interval(re%omega, x%lower, x%upper, lower=.true., upper=.true., &
            logscale=x%logscale)

         omega(:) = cmplx(re%omega, x%eta, dp)
      end if

      if (x%gap .or. x%points .gt. 0) then
         do i = 1, x%bands
            call coefficients(im%omega(:nP - 1), &
               cmplx(im%Delta(:nP - 1, i), kind=dp))

            if (x%gap) then
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

               call coefficients(im%omega(:nP - 1), &
                  cmplx(im%Z(:nP - 1, i), kind=dp))

               !$omp parallel do
               do n = 1, x%points
                  re%Z(n, i) = continuation(omega(n))
               end do
               !$omp end parallel do

               call coefficients(im%omega(:nP - 1), &
                  cmplx(im%phi(:nP - 1, i), kind=dp))

               !$omp parallel do
               do n = 1, x%points
                  re%phi(n, i) = continuation(omega(n))
               end do
               !$omp end parallel do

               if (x%ldos) then
                  call coefficients(im%omega(:nP - 1), &
                     cmplx(im%chi(:nP - 1, i), kind=dp))

                  !$omp parallel do
                  do n = 1, x%points
                     re%chi(n, i) = continuation(omega(n))
                  end do
                  !$omp end parallel do
               end if

               call coefficients(im%omega(:nP - 1), im%Sigma(:nP - 1, i))

               !$omp parallel do
               do n = 1, x%points
                  re%Sigma(n, i) = continuation(omega(n))
               end do
               !$omp end parallel do
            end if
         end do
      end if
   end subroutine realize
end module real_axis
