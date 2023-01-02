! Copyright (C) 2016-2023 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module real_axis
   use global
   use real_axis_pade
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
         allocate(re%status(x%bands))
      end if

      if (x%resolution .gt. 0) then
         allocate(re%omega(x%resolution))
         allocate(re%Delta(x%resolution, x%bands))
         allocate(re%Z(x%resolution, x%bands))

         if (x%ldos) allocate(re%chi(x%resolution, x%bands))
      end if

      if (x%measurable .or. x%resolution .gt. 0) then
         do i = 1, x%bands
            call coefficients(im%omega, im%Delta(:, i))

            if (x%measurable) then
               re%Delta0(i) = 1
               re%status(i) = -1

               do n = 1, x%limit
                  Delta0 = real(continuation(cmplx(re%Delta0(i), kind=dp)))

                  if (re%Delta0(i) .ap. Delta0) re%status(i) = n

                  re%Delta0(i) = Delta0

                  if (n .eq. re%status(i)) exit
               end do
            end if

            if (x%resolution .gt. 0) then
               call interval(re%omega, x%lower, x%upper, &
                  lower=.true., upper=.true.)

               allocate(omega(x%resolution))

               omega(:) = cmplx(re%omega, x%eta, dp)

               re%Delta(:, i) = continuation(omega)

               call coefficients(im%omega, im%Z(:, i))
               re%Z(:, i) = continuation(omega)

               if (x%ldos) then
                  call coefficients(im%omega, im%chi(:, i))
                  re%chi(:, i) = continuation(omega)
               end if
            end if
         end do
      end if
   end subroutine realize
end module real_axis
