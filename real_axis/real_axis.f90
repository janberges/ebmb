module real_axis
   use global
   use real_axis_pade
   implicit none

   private
   public :: realize

contains

   subroutine realize(x, im, re)
      type(universal), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(out) :: re

      integer :: i, n
      real(dp) :: Delta0

      if (x%measurable) then
         allocate(re%Delta0(x%bands))
         allocate(re%status(x%bands))
      end if

      if (x%resolution .gt. 0) then
         allocate(re%omega(x%resolution))
         allocate(re%Delta(x%resolution, x%bands))
         allocate(re%Z(x%resolution, x%bands))

         if (x%chi) allocate(re%chi(x%resolution, x%bands))
      end if

      if (x%measurable .or. x%resolution .gt. 0) then
         do i = 1, x%bands
            call coefficients(im%omega, im%Delta(:, i))

            if (x%measurable) then
               re%Delta0(i) = 1
               re%status(i) = -1

               do n = 1, x%limit
                  Delta0 = real(continuation(re%Delta0(i)), dp)

                  if (re%Delta0(i) .ap. Delta0) re%status(i) = n

                  re%Delta0(i) = Delta0

                  if (n .eq. re%status(i)) exit
               end do
            end if

            if (x%resolution .gt. 0) then
               call interval(re%omega, 0.0_dp, x%upper, &
                  lower=.true., upper=.false.)

               do n = 1, x%resolution
                  re%Delta(n, i) = continuation(re%omega(n))
               end do

               call coefficients(im%omega, im%Z(:, i))

               do n = 1, x%resolution
                  re%Z(n, i) = continuation(re%omega(n))
               end do

               if (x%chi) then
                  call coefficients(im%omega, im%chi(:, i))

                  do n = 1, x%resolution
                     re%chi(n, i) = continuation(re%omega(n))
                  end do
               end if
            end if
         end do
      end if
   end subroutine realize

   subroutine interval(x, a, b, lower, upper)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: a, b
      logical, intent(in), optional :: lower, upper

      integer :: i, j, k

      i = size(x)
      j = 1

      if (present(lower)) then
         if (lower) j = j - 1
      end if

      if (present(upper)) then
         if (upper) i = i - 1
      end if

      do k = 1, size(x)
         x(k) = i * a + j * b
         i = i - 1
         j = j + 1
      end do

      x = x / (i + j)
   end subroutine interval
end module real_axis
