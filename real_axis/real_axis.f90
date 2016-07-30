module real_axis
   use global
   use real_axis_pade
   implicit none

contains

   subroutine realize(i, im, re)
      type(universal), intent(in) :: i
      type(matsubara), intent(in) :: im
      type(continued), intent(out) :: re

      integer :: p, n
      real(dp) :: Delta0

      if (i%measurable) then
         allocate(re%Delta0(i%bands))
         allocate(re%status(i%bands))
      end if

      if (i%resolution .gt. 0) then
         allocate(re%omega(i%resolution))
         allocate(re%Delta(i%resolution, i%bands))
         allocate(re%Z(i%resolution, i%bands))

         if (i%chi) allocate(re%chi(i%resolution, i%bands))
      end if

      if (i%measurable .or. i%resolution .gt. 0) then
         do p = 1, i%bands
            call coefficients(im%omega, im%Delta(:, p))

            if (i%measurable) then
               re%Delta0(p) = 1
               re%status(p) = -1

               do n = 1, i%limit
                  Delta0 = real(continuation(re%Delta0(p)), dp)

                  if (re%Delta0(p) .ap. Delta0) re%status(p) = n

                  re%Delta0(p) = Delta0

                  if (n .eq. re%status(p)) exit
               end do
            end if

            if (i%resolution .gt. 0) then
               call interval(re%omega, 0.0_dp, i%upper, &
                  lower=.true., upper=.false.)

               do n = 1, i%resolution
                  re%Delta(n, p) = continuation(re%omega(n))
               end do

               call coefficients(im%omega, im%Z(:, p))

               do n = 1, i%resolution
                  re%Z(n, p) = continuation(re%omega(n))
               end do

               if (i%chi) then
                  call coefficients(im%omega, im%chi(:, p))

                  do n = 1, i%resolution
                     re%chi(n, p) = continuation(re%omega(n))
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
