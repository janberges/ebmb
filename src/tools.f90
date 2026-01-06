! Copyright (C) 2016-2026 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module tools
   use globals
   implicit none

   private
   public :: argument, bound, differential, interval, matches

contains

   function argument(n)
      character(:), allocatable :: argument
      integer, intent(in) :: n

      integer :: length

      call get_command_argument(n, length=length)

      allocate(character(length) :: argument)

      call get_command_argument(n, value=argument)
   end function argument

   real(dp) function bound(matrix)
      real(dp), intent(in) :: matrix(:, :)

      real(dp) :: R, C, S

      integer :: i

      R = 0.0_dp
      do i = 1, size(matrix, 1)
         S = sum(abs(matrix(i, :)))
         if (S .gt. R) R = S
      end do

      C = 0.0_dp
      do i = 1, size(matrix, 2)
         S = sum(abs(matrix(:, i)))
         if (S .gt. C) C = S
      end do

      bound = min(R, C)
   end function bound

   subroutine differential(x, dx)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: dx(:)

      integer :: n
      n = size(x)

      dx(1) = x(2) - x(1)
      dx(2:n - 1) = x(3:n) - x(1:n - 2)
      dx(n) = x(n) - x(n - 1)

      dx(:) = 0.5_dp * dx
   end subroutine differential

   subroutine interval(x, a, b, lower, upper, logscale)
      real(dp), intent(out) :: x(:)
      real(dp), intent(in) :: a, b
      logical, intent(in), optional :: lower, upper
      real(dp), intent(in), optional :: logscale

      integer :: i, j, k
      real(dp) :: l, m

      l = a
      m = b

      if (present(logscale)) then
         if (logscale .gt. 0.0_dp) then
            l = loga(l)
            m = loga(m)
         end if
      end if

      i = size(x)
      j = 1

      if (present(lower)) then
         if (lower) j = j - 1
      end if

      if (present(upper)) then
         if (upper) i = i - 1
      end if

      if (i + j .eq. 0) then
         i = 1
         j = 1
      end if

      do k = 1, size(x)
         x(k) = i * l + j * m
         i = i - 1
         j = j + 1
      end do

      x = x / (i + j)

      if (present(logscale)) then
         if (logscale .gt. 0.0_dp) x = expo(x)
      end if

   contains

      elemental function expo(x)
         real(dp) :: expo
         real(dp), intent(in) :: x

         expo = sign((exp(abs(x)) - 1.0_dp) / logscale, x)
      end function expo

      elemental function loga(x)
         real(dp) :: loga
         real(dp), intent(in) :: x

         loga = sign(log(abs(x * logscale) + 1.0_dp), x)
      end function loga
   end subroutine interval

   integer function matches(str, chr)
      character(*), intent(in) :: str
      character(1), intent(in) :: chr

      integer :: c

      matches = 0

      do c = 1, len(str)
         if (str(c:c) .eq. chr) matches = matches + 1
      end do
   end function matches
end module tools
