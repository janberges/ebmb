module tools
   use global
   implicit none

   private
   public :: argument, bound, differential, interval, matches

contains

   function argument(n)
      character(:), allocatable :: argument
      integer, intent(in) :: n

      integer :: size

      call get_command_argument(n, length=size)

      allocate(character(size) :: argument)

      call get_command_argument(n, value=argument)
   end function argument

   real(dp) function bound(matrix)
      real(dp), intent(in) :: matrix(:, :)

      real(dp) :: R, C, S

      integer :: i

      R = 0
      do i = 1, size(matrix, 1)
         S = sum(abs(matrix(i, :)))
         if (S .gt. R) R = S
      end do

      C = 0
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

      dx(:) = dx / 2
   end subroutine differential

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

   integer function matches(str, char)
      character(*), intent(in) :: str
      character(1), intent(in) :: char

      integer :: c

      matches = 0

      do c = 1, len(str)
         if (str(c:c) .eq. char) matches = matches + 1
      end do
   end function matches
end module tools
