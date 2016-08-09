module tools
   use global
   implicit none

   private
   public :: argument, interval, matches

contains

   function argument(n)
      character(:), allocatable :: argument
      integer, intent(in) :: n

      integer :: size

      call get_command_argument(n, length=size)

      allocate(character(size) :: argument)

      call get_command_argument(n, value=argument)
   end function argument

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
