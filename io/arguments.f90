module arguments
   implicit none

contains

   function argument(n)
      character(:), allocatable :: argument
      integer, intent(in) :: n

      integer :: size

      call get_command_argument(n, length=size)

      allocate(character(size) :: argument)

      call get_command_argument(n, value=argument)
   end function argument
end module arguments
