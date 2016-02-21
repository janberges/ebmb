program eb
   use arguments
   use eliashberg
   use global
   use io
   use realaxis
   implicit none

   integer :: n
   type(info) :: i

   do n = 1, command_argument_count()
      call load(argument(n), i)

      call solve(i)
      if (i%continue) call realize(i)

      call store(i)
   end do
end program eb
