program eb
   use arguments
   use eliashberg
   use global
   use io
   use mcmillan
   use realaxis
   implicit none

   integer :: n
   type(info) :: i

   do n = 1, command_argument_count()
      call load(argument(n), i)

      call estimate(i)
      call solve(i)
      if (i%continue) call realize(i)

      if (i%form .eq. 'text' .or. i%form .eq. 'both') call save_text(i)
      if (i%form .eq. 'data' .or. i%form .eq. 'both') call save_data(i)
   end do
end program eb
