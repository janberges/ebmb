program eb
   use arguments
   use eliashberg
   use global
   use io
   use realaxis
   use tc
   implicit none

   integer :: n
   type(info) :: i

   do n = 1, command_argument_count()
      call load(argument(n), i)

      call estimate(i)

      if (i%critical) then
         call bisection(i)
      else
         call solve(i)
         call realize(i)
      end if

      if (i%form .eq. 'text' .or. i%form .eq. 'both') call save_text(i)
      if (i%form .eq. 'data' .or. i%form .eq. 'both') call save_data(i)
   end do
end program eb
