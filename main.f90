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
   type(matsubara) :: im
   type(continued) :: re

   do n = 1, command_argument_count()
      call load(argument(n), i)

      call estimate(i)

      if (i%critical) then
         call bisection(i, im)
      else
         call solve(i, im)
         call realize(i, im, re)
      end if

      if (i%form .eq. 'text' .or. i%form .eq. 'both') call save_text(i, im, re)
      if (i%form .eq. 'data' .or. i%form .eq. 'both') call save_data(i, im, re)
   end do
end program eb
