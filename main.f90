program eb_local
   use arguments
   use eliashberg
   use global
   use io
   use realaxis
   use tc
   implicit none

   integer :: n

   type(universal) :: i
   type(matsubara) :: im
   type(continued) :: re

   do n = 1, command_argument_count()
      call load(argument(n), i)

      call estimate(i)

      if (i%critical) then
         call bisection(i, im)
      else
         call solve(i, im)
      end if

      call realize(i, im, re)

      call store(i, im, re)
   end do
end program eb_local
