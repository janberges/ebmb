program ebmb
   use eliashberg_constant_dos
   use eliashberg_variable_dos
   use global
   use io_load
   use io_tell
   use io_store
   use real_axis
   implicit none

   type(universal) :: x
   type(matsubara) :: im
   type(continued) :: re

   call load(x)

   if (x%chi) then
      call solve_variable_dos(x, im)
   else
      call solve_constant_dos(x, im)
   end if

   call realize(x, im, re)

   if (x%file .ne. 'none') call store(x, im, re)

   if (x%tell) call tell(x, im, re)
end program ebmb
