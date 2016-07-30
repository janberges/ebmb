program eb_local
   use eliashberg_constant_dos
   use eliashberg_variable_dos
   use global
   use io_load
   use io_store_data
   use io_store_text
   use real_axis
   use tc
   implicit none

   type(universal) :: x
   type(matsubara) :: im
   type(continued) :: re

   call load(x)

   call estimate(x)

   if (x%critical) then
      call bisection(x, im)
   else
      if (x%chi) then
         call solve_variable_dos(x, im)
      else
         call solve_constant_dos(x, im)
      end if
   end if

   call realize(x, im, re)

   if (x%form .eq. 'text' .or. x%form .eq. 'both') call store_text(x, im, re)
   if (x%form .eq. 'data' .or. x%form .eq. 'both') call store_data(x, im, re)
end program eb_local
