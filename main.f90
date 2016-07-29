program eb_local
   use arguments
   use eliashberg_constant_dos
   use eliashberg_variable_dos
   use global
   use io_load
   use io_store_data
   use io_store_text
   use realaxis
   use tc
   implicit none

   type(universal) :: i
   type(matsubara) :: im
   type(continued) :: re

   call load(i)

   call estimate(i)

   if (i%critical) then
      call bisection(i, im)
   else
      if (i%DOS) then
         call solve_variable_dos(i, im)
      else
         call solve_constant_dos(i, im)
      end if
   end if

   call realize(i, im, re)

   if (i%form .eq. 'text' .or. i%form .eq. 'both') call store_text(i, im, re)
   if (i%form .eq. 'data' .or. i%form .eq. 'both') call store_data(i, im, re)
end program eb_local
