program eb_local
   use bisection_critical
   use bisection_tc
   use eliashberg_constant_dos
   use eliashberg_variable_dos
   use global
   use io_load
   use io_store_data
   use io_store_text
   use real_axis
   implicit none

   type(universal) :: x
   type(matsubara) :: im
   type(continued) :: re

   call load(x)

   call estimate(x)

   select case (x%mode)
      case ('self-energy')
         if (x%chi) then
            call solve_variable_dos(x, im)
         else
            call solve_constant_dos(x, im)
         end if

         call realize(x, im, re)

      case ('critical')
         call critical(x)

         print '(F0.15)', x%variable

      case ('tc')
         call tc(x, im)
   end select

   if (x%form .eq. 'text' .or. x%form .eq. 'both') call store_text(x, im, re)
   if (x%form .eq. 'data' .or. x%form .eq. 'both') call store_data(x, im, re)
end program eb_local
