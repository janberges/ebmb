program ebmb
   use eliashberg_self_energy
   use eliashberg_self_energy_cdos
   use global
   use io_load
   use io_store
   use io_tell
   use real_axis
   implicit none

   type(parameters) :: x
   type(matsubara) :: im
   type(continued) :: re

   call load(x)

   if (x%chi) then
      call self_energy(x, im)
   else
      call self_energy_cdos(x, im)
   end if

   call realize(x, im, re)

   if (x%file .ne. 'none') call store(x, im, re)

   if (x%tell) call tell(x, im, re)
end program ebmb
