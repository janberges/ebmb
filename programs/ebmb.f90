! Copyright (C) 2021 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

program ebmb
   use dos
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
   type(occupancy) :: oc

   call load(x)

   if (x%ldos) then
      call self_energy(x, im, oc)
   else
      call self_energy_cdos(x, im)
   end if

   call realize(x, im, re)

   if (x%ldos) then
      call density_of_states(x, re, oc)
   end if

   if (x%file .ne. 'none') call store(x, im, re, oc)

   if (x%tell) call tell(x, im, re, oc)
end program ebmb
