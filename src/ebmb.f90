! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

program ebmb
   use dos
   use eliashberg_self_energy
   use eliashberg_self_energy_cdos
   use eliashberg_self_energy_real_axis
   use globals
   use io_load
   use io_store
   use io_tell
   use real_axis
   implicit none

   type(parameters) :: x
   type(matsubara) :: im
   type(continued) :: re
   type(occupancy) :: oc

   integer :: omp_num_threads

   call load(x)

   if (x%tell) then
      print "('This is  _   v2.0.0  __')"
      print "('    ___ | |_  __ __ ( (_')"
      print "('   / __)| _ \/  Y  \| _ \')"
      print "('   \___,|___/\  |  /|___/.')"

      omp_num_threads = 0
      !$omp parallel reduction(+:omp_num_threads)
      omp_num_threads = 1
      !$omp end parallel

      if (omp_num_threads .eq. 1) then
         print "(/, 'Running serially.')"
      else
         print "(/, 'Running on ', I0, ' threads.')", omp_num_threads
      end if
   end if

   if (x%realgw) then
      call self_energy_real_axis(x, im, re, oc)
   else if (x%ldos) then
      call self_energy(x, im, oc)
   else
      call self_energy_cdos(x, im)
   end if

   if (x%Sigma) call combine_self_energy_components(x, im)

   if (.not. x%realgw) then
      call realize(x, im, re)

      if (x%ldos .and. x%points .gt. 0) then
         if (x%stable) then
            call density_of_states_stable(x, im, re, oc)
         else
            call density_of_states(x, re, oc)
         end if
      end if
   end if

   if (x%output .ne. 'none') call store(x, im, re, oc)

   if (x%tell) call tell(x, im, re, oc)
end program ebmb
