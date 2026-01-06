! Copyright (C) 2016-2026 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

program tc
   use eliashberg_self_energy
   use eliashberg_self_energy_cdos
   use formatting
   use globals
   use io_load
   implicit none

   type(parameters) :: x
   type(matsubara) :: im
   type(occupancy) :: oc

   integer :: i, j ! band indices

   real(dp), allocatable :: upper(:), lower(:), T(:) ! bounds and Tc's

   character(:), allocatable :: head, body ! edit descriptors

   call load(x)

   if (x%tell) then
      call measure(x%flomat)

      head = edit('(Aw)')
      body = edit('(x)')

      print head, 'T/K'
      print rule(1)
   end if

   allocate(T(x%bands))

   allocate(upper(x%bands))
   allocate(lower(x%bands))

   lower(:) = -1.0_dp
   upper(:) = -1.0_dp

   call bounds

   bands: do i = 1, x%bands
      x%T = upper(i)

      do while (lower(i) .lt. 0.0_dp)
         if (x%T .le. x%error) then
            T(i) = 0.0_dp
            cycle bands
         end if

         x%T = x%T * (1.0_dp - x%rate)
         call bounds
      end do

      x%T = lower(i)

      do while (upper(i) .lt. 0.0_dp)
         x%T = x%T * (1.0_dp + x%rate)
         call bounds
      end do

      do
         x%T = 0.5_dp * (lower(i) + upper(i))

         if (upper(i) - lower(i) .le. 2.0_dp * x%error) then
            T(i) = x%T
            cycle bands
         end if

         call bounds
      end do
   end do bands

   if (x%tell) then
      print *
      print head, 'Tc/K'
      print rule(1)
      print body, T
   end if

   if (x%output .ne. 'none') then
      open (fun, &
         file=x%output, action='write', status='replace', access='stream')
      write (fun) T
      close (fun)
   end if

contains

   subroutine bounds
      if (x%tell) print body, x%T

      if (x%ldos) then
         call self_energy(x, im, oc)
      else
         call self_energy_cdos(x, im)
      end if

      do j = 1, x%bands
         if (abs(im%Delta(0, j)) .le. x%zero) then
            if (upper(j) .gt. x%T .or. upper(j) .lt. 0.0_dp) upper(j) = x%T
         else
            if (lower(j) .lt. x%T .or. lower(j) .lt. 0.0_dp) lower(j) = x%T
         end if
      end do
   end subroutine bounds
end program tc
