program tc
   use eliashberg_self_energy
   use eliashberg_self_energy_cdos
   use formatting
   use global
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
      call measure(x%form)

      head = edit('(Aw)')
      body = edit('(x)')

      print head, 'T/K'
      print rule(1)
   end if

   allocate(T(x%bands))

   allocate(upper(x%bands))
   allocate(lower(x%bands))

   lower(:) = -1
   upper(:) = -1

   call bounds

   BANDS: do i = 1, x%bands
      x%T = upper(i)

      do while (lower(i) .lt. 0)
         if (x%T .le. x%error) then
            T(i) = 0
            cycle BANDS
         end if

         x%T = x%T * (1 - x%rate)
         call bounds
      end do

      x%T = lower(i)

      do while (upper(i) .lt. 0)
         x%T = x%T * (1 + x%rate)
         call bounds
      end do

      do
         x%T = (lower(i) + upper(i)) / 2

         if (upper(i) - lower(i) .le. 2 * x%error) then
            T(i) = x%T
            cycle BANDS
         end if

         call bounds
      end do
   end do BANDS

   if (x%tell) then
      print *
      print head, 'Tc/K'
      print rule(1)
      print body, T
   end if

   if (x%file .ne. 'none') then
      open (unit, &
         file=x%file, action='write', status='replace', access='stream')
      write (unit) T
      close (unit)
   end if

contains

   subroutine bounds
      if (x%tell) print body, x%T

      if (x%chi) then
         call self_energy(x, im, oc)
      else
         call self_energy_cDOS(x, im)
      end if

      do j = 1, x%bands
         if (abs(im%Delta(0, j)) .le. x%zero) then
            if (upper(j) .gt. x%T .or. upper(j) .lt. 0) upper(j) = x%T
         else
            if (lower(j) .lt. x%T .or. lower(j) .lt. 0) lower(j) = x%T
         end if
      end do
   end subroutine bounds
end program tc
