! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

program critical
   use eliashberg_eigenvalue
   use eliashberg_eigenvalue_cdos
   use global
   use io_load
   implicit none

   type(parameters), target :: x

   real(dp), pointer :: variable => null() ! parameter to be optimized

   procedure(eigenvalue), pointer :: solver => null() ! solver to be used

   real(dp) :: bound(2) ! bisection bounds

   real(dp) :: status  ! greatest eigenvalue
   real(dp) :: status0 ! ... in previous step

   logical :: sc1 ! bound(1) within superconducting phase?
   logical :: try ! still trying out direction?

   integer :: i, j ! band indices
   integer :: error ! I/O status

   call load(x)

   variable => x%T

   if (x%T .lt. 0.0_dp) then
      variable => x%T
      variable = -variable
   end if

   if (x%omegaE .lt. 0.0_dp) then
      variable => x%omegaE
      variable = -variable
   end if

   do i = 1, x%bands
      do j = 1, x%bands
         if (x%lambda(j, i) .lt. 0.0_dp) then
            variable => x%lambda(j, i)
            variable = -variable
         end if

         if (x%muStar(j, i) .lt. 0.0_dp) then
            variable => x%muStar(j, i)
            variable = -variable
         end if
      end do
   end do

   if (x%ldos) then
      solver => eigenvalue
   else
      solver => eigenvalue_cdos
   end if

   call solver(status, x)

   status0 = status

   sc1 = status .ge. 1.0_dp
   try = .true.

   do
      bound(1) = variable
      variable = variable * (1.0_dp + x%rate)

      call solver(status, x)

      if (status .ap. status0) then
         print "('Error: Stationary point')"
         stop 1
      end if

      if (sc1 .neqv. status .ge. 1.0_dp) exit

      if (sc1 .eqv. status .gt. status0) then
         if (try) then
            variable = bound(1)
            x%rate = -x%rate
            try = .false.
            cycle
         end if

         print "('Error: Local extremum')"
         stop 1
      end if

      status0 = status
   end do

   bound(2) = variable

   do
      variable = 0.5_dp * sum(bound)

      if (abs(variable - bound(1)) .le. x%error) exit

      call solver(status, x)

      if (sc1 .eqv. status .ge. 1.0_dp) then
         bound(1) = variable
      else
         bound(2) = variable
      end if
   end do

   if (x%tell) print '(' // trim(x%form) // ')', variable

   if (x%file .ne. 'none') then
      open (unit, file=x%file, action='write', status='replace', &
         access='stream', iostat=error)

      if (error .ne. 0) then
         print "('Error: Cannot open output file ""', A, '""')", trim(x%file)
         stop 1
      end if

      write (unit) variable
      close (unit)
   end if
end program critical
