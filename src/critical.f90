! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

program critical
   use eliashberg_eigenvalue
   use eliashberg_eigenvalue_cdos
   use globals
   use io_load
   implicit none

   type(parameters), target :: x

   real(dp), pointer :: variable => null() ! parameter to be optimized

   procedure(eigenvalue), pointer :: solver => null() ! solver to be used

   real(dp) :: bound(2) ! bisection bounds

   real(dp) :: ev  ! greatest eigenvalue
   real(dp) :: ev0 ! ... in previous step

   logical :: sc1 ! bound(1) within superconducting phase?
   logical :: try ! still trying out direction?

   integer :: i, j ! band indices
   integer :: error ! I/O status

   integer :: argmin(2) ! indices of minimum matrix element

   call load(x)

   if (x%T .lt. 0.0_dp) then
      variable => x%T
   else if (x%omegaE .lt. 0.0_dp) then
      variable => x%omegaE
   else if (any(x%lambda .lt. 0.0_dp)) then
      argmin(:) = minloc(x%lambda)
      variable => x%lambda(argmin(1), argmin(2))
   else if (any(x%muStar .lt. 0.0_dp)) then
      argmin(:) = minloc(x%muStar)
      variable => x%muStar(argmin(1), argmin(2))
   else
      variable => x%T
      x%T = -x%T
   end if

   variable = -variable

   if (x%ldos) then
      solver => eigenvalue
   else
      solver => eigenvalue_cdos
   end if

   call solver(ev, x)

   ev0 = ev

   sc1 = ev .ge. 1.0_dp
   try = .true.

   do
      bound(1) = variable
      variable = variable * (1.0_dp + x%rate)

      call solver(ev, x)

      if (ev .ap. ev0) then
         print "('Error: Stationary point')"
         stop 1
      end if

      if (sc1 .neqv. ev .ge. 1.0_dp) exit

      if (sc1 .eqv. ev .gt. ev0) then
         if (try) then
            variable = bound(1)
            x%rate = -x%rate
            try = .false.
            cycle
         end if

         print "('Error: Local extremum')"
         stop 1
      end if

      ev0 = ev
   end do

   bound(2) = variable

   do
      variable = 0.5_dp * sum(bound)

      if (abs(variable - bound(1)) .le. x%error) exit

      call solver(ev, x)

      if (sc1 .eqv. ev .ge. 1.0_dp) then
         bound(1) = variable
      else
         bound(2) = variable
      end if
   end do

   if (x%tell) print '(' // trim(x%flomat) // ')', variable

   if (x%output .ne. 'none') then
      open (fun, file=x%output, action='write', status='replace', &
         access='stream', iostat=error)

      if (error .ne. 0) then
         print "('Error: Cannot open output file ""', A, '""')", trim(x%output)
         stop 1
      end if

      write (fun) variable
      close (fun)
   end if
end program critical
