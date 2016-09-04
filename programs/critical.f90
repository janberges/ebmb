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

   call load(x)

   variable => x%T

   if (x%T .lt. 0) then
      variable => x%T
      variable = -variable
   end if

   if (x%omegaE .lt. 0) then
      variable => x%omegaE
      variable = -variable
   end if

   do i = 1, x%bands
      do j = 1, x%bands
         if (x%lambda(j, i) .lt. 0) then
            variable => x%lambda(j, i)
            variable = -variable
         end if

         if (x%muStar(j, i) .lt. 0) then
            variable => x%muStar(j, i)
            variable = -variable
         end if
      end do
   end do

   if (x%chi) then
      solver => eigenvalue
   else
      solver => eigenvalue_cdos
   end if

   call solver(status, x)

   status0 = status

   sc1 = status .ge. 1
   try = .true.

   do
      bound(1) = variable
      variable = variable * (1 + x%rate)

      call solver(status, x)

      if (status .eq. status0) stop 'stationary point'

      if (sc1 .neqv. status .ge. 1) exit

      if (sc1 .eqv. status .gt. status0) then
         if (try) then
            variable = bound(1)
            x%rate = -x%rate
            try = .false.
            cycle
         end if

         stop 'local extremum'
      end if

      status0 = status
   end do

   bound(2) = variable

   do
      variable = sum(bound) / 2

      if (abs(variable - bound(1)) .le. x%error) exit

      call solver(status, x)

      if (sc1 .eqv. status .ge. 1) then
         bound(1) = variable
      else
         bound(2) = variable
      end if
   end do

   if (x%tell) print '(' // trim(x%form) // ')', variable

   if (x%file .ne. 'none') then
      open (unit, &
         file=x%file, action='write', status='replace', access='stream')
      write (unit) variable
      close (unit)
   end if
end program critical
