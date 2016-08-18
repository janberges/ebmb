program critical
   use eliashberg_eigenvalue
   use global
   use io_load
   implicit none

   type(parameters), target :: x

   real(dp), pointer :: variable => null() ! parameter to be optimized

   real(dp) :: inner ! bound within superconducting phase
   real(dp) :: outer ! bound beyond superconducting phase

   real(dp) :: status  ! greatest eigenvalue
   real(dp) :: status0 ! ... in previous step

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

   try = .true.

   call eigenvalue(status, x)
   status0 = status

   if (status .ge. 1) then
      do while (status .ge. 1)
         inner = variable

         variable = variable * (1 + x%rate)

         call eigenvalue(status, x)

         if (status .gt. status0) then
            if (try) then
               variable = variable / (1 + x%rate)
               x%rate = -x%rate
               try = .false.
            else
               stop 'stuck at local minimum'
            end if
         else
            status0 = status
         end if
      end do

      outer = variable
   else
      do while (status .lt. 1)
         outer = variable

         variable = variable * (1 - x%rate)

         call eigenvalue(status, x)

         if (status .lt. status0) then
            if (try) then
               variable = variable / (1 - x%rate)
               x%rate = -x%rate
               try = .false.
            else
               stop 'stuck at local maximum'
            end if
         else
            status0 = status
         end if
      end do

      inner = variable
   end if

   do
      variable = (inner + outer) / 2

      if (abs(outer - inner) .le. 2 * x%error) exit

      call eigenvalue(status, x)

      if (status .ge. 1) then
         inner = variable
      else
         outer = variable
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
