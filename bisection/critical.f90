module bisection_critical
   use eliashberg_eigenvalue
   use global
   implicit none

   private
   public :: critical

contains

   subroutine critical(x)
      type(universal), target, intent(inout) :: x

      real(dp) :: inner ! bound within superconducting phase
      real(dp) :: outer ! bound beyond superconducting phase

      real(dp) :: status  ! greatest eigenvalue
      real(dp) :: status0 ! ... in previous step

      logical :: try ! still trying out direction?

      integer :: i, j ! band indices

      x%variable => x%T

      if (x%T .lt. 0) then
         x%variable => x%T
         x%variable = -x%variable
      end if

      if (x%omegaE .lt. 0) then
         x%variable => x%omegaE
         x%variable = -x%variable
      end if

      do i = 1, x%bands
         do j = 1, x%bands
            if (x%lambda(j, i) .lt. 0) then
               x%variable => x%lambda(j, i)
               x%variable = -x%variable
            end if

            if (x%muStar(j, i) .lt. 0) then
               x%variable => x%muStar(j, i)
               x%variable = -x%variable
            end if
         end do
      end do

      try = .true.

      call greatest_eigenvalue(status, x)
      status0 = status

      if (status .ge. 1) then
         do while (status .ge. 1)
            inner = x%variable

            x%variable = x%variable * (1 + x%rate)

            call greatest_eigenvalue(status, x)

            if (status .gt. status0) then
               if (try) then
                  x%variable = x%variable / (1 + x%rate)
                  x%rate = -x%rate
                  try = .false.
               else
                  stop 'stuck at local minimum'
               end if
            else
               status0 = status
            end if
         end do

         outer = x%variable
      else
         do while (status .lt. 1)
            outer = x%variable

            x%variable = x%variable * (1 - x%rate)

            call greatest_eigenvalue(status, x)

            if (status .lt. status0) then
               if (try) then
                  x%variable = x%variable / (1 - x%rate)
                  x%rate = -x%rate
                  try = .false.
               else
                  stop 'stuck at local maximum'
               end if
            else
               status0 = status
            end if
         end do

         inner = x%variable
      end if

      do
         x%variable = (inner + outer) / 2

         if (abs(outer - inner) .le. 2 * x%error) exit

         call greatest_eigenvalue(status, x)

         if (status .ge. 1) then
            inner = x%variable
         else
            outer = x%variable
         end if
      end do
   end subroutine critical
end module bisection_critical
