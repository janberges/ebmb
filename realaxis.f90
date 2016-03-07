module realaxis
   use global
   use intervals
   use pade
   implicit none

contains

   subroutine realize(i)
      type(info), intent(inout) :: i

      integer :: n
      real(dp) :: Delta0

      if (i%measurable .or. i%resolution .gt. 0) then
         call coefficients(i%omega, i%Delta)
      end if

      if (i%measurable) then
         i%Delta0 = 1
         i%statusDelta0 = -1

         do n = 1, i%limit
            Delta0 = real(continuation(i%Delta0), dp)

            if (i%Delta0 .ap. Delta0) i%statusDelta0 = n

            i%Delta0 = Delta0

            if (n .eq. i%statusDelta0) exit
         end do
      end if

      if (i%resolution .gt. 0) then
         if (allocated(i%omega_)) deallocate(i%omega_)
         allocate(i%omega_(i%resolution))

         call interval(i%omega_, 0.0_dp, i%upper, lower=.true., upper=.false.)

         if (allocated(i%Delta_)) deallocate(i%Delta_)
         allocate(i%Delta_(i%resolution))

         do n = 1, i%resolution
            i%Delta_(n) = continuation(i%omega_(n))
         end do

         call coefficients(i%omega, i%Z)

         if (allocated(i%Z_)) deallocate(i%Z_)
         allocate(i%Z_(i%resolution))

         do n = 1, i%resolution
            i%Z_(n) = continuation(i%omega_(n))
         end do
      end if
   end subroutine realize
end module realaxis
