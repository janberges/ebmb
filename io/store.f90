module io_store
   use global
   implicit none

   private
   public :: store

contains

   subroutine store(x, im, re, oc)
      type(parameters), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re
      type(occupancy), intent(in) :: oc

      open (unit, &
         file=x%file, action='write', status='replace', access='stream')

      write (unit) 'INT:DIM:', 0_i4
      write (unit) 'status:', im%status

      write (unit) 'REAL:DIM:', 1_i4, size(im%omega, kind=i4)
      write (unit) 'iomega:', im%omega

      if (x%bands .gt. 1) &
         write (unit) 'DIM:', 2_i4, x%bands, size(im%omega, kind=i4)

      write (unit) 'Z:', im%Z
      write (unit) 'Delta:', im%Delta

      if (x%chi) then
         write (unit) 'chi:', im%chi

         write (unit) 'DIM:', 0_i4

         write (unit) 'n0:', oc%n0
         write (unit) "n:", oc%n

         write (unit) 'mu0:', oc%mu0
         write (unit) "mu:", oc%mu
      end if

      write (unit) 'DIM:'

      if (x%bands .gt. 1) then
         write (unit) 1_i4, x%bands
      else
         write (unit) 0_i4
      end if

      write (unit) 'phiC:', im%phiC

      if (x%measurable) then
         write (unit) 'INT:status0:', re%status
         write (unit) 'REAL:Delta0:', re%Delta0
      end if

      if (x%resolution .gt. 0) then
         write (unit) 'DIM:', 1_i4, x%resolution

         write (unit) 'omega:', re%omega

         if (x%bands .gt. 1) write (unit) 'DIM:', 2_i4, x%bands, x%resolution

         write (unit) 'Re[Z]:', real(re%Z)
         write (unit) 'Im[Z]:', aimag(re%Z)

         write (unit) 'Re[Delta]:', real(re%Delta)
         write (unit) 'Im[Delta]:', aimag(re%Delta)

         if (x%chi) then
            write (unit) 'Re[chi]:', real(re%chi)
            write (unit) 'Im[chi]:', aimag(re%chi)
         end if
      end if

      close (unit)
   end subroutine store
end module io_store
