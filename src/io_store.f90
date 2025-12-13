! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module io_store
   use globals
   implicit none

   private
   public :: store

contains

   subroutine store(x, im, re, oc)
      type(parameters), intent(in) :: x
      type(matsubara), intent(in) :: im
      type(continued), intent(in) :: re
      type(occupancy), intent(in) :: oc

      integer :: error

      open (fun, file=x%output, action='write', status='replace', &
         access='stream', iostat=error)

      if (error .ne. 0) then
         print "('Error: Cannot open output file ""', A, '""')", trim(x%output)
         stop 1
      end if

      write (fun) 'INT:DIM:', 0_i4
      write (fun) 'status:', im%steps

      write (fun) 'REAL:DIM:', 1_i4, size(im%omega, kind=i4)
      write (fun) 'iomega:', im%omega

      if (x%bands .gt. 1) &
         write (fun) 'DIM:', 2_i4, x%bands, size(im%omega, kind=i4)

      write (fun) 'domega:', aimag(im%Sigma)
      write (fun) 'Z:', im%Z
      write (fun) 'Delta:', im%Delta

      if (x%ldos) then
         write (fun) 'chi:', im%chi

         write (fun) 'DIM:', 0_i4

         write (fun) 'states:', oc%states

         if (x%points .gt. 0) write (fun) 'inspect:', oc%inspect

         write (fun) 'n0:', oc%n0
         write (fun) "n:", oc%n

         write (fun) 'mu0:', oc%mu0
         write (fun) "mu:", oc%mu
      end if

      if (x%la2F) then
         write (fun) 'DIM:'

         if (x%bands .gt. 1) then
            write (fun) 2_i4, x%bands, x%bands
         else
            write (fun) 0_i4
         end if

         write (fun) 'lambda:', x%lambda

         write (fun) 'DIM:', 0_i4
         write (fun) 'omegaE:', x%omegaE
         write (fun) 'omegaLog:', x%omegaLog
         write (fun) 'omega2nd:', x%omega2nd
      end if

      write (fun) 'DIM:'

      if (x%bands .gt. 1) then
         write (fun) 1_i4, x%bands
      else
         write (fun) 0_i4
      end if

      write (fun) 'phiC:', im%phiC

      if (x%ldos) write (fun) 'chiC:', im%chiC

      if (x%measurable) then
         write (fun) 'INT:status0:', re%steps
         write (fun) 'REAL:Delta0:', re%Delta0
      end if

      if (x%points .gt. 0) then
         write (fun) 'DIM:', 1_i4, x%points

         write (fun) 'omega:', re%omega

         if (x%bands .gt. 1) write (fun) 'DIM:', 2_i4, x%bands, x%points

         write (fun) 'Re[Z]:', real(re%Z)
         write (fun) 'Im[Z]:', aimag(re%Z)

         write (fun) 'Re[Delta]:', real(re%Delta)
         write (fun) 'Im[Delta]:', aimag(re%Delta)

         if (x%ldos) then
            write (fun) 'Re[chi]:', real(re%chi)
            write (fun) 'Im[chi]:', aimag(re%chi)

            write (fun) 'DOS:', re%dos
         end if

         write (fun) 'Re[Sigma]:', real(re%Sigma)
         write (fun) 'Im[Sigma]:', aimag(re%Sigma)
      end if

      close (fun)
   end subroutine store
end module io_store
