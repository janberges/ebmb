! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module eliashberg_eigenvalue
   use eigenvalues
   use eliashberg_self_energy
   use global
   implicit none

   private
   public :: eigenvalue

contains

   subroutine eigenvalue(status, x)
      type(parameters), intent(in) :: x

      real(dp), intent(out) :: status ! greatest eigenvalue

      real(dp), allocatable :: kernel(:, :) ! Eliashberg kernel
      real(dp), allocatable, save :: phi(:) ! order parameter

      type(matsubara) :: im
      type(occupancy) :: oc

      call self_energy(x, im, oc, kernel)

      if (x%power .and. x%bands .eq. 1) then
         if (allocated(phi)) then
            if (size(phi) .ne. size(kernel, 2)) deallocate(phi)
         end if

         if (.not. allocated(phi)) then
            allocate(phi(size(kernel, 2)))

            phi(:) = 0
            phi(1) = 1
         end if

         call power_method(kernel, phi, status)
      else
         status = maxval(real(spectrum(kernel), dp))
      end if
   end subroutine eigenvalue
end module eliashberg_eigenvalue
