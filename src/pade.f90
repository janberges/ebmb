! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module pade
   use globals
   implicit none

   private
   public :: coefficients, continuation

   integer :: n
   complex(qp), allocatable :: c(:, :)

contains

   subroutine coefficients(z, u)
      real(dp), intent(in) :: z(:)
      complex(dp), intent(in) :: u(:)

      complex(dp), parameter :: i = (0.0_dp, 1.0_dp)
      integer :: p

      if (all(u .ap. u(1))) then
         n = 1
      else
         n = size(z)
      end if

      if (allocated(c)) deallocate(c)
      allocate(c(n, n))

      c(1, :) = u(:n)

      do p = 2, n
         c(p, p:) = (c(p - 1, p - 1) - c(p - 1, p:)) &
            / (i * (z(p:) - z(p - 1)) * c(p - 1, p:))

         c(p, p - 1) = -i * z(p - 1) * c(p, p)
      end do
   end subroutine coefficients

   elemental function continuation(x)
      complex(dp) :: continuation
      complex(dp), intent(in) :: x

      complex(qp) :: frac
      integer :: p

      frac = (1.0_qp, 0.0_qp)

      do p = n, 2, -1
         frac = (1.0_dp, 0.0_dp) + (c(p, p) * x + c(p, p - 1)) / frac
      end do

      frac = c(1, 1) / frac

      continuation = cmplx(frac, kind=dp)
   end function continuation
end module pade
