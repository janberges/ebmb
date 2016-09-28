module real_axis_pade
   use global
   implicit none

   private
   public :: coefficients, continuation

   integer :: n
   complex(qp), allocatable :: c(:, :)

contains

   subroutine coefficients(z, u)
      real(dp), intent(in) :: z(:), u(:)

      complex(dp), parameter :: i = (0, 1)
      integer :: p

      n = size(z)

      if (allocated(c)) deallocate(c)
      allocate(c(n, n))

      if (all(u .ap. 0.0_dp)) then
         c(:, :) = 0
         return
      end if

      c(1, :) = u

      do p = 2, n
         c(p, p:) = (c(p - 1, p - 1) - c(p - 1, p:)) &
            / (i * (z(p:) - z(p - 1)) * c(p - 1, p:))

         c(p, p - 1) = -i * z(p - 1) * c(p, p)
      end do
   end subroutine coefficients

   elemental function continuation(x)
      complex(dp) :: continuation
      real(dp), intent(in) :: x

      complex(qp) :: frac
      integer :: p

      frac = 1

      do p = n, 2, -1
         frac = 1 + (c(p, p) * x + c(p, p - 1)) / frac
      end do

      frac = c(1, 1) / frac

      continuation = cmplx(frac, kind=dp)
   end function continuation
end module real_axis_pade
