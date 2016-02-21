module pade
   use global
   implicit none

   private
   public :: coefficients, continuation

   integer :: n, p, q

   complex(dp), parameter :: i = (0, 1)

   complex(qp), allocatable :: c(:, :)
   complex(qp), allocatable :: d(:, :)

contains

   subroutine coefficients(z, u)
      real(dp), intent(in) :: z(:), u(:)

      n = size(z)

      allocate(c(n, n))

      c(1, :) = u

      do p = 2, n
         do q = p, n
            c(p, q) = (c(p - 1, p - 1) - c(p - 1, q)) &
               / (i * (z(q) - z(p - 1)) * c(p - 1, q))
         end do

         c(p, p - 1) = -i * z(p - 1) * c(p, p)
      end do

      allocate(d(2, 0:n))

      d(1, 0) = 0
      d(1, 1) = c(1, 1)
      d(2, :) = 1
   end subroutine coefficients

   function continuation(x)
      complex(dp) :: continuation
      real(dp), intent(in) :: x

      do p = 2, n
         d(:, p) = d(:, p - 1) + d(:, p - 2) * (c(p, p - 1) + c(p, p) * x)
      end do

      continuation = cmplx(d(1, n) / d(2, n), kind=dp)
   end function continuation
end module pade
