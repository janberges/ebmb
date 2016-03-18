module integration
   use global
   implicit none

contains

   subroutine differential(x, dx)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: dx(:)

      integer :: n
      n = size(x)

      dx(1) = x(2) - x(1)
      dx(2:n - 1) = x(3:n) - x(1:n - 2)
      dx(n) = x(n) - x(n - 1)

      dx(:) = dx / 2
   end subroutine differential
end module integration
