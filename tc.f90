module tc
   use eliashberg_constant_dos
   use eliashberg_variable_dos
   use global
   implicit none

contains

   subroutine estimate(x)
      type(universal), intent(inout) :: x

      real(dp) :: lambda, muStar

      lambda = sqrt(sum(x%lambda ** 2))
      muStar = sqrt(sum(x%muStar ** 2))

      x%TcMD = x%omegaE / (1.2_dp * kB) * exp(-1.04_dp * (1 + lambda) &
         / (lambda - 0.62_dp * lambda * muStar - muStar))
   end subroutine estimate

   subroutine bisection(x, im)
      type(universal), intent(inout) :: x
      type(matsubara), intent(out) :: im

      real(dp), parameter :: ratio = 0.1_dp

      integer :: i, j
      real(dp) :: lower(x%bands), upper(x%bands)

      print '(A13)', 'T/K'

      allocate(x%TcEB(x%bands))

      lower(:) = -1
      upper(:) = -1

      x%T = max(x%TcMD, x%bound)
      call tell
      call bounds

      BANDS: do i = 1, x%bands
         x%T = upper(i)

         do while (lower(i) .lt. 0)
            if (x%T .le. x%error) then
               call critical
               cycle BANDS
            end if

            if (x%T .eq. x%bound) then
               print '(A13)', '(bound)'

               x%TcEB(i) = -1
               cycle BANDS
            end if

            x%T = x%T * (1 - ratio)
            x%T = max(x%T, x%bound)
            call tell
            call bounds
         end do

         x%T = lower(i)

         do while (upper(i) .lt. 0)
            x%T = x%T * (1 + ratio)
            call tell
            call bounds
         end do

         do
            x%T = (lower(i) + upper(i)) / 2

            call tell

            if (upper(i) - lower(i) .le. 2 * x%error) then
               call critical
               cycle BANDS
            end if

            call bounds
         end do
      end do BANDS

   contains

      subroutine tell
         print '(F13.9)', x%T
      end subroutine tell

      subroutine critical
         x%TcEB(i) = x%T

         print '(A13)', '(critical)'
      end subroutine critical

      subroutine bounds
         if (x%chi) then
            call solve_variable_dos(x, im)
         else
            call solve_constant_dos(x, im)
         end if

         do j = 1, x%bands
            if (abs(im%Delta(0, j)) .le. x%small) then
               if (upper(j) .gt. x%T .or. upper(j) .lt. 0) upper(j) = x%T
            else
               if (lower(j) .lt. x%T .or. lower(j) .lt. 0) lower(j) = x%T
            end if
         end do
      end subroutine bounds

   end subroutine bisection
end module tc
