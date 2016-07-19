module tc
   use eliashberg
   use global
   implicit none

contains

   subroutine estimate(i)
      type(universal), intent(inout) :: i

      real(dp) :: lambda, muStar

      lambda = sqrt(sum(i%lambda ** 2))
      muStar = sqrt(sum(i%muStar ** 2))

      i%TcMD = i%omegaE / 1.2_dp * exp(-1.04_dp * (1 + lambda) &
         / (lambda - 0.62_dp * lambda * muStar - muStar))
   end subroutine estimate

   subroutine bisection(i, im)
      type(universal), intent(inout) :: i
      type(matsubara), intent(out) :: im

      real(dp), parameter :: ratio = 0.1_dp

      integer :: p, q
      real(dp) :: lower(i%bands), upper(i%bands)

      print '(A13)', 'T/K'

      allocate(i%TcEB(i%bands))

      lower(:) = -1
      upper(:) = -1

      i%T = i%TcMD
      call bounds

      do p = 1, i%bands
         i%T = upper(p)

         do while (lower(p) .lt. 0)
            i%T = i%T * (1 - ratio)
            call bounds
         end do

         i%T = lower(p)

         do while (upper(p) .lt. 0)
            i%T = i%T * (1 + ratio)
            call bounds
         end do

         do
            i%T = (lower(p) + upper(p)) / 2
            call bounds

            if (abs(im%Delta(0, p)) .le. i%small .and. &
               (upper(p) - lower(p)) / 2 .le. i%error) then

               print '(A13)', '(ok)'

               i%TcEB(p) = i%T
               exit
            end if
         end do
      end do

   contains

      subroutine bounds
         print '(F13.9)', i%T / kB

         call solve(i, im)

         do q = 1, i%bands
            if (abs(im%Delta(0, q)) .le. i%small) then
               if (upper(q) .gt. i%T .or. upper(q) .lt. 0) upper(q) = i%T
            else
               if (lower(q) .lt. i%T .or. lower(q) .lt. 0) lower(q) = i%T
            end if
         end do
      end subroutine bounds

   end subroutine bisection
end module tc
