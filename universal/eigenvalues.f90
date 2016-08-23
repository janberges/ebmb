module eigenvalues
   use global
   use tools, only: bound
   implicit none
   private

   public :: spectrum, power_method

   interface
      subroutine dgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &
         work, lwork, info)

         use global

         character, intent(in) :: jobvl, jobvr

         integer, intent(in) :: n, lda, ldvl, ldvr, lwork
         integer, intent(out) :: info

         real(dp), intent(inout) :: a(lda, *)
         real(dp), intent(out) :: wr(*), wi(*), vl(lda, *), vr(lda, *), work(*)
      end subroutine dgeev
   end interface

contains

   function spectrum(matrix, error)
      real(dp), intent(in) :: matrix(:, :)
      integer, intent(out), optional :: error

      complex(dp) :: spectrum(size(matrix, 1))

      integer :: n, info

      real(dp) :: a(size(matrix, 1), size(matrix, 2))
      real(dp) :: wr(size(matrix, 1)), wi(size(matrix, 1))
      real(dp) :: v(1, 1), work(3 * size(matrix, 1))

      a(:, :) = matrix

      n = size(matrix, 1)

      call dgeev(           &
         & jobvl = 'N',     &
         & jobvr = 'N',     &
         &     n = n,       &
         &     a = a(1, 1), &
         &   lda = n,       &
         &    wr = wr(1),   &
         &    wi = wi(1),   &
         &    vl = v(1, 1), &
         &  ldvl = n,       &
         &    vr = v(1, 1), &
         &  ldvr = n,       &
         &  work = work(1), &
         & lwork = 3 * n,   &
         &  info = info     )

      spectrum = cmplx(wr, wi, dp)

      if (present(error)) error = info
   end function spectrum

   subroutine power_method(matrix, vector, value)
      real(dp), intent(inout) :: matrix(:, :), vector(:)
      real(dp), intent(out) :: value

      real(dp) :: shift, value0

      integer :: i

      shift = bound(matrix)

      do i = 1, size(matrix, 1)
         matrix(i, i) = matrix(i, i) + shift
      end do

      value0 = -1

      do
         vector(:) = matmul(matrix, vector)

         value = sqrt(sum(vector ** 2))
         vector(:) = vector / value

         if (value .ap. value0) exit

         value0 = value
      end do

      value = value - shift
   end subroutine power_method
end module eigenvalues
