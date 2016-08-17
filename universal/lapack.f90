module lapack
   use global
   implicit none
   private

   public :: eigenvalues

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

   subroutine eigenvalues(matrix, values, error)
      real(dp), intent(in) :: matrix(:, :)
      complex(dp), intent(out) :: values(:)
      integer, intent(out), optional :: error

      integer :: n, lwork, info

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

      values = cmplx(wr, wi, dp)

      if (present(error)) error = info
   end subroutine eigenvalues
end module lapack
