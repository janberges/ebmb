module eliashberg_eigenvalue
   use global
   implicit none

   private
   public :: eigenvalue

contains

   subroutine eigenvalue(status, x)
      type(parameters), intent(in) :: x

      real(dp), intent(out) :: status  ! greatest eigenvalue
      real(dp), save        :: status0 ! ... in previous step

      real(dp), allocatable, save :: &
         lambda(:, :, :), & ! frequency-dependent electron-phonon coupling
         renorm(:, :, :), & ! frequency-diagonal renormalization contribution
         muStar(:, :),    & ! rescaled Coulomb pseudo-potential
         matrix(:, :),    & ! Eliashberg matrix
         vector(:)          ! energy gap

      integer :: no ! index of overall cutoff frequency
      integer :: nC ! index of Coulomb cutoff frequency

      integer, save :: no0 = -1 ! 'no' from previous subroutine call

      integer :: i, j ! band indices
      integer :: n, m ! frequency indices
      integer :: p, q ! index offsets

      logical :: done ! eigenvalue converged?

      real(dp) :: nE ! 'index' defining omegaE as bosonic Matsubara frequency

      nE = x%omegaE / (2 * pi * kB * x%T)

      no = ceiling(x%cutoff  * nE - 0.5_dp)
      nC = ceiling(x%cutoffC * nE - 0.5_dp)

      if (no .ne. no0) then
         if (no0 .ne. -1) then
            deallocate(lambda)
            deallocate(renorm)
            deallocate(muStar)
            deallocate(matrix)
            deallocate(vector)
         end if

         allocate(lambda(0:2 * no - 1, 0:x%bands - 1, 0:x%bands - 1))
         allocate(renorm(0:    no - 1, 0:x%bands - 1, 0:x%bands - 1))

         allocate(muStar(0:x%bands - 1, 0:x%bands - 1))

         allocate(matrix(0:x%bands * no - 1, 0:x%bands * no - 1))
         allocate(vector(0:x%bands * no - 1))

         vector = (x%bands * no) ** (-0.5_dp)
         status0 = 1

         no0 = no
      end if

      do n = 0, 2 * no - 1
         lambda(n, :, :) = x%lambda / (1 + (n / nE) ** 2)
      end do

      if (x%imitate) then
         do n = 0, no - 1
            renorm(n, :, :) = 0

            do m = 0, no - 1
               renorm(n, :, :) = renorm(n, :, :) &
                  + lambda(abs(n - m), :, :) - lambda(n + m + 1, :, :)
            end do
         end do
      else
         renorm(0, :, :) = lambda(0, :, :)

         do n = 1, no - 1
            renorm(n, :, :) = renorm(n - 1, :, :) + 2 * lambda(n, :, :)
         end do
      end if

      if (x%rescale) then
         muStar(:, :) = x%muStar / (1 + x%muStar * log(nE / (nC + 0.5_dp)))
      else
         muStar(:, :) = x%muStar
      end if

      do i = 0, x%bands - 1
         p = i * no

         do j = 0, x%bands - 1
            q = j * no

            matrix(q     :q + nC - 1, p:p + no - 1) = -2 * muStar(j, i)
            matrix(q + nC:q + no - 1, p:p + no - 1) = 0

            do n = 0, no - 1
               matrix(q + n, p + n) = matrix(q + n, p + n) - renorm(n, j, i)

               do m = 0, no - 1
                  matrix(q + m, p + n) = (matrix(q + m, p + n) &
                     + lambda(abs(n - m), j, i) + lambda(n + m + 1, j, i))
               end do
            end do
         end do
      end do

      do m = 0, no - 1
         matrix(m::no, :) = matrix(m::no, :) / (2 * m + 1)
      end do

      do p = 0, x%bands * no - 1
         matrix(p, p) = matrix(p, p) + x%shift
      end do

      done = .false.

      do while (.not. done)
         vector = matmul(matrix, vector)
         status = sqrt(dot_product(vector, vector))

         if (status .ap. status0) done = .true.

         vector = vector / status
         status0 = status
      end do

      status = status - x%shift
   end subroutine eigenvalue
end module eliashberg_eigenvalue
