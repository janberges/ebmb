module eliashberg_eigenvalue_cdos
   use eigenvalues
   use global
   implicit none

   private
   public :: eigenvalue_cdos

contains

   subroutine eigenvalue_cdos(status, x)
      type(parameters), intent(in) :: x

      real(dp), intent(out) :: status ! greatest eigenvalue

      real(dp), allocatable, save :: &
         lambda(:, :, :), & ! adjusted phonon Green function
         muStar(:, :),    & ! Coulomb pseudo-potential
         matrix(:, :),    & ! Eliashberg kernel
         vector(:),       & ! order parameter
         values(:),       & ! all eigenvalues
         diag  (:)          ! diagonal renormalization contribution

      integer :: no ! index of overall cutoff frequency
      integer :: nC ! index of Coulomb cutoff frequency

      integer, save :: no0 = -1 ! 'no' from previous subroutine call

      integer :: i, j ! band indices
      integer :: n, m ! frequency indices
      integer :: p, q ! index offsets

      real(dp) :: nE ! 'index' defining omegaE as bosonic Matsubara frequency

      nE = x%omegaE / (2 * pi * kB * x%T)

      no = ceiling(x%cutoff  * nE - 0.5_dp)
      nC = ceiling(x%cutoffC * nE - 0.5_dp)

      if (no .ne. no0) then
         if (no0 .ne. -1) then
            deallocate(lambda)
            deallocate(muStar)
            deallocate(matrix)
            deallocate(vector)
            deallocate(values)
            deallocate(diag)
         end if

         allocate(lambda(1 - no:2 * no - 1, 0:x%bands - 1, 0:x%bands - 1))
         allocate(muStar(                   0:x%bands - 1, 0:x%bands - 1))

         allocate(matrix(0:x%bands * no - 1, 0:x%bands * no - 1))
         allocate(vector(0:x%bands * no - 1))
         allocate(values(0:x%bands * no - 1))

         allocate(diag(0:x%bands * no - 1))

         vector(:) = 0
         vector(0) = 1

         no0 = no
      end if

      do n = 1 - no, 2 * no - 1
         lambda(n, :, :) = x%lambda / (1 + (n / nE) ** 2)
      end do

      if (x%rescale) then
         muStar(:, :) = x%muStar / (1 + x%muStar * log(nE / (nC + 0.5_dp)))
      else
         muStar(:, :) = x%muStar
      end if

      do i = 0, x%bands - 1
         p = i * no

         do j = 0, x%bands - 1
            q = j * no

            do n = 0, no - 1
               do m = 0, no - 1
                  matrix(q + m, p + n) &
                     = lambda(n - m, j, i) + lambda(n + m + 1, j, i)
               end do
            end do

            matrix(q:q + nC - 1, p:p + no - 1) = &
            matrix(q:q + nC - 1, p:p + no - 1) - 2 * muStar(j, i)
         end do
      end do

      do i = 0, x%bands - 1
         p = i * no

         if (x%imitate) then
            do n = 0, no - 1
               diag(p + n) = sum &
                  (lambda(n:n - no + 1:-1, :, i) - lambda(n + 1:n + no, :, i))
            end do
         else
            diag(p) = sum(lambda(0, :, i))

            do n = 1, no - 1
               diag(p + n) = diag(p + n - 1) + 2 * sum(lambda(n, :, i))
            end do
         end if
      end do

      do i = 0, x%bands * no - 1
         matrix(i, i) = matrix(i, i) - diag(i)
      end do

      do m = 0, no - 1
         matrix(m::no, :) = matrix(m::no, :) / (2 * m + 1)
      end do

      if (x%power .and. x%bands .eq. 1) then
         call power_method(matrix, vector, status)
      else
         values(:) = real(spectrum(matrix), dp)
         status = maxval(values)
      end if
   end subroutine eigenvalue_cdos
end module eliashberg_eigenvalue_cdos
