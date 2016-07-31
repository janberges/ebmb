module eliashberg_eigenvalue
   use global
   implicit none

contains

   subroutine greatest_eigenvalue(status, x)
      type(universal), intent(in) :: x

      real(dp), intent(out) :: status  ! greatest eigenvalue
      real(dp), save        :: status0 ! ... in previous step

      real(dp), allocatable, save :: &
         lambda(:, :, :), & ! frequency-dependent electron-phonon coupling
         renorm(:, :, :), & ! frequency-diagonal renormalization contribution
         muStar(:, :),    & ! rescaled Coulomb pseudo-potential
         matrix(:, :),    & ! Eliashberg matrix
         vector(:)          ! energy gap

      integer       :: u, l    ! number of Matsubara frequencies (with mu*)
      integer, save :: u0 = -1 ! ... in previous subroutine call

      integer :: i, j ! band indices
      integer :: n, m ! frequency indices
      integer :: p, q ! index offsets

      logical :: done ! eigenvalue converged?

      real(dp) :: nE ! 'index' defining omegaE as bosonic Matsubara frequency

      nE = x%omegaE / (2 * pi * kB * x%T)

      u = ceiling(x%upper * nE - 0.5_dp)
      l = ceiling(x%lower * nE - 0.5_dp)

      if (u .ne. u0) then
         if (u0 .ne. -1) then
            deallocate(lambda)
            deallocate(renorm)
            deallocate(muStar)
            deallocate(matrix)
            deallocate(vector)
         end if

         allocate(lambda(0:2 * u - 1, 0:x%bands - 1, 0:x%bands - 1))
         allocate(renorm(0:    u - 1, 0:x%bands - 1, 0:x%bands - 1))

         allocate(muStar(0:x%bands - 1, 0:x%bands - 1))

         allocate(matrix(0:x%bands * u - 1, 0:x%bands * u - 1))
         allocate(vector(0:x%bands * u - 1))

         vector = (x%bands * u) ** (-0.5_dp)
         status0 = 1

         u0 = u
      end if

      do n = 0, 2 * u - 1
         lambda(n, :, :) = x%lambda / (1 + (n / nE) ** 2)
      end do

      if (x%cutoffZ) then
         do n = 0, u - 1
            renorm(n, :, :) = 0

            do m = 0, u - 1
               renorm(n, :, :) = renorm(n, :, :) &
                  + lambda(abs(n - m), :, :) - lambda(n + m + 1, :, :)
            end do
         end do
      else
         renorm(0, :, :) = lambda(0, :, :)

         do n = 1, u - 1
            renorm(n, :, :) = renorm(n - 1, :, :) + 2 * lambda(n, :, :)
         end do
      end if

      if (x%rescale) then
         muStar(:, :) = x%muStar / (1 + x%muStar * log(nE / (l + 0.5_dp)))
      else
         muStar(:, :) = x%muStar
      end if

      do i = 0, x%bands - 1
         p = i * u

         do j = 0, x%bands - 1
            q = j * u

            matrix(q    :q + l - 1, p:p + u - 1) = -2 * muStar(j, i)
            matrix(q + l:q + u - 1, p:p + u - 1) = 0

            do n = 0, u - 1
               matrix(q + n, p + n) = matrix(q + n, p + n) - renorm(n, j, i)

               do m = 0, u - 1
                  matrix(q + m, p + n) = (matrix(q + m, p + n) &
                     + lambda(abs(n - m), j, i) + lambda(n + m + 1, j, i))
               end do
            end do
         end do
      end do

      do m = 0, u - 1
         matrix(m::u, :) = matrix(m::u, :) / (2 * m + 1)
      end do

      do p = 0, x%bands * u - 1
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
   end subroutine greatest_eigenvalue
end module eliashberg_eigenvalue
