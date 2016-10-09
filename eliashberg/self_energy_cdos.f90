module eliashberg_self_energy_cdos
   use global
   implicit none

   private
   public :: self_energy_cdos

contains

   subroutine self_energy_cdos(x, im)
      type(parameters), intent(in) :: x
      type(matsubara), intent(out) :: im

      real(dp) :: nE, Z, Delta

      real(dp), allocatable :: lambda(:, :, :), mu(:, :, :)
      real(dp), allocatable :: muStar(:, :), A(:, :)

      integer :: step, i, j, n, m, no, nC
      logical :: done

      nE = x%omegaE / (2 * pi * kB * x%T)

      no = ceiling(x%cutoff  * nE - 0.5_dp)
      nC = ceiling(x%cutoffC * nE - 0.5_dp)

      allocate(im%omega(0:no - 1))

      do n = 0, no - 1
         im%omega(n) = (2 * n + 1) * pi * kB * x%T
      end do

      allocate(lambda(1 - no:2 * no - 1, x%bands, x%bands))

      do n = 1 - no, 2 * no - 1
         lambda(n, :, :) = x%lambda / (1 + (n / nE) ** 2)
      end do

      allocate(muStar(x%bands, x%bands))

      if (x%rescale) then
         muStar(:, :) = x%muStar / (1 + x%muStar * log(nE / (nC + 0.5_dp)))
      else
         muStar(:, :) = x%muStar
      end if

      allocate(mu(0:no - 1, x%bands, x%bands))

      do n = 0, nC - 1
         mu(n, :, :) = -2 * muStar
      end do

      mu(nC:, :, :) = 0

      allocate(im%Z(0:no - 1, x%bands))

      im%Z(:, :) = 1

      allocate(im%Delta(0:no - 1, x%bands))

      im%Delta(:, :) = 0

      if (.not. x%normal) im%Delta(0, :) = 1

      allocate(A(0:no - 1, x%bands))

      do i = 1, x%bands
         A(:, i) = 1 / sqrt(im%omega ** 2 + im%Delta(:, i) ** 2)
      end do

      im%status = -1

      do step = 1, x%limit
         done = .true.

         do i = 1, x%bands
            do n = 0, no - 1
               Z = 0
               Delta = 0

               do j = 1, x%bands
                  do m = 0, no - 1
                     Z = Z + im%omega(m) * A(m, j) &
                        * (lambda(n - m, j, i) - lambda(n + m + 1, j, i))

                     Delta = Delta + im%Delta(m, j) * A(m, j) * (mu(m, j, i) &
                        +  lambda(n - m, j, i) + lambda(n + m + 1, j, i))
                  end do
               end do

               Z = 1 + pi * kB * x%T * Z / im%omega(n)
               Delta = pi * kB * x%T * Delta / Z

               done = done &
                  .and. (im%Z    (n, i) .ap. Z) &
                  .and. (im%Delta(n, i) .ap. Delta)

               im%Z    (n, i) = Z
               im%Delta(n, i) = Delta

               A(n, i) = 1 / sqrt(im%omega(n) ** 2 + Delta ** 2)
            end do
         end do

         if (done) then
            im%status = step
            exit
         end if
      end do

      allocate(im%phiC(x%bands))

      do i = 1, x%bands
         im%phiC(i) = pi * kB * x%T * sum(im%Delta * A * mu(:, :, i))
      end do
   end subroutine self_energy_cdos
end module eliashberg_self_energy_cdos
