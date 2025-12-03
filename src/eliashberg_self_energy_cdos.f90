! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module eliashberg_self_energy_cdos
   use eliashberg_spectral_function
   use globals
   implicit none

   private
   public :: self_energy_cdos

contains

   subroutine self_energy_cdos(x, im)
      type(parameters), intent(in) :: x
      type(matsubara), intent(out) :: im

      real(dp) :: nE

      real(dp), allocatable :: lambda(:, :, :), mu(:, :, :)
      real(dp), allocatable :: muStar(:, :), A(:, :)
      real(dp), allocatable :: Z(:, :), Delta(:, :)

      integer :: step, i, j, n, m, no, nC

      if (.not. x%normal .and. x%steps .le. 10) then
         print "('Warning: Superconducting solution should be self-consistent')"
      end if

      nE = x%omegaE / (2.0_dp * pi * kB * x%T)

      no = ceiling(x%cutoff  * nE - 0.5_dp)
      nC = ceiling(x%cutoffC * nE - 0.5_dp)

      if (no .lt. 1)  no = 1

      allocate(im%omega(0:no - 1))

      do n = 0, no - 1
         im%omega(n) = (2 * n + 1) * pi * kB * x%T
      end do

      allocate(lambda(1 - no:2 * no - 1, x%bands, x%bands))

      do n = 1 - no, 2 * no - 1
         if (x%la2F) then
            call lambda_from_a2F(x, lambda(n, :, :), n)
         else
            lambda(n, :, :) = x%lambda / (1.0_dp + (n / nE) ** 2)
         end if
      end do

      allocate(muStar(x%bands, x%bands))

      if (x%rescale) then
         muStar(:, :) = x%muStar / (1.0_dp + x%muStar * log(nE / (nC + 0.5_dp)))
      else
         muStar(:, :) = x%muStar
      end if

      allocate(mu(0:no - 1, x%bands, x%bands))

      do n = 0, nC - 1
         mu(n, :, :) = -2.0_dp * muStar
      end do

      mu(nC:, :, :) = 0.0_dp

      allocate(im%Z(0:no - 1, x%bands))

      im%Z(:, :) = 1.0_dp

      allocate(im%Delta(0:no - 1, x%bands))

      im%Delta(:, :) = 0.0_dp

      if (.not. x%normal) im%Delta(0, :) = 1.0_dp

      allocate(A(0:no - 1, x%bands))

      allocate(Z(0:no - 1, x%bands))
      allocate(Delta(0:no - 1, x%bands))

      im%steps = -1

      do step = 1, x%steps
         do i = 1, x%bands
            A(:, i) = 1.0_dp / sqrt(im%omega ** 2 + im%Delta(:, i) ** 2)
         end do

         Z(:, :) = im%Z
         Delta(:, :) = im%Delta

         im%Z(:, :) = 0.0_dp
         im%Delta(:, :) = 0.0_dp

         do i = 1, x%bands
            !$omp parallel do
            do n = 0, no - 1
               do j = 1, x%bands
                  if (x%diag .and. i .ne. j) cycle

                  do m = 0, no - 1
                     im%Z(n, i) = im%Z(n, i) + im%omega(m) * A(m, j) &
                        * (lambda(n - m, j, i) - lambda(n + m + 1, j, i))

                     im%Delta(n, i) = im%Delta(n, i) + Delta(m, j) * A(m, j) &
                        * (lambda(n - m, j, i) + lambda(n + m + 1, j, i) &
                           + mu(m, j, i))
                  end do
               end do

            end do
            !$omp end parallel do

            im%Z(:, i) = 1.0_dp + pi * kB * x%T * im%Z(:, i) / im%omega
         end do

         im%Delta(:, :) = pi * kB * x%T * im%Delta / Z

         if (all(im%Z .ap. Z) .and. all(im%Delta .ap. Delta)) then
            im%steps = step
            exit
         end if
      end do

      allocate(im%phiC(x%bands))

      do i = 1, x%bands
         im%phiC(i) = pi * kB * x%T * sum(im%Delta * A * mu(:, :, i))
      end do

      if (x%Sigma) then
         allocate(im%Sigma(0:no - 1, x%bands))

         do i = 1, x%bands
            im%Sigma(:, i) = cmplx(im%phi(:, i), &
               im%omega * (1.0_dp - im%Z(:, i)), dp)
         end do
      end if
   end subroutine self_energy_cdos
end module eliashberg_self_energy_cdos
