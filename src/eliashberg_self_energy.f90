! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module eliashberg_self_energy
   use eliashberg_spectral_function, only: lambda_from_a2F
   use globals
   use tools, only: differential
   implicit none

   private
   public :: self_energy, initialize, weight

   logical :: initial = .true.

   real(dp), allocatable :: weight(:, :), matsum(:)

contains

   subroutine self_energy(x, im, oc, kernel)
      type(parameters), intent(in) :: x
      type(matsubara), intent(out) :: im
      type(occupancy), intent(out) :: oc

      real(dp), allocatable, intent(out), optional :: kernel(:, :)

      real(dp) :: nE, domega, residue

      real(dp), allocatable :: g(:, :, :), U(:, :, :), dosef(:), residues(:)
      real(dp), allocatable :: muC(:, :), muStar(:, :), A(:, :), B(:, :)
      real(dp), allocatable :: Z(:, :), phi(:, :), chi(:, :)

      real(dp), allocatable :: integral_Z  (:, :)
      real(dp), allocatable :: integral_phi(:, :)
      real(dp), allocatable :: integral_chi(:, :)

      integer :: step, i, j, n, m, p, q, no, nC

      if (.not. x%normal .and. x%steps .le. 10) then
         print "('Warning: Superconducting solution should be self-consistent')"
      end if

      if (initial) call initialize(x, oc)

      domega = 2.0_dp * pi * kB * x%T

      nE = x%omegaE / domega

      no = ceiling(x%cutoff  * nE - 0.5_dp)
      nC = ceiling(x%cutoffC * nE - 0.5_dp)

      if (no .lt. 1)  no = 1

      allocate(im%omega(0:no - 1))

      do n = 0, no - 1
         im%omega(n) = domega * (n + 0.5_dp)
      end do

      allocate(im%Z(0:no - 1, x%bands))

      im%Z(:, :) = 1.0_dp

      allocate(im%phi(0:no - 1, x%bands))

      im%phi(:, :) = 0.0_dp

      allocate(im%chi(0:no - 1, x%bands))

      im%chi(:, :) = 0.0_dp

      allocate(im%chiC(x%bands))

      im%chiC(:) = 0.0_dp

      allocate(Z(0:no - 1, x%bands))
      allocate(phi(0:no - 1, x%bands))
      allocate(chi(0:no - 1, x%bands))

      allocate(A(0:no - 1, x%bands))
      allocate(B(0:no - 1, x%bands))

      allocate(integral_Z  (0:no - 1, x%bands))
      allocate(integral_phi(0:no - 1, x%bands))
      allocate(integral_chi(0:no - 1, x%bands))

      allocate(residues(x%bands))

      if (x%n .ge. 0.0_dp) then
         oc%mu = (x%energy(1) * (2.0_dp * oc%states - x%n) &
            + x%energy(size(x%energy)) * x%n) / (2.0_dp * oc%states)
      else
         oc%mu = x%mu
      end if

      call dos(x%n, x%n .ge. 0, .true.)

      oc%n0 = oc%n
      oc%mu0 = oc%mu

      allocate(dosef(x%bands))

      if (x%divdos) then
         dosef(:) = x%dos(minloc(abs(x%energy - oc%mu), 1), :)
      else
         dosef(:) = 1.0_dp
      end if

      allocate(g(1 - no:2 * no - 1, x%bands, x%bands))

      !$omp parallel do
      do n = 1 - no, 2 * no - 1
         if (x%la2F) then
            call lambda_from_a2F(x, g(n, :, :), n)
         else
            g(n, :, :) = x%lambda / (1.0_dp + (n / nE) ** 2)
         end if

         do i = 1, x%bands
            g(n, :, i) = g(n, :, i) / dosef
         end do
      end do
      !$omp end parallel do

      allocate(muC(x%bands, x%bands))

      if (x%unscale) then
         where (x%energy .ap. oc%mu)
            matsum = 1.0_dp / x%omegaE
         elsewhere
            matsum = x%energy - oc%mu
            matsum = atan2(matsum, x%omegaE) / matsum
         end where

         do i = 1, x%bands
            residue = sum(weight(:, i) / dosef(i) * matsum) / pi

            muC(i, :) = x%muStar(i, :) / (1.0_dp - x%muStar(i, :) * residue)
         end do
      else
         muC(:, :) = x%muStar
      end if

      allocate(muStar(x%bands, x%bands))

      if (x%rescale) then
         where (x%energy .ap. oc%mu)
            matsum = 1.0_dp / (domega * (nC + 0.5_dp))
         elsewhere
            matsum = x%energy - oc%mu
            matsum = atan2(matsum, domega * (nC + 0.5_dp)) / matsum
         end where

         do i = 1, x%bands
            residue = sum(weight(:, i) / dosef(i) * matsum) / pi

            muStar(i, :) = muC(i, :) / (1.0_dp + muC(i, :) * residue)
         end do
      else
         muStar(:, :) = muC
      end if

      allocate(U(0:no - 1, x%bands, x%bands))

      do n = 0, nC - 1
         U(n, :, :) = -2.0_dp * muStar

         do i = 1, x%bands
            U(n, :, i) = U(n, :, i) / dosef
         end do
      end do

      U(nC:, :, :) = 0.0_dp

      if (.not. (x%normal .or. present(kernel))) integral_phi(0, :) = 1.0_dp

      im%steps = -1

      do step = 1, x%steps
         Z(:, :) = im%Z
         phi(:, :) = im%phi
         chi(:, :) = im%chi

         im%Z(:, :) = 0.0_dp
         im%phi(:, :) = 0.0_dp
         im%chi(:, :) = 0.0_dp

         do i = 1, x%bands
            !$omp parallel do
            do n = 0, no - 1
               do j = 1, x%bands
                  if (x%diag .and. i .ne. j) cycle

                  do m = 0, no - 1
                     im%Z(n, i) = im%Z(n, i) + integral_Z(m, j) &
                        * (g(n - m, j, i) - g(n + m + 1, j, i))

                     im%phi(n, i) = im%phi(n, i) + integral_phi(m, j) &
                        * (g(n - m, j, i) + g(n + m + 1, j, i) + U(m, j, i))

                     if (x%chi) then
                        im%chi(n, i) = im%chi(n, i) - integral_chi(m, j) &
                           * (g(n - m, j, i) + g(n + m + 1, j, i))
                     end if
                  end do
               end do
            end do
            !$omp end parallel do

            im%Z(:, i) = 1.0_dp + im%Z(:, i) * kB * x%T / im%omega
         end do

         im%phi(:, :) = im%phi * kB * x%T
         im%chi(:, :) = im%chi * kB * x%T

         if (x%chi .and. x%chiC) then
            call calculate_residue(nC, .false.)

            do i = 1, x%bands
               im%chiC(i) = sum( &
                  (2.0_dp * kB * x%T * sum(integral_chi(:nC - 1, :), 1) &
                     + residues / 2.0_dp) * muC(:, i) / dosef)

               im%chi(:, i) = im%chi(:, i) + im%chiC(i)
            end do
         end if

         call dos(oc%n0, x%conserve, .false.)

         if (all(im%Z .ap. Z) .and. all(im%phi .ap. phi) &
               .and. all(im%chi .ap. chi)) then
            im%steps = step
            exit
         end if
      end do

      allocate(im%Delta(0:no - 1, x%bands))

      im%Delta(:, :) = im%phi / im%Z

      allocate(im%phiC(x%bands))

      do i = 1, x%bands
         im%phiC(i) = kB * x%T * sum(integral_phi * U(:, :, i))
      end do

      if (x%Sigma) then
         allocate(im%Sigma(0:no - 1, x%bands))

         do i = 1, x%bands
            im%Sigma(:, i) = cmplx(im%phi(:, i) + im%chi(:, i), &
               im%omega * (1.0_dp - im%Z(:, i)), dp)
         end do
      end if

      if (present(kernel)) then
         allocate(kernel(x%bands * no, x%bands * no))

         do i = 1, x%bands
            p = i * no
            do j = 1, x%bands
               q = j * no
               !$omp parallel do
               do n = 0, no - 1
                  do m = 0, no - 1
                     kernel(q - m, p - n) = kB * x%T * A(m, j) &
                        * (g(n - m, j, i) + g(n + m + 1, j, i) + U(m, j, i))
                  end do
               end do
               !$omp end parallel do
            end do
         end do
      end if

      if (x%readjust) call dos(oc%n0, .true., .false.)

   contains

      subroutine dos(ntarget, optimize, exact)
         real(dp), intent(in) :: ntarget
         logical, intent(in) :: optimize, exact

         do
            do i = 1, x%bands
               !$omp parallel do
               do n = 0, no - 1
                  call integrate(n, i)
               end do
               !$omp end parallel do
            end do

            call calculate_residue(no, exact)

            oc%n = oc%states - 4.0_dp * kB * x%T * sum(integral_chi) - residue

            if (abs(oc%n - ntarget) .le. x%toln .or. .not. optimize) exit

            oc%mu = ((ntarget - oc%states + residue) / (4.0_dp * kB * x%T) &
               + sum(A * im%chi + B)) / sum(A)
         end do
      end subroutine dos

      subroutine integrate(n, i)
         integer, intent(in) :: n, i
         real(dp) :: trapezia(size(x%energy))

         trapezia(:) = weight(:, i) / ((im%omega(n) * im%Z(n, i)) ** 2 &
            + (x%energy - oc%mu + im%chi(n, i)) ** 2 + im%phi(n, i) ** 2)

         A(n, i) = sum(trapezia)
         B(n, i) = sum(trapezia * x%energy)

         integral_Z  (n, i) = A(n, i) *  im%Z  (n, i) * im%omega(n)
         integral_phi(n, i) = A(n, i) *  im%phi(n, i)
         integral_chi(n, i) = A(n, i) * (im%chi(n, i) - oc%mu) + B(n, i)
      end subroutine integrate

      subroutine calculate_residue(nM, exact)
         integer, intent(in) :: nM
         logical, intent(in) :: exact

         do i = 1, x%bands
            matsum(:) = x%energy - oc%mu + im%chiC(i)

            if (exact) then
               residues(i) = sum(weight(:, i) &
                  * tanh(matsum / (2.0_dp * kB * x%T)))

               do n = 0, nM - 1
                  residues(i) = residues(i) &
                     - 4.0_dp * kB * x%T * sum(weight(:, i) &
                        * matsum / (im%omega(n) ** 2 + matsum ** 2))
               end do
            else
               residues(i) = 2.0_dp / pi * sum(weight(:, i) &
                  * atan2(matsum, domega * (nM + 0.5_dp)))
            end if
         end do

         residue = sum(residues)
      end subroutine calculate_residue

   end subroutine self_energy

   subroutine initialize(x, oc)
      type(parameters), intent(in) :: x
      type(occupancy), intent(out) :: oc

      integer :: i

      initial = .false.

      if (allocated(weight)) deallocate(weight)
      allocate(weight(size(x%energy), x%bands))

      if (allocated(matsum)) deallocate(matsum)
      allocate(matsum(size(x%energy)))

      call differential(x%energy, weight(:, 1))

      do i = 2, x%bands
         weight(:, i) = weight(:, 1)
      end do

      weight(:, :) = weight * x%dos

      oc%states = sum(weight)
   end subroutine initialize
end module eliashberg_self_energy
