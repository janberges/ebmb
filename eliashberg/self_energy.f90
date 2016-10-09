module eliashberg_self_energy
   use global
   use tools, only: differential
   implicit none

   private
   public :: self_energy, initialize

   logical :: initial = .true.

   real(dp) :: states
   real(dp), allocatable :: weight(:, :), trapezia(:), matsum(:)

contains

   subroutine self_energy(x, im, oc, kernel)
      type(parameters), intent(in) :: x
      type(matsubara), intent(out) :: im
      type(occupancy), intent(out) :: oc

      real(dp), allocatable, intent(out), optional :: kernel(:, :)

      real(dp) :: nE, Z, phi, chi, mu, domega, A0, B0, residue

      real(dp), allocatable :: g(:, :, :), U(:, :, :)
      real(dp), allocatable :: muStar(:, :), A(:, :), B(:, :)

      real(dp), allocatable :: integral_Z  (:, :)
      real(dp), allocatable :: integral_phi(:, :)
      real(dp), allocatable :: integral_chi(:, :)

      integer :: step, i, j, n, m, p, q, no, nC, f
      logical :: done

      if (initial) call initialize(x)

      if (0 .lt. x%n .and. x%n .lt. 2) then
         oc%n = x%n

         oc%mu &
            = (x%energy(1) * (2 - oc%n) + x%energy(size(x%energy)) * oc%n) / 2

         done = .false.

         do while (.not. done)
            where (x%energy .ap. oc%mu)
               matsum = 1 / (2 * kB * x%T)
            elsewhere
               matsum = x%energy - oc%mu
               matsum = tanh(matsum / (2 * kB * x%T)) / matsum
            end where

            A0 = 0
            B0 = 0

            do i = 1, x%bands
               trapezia(:) = weight(:, i) * matsum

               A0 = A0 + sum(trapezia)
               B0 = B0 + sum(trapezia * x%energy)
            end do

            mu = (oc%n - 1 + B0) / A0

            if (oc%mu .ap. mu) done = .true.

            oc%mu = mu
         end do
      else
         oc%mu = x%mu

         oc%n = 1

         matsum(:) = tanh((x%energy - x%mu) / (2 * kB * x%T))

         do i = 1, x%bands
            oc%n = oc%n - sum(weight(:, i) * matsum)
         end do
      end if

      oc%n0  = oc%n
      oc%mu0 = oc%mu

      f = minloc(abs(x%energy - oc%mu), 1)

      domega = 2 * pi * kB * x%T

      nE = x%omegaE / domega

      no = ceiling(x%cutoff  * nE - 0.5_dp)
      nC = ceiling(x%cutoffC * nE - 0.5_dp)

      allocate(im%omega(0:no - 1))

      do n = 0, no - 1
         im%omega(n) = domega * (n + 0.5_dp)
      end do

      allocate(g(1 - no:2 * no - 1, x%bands, x%bands))

      do n = 1 - no, 2 * no - 1
         g(n, :, :) = x%lambda / (1 + (n / nE) ** 2)

         do i = 1, x%bands
            g(n, :, i) = g(n, :, i) * states / x%dos(f, :)
         end do
      end do

      allocate(muStar(x%bands, x%bands))

      muStar(:, :) = x%muStar / (1 + x%muStar &
         * log(2 * x%omegaE / (x%energy(size(x%energy)) - x%energy(1))))

      if (x%rescale) then
         where (x%energy .ap. oc%mu)
            matsum = 1 / (domega * (nC + 0.5_dp))
         elsewhere
            matsum = x%energy - oc%mu
            matsum = atan(matsum / (domega * (nC + 0.5_dp))) / matsum
         end where

         residue = 0

         do i = 1, x%bands
            residue = residue + sum(weight(:, i) / x%dos(f, i) * matsum)
         end do

         residue = residue * states / pi

         muStar(:, :) = muStar / (1 + muStar * residue)
      end if

      allocate(U(0:no - 1, x%bands, x%bands))

      do n = 0, nC - 1
         U(n, :, :) = -2 * muStar

         do i = 1, x%bands
            U(n, :, i) = U(n, :, i) * states / x%dos(f, :)
         end do
      end do

      U(nC:, :, :) = 0

      allocate(im%Z(0:no - 1, x%bands))

      im%Z(:, :) = 1

      allocate(im%phi(0:no - 1, x%bands))

      im%phi(:, :) = 0

      if (.not. (x%normal .or. present(kernel))) im%phi(0, :) = 1

      allocate(im%chi(0:no - 1, x%bands))

      im%chi(:, :) = 0

      allocate(A(0:no - 1, x%bands))
      allocate(B(0:no - 1, x%bands))

      allocate(integral_Z  (0:no - 1, x%bands))
      allocate(integral_phi(0:no - 1, x%bands))
      allocate(integral_chi(0:no - 1, x%bands))

      do i = 1, x%bands
         do n = 0, no - 1
            call integrate(n, i)
         end do
      end do

      im%status = -1

      do step = 1, x%limit
         done = .true.

         do i = 1, x%bands
            do n = 0, no - 1
               Z = 0
               phi = 0
               chi = 0

               do j = 1, x%bands
                  do m = 0, no - 1
                     Z = Z + integral_Z(m, j) &
                        * (g(n - m, j, i) - g(n + m + 1, j, i))

                     phi = phi + integral_phi(m, j) &
                        * (g(n - m, j, i) + g(n + m + 1, j, i) + U(m, j, i))

                     chi = chi - integral_chi(m, j) &
                        * (g(n - m, j, i) + g(n + m + 1, j, i))
                  end do
               end do

               Z = 1 + Z * kB * x%T / im%omega(n)
               phi = phi * kB * x%T
               chi = chi * kB * x%T

               done = done &
                  .and. (im%Z  (n, i) .ap. Z) &
                  .and. (im%phi(n, i) .ap. phi) &
                  .and. (im%chi(n, i) .ap. chi)

               im%Z  (n, i) = Z
               im%phi(n, i) = phi
               im%chi(n, i) = chi

               call integrate(n, i)
            end do
         end do

         if (x%conserve) then
            matsum(:) = atan((x%energy - oc%mu) / (domega * (no + 0.5_dp))) &
               / domega

            residue = 0

            do i = 1, x%bands
               residue = residue + sum(weight(:, i) * matsum)
            end do

            mu = ((oc%n - 1) / (4 * kB * x%T) + sum(A * im%chi + B) + residue) &
               / sum(A)

            done = done .and. (oc%mu .ap. mu)

            oc%mu = mu
         end if

         if (done) then
            im%status = step
            exit
         end if
      end do

      allocate(im%Delta(0:no - 1, x%bands))

      im%Delta(:, :) = im%phi / im%Z

      allocate(im%phiC(x%bands))

      do i = 1, x%bands
         im%phiC(i) = kB * x%T * sum(integral_phi * U(:, :, i))
      end do

      oc%n = 1 - 4 * kB * x%T * (sum(integral_chi) + residue)

      if (present(kernel)) then
         allocate(kernel(x%bands * no, x%bands * no))

         do i = 1, x%bands
            p = i * no
            do j = 1, x%bands
               q = j * no
               do n = 0, no - 1
                  do m = 0, no - 1
                     kernel(q - m, p - n) = kB * x%T * A(m, j) &
                        * (g(n - m, j, i) + g(n + m + 1, j, i) + U(m, j, i))
                  end do
               end do
            end do
         end do
      end if

   contains

      subroutine integrate(n, i)
         integer, intent(in) :: n, i

         trapezia(:) = weight(:, i) / ((im%omega(n) * im%Z(n, i)) ** 2 &
            + (x%energy - oc%mu + im%chi(n, i)) ** 2 + im%phi(n, i) ** 2)

         A(n, i) = sum(trapezia)
         B(n, i) = sum(trapezia * x%energy)

         integral_Z  (n, i) = A(n, i) *  im%Z  (n, i) * im%omega(n)
         integral_phi(n, i) = A(n, i) *  im%phi(n, i)
         integral_chi(n, i) = A(n, i) * (im%chi(n, i) - oc%mu) + B(n, i)
      end subroutine integrate

   end subroutine self_energy

   subroutine initialize(x)
      type(parameters), intent(in) :: x

      integer :: i

      initial = .false.

      if (allocated(weight)) deallocate(weight)
      allocate(weight(size(x%energy), x%bands))

      if (allocated(trapezia)) deallocate(trapezia)
      allocate(trapezia(size(x%energy)))

      if (allocated(matsum)) deallocate(matsum)
      allocate(matsum(size(x%energy)))

      call differential(x%energy, weight(:, 1))

      do i = 2, x%bands
         weight(:, i) = weight(:, 1)
      end do

      weight(:, :) = weight * x%dos

      states = sum(weight)

      weight(:, :) = weight / states
   end subroutine initialize
end module eliashberg_self_energy
