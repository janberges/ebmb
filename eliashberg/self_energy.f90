module eliashberg_self_energy
   use global
   use tools, only: differential
   implicit none

   private
   public :: self_energy, initialize

   logical :: initial = .true.

   real(dp), allocatable :: weight(:, :), trapezia(:), matsum(:)

contains

   subroutine self_energy(x, im)
      type(parameters), intent(in) :: x
      type(matsubara), intent(out) :: im

      real(dp) :: nE, Z, phi, chi, mu, domega, occupation, A0, B0

      real(dp), allocatable :: g(:, :, :), U(:, :, :)
      real(dp), allocatable :: muStar(:, :), A(:, :), B(:, :)

      real(dp), allocatable :: integral_Z  (:, :)
      real(dp), allocatable :: integral_phi(:, :)
      real(dp), allocatable :: integral_chi(:, :)

      integer :: step, i, j, n, m, no, nC, f
      logical :: done

      if (initial) call initialize(x)

      matsum(:) = 1 - tanh((x%energy - x%mu) / (2 * kB * x%T))

      occupation = 0

      do i = 1, x%bands
         occupation = occupation + sum(weight(:, i) * matsum)
      end do

      f = minloc(abs(x%energy - x%mu), 1)

      im%mu = x%energy(f)

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
            g(n, :, i) = g(n, :, i) / x%dos(f, :)
         end do
      end do

      allocate(muStar(x%bands, x%bands))

      if (x%rescale) then
         muStar(:, :) = x%muStar / (1 + x%muStar * log(nE / (nC + 0.5_dp)))
      else
         muStar(:, :) = x%muStar
      end if

      allocate(U(0:no - 1, x%bands, x%bands))

      do n = 0, nC - 1
         U(n, :, :) = -2 * muStar

         do i = 1, x%bands
            U(n, :, i) = U(n, :, i) / x%dos(f, :)
         end do
      end do

      U(nC:, :, :) = 0

      allocate(im%Z(0:no - 1, x%bands))

      im%Z(:, :) = 1

      allocate(im%phi(0:no - 1, x%bands))

      im%phi(:, :) = 0
      im%phi(0, :) = 1

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

         where (x%energy .ap. im%mu)
            matsum = 1 / ((no + 0.5_dp) * domega ** 2)
         elsewhere
            matsum = x%energy - im%mu
            matsum = atan(matsum / (domega * (no + 0.5_dp))) / (domega * matsum)
         end where

         A0 = 0
         B0 = 0

         do i = 1, x%bands
            trapezia(:) = weight(:, i) * matsum

            A0 = A0 + sum(trapezia)
            B0 = B0 + sum(trapezia * x%energy)
         end do

         mu = ((occupation - 1) / (4 * kB * x%T) + sum(A * im%chi + B) + B0) &
            / (sum(A) + A0)

         done = done .and. (im%mu .ap. mu)

         im%mu = mu

         if (done) then
            im%status = step
            exit
         end if
      end do

      allocate(im%Delta(0:no - 1, x%bands))

      im%Delta(:, :) = im%phi / im%Z

      allocate(im%phiC(x%bands))

      do i = 1, x%bands
         im%phiC(i) = kB * x%T * sum(im%phi * integral_phi * U(:, :, i))
      end do

   contains

      subroutine integrate(n, i)
         integer, intent(in) :: n, i

         trapezia(:) = weight(:, i) / ((im%omega(n) * im%Z(n, i)) ** 2 &
            + (x%energy - im%mu + im%chi(n, i)) ** 2 + im%phi(n, i) ** 2)

         A(n, i) = sum(trapezia)
         B(n, i) = sum(trapezia * x%energy)

         integral_Z  (n, i) = A(n, i) *  im%Z  (n, i) * im%omega(n)
         integral_phi(n, i) = A(n, i) *  im%phi(n, i)
         integral_chi(n, i) = A(n, i) * (im%chi(n, i) - im%mu) + B(n, i)
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
   end subroutine initialize
end module eliashberg_self_energy
