! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module eliashberg_self_energy_real_axis
   use eliashberg_self_energy, &
      only: initialize_dos => initialize, weight_dos => weight
   use eliashberg_spectral_function, &
      only: initialize_a2F => initialize, weight_a2F => weight
   use global
   use tools, only: differential, interval
   implicit none

   private
   public :: self_energy_real_axis

contains

   subroutine self_energy_real_axis(x, im, re, oc)
      type(parameters), intent(inout) :: x
      type(matsubara), intent(out) :: im
      type(continued), intent(out) :: re
      type(occupancy), intent(out) :: oc

      integer :: step, i, j, n, m, no
      real(dp), parameter :: xmax = log(huge(1.0_dp) / 2.0_dp - 1.0_dp)
      real(dp) :: beta, fermi, domega, const
      real(dp), allocatable :: weight(:), bose(:), dosef(:)
      real(dp), allocatable :: w1(:), w2(:), n1(:, :), n2(:, :), r1(:), r2(:)
      complex(dp), allocatable :: omega(:), G0(:, :), G(:, :), Sigma(:, :)
      complex(dp), allocatable :: c1(:), c2(:)

      character(:), allocatable :: absent

      if (.not. x%normal) then
         absent = 'normal state'
      else if (x%points .le. 0) then
         absent = 'frequency points'
      else if (.not. x%ldos) then
         absent = 'density of states'
      else if (.not. x%la2F) then
         absent = 'Eliashberg spectral function'
      else if (.not. x%chi) then
         absent = 'considering energy shift'
      else
         absent = 'none'
      end if

      if (absent .ne. 'none') then
         print "('Error: Real-axis GW requires ', A)", absent
         stop 1
      end if

      call initialize_dos(x, oc)
      call initialize_a2F(x)

      domega = 2 * pi * kB * x%T

      no = ceiling(x%cutoff * x%omegaE / domega - 0.5_dp)

      allocate(im%omega(0:no - 1))
      allocate(im%Z(0:no - 1, x%bands))
      allocate(im%phi(0:no - 1, x%bands))
      allocate(im%chi(0:no - 1, x%bands))
      allocate(im%Delta(0:no - 1, x%bands))
      allocate(im%chiC(x%bands))
      allocate(im%phiC(x%bands))

      allocate(re%omega(x%points))
      allocate(re%Z(x%points, x%bands))
      allocate(re%chi(x%points, x%bands))
      allocate(re%Delta(x%points, x%bands))
      allocate(re%dos(x%points, x%bands))

      allocate(weight(x%points))
      allocate(omega(x%points))
      allocate(G(x%points, x%bands))
      allocate(Sigma(x%points, x%bands))
      allocate(n1(size(x%omega), x%bands))
      allocate(n2(size(x%omega), x%bands))
      allocate(dosef(x%bands))

      beta = 1.0_dp / (kB * x%T)

      bose = bose_fun(x%omega)

      do n = 0, no - 1
         im%omega(n) = domega * (n + 0.5_dp)
      end do

      call interval(re%omega, x%lower, x%upper, lower=.true., upper=.true., &
         logscale=x%logscale)

      call differential(re%omega, weight)

      weight = -weight / pi

      omega = cmplx(re%omega, x%eta, dp)

      if (x%n .ge. 0) then
         oc%mu = (x%energy(1) * (2 * oc%states - x%n) &
            + x%energy(size(x%energy)) * x%n) / (2 * oc%states)
      else
         oc%mu = x%mu
      end if

      Sigma(:, :) = (0.0_dp, 0.0_dp)

      call dos(x%n, x%n .ge. 0)

      oc%n0 = oc%n
      oc%mu0 = oc%mu

      if (x%divdos) then
         dosef(:) = x%dos(minloc(abs(x%energy - oc%mu), 1), :)
      else
         dosef(:) = 1.0_dp
      end if

      do step = 1, x%limit
         if (x%tell) print "('GW iteration ', I0)", step

         Sigma(:, :) = (0.0_dp, 0.0_dp)

         do j = 1, x%bands
            do m = 1, x%points
               call prepare(m, j)

               do i = 1, x%bands
                  do n = 1, x%points
                     Sigma(n, i) = Sigma(n, i) + sum( &
                        n1(:, i) / (omega(n) + w1) + n2(:, i) / (omega(n) + w2))
                  end do

                  Sigma(:, i) = Sigma(:, i) + weight(m) * aimag(G(m, j)) &
                     * (0.5_dp - fermi) * x%muStar(j, i) / dosef(j)
               end do
            end do
         end do

         G0 = G

         call dos(oc%n0, x%conserve)

         if (all(abs(G0 - G) .ap. 0.0_dp)) exit
      end do

      im%Z(:, :) = 0.0_dp
      im%phi(:, :) = 0.0_dp
      im%chi(:, :) = 0.0_dp
      im%Delta(:, :) = 0.0_dp
      im%chiC(:) = 0.0_dp
      im%phiC(:) = 0.0_dp

      re%Z(:, :) = (0.0_dp, 0.0_dp)
      re%chi(:, :) = (0.0_dp, 0.0_dp)
      re%Delta(:, :) = (0.0_dp, 0.0_dp)

      do j = 1, x%bands
         do m = 1, x%points
            call prepare(m, j)

            do i = 1, x%bands
               do n = 0, no - 1
                  r1 = n1(:, i) / (w1 ** 2 + im%omega(n) ** 2)
                  r2 = n2(:, i) / (w2 ** 2 + im%omega(n) ** 2)

                  im%Z(n, i) = im%Z(n, i) - sum(im%omega(n) * (r1 + r2))
                  im%chi(n, i) = im%chi(n, i) + sum(w1 * r1 + w2 * r2)
               end do

               do n = 1, x%points
                  c1 = n1(:, i) / (w1 ** 2 - omega(n) ** 2)
                  c2 = n2(:, i) / (w2 ** 2 - omega(n) ** 2)

                  re%Z(n, i) = re%Z(n, i) - sum(omega(n) * (c1 + c2))
                  re%chi(n, i) = re%chi(n, i) + sum(w1 * c1 + w2 * c2)
               end do

               const = weight(m) * aimag(G(m, j)) &
                  * (0.5_dp - fermi) * x%muStar(j, i) / dosef(j)

               im%chi(:, i) = im%chi(:, i) + const
               re%chi(:, i) = re%chi(:, i) + const

               im%chiC(i) = im%chiC(i) + const
            end do
         end do
      end do

      do i = 1, x%bands
         im%Z(:, i) = 1.0_dp - im%Z(:, i) / im%omega
         re%Z(:, i) = 1.0_dp - re%Z(:, i) / omega
      end do

      re%dos(:, :) = -aimag(G) / pi

   contains

      subroutine dos(ntarget, optimize)
         real(dp), intent(in) :: ntarget
         logical, intent(in) :: optimize

         real(dp) :: bell(x%points), w(x%points)

         where (re%omega .ap. 0.0_dp)
            bell = 0.5_dp * beta
         elsewhere
            bell = tanh(0.5_dp * re%omega * beta) / re%omega
         end where

         do
            do i = 1, x%bands
               do n = 1, x%points
                  G(n, i) = sum(weight_dos(:, i) &
                     / (omega(n) - x%energy + oc%mu - Sigma(n, i)))
               end do
            end do

            w(:) = weight * sum(aimag(G), 2)

            oc%n = 2 * sum(w * fermi_fun(re%omega))
            oc%inspect = sum(w)

            if (oc%inspect .ap. 0.0_dp) then
               print "('Error: Too small energy window has drifted away')"
               stop 1
            end if

            if ((oc%n .ap. ntarget) .or. .not. optimize) exit

            w(:) = w * bell

            oc%mu = oc%mu + (ntarget - oc%inspect + sum(w * re%omega)) / sum(w)
         end do

         if (abs(oc%inspect - oc%states) .gt. 0.1_dp) then
            print "('Warning: Spectral function breaks sum rule')"
         end if
      end subroutine dos

      subroutine prepare(m, j)
         integer, intent(in) :: m, j
         real(dp) :: spec(size(x%omega), x%bands)

         fermi = fermi_fun(re%omega(m))

         spec = weight(m) * weight_a2F(:, j, :) * aimag(G(m, j)) / dosef(j)

         w1 = -re%omega(m) - x%omega
         w2 = -re%omega(m) + x%omega

         do i = 1, x%bands
            n1(:, i) = spec(:, i) * (1.0_dp - fermi + bose)
            n2(:, i) = spec(:, i) * (fermi + bose)
         end do
      end subroutine prepare

      elemental function fermi_fun(omega)
         real(dp) :: fermi_fun
         real(dp), intent(in) :: omega

         fermi_fun = 1.0_dp / (exp(min(omega * beta, xmax)) + 1.0_dp)
      end function fermi_fun

      elemental function bose_fun(omega)
         real(dp) :: bose_fun
         real(dp), intent(in) :: omega

         bose_fun = 1.0_dp / (exp(min(omega * beta, xmax)) - 1.0_dp)
      end function bose_fun
   end subroutine self_energy_real_axis
end module eliashberg_self_energy_real_axis
