! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module eliashberg_self_energy_real_axis
   use eliashberg_self_energy, only: initialize, weight
   use eliashberg_spectral_function, &
      only: initialize_a2F => initialize, weight_a2F => weight
   use global
   use tools, only: interval
   implicit none

   private
   public :: self_energy_real_axis

contains

   subroutine self_energy_real_axis(x, im, re, oc)
      type(parameters), intent(inout) :: x
      type(matsubara), intent(out) :: im
      type(continued), intent(out) :: re
      type(occupancy), intent(out) :: oc

      integer :: step, n, m
      real(dp), parameter :: xmax = log(huge(1.0_dp) / 2.0_dp - 1.0_dp)
      real(dp) :: prefactor, beta, fermi, dosef
      real(dp), allocatable :: bose(:), spec(:), w1(:), w2(:), n1(:), n2(:)
      complex(dp), allocatable :: omega(:), G0(:), G(:), Sigma(:), c1(:), c2(:)

      character(:), allocatable :: absent

      if (.not. x%normal) then
         absent = 'normal state'
      else if (x%bands .ne. 1) then
         absent = 'single band (for now)'
      else if (x%resolution .le. 0) then
         absent = 'finite resolution'
      else if (.not. x%ldos) then
         absent = 'density of states'
      else if (.not. x%la2F) then
         absent = 'Eliashberg spectral function'
      else if (x%n .ge. 0.0_dp .or. (x%mu .na. 0.0_dp)) then
         absent = 'constant Fermi level at zero (for now)'
      else if (any(x%muStar .na. 0.0_dp)) then
         absent = 'Coulomb pseudo-potential of zero (for now)'
      else
         absent = 'none'
      end if

      if (absent .ne. 'none') then
         print "('Error: Real-axis GW requires ', A)", absent
         stop 1
      end if

      call initialize(x, oc)
      call initialize_a2F(x)

      allocate(im%omega(0))
      allocate(im%Z(0, x%bands))
      allocate(im%phi(0, x%bands))
      allocate(im%chi(0, x%bands))
      allocate(im%Delta(0, x%bands))
      allocate(im%phiC(x%bands))

      allocate(re%omega(x%resolution))
      allocate(re%Z(x%resolution, x%bands))
      allocate(re%chi(x%resolution, x%bands))
      allocate(re%Delta(x%resolution, x%bands))
      allocate(re%dos(x%resolution, x%bands))

      allocate(omega(x%resolution))
      allocate(G0(x%resolution))
      allocate(G(x%resolution))
      allocate(Sigma(x%resolution))

      beta = 1.0_dp / (kB * x%T)

      if (x%divdos) then
         dosef = x%dos(minloc(abs(x%energy), 1), 1)
      else
         dosef = 1.0_dp
      end if

      bose = bose_fun(x%omega)

      call interval(re%omega, x%lower, x%upper, lower=.true., upper=.true.)

      prefactor = -(re%omega(2) - re%omega(1)) / pi

      omega = cmplx(re%omega, x%eta, dp)

      do n = 1, x%resolution
         G0(n) = sum(weight(:, 1) / (omega(n) - x%energy))
      end do

      oc%mu0 = x%mu
      oc%n0 = 2 * prefactor * sum(aimag(G0) * fermi_fun(re%omega))

      do step = 1, x%limit
         if (x%tell) print "('GW iteration ', I0)", step

         Sigma(:) = (0.0_dp, 0.0_dp)

         do m = 1, x%resolution
            call prepare(m)

            do n = 1, x%resolution
               Sigma(n) = Sigma(n) &
                  + sum(n1 / (omega(n) + w1) + n2 / (omega(n) + w2))
            end do
         end do

         Sigma(:) = prefactor / dosef * Sigma

         do n = 1, x%resolution
            G(n) = sum(weight(:, 1) / (omega(n) - x%energy - Sigma(n)))
         end do

         if (all(abs(G0 - G) .ap. 0.0_dp)) exit

         G0(:) = G
      end do

      oc%mu = x%mu
      oc%n = 2 * prefactor * sum(aimag(G) * fermi_fun(re%omega))

      im%phiC(:) = 0.0_dp

      re%Z(:, 1) = (0.0_dp, 0.0_dp)
      re%chi(:, 1) = (0.0_dp, 0.0_dp)
      re%Delta(:, 1) = (0.0_dp, 0.0_dp)

      do m = 1, x%resolution
         call prepare(m)

         do n = 1, x%resolution
            c1 = n1 / (w1 ** 2 - omega(n) ** 2)
            c2 = n2 / (w2 ** 2 - omega(n) ** 2)

            re%Z(n, 1) = re%Z(n, 1) - sum(omega(n) * (c1 + c2))
            re%chi(n, 1) = re%chi(n, 1) + sum(w1 * c1 + w2 * c2)
         end do
      end do

      re%Z(:, 1) = prefactor / dosef * re%Z(:, 1)
      re%chi(:, 1) = prefactor / dosef * re%chi(:, 1)

      re%Z(:, 1) = 1.0_dp - re%Z(:, 1) / omega

      re%dos(:, 1) = -aimag(G) / pi

   contains

      subroutine prepare(m)
         integer, intent(in) :: m

         fermi = fermi_fun(re%omega(m))

         spec = weight_a2F(:, 1, 1) * aimag(G0(m))

         w1 = -re%omega(m) - x%omega
         w2 = -re%omega(m) + x%omega

         n1 = spec * (1.0_dp - fermi + bose)
         n2 = spec * (fermi + bose)
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
