module eliashberg_variable_dos
   use global
   implicit none

contains

   subroutine solve_variable_dos(x, im)
      type(universal), intent(in) :: x
      type(matsubara), intent(inout) :: im

      real(dp) :: nE, Z, phi, chi

      real(dp), allocatable :: lambda(:, :, :), mu(:, :, :)
      real(dp), allocatable :: muStar(:, :), A(:, :), B(:, :)

      real(dp), allocatable, save :: weight(:, :), trapezoids(:)

      integer, save :: u0 = -1

      integer :: step, i, j, n, m, u, l
      logical :: done

      nE = x%omegaE / (2 * pi * kB * x%T)

      u = ceiling(x%upper * nE - 0.5_dp)
      l = ceiling(x%lower * nE - 0.5_dp)

      if (u .ne. u0) then
         if (u0 .eq. -1) then
            allocate(weight(size(x%energy), x%bands))
            allocate(trapezoids(size(x%energy)))

            call differential(x%energy, weight(:, 1))

            do i = 2, x%bands
               weight(:, i) = weight(:, 1)
            end do

            n = minloc(abs(x%energy), 1)

            do i = 1, x%bands
               weight(:, i) = weight(:, i) * x%dos(:, i) / x%dos(n, i)
            end do
         else
            deallocate(im%omega)
            deallocate(im%Z)
            deallocate(im%phi)
            deallocate(im%chi)
            deallocate(im%Delta)
            deallocate(im%phiC)
         end if

         allocate(im%omega(0:u - 1))
         allocate(im%Z    (0:u - 1, x%bands))
         allocate(im%phi  (0:u - 1, x%bands))
         allocate(im%chi  (0:u - 1, x%bands))
         allocate(im%Delta(0:u - 1, x%bands))

         allocate(im%phiC(x%bands))

         im%Z(:, :) = 1

         im%phi(:, :) = 0
         im%phi(0, :) = 1
         im%chi(:, :) = 0

         im%Delta(:, :) = 0
         im%Delta(0, :) = 1

         u0 = u
      end if

      do n = 0, u - 1
         im%omega(n) = (2 * n + 1) * pi * kB * x%T
      end do

      allocate(lambda(1 - u:2 * u - 1, x%bands, x%bands))

      do n = 1 - u, 2 * u - 1
         lambda(n, :, :) = x%lambda / (1 + (n / nE) ** 2)
      end do

      allocate(muStar(x%bands, x%bands))

      if (x%rescale) then
         muStar = x%muStar / (1 + x%muStar * log(nE / (l + 0.5_dp)))
      else
         muStar = x%muStar
      end if

      allocate(mu(0:u - 1, x%bands, x%bands))

      do n = 0, l - 1
         mu(n, :, :) = -2 * muStar
      end do

      mu(l:, :, :) = 0

      allocate(A(0:u - 1, x%bands))
      allocate(B(0:u - 1, x%bands))

      do i = 1, x%bands
         do n = 0, u - 1
            call integrate(n, i)
         end do
      end do

      im%status = -1

      do step = 1, x%limit
         done = .true.

         do i = 1, x%bands
            do n = 0, u - 1
               Z = 0
               phi = 0
               chi = 0

               do j = 1, x%bands
                  do m = 0, u - 1
                     Z = Z + im%omega(m) * im%Z(m, j) * A(m, j) &
                        * (lambda(n - m, j, i) - lambda(n + m + 1, j, i))

                     phi = phi + im%phi(m, j) * A(m, j) * (mu(m, j, i) &
                        +  lambda(n - m, j, i) + lambda(n + m + 1, j, i))

                     chi = chi - (im%chi(m, j) * A(m, j) + B(m, j)) &
                        * (lambda(n - m, j, i) + lambda(n + m + 1, j, i))
                  end do
               end do

               Z = 1 + Z * kB * x%T / im%omega(n)
               phi = phi * kB * x%T
               chi = chi * kB * x%T

               done = done &
                  .and. (im%Z(n, i) .ap. Z) &
                  .and. (im%phi(n, i) .ap. phi) &
                  .and. (im%chi(n, i) .ap. chi)

               im%Z(n, i) = Z
               im%phi(n, i) = phi
               im%chi(n, i) = chi

               call integrate(n, i)
            end do
         end do

         if (done) then
            im%status = step
            exit
         end if
      end do

      im%Delta(:, :) = im%phi / im%Z

      do i = 1, x%bands
         im%phiC(i) = kB * x%T * sum(im%phi * A * mu(:, :, i))
      end do

   contains

      subroutine integrate(n, i)
         integer, intent(in) :: n, i

         trapezoids(:) = weight(:, i) / ((im%omega(n) * im%Z(n, i)) ** 2 &
            + (x%energy + im%chi(n, i)) ** 2 + im%phi(n, i) ** 2)

         A(n, i) = sum(trapezoids)
         B(n, i) = sum(trapezoids * x%energy)
      end subroutine integrate

   end subroutine solve_variable_dos

   subroutine differential(x, dx)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: dx(:)

      integer :: n
      n = size(x)

      dx(1) = x(2) - x(1)
      dx(2:n - 1) = x(3:n) - x(1:n - 2)
      dx(n) = x(n) - x(n - 1)

      dx(:) = dx / 2
   end subroutine differential
end module eliashberg_variable_dos
