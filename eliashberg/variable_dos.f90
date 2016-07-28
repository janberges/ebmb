module eliashberg_variable_dos
   use global
   implicit none

contains

   subroutine solve_variable_dos(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(out) :: im

      real(dp) :: nE, Z, phi, chi

      real(dp), allocatable :: lambda(:, :, :), mu(:, :, :)
      real(dp), allocatable :: muStar(:, :), A(:, :), B(:, :), trapezoids(:)

      integer :: step, p, q, n, m, u, l
      logical :: done

      nE = i%omegaE / (2 * pi * i%T)

      u = ceiling(i%upper * nE - 0.5_dp)
      l = ceiling(i%lower * nE - 0.5_dp)

      allocate(im%omega(0:u - 1))

      do n = 0, u - 1
         im%omega(n) = (2 * n + 1) * pi * i%T
      end do

      allocate(lambda(1 - u:2 * u - 1, i%bands, i%bands))

      do n = 1 - u, 2 * u - 1
         lambda(n, :, :) = i%lambda / (1 + (n / nE) ** 2)
      end do

      allocate(muStar(i%bands, i%bands))

      if (i%rescale) then
         muStar = i%muStar / (1 + i%muStar * log(nE / (l + 0.5_dp)))
      else
         muStar = i%muStar
      end if

      allocate(mu(0:u - 1, i%bands, i%bands))

      do n = 0, l - 1
         mu(n, :, :) = -2 * muStar
      end do

      mu(l:, :, :) = 0

      allocate(im%Z(0:u - 1, i%bands))

      im%Z(:, :) = 1

      allocate(im%phi(0:u - 1, i%bands))
      allocate(im%chi(0:u - 1, i%bands))

      im%phi(:, :) = 0
      im%phi(0, :) = 1
      im%chi(:, :) = 0

      allocate(A(0:u - 1, i%bands))
      allocate(B(0:u - 1, i%bands))

      allocate(trapezoids(size(i%energy)))

      do p = 1, i%bands
         do n = 0, u - 1
            call integrate(n, p)
         end do
      end do

      im%status = -1

      do step = 1, i%limit
         done = .true.

         do p = 1, i%bands
            do n = 0, u - 1
               Z = 0
               phi = 0
               chi = 0

               do q = 1, i%bands
                  do m = 0, u - 1
                     Z = Z + im%omega(m) * im%Z(m, q) * A(m, q) &
                        * (lambda(n - m, q, p) - lambda(n + m + 1, q, p))

                     phi = phi + im%phi(m, q) * A(m, q) * (mu(m, q, p) &
                        +  lambda(n - m, q, p) + lambda(n + m + 1, q, p))

                     chi = chi - (im%chi(m, q) * A(m, q) + B(m, q)) &
                        * (lambda(n - m, q, p) + lambda(n + m + 1, q, p))
                  end do
               end do

               Z = 1 + Z * i%T / im%omega(n)
               phi = phi * i%T
               chi = chi * i%T

               done = done &
                  .and. (im%Z(n, p) .ap. Z) &
                  .and. (im%phi(n, p) .ap. phi) &
                  .and. (im%chi(n, p) .ap. chi)

               im%Z(n, p) = Z
               im%phi(n, p) = phi
               im%chi(n, p) = chi

               call integrate(n, p)
            end do
         end do

         if (done) then
            im%status = step
            exit
         end if
      end do

      allocate(im%Delta(0:u - 1, i%bands))

      im%Delta(:, :) = im%phi / im%Z

      allocate(im%phiC(i%bands))

      do p = 1, i%bands
         im%phiC(p) = i%T * sum(im%phi * A * mu(:, :, p))
      end do

   contains

      subroutine integrate(n, p)
         integer, intent(in) :: n, p

         trapezoids(:) = i%weight(:, p) / ((im%omega(n) * im%Z(n, p)) ** 2 &
            + (i%energy + im%chi(n, p)) ** 2 + im%phi(n, p) ** 2)

         A(n, p) = sum(trapezoids)
         B(n, p) = sum(trapezoids * i%energy)
      end subroutine integrate

   end subroutine solve_variable_dos
end module eliashberg_variable_dos
