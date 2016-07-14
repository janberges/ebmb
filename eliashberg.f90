module eliashberg
   use global
   implicit none

contains

   subroutine solve(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(out) :: im

      real(dp), allocatable :: lambda(:, :, :), mu(:, :, :)
      real(dp), allocatable :: A(:, :), B(:, :), trapezoids(:)

      real(dp) :: omegaC, Z, Delta, phi, chi

      integer :: step, p, q, n, m
      logical :: done

      im%u = ceiling((i%upper / (pi * i%T) - 1) / 2)
      im%l = ceiling((i%lower / (pi * i%T) - 1) / 2)

      allocate(im%omega(0:im%u - 1))

      do n = 0, im%u - 1
         im%omega(n) = (2 * n + 1) * pi * i%T
      end do

      allocate(lambda(1 - im%u:2 * im%u - 1, i%bands, i%bands))

      do n = 1 - im%u, 2 * im%u - 1
         lambda(n, :, :) = i%lambda / (1 + (2 * n * pi * i%T / i%omegaE) ** 2)
      end do

      allocate(mu(0:im%u - 1, i%bands, i%bands))

      omegaC = (2 * im%l + 1) * pi * i%T

      if (i%rescale) then
         im%muStar = i%muStar / (1 + i%muStar * log(i%omegaE / omegaC))
      else
         im%muStar = i%muStar
      end if

      do n = 0, im%l - 1
         mu(n, :, :) = -2 * im%muStar
      end do

      mu(im%l:, :, :) = 0

      allocate(im%Z(0:im%u - 1, i%bands))
      allocate(im%Delta(0:im%u - 1, i%bands))
      allocate(im%phiC(i%bands))

      im%Z(:, :) = 1

      allocate(A(0:im%u - 1, i%bands))

      im%status = -1

      if (i%DOS) then
         allocate(im%phi(0:im%u - 1, i%bands))
         allocate(im%chi(0:im%u - 1, i%bands))

         im%phi(:, :) = 1
         im%chi(:, :) = 0

         allocate(B(0:im%u - 1, i%bands))

         allocate(trapezoids(size(i%energy)))

         do p = 1, i%bands
            do n = 0, im%u - 1
               call integrate(n, p)
            end do
         end do

         do step = 1, i%limit
            done = .true.

            do p = 1, i%bands
               do n = 0, im%u - 1
                  Z = 0
                  phi = 0
                  chi = 0

                  do q = 1, i%bands
                     do m = 0, im%u - 1
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

         if (all(im%phi(0, :) .lt. 0)) im%phi(:, :) = -im%phi

         im%Delta(:, :) = im%phi / im%Z

         do p = 1, i%bands
            im%phiC(p) = i%T * sum(im%phi * A * mu(:, :, p))
         end do
      else
         im%Delta(:, :) = 1

         do p = 1, i%bands
            A(:, p) = 1 / sqrt(im%omega ** 2 + im%Delta(:, p) ** 2)
         end do

         do step = 1, i%limit
            done = .true.

            do p = 1, i%bands
               do n = 0, im%u - 1
                  Z = 0
                  Delta = 0

                  do q = 1, i%bands
                     do m = 0, im%u - 1
                        Z = Z + im%omega(m) * A(m, q) &
                           * (lambda(n - m, q, p) - lambda(n + m + 1, q, p))

                        Delta = Delta + im%Delta(m, q) * A(m, q) * (mu(m, q, p)&
                           +  lambda(n - m, q, p) + lambda(n + m + 1, q, p))
                     end do
                  end do

                  Z = 1 + pi * i%T * Z / im%omega(n)
                  Delta = pi * i%T * Delta / Z

                  done = done &
                     .and. (im%Z(n, p) .ap. Z) &
                     .and. (im%Delta(n, p) .ap. Delta)

                  im%Z(n, p) = Z
                  im%Delta(n, p) = Delta

                  A(n, p) = 1 / sqrt(im%omega(n) ** 2 + Delta ** 2)
               end do
            end do

            if (done) then
               im%status = step
               exit
            end if
         end do

         if (all(im%Delta(0, :) .lt. 0)) im%Delta(:, :) = -im%Delta

         do p = 1, i%bands
            im%phiC(p) = pi * i%T * sum(im%Delta * A * mu(:, :, p))
         end do
      end if

   contains

      subroutine integrate(n, p)
         integer, intent(in) :: n, p

         trapezoids(:) = i%weight(:, p) / ((im%omega(n) * im%Z(n, p)) ** 2 &
            + (i%energy + im%chi(n, p)) ** 2 + im%phi(n, p) ** 2)

         A(n, p) = sum(trapezoids)
         B(n, p) = sum(trapezoids * i%energy)
      end subroutine integrate

   end subroutine solve
end module eliashberg
