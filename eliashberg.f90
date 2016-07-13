module eliashberg
   use global
   implicit none

contains

   subroutine solve(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(out) :: im

      real(dp), allocatable :: lambda(:), mu(:), A(:), B(:), trapezoids(:)
      real(dp) :: omegaC, Z, Delta, phi, chi

      integer :: step, n, m
      logical :: done

      im%u = ceiling((i%upper / (pi * i%T) - 1) / 2)
      im%l = ceiling((i%lower / (pi * i%T) - 1) / 2)

      allocate(im%omega(0:im%u - 1))

      do n = 0, im%u - 1
         im%omega(n) = (2 * n + 1) * pi * i%T
      end do

      allocate(lambda(1 - im%u:2 * im%u - 1))

      do n = 1 - im%u, 2 * im%u - 1
         lambda(n) = i%lambda / (1 + (2 * n * pi * i%T / i%omegaE) ** 2)
      end do

      allocate(mu(0:im%u - 1))

      omegaC = (2 * im%l + 1) * pi * i%T

      if (i%rescale) then
         im%muStar = i%muStar / (1 + i%muStar * log(i%omegaE / omegaC))
      else
         im%muStar = i%muStar
      end if

      mu(:im%l - 1) = -2 * im%muStar
      mu(im%l:) = 0

      allocate(im%Z(0:im%u - 1))
      allocate(im%Delta(0:im%u - 1))

      im%Z(:) = 1

      allocate(A(0:im%u - 1))

      im%status = -1

      if (i%DOS) then
         allocate(im%phi(0:im%u - 1))
         allocate(im%chi(0:im%u - 1))

         im%phi(:) = 1
         im%chi(:) = 0

         allocate(B(0:im%u - 1))

         allocate(trapezoids(size(i%energy)))

         do n = 0, im%u - 1
            trapezoids(:) = i%weight / (im%omega(n) ** 2 + i%energy ** 2 + 1)

            A(n) = sum(trapezoids)
            B(n) = sum(trapezoids * i%energy)
         end do

         do step = 1, i%limit
            done = .true.

            do n = 0, im%u - 1
               Z = 0
               phi = 0
               chi = 0

               do m = 0, im%u - 1
                  Z = Z + im%omega(m) * im%Z(m) * A(m) &
                     * (lambda(n - m) - lambda(n + m + 1))

                  phi = phi + im%phi(m) * A(m) &
                     * (lambda(n - m) + lambda(n + m + 1) + mu(m))

                  chi = chi - (im%chi(m) * A(m) + B(m)) &
                     * (lambda(n - m) + lambda(n + m + 1))
               end do

               Z = 1 + Z * i%T / im%omega(n)
               phi = phi * i%T
               chi = chi * i%T

               done = done &
                  .and. (im%Z(n) .ap. Z) &
                  .and. (im%phi(n) .ap. phi) &
                  .and. (im%chi(n) .ap. chi)

               im%Z(n) = Z
               im%phi(n) = phi
               im%chi(n) = chi

               trapezoids(:) = i%weight &
                  / ((im%omega(n) * Z) ** 2 + (i%energy + chi) ** 2 + phi ** 2)

               A(n) = sum(trapezoids)
               B(n) = sum(trapezoids * i%energy)
            end do

            if (done) then
               im%status = step
               exit
            end if
         end do

         if (im%phi(0) .lt. 0) im%phi = -im%phi

         im%Delta(:) = im%phi / im%Z

         im%phiC = i%T * sum(im%phi * A * mu)
      else
         im%Delta(:) = 1

         A(:) = 1 / sqrt(1 + im%omega ** 2)

         do step = 1, i%limit
            done = .true.

            do n = 0, im%u - 1
               Z = 0
               Delta = 0

               do m = 0, im%u - 1
                  Z = Z + im%omega(m) * A(m) &
                     * (lambda(n - m) - lambda(n + m + 1))

                  Delta = Delta + im%Delta(m) * A(m) &
                     * (lambda(n - m) + lambda(n + m + 1) + mu(m))
               end do

               Z = 1 + pi * i%T * Z / im%omega(n)
               Delta = pi * i%T * Delta / Z

               done = done &
                  .and. (im%Z(n) .ap. Z) &
                  .and. (im%Delta(n) .ap. Delta)

               im%Z(n) = Z
               im%Delta(n) = Delta

               A(n) = 1 / sqrt(im%omega(n) ** 2 + Delta ** 2)
            end do

            if (done) then
               im%status = step
               exit
            end if
         end do

         if (im%Delta(0) .lt. 0) im%Delta = -im%Delta

         im%phiC = pi * i%T * sum(im%Delta * A * mu)
      end if
   end subroutine solve
end module eliashberg
