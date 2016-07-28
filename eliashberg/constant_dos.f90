module eliashberg_constant_dos
   use global
   implicit none

contains

   subroutine solve_constant_dos(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(out) :: im

      real(dp) :: nE, Z, Delta

      real(dp), allocatable :: lambda(:, :, :), mu(:, :, :)
      real(dp), allocatable :: muStar(:, :), A(:, :)

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

      allocate(im%Delta(0:u - 1, i%bands))

      im%Delta(:, :) = 0
      im%Delta(0, :) = 1

      allocate(A(0:u - 1, i%bands))

      do p = 1, i%bands
         A(:, p) = 1 / sqrt(im%omega ** 2 + im%Delta(:, p) ** 2)
      end do

      im%status = -1

      do step = 1, i%limit
         done = .true.

         do p = 1, i%bands
            do n = 0, u - 1
               Z = 0
               Delta = 0

               do q = 1, i%bands
                  do m = 0, u - 1
                     Z = Z + im%omega(m) * A(m, q) &
                        * (lambda(n - m, q, p) - lambda(n + m + 1, q, p))

                     Delta = Delta + im%Delta(m, q) * A(m, q) * (mu(m, q, p) &
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

      allocate(im%phiC(i%bands))

      do p = 1, i%bands
         im%phiC(p) = pi * i%T * sum(im%Delta * A * mu(:, :, p))
      end do
   end subroutine solve_constant_dos
end module eliashberg_constant_dos
