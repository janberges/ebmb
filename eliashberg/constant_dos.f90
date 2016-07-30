module eliashberg_constant_dos
   use global
   implicit none

contains

   subroutine solve_constant_dos(x, im)
      type(universal), intent(in) :: x
      type(matsubara), intent(inout) :: im

      real(dp) :: nE, Z, Delta

      real(dp), allocatable :: lambda(:, :, :), mu(:, :, :)
      real(dp), allocatable :: muStar(:, :), A(:, :)

      integer, save :: u0 = -1

      integer :: step, i, j, n, m, u, l
      logical :: done

      nE = x%omegaE / (2 * pi * kB * x%T)

      u = ceiling(x%upper * nE - 0.5_dp)
      l = ceiling(x%lower * nE - 0.5_dp)

      if (u .ne. u0) then
         if (u0 .ne. -1) then
            deallocate(im%omega)
            deallocate(im%Z)
            deallocate(im%Delta)
            deallocate(im%phiC)
         end if

         allocate(im%omega(0:u - 1))
         allocate(im%Z    (0:u - 1, x%bands))
         allocate(im%Delta(0:u - 1, x%bands))

         allocate(im%phiC(x%bands))

         im%Z(:, :) = 1

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

      do i = 1, x%bands
         A(:, i) = 1 / sqrt(im%omega ** 2 + im%Delta(:, i) ** 2)
      end do

      im%status = -1

      do step = 1, x%limit
         done = .true.

         do i = 1, x%bands
            do n = 0, u - 1
               Z = 0
               Delta = 0

               do j = 1, x%bands
                  do m = 0, u - 1
                     Z = Z + im%omega(m) * A(m, j) &
                        * (lambda(n - m, j, i) - lambda(n + m + 1, j, i))

                     Delta = Delta + im%Delta(m, j) * A(m, j) * (mu(m, j, i) &
                        +  lambda(n - m, j, i) + lambda(n + m + 1, j, i))
                  end do
               end do

               Z = 1 + pi * kB * x%T * Z / im%omega(n)
               Delta = pi * kB * x%T * Delta / Z

               done = done &
                  .and. (im%Z(n, i) .ap. Z) &
                  .and. (im%Delta(n, i) .ap. Delta)

               im%Z(n, i) = Z
               im%Delta(n, i) = Delta

               A(n, i) = 1 / sqrt(im%omega(n) ** 2 + Delta ** 2)
            end do
         end do

         if (done) then
            im%status = step
            exit
         end if
      end do

      do i = 1, x%bands
         im%phiC(i) = pi * kB * x%T * sum(im%Delta * A * mu(:, :, i))
      end do
   end subroutine solve_constant_dos
end module eliashberg_constant_dos
