module eliashberg
   use global
   implicit none

contains

   subroutine solve(i, im)
      type(info), intent(inout) :: i
      type(matsubara), intent(out) :: im

      real(dp), allocatable :: E(:), lambda(:), mu(:)
      real(dp) :: omegaC, Z, Delta

      integer :: step, n, m, u, l
      logical :: done

      u = nint((i%upper / (pi * i%kT) - 1) / 2)
      l = nint((i%lower / (pi * i%kT) - 1) / 2)

      allocate(im%omega(0:u - 1))

      do n = 0, u - 1
         im%omega(n) = (2 * n + 1) * pi * i%kT
      end do

      allocate(lambda(1 - u:2 * u - 1))

      do n = 1 - u, 2 * u - 1
         lambda(n) = i%lambda / (1 + (2 * n * pi * i%kT / i%omegaE) ** 2)
      end do

      allocate(mu(0:u - 1))

      omegaC = (2 * l + 1) * pi * i%kT

      im%muStar = i%muStar / (1 + i%muStar * log(i%omegaE / omegaC))

      mu(:l - 1) = -2 * im%muStar
      mu(l:) = 0

      allocate(im%Z(0:u - 1))
      allocate(im%Delta(0:u - 1))

      im%Z(:) = 1
      im%Delta(:) = 1

      allocate(E(0:u - 1))

      E(:) = sqrt(1 + im%omega ** 2)

      im%status = -1

      do step = 1, i%limit
         done = .true.

         do n = 0, u - 1
            Z = 0
            Delta = 0

            do m = 0, u - 1
               Z = Z + im%omega(m) / E(m) &
                  * (lambda(n - m) - lambda(n + m + 1))

               Delta = Delta + im%Delta(m) / E(m) &
                  * (lambda(n - m) + lambda(n + m + 1) + mu(m))
            end do

            Z = 1 + pi * i%kT * Z / im%omega(n)
            Delta = pi * i%kT * Delta / Z

            if ((im%Z(n) .na. Z) .or. (im%Delta(n) .na. Delta)) done = .false.

            im%Z(n) = Z
            im%Delta(n) = Delta

            E(n) = sqrt(im%omega(n) ** 2 + Delta ** 2)
         end do

         if (done) then
            im%status = step
            exit
         end if
      end do

      if (im%Delta(0) .lt. 0) im%Delta = -im%Delta

      im%phiC = pi * i%kT * sum(im%Delta / E * mu)
   end subroutine solve
end module eliashberg
