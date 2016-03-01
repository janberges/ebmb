module eliashberg
   use global
   implicit none

contains

   subroutine solve(i)
      type(info), intent(inout) :: i

      real(dp), allocatable :: E(:), lambda(:), mu(:)
      real(dp) :: Z, Delta

      integer :: step, n, m, u, l
      logical :: done

      if (i%lower .lt. 0) i%lower = i%upper

      u = nint((i%upper / (pi * i%kT) - 1) / 2)
      l = nint((i%lower / (pi * i%kT) - 1) / 2)

      i%upper = (2 * u + 1) * pi * i%kT
      i%lower = (2 * l + 1) * pi * i%kT

      allocate(i%omega(0:u - 1))

      do n = 0, u - 1
         i%omega(n) = (2 * n + 1) * pi * i%kT
      end do

      allocate(lambda(1 - u:2 * u - 1))

      do n = 1 - u, 2 * u - 1
         lambda(n) = i%lambda / (1 + (2 * n * pi * i%kT / i%omegaE) ** 2)
      end do

      allocate(mu(0:u - 1))

      mu(:l - 1) = -2 * i%muStar / (1 + i%muStar * log(i%omegaE / i%lower))
      mu(l:) = 0

      allocate(i%Z(0:u - 1))
      allocate(i%Delta(0:u - 1))

      i%Z(:) = 1
      i%Delta(:) = 1

      allocate(E(0:u - 1))

      E(:) = sqrt(1 + i%omega ** 2)

      i%status = -1

      do step = 1, i%limit
         done = .true.

         do n = 0, u - 1
            Z = 0
            Delta = 0

            do m = 0, u - 1
               Z = Z + i%omega(m) / E(m) &
                  * (lambda(n - m) - lambda(n + m + 1))

               Delta = Delta + i%Delta(m) / E(m) &
                  * (lambda(n - m) + lambda(n + m + 1) + mu(m))
            end do

            Z = 1 + pi * i%kT * Z / i%omega(n)
            Delta = pi * i%kT * Delta / Z

            if (abs(i%Delta(n) - Delta) .gt. i%tiny &
               .or. abs(i%Z(n) - Z) .gt. i%tiny) done = .false.

            i%Z(n) = Z
            i%Delta(n) = Delta

            E(n) = sqrt(i%omega(n) ** 2 + Delta ** 2)
         end do

         if (done) then
            i%status = step
            exit
         end if
      end do

      if (i%Delta(0) .lt. 0) i%Delta = -i%Delta

      i%phiC = pi * i%kT * sum(i%Delta / E * mu)
   end subroutine solve
end module eliashberg
