module eliashberg
   use global
   implicit none

   private
   public :: solve

   integer :: step, n, m
   logical :: done

contains

   subroutine solve(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(out) :: im

      real(dp) :: omegaC

      im%n = nint((i%upper / (pi * i%kT) - 1) / 2)
      im%m = nint((i%lower / (pi * i%kT) - 1) / 2)

      allocate(im%omega(0:im%n - 1))

      do n = 0, im%n - 1
         im%omega(n) = (2 * n + 1) * pi * i%kT
      end do

      allocate(im%lambda(1 - im%n:2 * im%n - 1))

      do n = 1 - im%n, 2 * im%n - 1
         im%lambda(n) = i%lambda / (1 + (2 * n * pi * i%kT / i%omegaE) ** 2)
      end do

      allocate(im%mu(0:im%n - 1))

      omegaC = (2 * im%m + 1) * pi * i%kT

      im%muStar = i%muStar / (1 + i%muStar * log(i%omegaE / omegaC))

      im%mu(:im%m - 1) = -2 * im%muStar
      im%mu(im%m:) = 0

      if (i%DOS) then
         call variableDOS(i, im)
      else
         call constantDOS(i, im)
      end if
   end subroutine solve

   subroutine constantDOS(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(inout) :: im

      real(dp), allocatable :: E(:)
      real(dp) :: Z, Delta

      allocate(im%Z(0:im%n - 1))
      allocate(im%Delta(0:im%n - 1))

      im%Z(:) = 1
      im%Delta(:) = 1

      allocate(E(0:im%n - 1))

      E(:) = sqrt(1 + im%omega ** 2)

      im%status = -1

      do step = 1, i%limit
         done = .true.

         do n = 0, im%n - 1
            Z = 0
            Delta = 0

            do m = 0, im%n - 1
               Z = Z + im%omega(m) / E(m) &
                  * (im%lambda(n - m) - im%lambda(n + m + 1))

               Delta = Delta + im%Delta(m) / E(m) &
                  * (im%lambda(n - m) + im%lambda(n + m + 1) + im%mu(m))
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

      im%phiC = pi * i%kT * sum(im%Delta / E * im%mu)
   end subroutine constantDOS

   subroutine variableDOS(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(inout) :: im

      real(dp), allocatable :: A(:), B(:), trapezoids(:)
      real(dp) :: Z, phi, chi

      allocate(im%Z(0:im%n - 1))
      allocate(im%phi(0:im%n - 1))
      allocate(im%chi(0:im%n - 1))

      im%Z(:) = 1
      im%phi(:) = 1
      im%chi(:) = 0

      allocate(A(0:im%n - 1))
      allocate(B(0:im%n - 1))

      allocate(trapezoids(size(i%energy)))

      do n = 0, im%n - 1
         trapezoids(:) = i%weight / (im%omega(n) ** 2 + i%energy ** 2 + 1)

         A(n) = sum(trapezoids)
         B(n) = sum(trapezoids * i%energy)
      end do

      im%status = -1

      do step = 1, i%limit
         done = .true.

         do n = 0, im%n - 1
            Z = 0
            phi = 0
            chi = 0

            do m = 0, im%n - 1
               Z = Z + im%omega(m) * im%Z(m) * A(m) &
                  * (im%lambda(n - m) - im%lambda(n + m + 1))

               phi = phi + im%phi(m) * A(m) &
                  * (im%lambda(n - m) + im%lambda(n + m + 1) + im%mu(m))

               chi = chi - (im%chi(m) * A(m) + B(m)) &
                  * (im%lambda(n - m) + im%lambda(n + m + 1))
            end do

            Z = 1 + Z * i%kT / im%omega(n)
            phi = phi * i%kT
            chi = chi * i%kT

            if ((im%Z(n) .na. Z) &
               .or. (im%phi(n) .na. phi) &
               .or. (im%chi(n) .na. chi)) &
               done = .false.

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

      im%Delta = im%phi
      im%Delta = im%Delta / im%Z

      if (im%Delta(0) .lt. 0) im%Delta = -im%Delta

      im%phiC = i%kT * sum(im%phi * A * im%mu)
   end subroutine variableDOS
end module eliashberg
