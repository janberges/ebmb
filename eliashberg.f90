module eliashberg
   use global
   implicit none

   private
   public :: solve

contains

   subroutine solve(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(out) :: im

      integer :: n, cutC
      real(dp) :: omegaC

      im%n = nint((i%upper / (pi * i%kT) - 1) / 2)
      cutC = nint((i%lower / (pi * i%kT) - 1) / 2)

      allocate(im%omega(0:im%n - 1))

      do n = 0, im%n - 1
         im%omega(n) = (2 * n + 1) * pi * i%kT
      end do

      allocate(im%lambda(1 - im%n:2 * im%n - 1))

      do n = 1 - im%n, 2 * im%n - 1
         im%lambda(n) = i%lambda / (1 + (2 * n * pi * i%kT / i%omegaE) ** 2)
      end do

      allocate(im%mu(0:im%n - 1))

      omegaC = (2 * cutC + 1) * pi * i%kT

      im%muStar = i%muStar / (1 + i%muStar * log(i%omegaE / omegaC))

      im%mu(:cutC - 1) = -2 * im%muStar
      im%mu(cutC:) = 0

      call constantDOS(i, im)
   end subroutine solve

   subroutine constantDOS(i, im)
      type(universal), intent(in) :: i
      type(matsubara), intent(inout) :: im

      real(dp), allocatable :: E(:)
      real(dp) :: Z, Delta

      integer :: step, n, m
      logical :: done

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
end module eliashberg
