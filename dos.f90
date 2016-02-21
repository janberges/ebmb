program dos
   use arguments
   use global
   use intervals
   implicit none

   real(dp) ::  E = 0.500_dp ! energy
   real(dp) :: dE = 0.007_dp ! half of energy-shell thickness
   real(dp) ::  t = 0.250_dp ! hopping parameter

   integer :: N = 1000 ! number of unit cells per dimension

   real(dp) :: epsilon
   real(dp), allocatable :: k(:)

   integer :: i, j, count

   character(:), allocatable :: setting, lhs, rhs

   do i = 1, command_argument_count()
      setting = argument(i)

      j = index(setting, '=')

      lhs = setting(:j - 1)
      rhs = setting(j + 1:)

      select case (lhs)
         case ( 'E'); read (rhs, *) E
         case ('dE'); read (rhs, *) dE
         case ( 't'); read (rhs, *) t
         case ( 'N'); read (rhs, *) N
      end select
   end do

   allocate(k(N))

   call interval(k, -pi, pi, lower=.true., upper=.false.)

   count = 0
   do i = 1, N
      do j = 1, N
         epsilon = -2 * t * (cos(k(i)) + cos(k(j)) - 2)
         if (E - dE .le. epsilon .and. epsilon .lt. E + dE) count = count + 1
      end do
   end do

   write (*, '(ES23.14E3)') count / (2 * dE * N ** 2)
end program dos
