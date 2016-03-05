program tc
   use arguments
   use eliashberg
   use global
   use io
   use mcmillan
   implicit none

   type(info) :: i

   real(dp) :: small = 1e-10_dp
   real(dp) :: error = 1e-10_dp
   real(dp) :: ratio = 0.1_dp

   character(:), allocatable :: arg, lhs, rhs

   logical :: exists

   real(dp) :: lower, upper

   integer :: n, m

   do n = 1, command_argument_count()
      arg = argument(n)

      inquire(file=arg, exist=exists)

      if (exists) then
         call load(arg, i)

         call estimate(i)

         i%kT = i%Tc * kB / qe

         call solve(i)

         if (i%Delta(0) .gt. small) then
            lower = i%kT

            do while (i%Delta(0) .gt. small)
               i%kT = i%kT * (1 + ratio)
               call solve(i)
            end do

            upper = i%kT
         else
            upper = i%kT

            do while (i%Delta(0) .le. small)
               i%kT = i%kT * (1 - ratio)
               call solve(i)
            end do

            lower = i%kT
         end if

         do
            i%kT = (lower + upper) / 2

            if ((upper - lower) / 2 .le. error) exit

            call solve(i)

            if (i%Delta(0) .gt. small) then
               lower = i%kT
            else
               upper = i%kT
            end if
         end do

         write (*, "(ES22.14E3, ' (McMillan)')") i%Tc
         write (*, "(ES22.14E3, ' (Eliashberg)')") i%kT * qe / kB
      else
         m = index(arg, '=')

         if (m .ne. 0) then
            lhs = arg(:m - 1)
            rhs = arg(m + 1:)

            select case (lhs)
               case ('small'); read (rhs, *) small
               case ('error'); read (rhs, *) error
               case ('ratio'); read (rhs, *) ratio
            end select
         end if
      end if
   end do
end program tc
