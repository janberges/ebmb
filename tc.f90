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
   real(dp) :: ratio = 2

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

         lower = i%kT / ratio
         upper = ratio * i%kT

         do
            call solve(i)

            if (i%Delta(0) .le. small) then
               upper = i%kT
            else
               lower = i%kT
            end if

            i%kT = (lower + upper) / 2

            if ((upper - lower) / 2 .le. error) exit
         end do

         write (*, "(ES21.14E3, ' (McMillan)')") i%Tc
         write (*, "(ES21.14E3, ' (Eliashberg)')") i%kT * qe / kB
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
