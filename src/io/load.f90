! Copyright (C) 2016-2024 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module io_load
   use eliashberg_spectral_function
   use global
   use tools, only: argument, matches
   implicit none

   private
   public :: load

contains

   subroutine load(x)
      type(parameters), intent(out) :: x

      character(:), allocatable :: setting ! command-line argument
      character(:), allocatable :: lhs, rhs ! left- and right-hand side

      character(:), allocatable :: lambda ! string defining lambda
      character(:), allocatable :: muStar ! string defining muStar
      character(:), allocatable :: muC    ! string defining muC

      character(:), allocatable :: dos_file ! file with density of states
      character(:), allocatable :: a2F_file ! file with Eliashberg function

      integer :: equals ! position of '='

      integer :: i ! band index
      integer :: n ! argument number

      integer :: error ! I/O status

      real(dp) :: elements ! number of elements in lambda and muStar

      real(dp) :: clip = 15.0_dp ! maximum real-axis frequency (omegaE)

      elements = x%bands ** 2
      x%bands = -1

      do n = 1, command_argument_count()
         setting = argument(n)

         equals = index(setting, '=')

         if (equals .lt. 2 .or. equals .eq. len(setting)) then
            print "('Error: Invalid argument ""', A, '""')", setting
            stop 1
         end if

         lhs = setting(:equals - 1)
         rhs = setting(equals + 1:)

         error = 0

         select case (lhs)
            case ('file'); x%file = rhs
            case ('form'); x%form = rhs

            case ('tell'); read (rhs, *, iostat=error) x%tell

            case ('T'); read (rhs, *, iostat=error) x%T

            case ('omegaE');  read (rhs, *, iostat=error) x%omegaE
            case ('cutoff');  read (rhs, *, iostat=error) x%cutoff
            case ('cutoffC'); read (rhs, *, iostat=error) x%cutoffC

            case ('lambda', 'lamda')
               lambda = rhs
               elements = matches(rhs, ',') + 1

            case ('muStar', 'mu*')
               muStar = rhs
               elements = matches(rhs, ',') + 1

            case ('muC')
               muC = rhs
               elements = matches(rhs, ',') + 1

            case ('bands'); read (rhs, *, iostat=error) x%bands

            case ('dos', 'DOS'); dos_file = rhs
            case ('a2f', 'a2F'); a2F_file = rhs

            case ('n');  read (rhs, *, iostat=error) x%n
            case ('mu'); read (rhs, *, iostat=error) x%mu

            case ('conserve'); read (rhs, *, iostat=error) x%conserve
            case ('chi');      read (rhs, *, iostat=error) x%chi

            case ('limit'); read (rhs, *, iostat=error) x%limit

            case ('epsilon'); read (rhs, *, iostat=error) eps
            case ('error');   read (rhs, *, iostat=error) x%error
            case ('zero');    read (rhs, *, iostat=error) x%zero
            case ('rate');    read (rhs, *, iostat=error) x%rate

            case ('lower'); read (rhs, *, iostat=error) x%lower
            case ('upper'); read (rhs, *, iostat=error) x%upper
            case ('clip');  read (rhs, *, iostat=error) clip

            case ('eta', '0+'); read (rhs, *, iostat=error) x%eta

            case ('resolution'); read (rhs, *, iostat=error) x%resolution
            case ('measurable'); read (rhs, *, iostat=error) x%measurable

            case ('unscale'); read (rhs, *, iostat=error) x%unscale
            case ('rescale'); read (rhs, *, iostat=error) x%rescale
            case ('imitate'); read (rhs, *, iostat=error) x%imitate

            case ('normal'); read (rhs, *, iostat=error) x%normal

            case ('power'); read (rhs, *, iostat=error) x%power

            case default
               print "('Error: Unknown parameter ""', A, '""')", lhs
               stop 1
         end select

         if (error .ne. 0) then
            print "('Error: Invalid value for parameter ""', A, '""')", lhs
            stop 1
         end if
      end do

      if (x%bands .eq. -1) x%bands = int(sqrt(elements))

      allocate(x%lambda(x%bands, x%bands))
      allocate(x%muStar(x%bands, x%bands))

      if (allocated(lambda)) then
         read (lambda, *, iostat=error) x%lambda

         if (error .ne. 0) then
            print "('Error: ""lambda"" needs ', I0, ' numbers')", size(x%lambda)
            stop 1
         end if
      else
         x%lambda(:, :) = 0

         do i = 1, x%bands
            x%lambda(i, i) = 1
         end do
      end if

      if (allocated(muC)) then
         read (muC, *, iostat=error) x%muStar

         if (error .ne. 0) then
            print "('""muC"" needs ', I0, ' numbers')", size(x%muStar)
            stop 1
         end if

         x%unscale = .false.
      else if (allocated(muStar)) then
         read (muStar, *, iostat=error) x%muStar

         if (error .ne. 0) then
            print "('""muStar"" needs ', I0, ' numbers')", size(x%muStar)
            stop 1
         end if
      else
         x%muStar(:, :) = 0
      end if

      if (allocated(dos_file)) then
         x%ldos = .true.
         call load_dos(dos_file, x)
      end if

      if (allocated(a2F_file)) then
         x%la2F = .true.
         call load_a2F(a2F_file, x)
         call integrate_a2F(x)
      end if

      if (x%cutoffC .lt. 0) x%cutoffC = x%cutoff

      if (x%upper .lt. x%lower) x%upper = clip * x%omegaE
   end subroutine load

   subroutine load_dos(file, x)
      character(*), intent(in) :: file
      type(parameters), intent(inout) :: x

      integer :: n, m

      real(dp) :: test
      integer :: error

      open (unit, file=file, action='read', status='old', iostat=error)

      if (error .ne. 0) then
         print "('Error: Cannot open DOS file ""', A, '""')", file
         stop 1
      end if

      n = 0 ! number of sample points

      do
         read (unit, *, iostat=error) test
         if (error .ne. 0) exit
         n = n + 1
      end do

      rewind unit

      allocate(x%energy(n)) ! free-electron energy (eV)
      allocate(x%dos(n, x%bands)) ! density of states (a.u.)

      do m = 1, n
         read (unit, *, iostat=error) x%energy(m), x%dos(m, :)

         if (error .ne. 0) then
            print "('Error: DOS file needs ', I0, ' numbers per line')", &
               x%bands + 1
            stop 1
         end if
      end do

      close (unit)
   end subroutine load_dos

   subroutine load_a2F(file, x)
      character(*), intent(in) :: file
      type(parameters), intent(inout) :: x

      integer :: n, m

      real(dp) :: test
      integer :: error

      open (unit, file=file, action='read', status='old', iostat=error)

      if (error .ne. 0) then
         print "('Error: Cannot open a2F file ""', A, '""')", file
         stop 1
      end if

      n = 0 ! number of sample points

      do
         read (unit, *, iostat=error) test
         if (error .ne. 0) exit
         if (test .na. 0.0_dp) n = n + 1
      end do

      rewind unit

      allocate(x%omega(n)) ! phonon energy (frequency argument)
      allocate(x%a2F(n, x%bands, x%bands)) ! Eliashberg spectral function

      x%omega = 0.0_dp

      do m = 1, n
         do while (x%omega(m) .ap. 0.0_dp)
            read (unit, *, iostat=error) x%omega(m), x%a2F(m, :, :)

            if (error .ne. 0) then
               print "('Error: a2F file needs ', I0, ' numbers per line')", &
                  x%bands ** 2 + 1
               stop 1
            end if
         end do
      end do

      close (unit)
   end subroutine load_a2F
end module io_load
