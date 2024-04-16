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

      integer :: equals ! position of '='

      integer :: i ! band index
      integer :: n ! argument number

      character(99) :: dos_file = 'none' ! file with density of states
      character(99) :: a2F_file = 'none' ! file with Eliashberg spectral function

      real(dp) :: elements ! number of elements in lambda and muStar

      real(dp) :: clip = 15.0_dp ! maximum real-axis frequency (omegaE)

      elements = x%bands ** 2
      x%bands = -1

      do n = 1, command_argument_count()
         setting = argument(n)

         equals = index(setting, '=')

         lhs = setting(:equals - 1)
         rhs = setting(equals + 1:)

         select case (lhs)
            case ('file'); read (rhs, '(A)') x%file
            case ('form'); read (rhs, '(A)') x%form

            case ('tell'); read (rhs, *) x%tell

            case ('T'); read (rhs, *) x%T

            case ('omegaE');  read (rhs, *) x%omegaE
            case ('cutoff');  read (rhs, *) x%cutoff
            case ('cutoffC'); read (rhs, *) x%cutoffC

            case ('lambda', 'lamda')
               lambda = rhs
               elements = matches(rhs, ',') + 1

            case ('muStar', 'mu*')
               muStar = rhs
               elements = matches(rhs, ',') + 1

            case ('muC')
               muC = rhs
               elements = matches(rhs, ',') + 1

            case ('bands'); read (rhs, *) x%bands

            case ('dos'); read (rhs, '(A)') dos_file
            case ('a2F'); read (rhs, '(A)') a2F_file

            case ('n');  read (rhs, *) x%n
            case ('mu'); read (rhs, *) x%mu

            case ('conserve'); read (rhs, *) x%conserve
            case ('chi');      read (rhs, *) x%chi

            case ('limit'); read (rhs, *) x%limit

            case ('epsilon'); read (rhs, *) epsilon
            case ('error');   read (rhs, *) x%error
            case ('zero');    read (rhs, *) x%zero
            case ('rate');    read (rhs, *) x%rate

            case ('lower'); read (rhs, *) x%lower
            case ('upper'); read (rhs, *) x%upper
            case ('clip');  read (rhs, *) clip

            case ('eta', '0+'); read (rhs, *) x%eta

            case ('resolution'); read (rhs, *) x%resolution
            case ('measurable'); read (rhs, *) x%measurable

            case ('unscale'); read (rhs, *) x%unscale
            case ('rescale'); read (rhs, *) x%rescale
            case ('imitate'); read (rhs, *) x%imitate

            case ('normal'); read (rhs, *) x%normal

            case ('power'); read (rhs, *) x%power

            case default
               print "('Ignored unknown parameter ''', A, '''')", lhs
         end select
      end do

      if (x%bands .eq. -1) x%bands = nint(sqrt(elements))

      allocate(x%lambda(x%bands, x%bands))
      allocate(x%muStar(x%bands, x%bands))

      if (allocated(lambda)) then
         read (lambda, *) x%lambda
      else
         x%lambda(:, :) = 0

         do i = 1, x%bands
            x%lambda(i, i) = 1
         end do
      end if

      if (allocated(muC)) then
         read (muC, *) x%muStar
         x%unscale = .false.
      else if (allocated(muStar)) then
         read (muStar, *) x%muStar
      else
         x%muStar(:, :) = 0
      end if

      if (dos_file .ne. 'none') then
         x%ldos = .true.
         call load_dos(dos_file, x)
      end if

      if (a2F_file .ne. 'none') then
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

      open (unit, file=file, action='read', status='old')

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
         read (unit, *) x%energy(m), x%dos(m, :)
      end do

      close (unit)
   end subroutine load_dos

   subroutine load_a2F(file, x)
      character(*), intent(in) :: file
      type(parameters), intent(inout) :: x

      integer :: n, m

      real(dp) :: test
      integer :: error

      open (unit, file=file, action='read', status='old')

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
            read (unit, *) x%omega(m), x%a2F(m, :, :)
         end do
      end do

      close (unit)
   end subroutine load_a2F
end module io_load
