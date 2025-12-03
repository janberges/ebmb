! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module io_load
   use eliashberg_self_energy, only: initialize_dos => initialize
   use eliashberg_spectral_function, only: initialize_a2F => initialize, &
      integrate_a2F
   use globals
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

      integer :: i, j ! band indices
      integer :: n ! argument number

      integer :: error ! I/O status

      integer :: elements ! number of elements in lambda and muStar

      if (x%diag) then
         elements = x%bands
      else
         elements = x%bands ** 2
      end if

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
            case ('file'); x%output = rhs
            case ('form'); x%flomat = rhs

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
            case ('diag');  read (rhs, *, iostat=error) x%diag

            case ('dos', 'DOS'); dos_file = rhs
            case ('a2f', 'a2F'); a2F_file = rhs

            case ('n');  read (rhs, *, iostat=error) x%n
            case ('mu'); read (rhs, *, iostat=error) x%mu

            case ('conserve'); read (rhs, *, iostat=error) x%conserve
            case ('readjust'); read (rhs, *, iostat=error) x%readjust
            case ('chi');      read (rhs, *, iostat=error) x%chi
            case ('chiC');     read (rhs, *, iostat=error) x%chiC
            case ('Sigma');    read (rhs, *, iostat=error) x%Sigma

            case ('steps', 'limit'); read (rhs, *, iostat=error) x%steps

            case ('epsilon'); read (rhs, *, iostat=error) eps
            case ('toln');    read (rhs, *, iostat=error) x%toln
            case ('error');   read (rhs, *, iostat=error) x%error
            case ('zero');    read (rhs, *, iostat=error) x%zero
            case ('rate');    read (rhs, *, iostat=error) x%rate

            case ('lower'); read (rhs, *, iostat=error) x%lower
            case ('upper'); read (rhs, *, iostat=error) x%upper

            case ('points', 'resolution'); read (rhs, *, iostat=error) x%points

            case ('logscale'); read (rhs, *, iostat=error) x%logscale

            case ('eta', '0+'); read (rhs, *, iostat=error) x%eta

            case ('measurable'); read (rhs, *, iostat=error) x%measurable

            case ('unscale'); read (rhs, *, iostat=error) x%unscale
            case ('rescale'); read (rhs, *, iostat=error) x%rescale
            case ('imitate'); read (rhs, *, iostat=error) x%imitate

            case ('divdos'); read (rhs, *, iostat=error) x%divdos
            case ('stable'); read (rhs, *, iostat=error) x%stable
            case ('normal'); read (rhs, *, iostat=error) x%normal
            case ('realgw'); read (rhs, *, iostat=error) x%realgw
            case ('eta0Im'); read (rhs, *, iostat=error) x%eta0Im

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

      if (x%bands .eq. -1) then
         if (x%diag) then
            x%bands = elements
         else
            x%bands = int(sqrt(real(elements, dp)))
         end if
      end if

      allocate(x%lambda(x%bands, x%bands))
      allocate(x%muStar(x%bands, x%bands))

      x%lambda(:, :) = 0.0_dp

      do i = 1, x%bands
         x%lambda(i, i) = 1.0_dp
      end do

      if (allocated(lambda)) then
         if (x%diag) then
            read (lambda, *, iostat=error) (x%lambda(i, i), i = 1, x%bands)
         else
            read (lambda, *, iostat=error) x%lambda
         end if

         if (error .ne. 0) then
            print "('Error: ""lambda"" needs ', I0, ' numbers')", size(x%lambda)
            stop 1
         end if
      end if

      x%muStar(:, :) = 0.0_dp

      if (allocated(muC)) then
         if (x%diag) then
            read (muC, *, iostat=error) (x%muStar(i, i), i = 1, x%bands)
         else
            read (muC, *, iostat=error) x%muStar
         end if

         if (error .ne. 0) then
            print "('""muC"" needs ', I0, ' numbers')", size(x%muStar)
            stop 1
         end if

         x%unscale = .false.
      else if (allocated(muStar)) then
         if (x%diag) then
            read (muStar, *, iostat=error) (x%muStar(i, i), i = 1, x%bands)
         else
            read (muStar, *, iostat=error) x%muStar
         end if

         if (error .ne. 0) then
            print "('""muStar"" needs ', I0, ' numbers')", size(x%muStar)
            stop 1
         end if
      end if

      if (allocated(dos_file)) then
         if (dos_file .ne. 'none') then
            x%ldos = .true.
            call load_dos(dos_file, x)
            call initialize_dos(x)
         end if
      end if

      if (allocated(a2F_file)) then
         if (a2F_file .ne. 'none') then
            x%la2F = .true.
            call load_a2F(a2F_file, x)
            call initialize_a2F(x)
            call integrate_a2F(x)
         end if
      end if

      if (.not. x%diag) then
         x%diag = .true.

         outer: do i = 1, x%bands
            do j = 1, x%bands
               if (i .eq. j) cycle

               if ((x%lambda(j, i) .na. 0.0_dp) .or. &
                   (x%muStar(j, i) .na. 0.0_dp)) then
                  x%diag = .false.
                  exit outer
               end if
            end do
         end do outer
      end if

      if (x%cutoffC .lt. 0.0_dp) x%cutoffC = x%cutoff

      if (x%upper .lt. x%lower) x%upper = x%cutoff * x%omegaE
   end subroutine load

   subroutine load_dos(filename, x)
      character(*), intent(in) :: filename
      type(parameters), intent(inout) :: x

      integer :: n, m

      real(dp) :: test
      integer :: error

      open (fun, file=filename, action='read', status='old', iostat=error)

      if (error .ne. 0) then
         print "('Error: Cannot open DOS file ""', A, '""')", filename
         stop 1
      end if

      n = 0 ! number of sample points

      do
         read (fun, *, iostat=error) test
         if (error .ne. 0) exit
         n = n + 1
      end do

      rewind fun

      allocate(x%energy(n)) ! free-electron energy (eV)
      allocate(x%dos(n, x%bands)) ! density of states (a.u.)

      do m = 1, n
         read (fun, *, iostat=error) x%energy(m), x%dos(m, :)

         if (error .ne. 0) then
            print "('Error: DOS file needs ', I0, ' numbers per line')", &
               x%bands + 1
            stop 1
         end if
      end do

      close (fun)
   end subroutine load_dos

   subroutine load_a2F(filename, x)
      character(*), intent(in) :: filename
      type(parameters), intent(inout) :: x

      integer :: i, n, m

      real(dp) :: test
      integer :: error

      open (fun, file=filename, action='read', status='old', iostat=error)

      if (error .ne. 0) then
         print "('Error: Cannot open a2F file ""', A, '""')", filename
         stop 1
      end if

      n = 0 ! number of sample points

      do
         read (fun, *, iostat=error) test
         if (error .ne. 0) exit
         if (test .na. 0.0_dp) n = n + 1
      end do

      rewind fun

      allocate(x%omega(n)) ! phonon energy (frequency argument)
      allocate(x%a2F(n, x%bands, x%bands)) ! Eliashberg spectral function

      x%omega = 0.0_dp

      do m = 1, n
         do while (x%omega(m) .ap. 0.0_dp)
            if (x%diag) then
               x%a2F(m, :, :) = 0.0_dp
               read (fun, *, iostat=error) x%omega(m), (x%a2F(m, i, i), &
                  i = 1, x%bands)
            else
               read (fun, *, iostat=error) x%omega(m), x%a2F(m, :, :)
            end if

            if (error .ne. 0) then
               print "('Error: a2F file needs ', I0, ' numbers per line')", &
                  x%bands ** 2 + 1
               stop 1
            end if
         end do
      end do

      close (fun)
   end subroutine load_a2F
end module io_load
