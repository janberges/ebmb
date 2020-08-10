module io_load
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

      real(dp) :: elements ! number of elements in lambda and muStar

      real(dp) :: clip = 15.0_dp ! maximum real-axis frequency (omegaE)

      elements = x%bands ** 2

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

            case ('dos'); read (rhs, *) dos_file

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

      x%bands = nint(sqrt(elements))

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

      n = 0 ! density-of-states resolution

      do
         read (unit, *, iostat=error) test
         if (error .ne. 0) exit
         n = n + 1
      end do

      rewind unit

      allocate(x%energy(n)) ! free-electron energy (eV)
      allocate(x%dos(n, x%bands)) ! density of states (a.u.)

      do m = 1, n
         read(unit, *) x%energy(m), x%dos(m, :)
      end do

      close (unit)
   end subroutine load_dos
end module io_load
