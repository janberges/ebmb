module io_load
   use global
   implicit none

   private
   public :: load

contains

   function argument(n)
      character(:), allocatable :: argument
      integer, intent(in) :: n

      integer :: size

      call get_command_argument(n, length=size)

      allocate(character(size) :: argument)

      call get_command_argument(n, value=argument)
   end function argument

   subroutine load(x)
      type(parameters), intent(out) :: x

      character(:), allocatable :: setting ! command-line argument
      character(:), allocatable :: lhs, rhs ! left- and right-hand side

      character(:), allocatable :: lambda ! string defining lambda
      character(:), allocatable :: muStar ! string defining muStar

      integer :: equals ! position of '='

      integer :: i ! band index
      integer :: n ! argument number

      character(:), allocatable :: dos_file

      real(dp) :: elements ! number of elements in lambda and muStar

      elements = x%bands ** 2

      do n = 1, command_argument_count()
         setting = argument(n)

         equals = index(setting, '=')

         lhs = setting(:equals - 1)
         rhs = setting(equals + 1:)

         select case (lhs)
            case ('file'); read (rhs, *) x%file
            case ('form'); read (rhs, *) x%form
            case ('tell'); read (rhs, *) x%tell

            case ('T'); read (rhs, *) x%T

            case ('omegaE');  read (rhs, *) x%omegaE
            case ('cutoff');  read (rhs, *) x%cutoff
            case ('cutoffC'); read (rhs, *) x%cutoffC

            case ('lambda', 'lamda')
               lambda = rhs
               elements = values(rhs)

            case ('muStar', 'mu*')
               muStar = rhs
               elements = values(rhs)

            case ('dos'); dos_file = rhs

            case ('limit'); read (rhs, *) x%limit

            case ('epsilon'); read (rhs, *) epsilon
            case ('error');   read (rhs, *) x%error
            case ('zero');    read (rhs, *) x%zero
            case ('rate');    read (rhs, *) x%rate
            case ('shift');   read (rhs, *) x%shift

            case ('resolution'); read (rhs, *) x%resolution
            case ('measurable'); read (rhs, *) x%measurable

            case ('rescale'); read (rhs, *) x%rescale
            case ('imitate'); read (rhs, *) x%imitate

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

      if (allocated(muStar)) then
         read (muStar, *) x%muStar
      else
         x%muStar(:, :) = 0
      end if

      if (allocated(dos_file)) then
         x%chi = .true.
         call load_dos(dos_file, x)
      end if

      if (x%cutoffC .lt. 0) x%cutoffC = x%cutoff
   end subroutine load

   integer function values(list)       ! number of ...
      character(*), intent(in) :: list ! comma-separated values

      integer :: c ! character position

      values = 1

      do c = 1, len(list)
         if (list(c:c) .eq. ',') values = values + 1
      end do
   end function values

   subroutine load_dos(file, x)
      character(*), intent(in) :: file
      type(parameters), intent(inout) :: x

      integer :: n, m
      integer, parameter :: unit = 11

      open (unit, file=file, action='read', status='old')

      read (unit, *) n ! density-of-states resolution

      allocate(x%energy(n)) ! free-electron energy (eV)
      allocate(x%dos(n, x%bands)) ! density of states (a.u.)

      do m = 1, n
         read(unit, *) x%energy(m), x%dos(m, :)
      end do

      close (unit)
   end subroutine load_dos
end module io_load
