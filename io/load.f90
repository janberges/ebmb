module io_load
   use arguments
   use global
   use integration
   implicit none

   private
   public :: load

contains

   subroutine load(i)
      type(universal), intent(out) :: i

      character(:), allocatable :: setting ! command-line argument
      character(:), allocatable :: lhs, rhs ! left- and right-hand side

      character(:), allocatable :: lambda ! string defining lambda
      character(:), allocatable :: muStar ! string defining muStar

      integer :: equals ! position of '='

      integer :: p ! band index
      integer :: n ! argument number

      character(:), allocatable :: DOSfile

      real(dp) :: elements ! number of elements in lambda and muStar

      elements = i%bands ** 2

      i%name = 'eb_local_untitled'

      do n = 1, command_argument_count()
         setting = argument(n)

         equals = index(setting, '=')

         lhs = setting(:equals - 1)
         rhs = setting(equals + 1:)

         select case (lhs)
            case ('name'); i%name = rhs

            case ('T'); read (rhs, *) i%T

            case ('error'); read (rhs, *) i%error
            case ('bound'); read (rhs, *) i%bound
            case ('small'); read (rhs, *) i%small

            case ('omegaE'); read (rhs, *) i%omegaE

            case ('lambda', 'lamda')
               lambda = rhs
               elements = values(rhs)

            case ('muStar', 'mu*')
               muStar = rhs
               elements = values(rhs)

            case ('DOSfile'); DOSfile = rhs

            case ('upper'); read (rhs, *) i%upper
            case ('lower'); read (rhs, *) i%lower

            case ('limit'); read (rhs, *) i%limit

            case ('measurable'); read (rhs, *) i%measurable
            case ('resolution'); read (rhs, *) i%resolution

            case ('form'); read (rhs, *) i%form
            case ('edit'); read (rhs, *) i%edit

            case ('standalone'); read (rhs, *) i%standalone
            case ('rescale');    read (rhs, *) i%rescale

            case ('epsilon'); read (rhs, *) negligible_difference
         end select
      end do

      if (i%T .lt. 0) then
         i%critical = .true.
         i%measurable = .false.
         i%resolution = 0
      end if

      i%T = kB * i%T

      i%error = kB * i%error
      i%bound = kB * i%bound

      i%bands = nint(sqrt(elements))

      allocate(i%lambda(i%bands, i%bands))
      allocate(i%muStar(i%bands, i%bands))

      if (allocated(lambda)) then
         read (lambda, *) i%lambda
      else
         i%lambda(:, :) = 0

         do p = 1, i%bands
            i%lambda(p, p) = 1
         end do
      end if

      if (allocated(muStar)) then
         read (muStar, *) i%muStar
      else
         i%muStar(:, :) = 0
      end if

      if (allocated(DOSfile)) then
         i%DOS = .true.
         call load_dos(DOSfile, i)
      end if

      if (i%lower .lt. 0) i%lower = i%upper
   end subroutine load

   integer function values(list)       ! number of ...
      character(*), intent(in) :: list ! comma-separated values

      integer :: c ! character position

      values = 1

      do c = 1, len(list)
         if (list(c:c) .eq. ',') values = values + 1
      end do
   end function values

   subroutine load_dos(DOSfile, i)
      character(*), intent(in) :: DOSfile
      type(universal), intent(inout) :: i

      integer :: n, m, p
      integer, parameter :: unit = 11

      open (unit, file=DOSfile, action='read', status='old')

      read (unit, *) n ! density-of-states resolution

      allocate(i%energy(n)) ! free-electron energy (eV)
      allocate(i%density(n, i%bands)) ! density of states (a.u.)
      allocate(i%weight(n, i%bands)) ! integration weight (eV)

      do m = 1, n
         read(unit, *) i%energy(m), i%density(m, :)
      end do

      close (unit)

      call differential(i%energy, i%weight(:, 1))

      do p = 2, i%bands
         i%weight(:, p) = i%weight(:, 1)
      end do

      n = minloc(abs(i%energy), 1) ! index of Fermi level

      do p = 1, i%bands
         i%weight(:, p) = i%weight(:, p) * i%density(:, p) / i%density(n, p)
      end do
   end subroutine load_dos
end module io_load
