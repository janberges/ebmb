! Copyright (C) 2016-2025 Jan Berges
! This program is free software under the terms of the GNU GPLv3 or later.

module formatting
   use globals
   implicit none

   private
   public :: measure, edit, rule

   integer :: width
   character(:), allocatable :: w, x

contains

   subroutine measure(flomat)
      character(*), intent(in) :: flomat

      character(100) :: test

      x = trim(flomat)

      write (test, "(" // x // ", '|')") pi
      width = index(test, '|') - 1

      write (test, '(I0)') width
      w = trim(test)
   end subroutine measure

   function edit(descriptor)
      character(:), allocatable :: edit
      character(*), intent(in) :: descriptor

      integer :: n

      edit = descriptor

      do
         n = scan(edit, 'wx')
         if (n .eq. 0) return

         select case (edit(n:n))
            case ('w'); edit = edit(:n - 1) // w // edit(n + 1:)
            case ('x'); edit = edit(:n - 1) // x // edit(n + 1:)
         end select
      end do
   end function edit

   function rule(n)
      character(:), allocatable :: rule
      integer, intent(in) :: n

      rule = "('" // repeat('_', n * width) // "')"
   end function rule
end module formatting
