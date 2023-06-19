module datline_h

  use const_kind_m

! public types

!> type storing one logical line of dat file input
  type, public :: datline_t
     integer(ki4) :: nlines=4 !< type variable
   !     character(len=80*nlines) :: line
     character(len=80*4) :: line !< type variable
     integer(ki4)  :: inline !< type variable
     character(len=5) :: type !< type variable
     integer(ki4)  :: n8 !< type variable
     character(len=8), dimension(:), allocatable :: flds8 !< type variable
     integer(ki4)  :: n16 !< type variable
     character(len=16), dimension(:), allocatable :: flds16 !< type variable
  end type datline_t

end module datline_h
