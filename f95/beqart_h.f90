!> @addtogroup groupname3
!> @{
module beqart_h
!> @}
  use const_kind_m
  use fmesh_h

  implicit none

!> data structure for 3D applied magnetic field
  type, public :: beqart_t

     real(kr8), dimension(:,:,:), allocatable:: bx !< field array
     real(kr8), dimension(:,:,:), allocatable:: by !< field array
     real(kr8), dimension(:,:,:), allocatable:: bz  !< field array

     type (fmesh_t) :: fmesh  !< fmesh
  end type beqart_t

end module beqart_h
