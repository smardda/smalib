!> @addtogroup groupname3
!> @{
module fmesh_h
!> @}
  use const_kind_m
  use position_h

! public types
  type, public :: fmesh_t

     character(len=80) :: formula !< mesh formula
     integer(ki4)  :: ndimf !< dimension of mesh array (1, 2 or 3)
     integer(ki4)  :: nstagf !< whether mesh staggered (1) or not (0)

     integer(ki4)  :: nxf !< first dimension of field array
     integer(ki4)  :: nyf !< second dimension of field array
     integer(ki4)  :: nzf !< third dimension of field array

     real(kr8), dimension (:), allocatable :: xf !< mesh-points for field in first direction
     real(kr8), dimension (:), allocatable :: yf !< mesh-points for field in second direction
     real(kr8), dimension (:), allocatable :: zf !< mesh-points for field in third direction

     real(kr8) :: dxf !< field mesh spacing size in first direction
     real(kr8) :: dyf !< field mesh spacing size in second direction
     real(kr8) :: dzf !< field mesh spacing size in third direction

     real(kr8) :: rdxf !< reciprocal of field mesh spacing size in first direction
     real(kr8) :: rdyf !< reciprocal of field mesh spacing size in second direction
     real(kr8) :: rdzf !< reciprocal of field mesh spacing size in third direction

     real(kr8), dimension (3) :: xa0f !< lower and upper co-ordinate bound of mesh as vectors
     real(kr8), dimension (3) :: xa1f !< lower and upper co-ordinate bound of mesh as vectors

     real(kr8) :: xlengthf !< physical extent of mesh  in first direction
     real(kr8) :: ylengthf !< physical extent of mesh  in second direction
     real(kr8) :: zlengthf !< physical extent of mesh  in third direction

     real(kr8) :: x0f !< origin of mesh in first direction
     real(kr8) :: y0f !< origin of mesh in second direction
     real(kr8) :: z0f !< origin of mesh in third direction

     real(kr8) :: x1f !< far point of mesh in first direction
     real(kr8) :: y1f !< far point of mesh in second direction
     real(kr8) :: z1f !< far point of mesh in third direction

     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:), allocatable :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable :: npar !< general integer parameters

     real(kr8), dimension (3) :: xL0 !< lower and upper co-ordinate bound of mesh as vectors
     real(kr8), dimension (3) :: xL1 !< lower and upper co-ordinate bound of mesh as vectors

     logical :: iltfm !< mesh is subject to transformation
     character*2 :: lunit !< units, either metres (me) or mm (mx,x/=e)
     character*2 :: punit !< plot units, either metres (me) or mm (mx,x/=e)
     type(tfmdata_t) :: tfmdata   !< describe mesh transform (apply to put in CAD frame)

   !     integer(ki4) :: nx_in !< input field mesh size
   !     integer(ki4) :: ny_in !< input field mesh size
   !     integer(ki4) :: nz_in !< input field mesh size
   !     real(kr8), dimension (:), allocatable:: x_in !< input field mesh
   !     real(kr8), dimension (:), allocatable::  y_in !< input field mesh
   !     real(kr8), dimension (:), allocatable::  z_in !< input field mesh
   !     real(kr8) :: dx_in !< input field mesh spacing size
   !     real(kr8) ::  dy_in !< input field mesh spacing size
   !     real(kr8) ::  dz_in !< input field mesh spacing size
   !
   !     real(kr8) :: length_x !< physical extent of mesh
   !     real(kr8) ::  length_y !< physical extent of mesh
   !     real(kr8) :: length_z !< physical extent of mesh
   !     real(kr8) :: x0 !< origin of mesh in first direction
   !     real(kr8) :: y0 !< origin of mesh in second direction
   !     real(kr8) :: z0 !< origin of mesh in third direction
   !
   !     integer(ki4)  :: fnx !< first dimension of field array (100)
   !     integer(ki4)  :: fny !< second dimension of field array (50)
   !     integer(ki4)  :: fnz !< third dimension of field array (200)
   !     real(kr8) :: fdx !< field mesh spacing sizes in each direction
   !     real(kr8) :: fdy !< field mesh spacing sizes in each direction
   !     real(kr8) :: fdz !< field mesh spacing sizes in each direction

  end type fmesh_t

end module fmesh_h
