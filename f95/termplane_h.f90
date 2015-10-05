module termplane_h

  use const_kind_m

! public types
!! coordinate directions \f$ (R,Z,\xi) \f$ correspond to resp. \f$ 1,2,3 \f$
  type, public :: termplane_t
     real(kr4), dimension(:), allocatable :: termplane !< termplane vector
     real(kr4), dimension(:,:), allocatable :: termstore !< termplane storage of vector
     integer(ki4), dimension(:,:), allocatable :: termplanedir !< termplane direction \f$ \pm 1-3 \f$
     integer(ki4)  :: ntermplane !< number of active termination planes
  end type termplane_t

end module termplane_h
