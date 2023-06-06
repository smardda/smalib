!> @addtogroup groupname3
!> @{
module pcle_h
!> @}
  use const_kind_m

! public types
  type, public :: posnode_t
     real(kr4), dimension(3) :: posvec !< pcle vector
     integer(ki4)  :: node !< corresponding node number
  end type posnode_t

end module pcle_h
