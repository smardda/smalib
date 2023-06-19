module geofil_h

  use const_kind_m
  use geobjlist_h
  use gfcontrol_h

! public types
!> geofil object
  type, public :: geofil_t
     type(geobjlist_t) :: geobjl !< geometrical objects
     type(gfnumerics_t)  :: n  !< numerical parameters
   !> integer scalar describing geometry such as triangles
   !! denoting e.g. array of bodies
     integer(ki4), dimension(:), allocatable :: scag !< .
     integer(ki4) :: nscag !< number of geometry scalars in array, e.g. number of body entries

  end type geofil_t

end module geofil_h
