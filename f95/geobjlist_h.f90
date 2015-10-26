module geobjlist_h

  use const_kind_m
  use position_h
  use geobj_m


! public types

!> type storing geobj coord data and useful bits
  type, public :: geobjlist_t
     integer(ki4)  :: ngtype !< type of geobj (if all the same)
     type(posveclis_t) :: posl !< list of positions
     integer(ki4), dimension(:), allocatable :: nodl !< list of nodes
     integer(ki4)  :: np !< number of position vectors
     type(geobj1_t), dimension(:), allocatable :: obj !< type variable
     type(geobj2_t), dimension(:), allocatable :: obj2 !< type variable
     integer(ki4)  :: nnod !< no of nodes
     integer(ki4)  :: ng !< no of geobjs
     integer(ki4)  :: ngunassigned !< no of geobjs not assigned
     integer(ki2)  :: nquant !< estimated then actual quantising number (log)
     integer(ki4)  :: minobj !< minimum number of geobj in a bin
     integer(ki4)  :: nparam=0 !< integer parameter from file header
     real(kr4)  :: rparam=0._kr4 !< real parameter from file header
     real(kr4), dimension(3,2) :: coordbb !< bb of geobj coord
     real(kr4), dimension(3,2) :: binbb !< bb for geobj binning
     real(kr4) :: tolerance !< max distance from face
     real(kr4) :: minmaxtolerance !< min max distance from face
     integer(ki4) :: nbdcub=0 !< number of extra bounding cubes
     real(kr4) :: dilen   !< inner cube separation from geometry
     real(kr4) :: dolen   !< outer cube separation from inner
     integer(ki4)  :: nwset !< nonzero if obj%weight set
     type(quantfm_t)    :: quantfm !< geobj \f$ x \f$ to mesh units  scaling
     type(tfmdata_t)    :: tfmdata !< position \f$ x \f$ to \f$ x \f$ scaling
  end type geobjlist_t

end module geobjlist_h
