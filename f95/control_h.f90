!> @addtogroup groupname3
!> @{
module control_h
!> @}
  use const_kind_m
  use position_h

! public types

  type, public :: numerics_t
   !! numerical run parameters
     integer(ki4) :: nsize !< size of binary tree
     integer(ki4) :: nsizep !< size of binary tree pter
     integer(ki4) :: nsizee !< size of exten array
     integer(ki4) :: nsizeh !< size of linked list hoc
     integer(ki4) :: nsizel !< size of linked list
     integer(ki2) :: ndepth !< max depth of binary tree
     integer(ki4) :: nttype !< type of binary tree
     integer(ki4) :: nttalg !< type of binary tree algorithm for top node
     integer(ki2), dimension(3) :: nxyz !< dimensions of top node (nttype=2)
     real(kr4), dimension(3) :: hxyz !< specified mesh sizes
     integer(ki4) :: ngeobj !< number of geobj records to be read
     integer(ki4) :: mtype !< type of margin for quantisation
     integer(ki4) :: geomtype !< geobj type
     integer(ki4) :: mingeobjinbin !< minimum number of geobj records in a bin
     integer(ki2) :: nquante !< estimated log number of mesh cells
     real(kr4) :: maxtolerance !< max max distance of geobj from its side
     real(kr4) :: mintolerance !< min max distance of geobj from its side
     integer(ki4) :: nbdcub !< number of extra bounding cubes
     integer(ki4) :: cornflag !< flags (=1) if corner vectors to be used
     real(kr4), dimension(3) :: lowcorner !< vector defining corner with lesser component values
     real(kr4), dimension(3) :: upcorner !< vector defining corner with greater component values
     real(kr4) :: dilen   !< inner cube separation from geometry
     real(kr4) :: dolen   !< outer cube separation from inner
     type(quantfm_t) :: geobj_coord_tfm !< geobj \f$ x \f$ to mesh units scaling
     type(tfmdata_t) :: position_coord_tfm !< position \f$ x \f$ to \f$ x \f$ scaling
     real(kr4), dimension(3,2) :: coordbb !< bounding box of geobj coord
     real(kr4), dimension(3,2) :: binbb !< bounding box for geobj binning
  end type numerics_t

  type, public :: mtnumerics_t
     type(numerics_t) :: n !< data for hdsgen
     integer(ki4) :: vcontent !< interpret vtk test file data
     real(kr8) :: mdt !< nominal timestep if velocity in vtk test file
     character(len=80) :: fldnam       !< fieldname
     character(3) :: qfilesuf !< suffix of mtest input file
  end type mtnumerics_t

! public types
  type, public :: files_t
   !! file names
     character(len=80) :: hdsgenout       !< output data
     character(len=80) :: dengenout       !< output data
     character(len=80) :: legtfmout       !< output data
     character(len=80) :: moveout       !< output data
     character(len=80) :: log            !< log file
     character(len=80) :: vtkdata         !< Objectg vtk format data file
     character(len=80) :: hdsdata         !< Objectg hds format data file
     character(len=80) :: qrydata         !< Query data file
     character(len=80) :: hdslist            !< HDS output
     character(len=80) :: hdsv           !< DUPLICATE vtk plot of HDS
     character(len=80) :: hdsm           !< vtk plot of HDS in mapped coords
     character(len=80) :: hdsq          !< vtk plot of HDS in mapped quantised coords
     character(len=80) :: geobjq     !< vtk plot of assigned geobj
     character(len=80) :: lostgeobj  !< vtk plot file of unassigned geobj
     character(len=80) :: allgeobj  !< vtk plot file of all geobj
     character(len=80) :: allgeobjq  !< DUPLICATE vtk plot file of all geobj quantised
     character(len=80) :: geoptq  !< vtk plot file of all geobj in mapped quantised coords
     character(len=80) :: densitygeobj  !< vtk plot file of density geobj
     character(len=80) :: den!< vtk plot file of density field
     character(len=80) :: denq!< vtk plot file of quantised density field
     character(len=80) :: mov!< vtk plot file of particle moves
     character(len=80) :: movq!< vtk plot file of quantised particle moves
     character(len=80) :: pts!< vtk plot file of particle points
     character(len=80) :: ptsq!< vtk plot file of quantised particle points
     character(len=80) :: tfm !< vtk plot file of transformed density field
     character(len=80) :: rtt !< rtt file of transformed density field
     character(len=80) :: gnu !< gnuplot file of transformed density field
     character(len=80) :: rest !< rtt file of restored density field
  end type files_t

  type, public :: numbers_t
   !! run parameters for transform control
     integer(ki4) :: nestart !< energy level control
     integer(ki4) :: nestop !< energy level control
     integer(ki4) :: degree !< expansion degree in angles
     real(kr4) :: minenergy !< energy level control
     real(kr4) :: maxenergy !< energy level control
     real(kr4), dimension(3) :: coord !< source coord
     real(kr4), dimension(3) :: uvec !< u vector
     real(kr4), dimension(3) :: vvec !< v vector
     real(kr4) :: strength !< source strength
     character(len=80) :: name  !< source name
  end type numbers_t

  type, public :: plots_t
   !! vtk plot output selectors
     logical :: hds   !< DUPLICATE vtk plot of hds
     logical :: hdsm   !< vtk plot of hds
     logical :: hdsbin     !< DUPLICATE vtk plot of HDS bbs
     logical :: hdsq   !< vtk plot of hds in mapped quantised coords
     logical :: geobj   !< vtk plot of geobj
     logical :: geobjq   !< vtk plot of geobjq quantised
     logical :: lostgeobj  !< vtk plot file of unassigned geobj
     logical :: allgeobj  !< DUPLICATE vtk plot file of all geobj
     logical :: allgeobjq  !< DUPLICATE vtk plot file of all geobj quantised
     logical :: geoptq  !< vtk plot file of all geobj in mapped quantised coords
     logical :: densitygeobj  !< vtk plot file of density geobj
     logical :: den   !< vtk plot of den
     logical :: denq   !< vtk plot of quantised den
     logical :: tfm   !< vtk plot of transformed den
     logical :: rest   !< vtk plot of restored den
     logical :: mov   !< vtk plot of mov
     logical :: movq   !< vtk plot of quantised mov
     logical :: pts   !< vtk plot of pts
     logical :: ptsq   !< vtk plot of quantised pts
  end type plots_t


end module control_h
