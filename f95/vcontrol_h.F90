module vcontrol_h

  use const_kind_m
  use position_h

! public types

  type, public :: vtfmdata_t
     integer(ki4) :: ntfmtyp !< number of types of position vector transform
     real(kr4), dimension(:,:), allocatable :: scale !< scaling vector of transform
     real(kr4), dimension(:,:), allocatable :: offset !< offset vector
     real(kr4), dimension(:,:,:), allocatable :: matrix  !< maxtrix of transform
     integer(ki4), dimension(:), allocatable :: ntfm !< type of position vector transform
     character(len=20), dimension(:), allocatable :: id  !< identifier for transform
  end type vtfmdata_t

  type, public :: vnumerics_t
   !! numerical run parameters
     logical :: split !<  split file by attribute if set
     logical :: same !< make attribute take same value
     integer(ki4) :: nvalue !<  value for attribute to take
     character(len=80) :: name !<  name of attribute to be split / homogenised
     integer(ki4) :: npans !< number of panels for which transform defined
     integer(ki4), dimension(:), allocatable :: pantfm !< number of transform to apply to body
     integer(ki4), dimension(:,:), allocatable :: panbod !< second index is no. of panel corresponding to body
     character(len=6) :: angles  !< units, either radian or degree
     type(vtfmdata_t) :: vptfm !< array of position \f$ x \f$ to \f$ x \f$ scaling
     type(vtfmdata_t) :: vpantfm !< for each panel, array of position \f$ x \f$ to \f$ x \f$ scaling
     logical :: paneltfm !< apply transform if .TRUE.
     integer(ki4) :: maxindx !< dimension of bods index array
     integer(ki4) :: maxbodsf !< used to generate unique bods numbers over many files
     logical :: preserve !< bods remain distinct
     logical :: extract !< extract objects according to criterion
     character(len=80) :: key !< key for extraction
     real(kr8), dimension(2) :: centre !< centre of discharge in \f$ (R,Z) \f$
     real(kr8) :: angmin !< minimum angle
     real(kr8) :: angmax !< maximum angle
  end type vnumerics_t


! public types
  type, public :: vfiles_t
   !! file names
     character(len=80)  :: vout       !< output data
     character(len=80)  :: log            !< log file
     character(len=80), dimension(:), allocatable :: vtkdata !< vtk input data files
     integer(ki4), dimension(:), allocatable :: vtkcopies !< number of copies of each vtkdata entry
     character(len=80), dimension(:), allocatable :: vtklabel !< describes each vtkdata entry
     integer(ki4) :: nvtkdata !< actual size of vtkdata
     character(len=80)  :: vtkout         !< vtk output data file
  end type vfiles_t


end module vcontrol_h
