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
  end type vtfmdata_t

  type, public :: vnumerics_t
   !! numerical run parameters
     integer(ki4) :: npans !< number of panels for which transform defined
     integer(ki4), dimension(:), allocatable :: pantfm !< number of transform to apply to body
     integer(ki4), dimension(:,:), allocatable :: panbod !< second index is no. of panel corresponding to body
     type(vtfmdata_t) :: vptfm !< array of position \f$ x \f$ to \f$ x \f$ scaling
     type(vtfmdata_t) :: vpantfm !< for each panel, array of position \f$ x \f$ to \f$ x \f$ scaling
  end type vnumerics_t


! public types
  type, public :: vfiles_t
   !! file names
     character(len=80)  :: vout       !< output data
     character(len=80)  :: log            !< log file
     character(len=80)  :: vtkdata         !< vtk input data file
     character(len=80)  :: vtkout         !< vtk output data file
  end type vfiles_t


end module vcontrol_h
