module gfcontrol_h

  use const_kind_m

!! public types

!> numerical run parameters
  type, public :: gfnumerics_t
     logical :: calcangle !< calculate angle
     integer(ki4), dimension(0:10) :: objadd !< number of objects of given description to add to geobjlist
     character(len=80) :: namekey !< used to identify scalar variable, usually identifies body
  end type gfnumerics_t

!> file names
  type, public :: gffiles_t
     character(len=80)  :: out       !< output data
     character(len=80)  :: log           !< log file
     character(len=80)  :: vtk !< vtk input data file
     character(len=80)  :: vtkout         !< vtk output data file
  end type gffiles_t

!> plot output selectors
  type, public :: gfplots_t
     logical  :: vtk   !< INERT vtk plot selector
  end type gfplots_t

end module gfcontrol_h
