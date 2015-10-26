module pcontrol_h

  use const_kind_m


! public types

  type, public :: pfiles_t
   !! file names
     character(len=80)  :: powcalout       !< output data
     character(len=80)  :: log            !< log file
     character(len=80)  :: geoq         !< geoq output format file name
     character(len=80)  :: vtkdata         !< vtk input data file
     character(len=80)  :: hdsdata         !< hds input data file
     character(len=80)  :: laudata         !< launch data
     character(len=80)  :: vtkres         !< vtk input res file
     character(len=80)  :: cartv   !< vtk plot file of power per big triangle in Cartesians
     character(len=80)  :: allcartv   !< vtk plot file of power for all triangles in Cartesians
     character(len=80)  :: wall   !< vtk plot file of power for shadowing wall triangles in Cartesians
     character(len=80)  :: ptzv    !< vtk plot file of all powers in \f$ (\psi,\theta,\zeta) \f$
     character(len=80)  :: flincart  !< plot file for field line following
     character(len=80)  :: flinends  !< gnuplot file for field line ends
     character(len=80)  :: flinptz  !< plot file for field line following in \f$ (\psi,\theta,\zeta) \f$
     character(len=80)  :: allptzv    !< vtk plot file of powers in \f$ (\psi,\theta,\zeta) \f$
     character(len=80)  :: geofld  !< vtk plot file of geometry and field
     character(len=80)  :: geofldq  !< vtk plot file of quantised geometry and field
     character(len=80)  :: gnu !< gnuplot file TBD
  end type pfiles_t

  type, public :: pplots_t
   !! vtk plot output selectors
     logical  :: cartv   !< power in Cartesians
     logical  :: allcartv   !< all power in Cartesians
     logical  :: wall   !< wall power in Cartesians
     logical  :: ptzv    !< vtk plot file of powers in \f$ (\psi,\theta,\zeta) \f$
     logical  :: flincart    !< plot file for field line following
     logical  :: flinends    !< gnuplot field line endpoints
     logical  :: flinptz    !< plot file for field line following in \f$ (\psi,\theta,\zeta) \f$
     logical  :: allptzv    !< vtk plot file of all powers in \f$ (\psi,\theta,\zeta) \f$
     logical  :: geofld  !< vtk plot file of geometry and field
     logical  :: geofldq  !< vtk plot file of quantised geometry and field
     logical  :: gnu !< gnuplot file TBD
  end type pplots_t

end module pcontrol_h
