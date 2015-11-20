module powcal_h

  use const_kind_m
  use position_h
  use geobjlist_h
  use geobj_m
  use beq_h
  use odes_h
  use powres_h
  use edgprof_h
  use edgprof_m
  use spl2d_m
  use spl3d_m
  use termplane_h


! public types

!> numerical control parameters for all power calculation
  type, public :: pnumerics_t
     real(kr8) :: f !< power split ion to electron direction
     real(kr8) :: lmid !< power decay length at outer midplane (m)
     real(kr8) :: ploss !< power crossing the separatrix (W)
     real(kr8) :: sigma !< diffusion length (m)
     integer(ki4)  :: nlevel !< refinement level
     integer(ki4)  :: shadow !< shadowing calculation if positive
     real(kr8) :: qpara0 !< parallel reference flux (W/sqm)
   ! additions for powcal
     character(len=80) :: caltype !< type of calculation ('afws','local', 'msus','global', 'msum','middle')
     integer(ki4)  :: mtrack !< maximum number of tracks to follow (<= hardwired parameter)
     integer(ki4)  :: ntrack !< if non-zero, number of tracks to follow
     integer(ki4), dimension(:), allocatable  :: trackno !< array containing track element numbers
     integer(ki4)  :: nanalau !< type of analytic definition of launch
     integer(ki4)  :: nlaupow !< calculate power deposited on launch geometry
     integer(ki4)  :: nshapow !< calculate power deposited on shadowing geometry
     logical  :: ltermplane !< termination planes present
     type(termplane_t) :: termp !< termination planes
  end type pnumerics_t

!> type describing launch calculation data
  type, public :: powlau_t
     integer(ki4)  :: ndummy  !< type variable
  end type powlau_t

!> type describing all power calculation data
  type, public :: powcal_t
     type(powres_t) :: powres !< powres
     type(edgprof_t) :: edgprof !< edgprof
     type(powlau_t) :: powlau !< powlau
     type(odes_t) :: odes !< odes
     type(pnumerics_t) :: n !< control  parameters
  end type powcal_t

end module powcal_h
