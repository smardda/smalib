module odes_h

  use const_kind_m
  use position_h


! public types

!> numerical control parameters
  type, public :: onumerics_t
     real(kr8) :: epsr !< relative error
     real(kr8) :: epsa !< absolute error
     real(kr8) :: dt0 !< initial timestep
     integer(ki4)  :: stepmax !< max number of steps
     real(kr8) :: tmax  !< maximum allowed value of dependent variable
     integer(ki4) :: npdim !< dimension of system
     integer(ki4)  :: nstartcon !< control start of ODE trajectory
     integer(ki4)  :: ntermcon !< control termination of ODE trajectory
     real(kr8), dimension(3)  :: termcon !< parameters to control termination of ODE trajectory
  end type onumerics_t

!> type storing ode solution data
  type, public :: odes_t
     real(kr8):: epsmach !< machine  \f$ \epsilon \f$
     real(kr8):: epsrmin !< min permitted value of epsr
     real(kr8):: epsbig !< no. significantly bigger than epsmach
     integer(ki4) :: nfacup !< magic factor for timestep control
     integer(ki4) :: nfacdn !< magic factor for timestep control
     real(kr8):: ffac !< magic factor for timestep control
     integer(ki4) :: npdim !< dimension of system
     integer(ki4) :: nsord !< max scheme order
     real(kr8) :: facu !< nfacup
     real(kr8) :: facut !< (ffac/facu)**nsord
     real(kr8) :: facd !< 1/nfacdn
     real(kr8) :: facdt !< (ffac/facd)**nsord
     real(kr8) :: rsord !< 1/nsord
     real(kr8) :: rsordb !< 1/(6*nsord)
     real(kr8) :: reps !< reciprocal of error epsr
     real(kr8) :: epsar !< type variable
     real(kr8) :: dtmin !< minimum acceptable timestep
     real(kr8) :: t  !< dependent variable, aka time \f$ t \f$ (not used in rjfunct)
     real(kr8) :: dt  !< timestep
     real(kr8) :: dto  !< previous timestep
     type(posveclis_t) :: vecp !< list of position vectors
     integer(ki4)  :: sstatus !< solver status
     real(kr8) :: rlength !< trajectory length, running total (use vecp-3 instead)
     real(kr8) :: dlength !< latest increment to trajectory length (use vecp-3 instead)
     real(kr8), dimension(:), allocatable :: length !< trajectory length (use vecp-3 instead)
     integer(ki4)  :: ndt !< number of entries in track (successful timesteps)
     real(kr8) :: glt !< generic trajectory limit, typically \f$ 2\pi \f$
     real(kr8), dimension(:), allocatable :: g3do !< difference between solutions on previous step
     real(kr8), dimension(:), allocatable :: poso !< trajectory old start position
     real(kr8), dimension(:), allocatable :: posa !< trajectory extremal position (minimum)
     real(kr8), dimension(:), allocatable :: posb !< trajectory extremal position (maximum)
     integer(ki4), dimension(:), allocatable :: posk !< trajectory sector
     type(onumerics_t) :: n !< control  parameters
  end type odes_t

end module odes_h
