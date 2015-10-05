module beq_h

  use const_kind_m
  use spl2d_m
  use spl3d_m

! public types

!> numerical control parameters
  type, public :: bnumerics_t
     integer(ki4) :: cenopt !< \f$ C_{opt} \f$, option for \f$ R_{cen}, Z_{cen} \f$
     integer(ki4) :: psiopt !< \f$ \psi_{opt} \f$, option for \f$ \psi{\min}, \psi_{\max}\f$
     integer(ki4) :: bdryopt !< \f$ m_{opt} \f$, option for limiting \f$ \psi \f$ (ltr or X-point)
     integer(ki4) :: nopt !< \f$ N_{opt} \f$, option for \f$ N_{\psi}, N_{\theta} \f$
     integer(ki4) :: thetaopt !< \f$ \theta_{opt} \f$, option for \f$ \theta_{\min}, \theta_{\max} \f$
     integer(ki4) :: zetaopt !< \f$ \zeta_{opt} \f$, option for \f$ \zeta_{\min}, \zeta_{\max} \f$
     integer(ki4) :: xiopt !< \f$ \xi_{opt} \f$, option for \f$ \xi_{\min}, \xi_{\max} \f$
     real(kr8) :: rcen !< \f$ R_{cen} \f$
     real(kr8) :: zcen !< \f$ Z_{cen} \f$
     real(kr8) :: psimin !< \f$ \psi_{\min} \f$
     real(kr8) :: psimax !< \f$ \psi_{\max} \f$
     integer(ki4) :: npsi !< \f$ N_{\psi} \f$
     real(kr8) :: delpsi !< \f$ \delta\psi \f$
     real(kr8) :: thetamin !< \f$ \theta_{\min} \f$
     real(kr8) :: thetamax !< \f$ \theta_{\max} \f$
     real(kr8) :: deltheta !< \f$ \delta\theta \f$
     real(kr8) :: zetamin !< \f$ \zeta_{\min} \f$
     real(kr8) :: zetamax !< \f$ \zeta_{\max} \f$
     real(kr8) :: ximax !< maximum \f$ \xi \f$ value
     real(kr8) :: ximin !< minimum \f$ \xi \f$ value
     real(kr8), dimension(3,2) :: xbb !< limits on quantised vectors
     real(kr8) :: rmove !< \f$ R_{mov} \f$
     real(kr8) :: zmove !< \f$ Z_{mov} \f$
     real(kr8) :: fscale !< Scale toroidal field \f$ f \f$
     real(kr8) :: psiref !< \f$ \psi_X \f$
     real(kr8) :: thetaref !< \f$ \theta_X \f$
     integer(ki4) :: ntheta !< \f$ N_{\theta} \f$
     integer(ki4) :: nzetp !< \f$ N_{\zeta P} \f$
     integer(ki4) :: fldspec !< specifies way equilibrium described (1 for J, 2 for gradpsi)
     character(len=80) :: eqtype !< type of equilibrium ('eqdsk' for EQDSK,'equ' for FIESTA, 'ana' for analytic)
     character(len=80) :: vacfile !< specifies file for spl3d format vacuum field
     integer(ki4) :: mrip !< \f$ N \f$ number of ripple coils
     real(kr8) :: irip !< unused parameter of ripple coils
     real(kr8) :: arip !< \f$ a \f$ for ripple coils
  end type bnumerics_t

!> data structure describing equilibrium field
  type, public :: beq_t
     real(kr8), dimension(:), allocatable :: f !< \f$ f(\psi) \f$ from EFIT
     real(kr8) :: rmin !< \f$ R_{\min} \f$
     real(kr8) :: rmax !< \f$ R_{\max} \f$
     integer(ki4) :: mr !< \f$ N_R \f$
     real(kr8) :: dr !< \f$ (R_{\max}-R_{\min})/N_R \f$
     real(kr8) :: zmin !< \f$ Z_{\min} \f$
     real(kr8) :: zmax !< \f$ Z_{\max} \f$
     integer(ki4) :: mz !< \f$ N_Z \f$
     real(kr8) :: dz !< \f$ (Z_{\max}-Z_{\min})/N_z \f$
     real(kr8), dimension(:), allocatable :: i !< \f$ I(\psi) \f$
     type(spl2d_t) :: psi !< \f$ \psi(R,Z) \f$ from input
     type(spl2d_t) :: dpsidr !< \f$ \frac{\partial\psi}{\partial R}(R,Z) \f$
     type(spl2d_t) :: dpsidz !< \f$ \frac{\partial\psi}{\partial Z}(R,Z) \f$
     type(spl2d_t) :: rispldr !< \f$ \frac{R}{I}\frac{\partial\psi}{\partial R}(R,Z) \f$
     type(spl2d_t) :: rispldz !< \f$ \frac{R}{I}\frac{\partial\psi}{\partial Z}(R,Z) \f$
     type(spl2d_t) :: r !< \f$ R(\psi,\theta) \f$
     type(spl2d_t) :: z !< \f$ Z(\psi,\theta) \f$
     type(spl2d_t) :: rjac !< \f$ R/J(\psi,\theta) \f$
     type(bnumerics_t) :: n !< control  parameters
     real(kr8) :: dpsi !< \f$ (\psi_{\max}-\psi_{\min})/N_\psi \f$
     real(kr8) :: psiaxis !< \f$ \psi \f$ on axis
     real(kr8) :: psibdry !< \f$ \psi \f$ at reference boundary (limiter, X-point or user)
     real(kr8) :: psiltr !< \f$ \psi \f$ at limiter
     real(kr8) :: psixpt !< \f$ \psi \f$ at X-point
     real(kr8) :: thetaxpt !< \f$ \theta \f$ at X-point
     real(kr8) :: psiotr !< Extreme value of \f$ \psi \f$ for geometry which is not limiter extremum
     real(kr8) :: psinorm !< normalisation of \f$ \psi \f$, i.e. a representative magnitude
     real(kr8) :: psiqbdry !< \f$ \psi \f$ at boundary set by EQDSK
     real(kr8) :: thetagmax !< maximum \f$ \theta \f$ value of geometry
     real(kr8) :: thetagmin !< minimum \f$ \theta \f$ value of geometry
     real(kr8) :: rqcen !< \f$ R \f$ at central extremum of \f$ \psi \f$
     real(kr8) :: zqcen !< \f$ Z \f$ at central extremum of \f$ \psi \f$
     real(kr8) :: rbdry !< \f$ R \f$ at reference boundary
     real(kr8) :: bpbdry !< \f$ B_p \f$ at reference boundary
     real(kr8) :: btotbdry !< Total \f$ B \f$ at reference boundary
     real(kr8) :: ximaxm !< maximum \f$ \xi \f$ value \f$ -\epsilon \f$
     real(kr8) :: ximinp !< minimum \f$ \xi \f$ value \f$ +\epsilon \f$
     real(kr8), dimension(3) :: domlen !< dimensions of domain (quantised)
     integer(ki4) :: nzets !< \f$ N_{\zeta s} \f$
     integer(ki4) :: ntmin !< minimum valid value of \f$ N_\theta \f$
     integer(ki4) :: ntmax !< maximum valid value of \f$ N_\theta \f$
     real(kr8) :: dtheta !< \f$ (\theta_{\max}-\theta_{\min})/N_\theta \f$
     real(kr8), dimension(:), allocatable :: srmin !< type variable
 &   !< \f$ r_{\min}(\theta) \f$, set of \f$ r \f$  corresponding to \f$ \psi_{\min} \f$
     real(kr8), dimension(:), allocatable :: srmax !< type variable
 &   !< \f$ r_{\max}(\theta) \f$, set of \f$ r \f$  corresponding to \f$ \psi_{\max} \f$
     real(kr8) :: ivac !< \f$ I \f$ for vacuum field
     type(spl3d_t)  :: vacfld  !< vacuum field structure
  end type beq_t

! public variables
  integer(ki4), public, parameter :: beq_spline_order=4 !< set to four for cubic splines (use of tc01a assumes cubic splines)
  logical, public :: beq_nobinq !< flag whether \f$ B \f$ data in EQDSK file
  logical, parameter :: BEQ_OVERRIDE_AFWS=.TRUE. !< use values required by AFWS project
  logical, parameter :: BEQ_OVERRIDE_FTU=.TRUE. !< use values required by FTU project

end module beq_h
