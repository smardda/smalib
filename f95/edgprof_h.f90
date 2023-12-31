module edgprof_h

  use const_kind_m


! public types

!> type storing edge profile data
  type, public :: edgprof_t
     character(len=80), dimension(:), allocatable  :: formula !< profile formula
     real(kr8) :: f !< power split
     real(kr8), dimension(:), allocatable :: lmid !< decay length
     real(kr8), dimension(:), allocatable :: ploss!< power loss
     real(kr8), dimension(:), allocatable :: sigma !< diffusion length
     real(kr8) :: qpara0 !< parallel power flux (W/sqm)
     real(kr8), dimension(:), allocatable :: lmidnr !< decay length near SOL
     real(kr8) :: rqpara0 !< ratio of parallel power flux near SOL
     integer(ki4) :: npos !< number of positions
     real(kr8), dimension(:), allocatable   :: pos !< positions
     real(kr8), dimension(:), allocatable   :: prof !< deposition profile
     real(kr8) :: fint !< integral over profile
     character(len=80) :: postype !< type of positions
     integer(ki4) :: nrpams !< number of real parameters
     integer(ki4) :: nipams !< number of integer parameters
     real(kr8), dimension(:), allocatable   :: rpar !< general real parameters
     integer(ki4), dimension(:), allocatable   :: npar !< general integer parameters

     real(kr8), dimension(:), allocatable :: rblfac !< factor \f$ \frac{1}{2 \pi \lambda_{mid} R_m B_{pm}} \f$
     real(kr8), dimension(:), allocatable :: fpfac !< factor \f$ \frac{F P_{loss}}}{2 \pi \lambda_{mid} R_m B_{pm}} \f$
     real(kr8), dimension(:), allocatable :: rblfacnr !< factor \f$ \frac{1}{2 \pi \lambda_{mid,near} R_m B_{pm}} \f$
     real(kr8), dimension(:), allocatable :: fpfacnr !< factor \f$ \frac{F P_{loss}}}{2 \pi \lambda_{mid,near} R_m B_{pm}} \f$
     real(kr8), dimension(:), allocatable :: slfac !< factor \f$ \frac{\sigma}{2 \lambda_{mid}} \f$

  end type edgprof_t

end module edgprof_h
