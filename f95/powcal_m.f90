module powcal_m

  use const_kind_m
  use const_numphys_h
  use control_h
  use log_m
  use position_h
  use spl2d_m
  use spl3d_m
  use powcal_h
  use powres_h
!     use powres_m
  use geobjlist_h
  use beq_h
  use vfile_m
  use position_m
  use beq_m
  use btree_m
  use powelt_m
  use geobjlist_m
  use nrsolve_m
  use pcontrol_h
  use termplane_h
  use termplane_m
  use edgprof_h
  use edgprof_m

  implicit none
  private

! public subroutines
  public :: &
  powcal_readv, &   !< read (vtk) powcal data
  powcal_init, &   !< initialise  powcal data structure
  powcal_readcon, &   !< namelist input of controls
  powcal_move, &   !< loop over entries in results
  powcal_calc, &   !< coordinate calculation of power deposition
  powcal_writev, &   !< write (vtk) powcal data
  powcal_refine, &   !< set refinement level
  powcal_write, &  !< output powcal data
  powcal_xferpow, & !< transfer power to shadow wall geometry
  powcal_delete  !< delete powcal data

!public variables

! private types

! private variables
  character(*), parameter :: m_name='powcal_m' !< module name
  integer   :: status   !< error status
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer(ki4) :: level=1   !< refinement level, dynamic if inlevel>1
  integer(ki4) :: inlevel=1   !< number of subelements at refinement level
  integer(ki4) :: imlevel=1   !< number of subelements at previous refinement level +1
  integer(ki4) :: nin   !< input channel for geobj data
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: idum !< dummy integer
  real(kr8) :: zdum !< dummy real
  logical :: iltest !< logical flag
  integer(ki4) :: infilelevel   !< max refinement level from input

  contains
!---------------------------------------------------------------------
!> read (vtk) powcal data
subroutine powcal_readv(self,vtkfile)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  character(*), intent(in)       :: vtkfile !< name of input file


  !! local
  character(*), parameter :: s_name='powcal_readv' !< subroutine name

  call geobjlist_read(self%powres%geobjl,vtkfile)
  idum=0
  call position_readveclis(self%powres%vecx,vtkfile,'Xcart',idum,1)
  self%powres%vecb%np=self%powres%vecx%np
  idum=0
  call position_readveclis(self%powres%vecb,vtkfile,'Bcart',idum,1)

end subroutine powcal_readv
!---------------------------------------------------------------------
!> initialise  powcal data structure
subroutine powcal_init(self,numerics,plot,gnumerics)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< object data structure
  type(pnumerics_t), intent(in) :: numerics   !< object control data structure
  type(pplots_t), intent(in) :: plot   !< vtk plot selectors
  type(numerics_t), intent(inout), optional :: gnumerics   !< object control data structure

  !! local
  character(*), parameter :: s_name='powcal_init' !< subroutine name
  real(kr8) :: zsign !< sign factor

  ! initialise control structure
  self%n=numerics

  ! set sign of psi coordinate transform
  call beq_rsigset(self%powres%beq)

  ! initialise powelt data
  call powelt_init

  ! check refinement level
  infilelevel=max(1,self%powres%geobjl%nparam)
  call log_value("data file has refinement level ",infilelevel)
  call log_value("input demands refinement level ",self%n%nlevel)
  if (self%n%nlevel>infilelevel) then
     call log_error(m_name,s_name,1,error_fatal,'Requested refinement level too high')
  end if
  ! set size of pow array
  allocate(self%powres%pow(self%powres%geobjl%ng), stat=status)
  call log_alloc_check(m_name,s_name,1,status)

  ! calculate number of level 1 triangles
  self%powres%npowe=self%powres%geobjl%ng/powelt_table(infilelevel,2)

  ! set size of power statistics arrays
  allocate(self%powres%pows(self%powres%npowe), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  allocate(self%powres%powa(self%powres%npowe), stat=status)
  call log_alloc_check(m_name,s_name,3,status)

  zsign=beq_rsig()
  call edgprof_factors(self%edgprof,&
 &self%powres%beq%rbdry,self%powres%beq%bpbdry,self%powres%beq%btotbdry,zsign)
  ! track controls
  self%powres%flincart=plot%flincart
  self%powres%flinends=plot%flinends
  self%powres%flinptz=plot%flinptz

  calcn_type: select case (self%n%caltype)
  case('msus')
     ! set size of psista array
     allocate(self%powres%psista(self%powres%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,4,status)
  case('afws')
     ! limits on psi-theta interpolation
     call spl2d_ptlimits(self%powres%beq%rjac,self%powres%psimin,self%powres%psimax,&
     self%powres%thetamin,self%powres%thetamax)
  case default
  end select calcn_type

  if (self%n%shadow/=0) then
     if (self%n%shadow==-1) then
        ! testing, make quantising transform the identity
        gnumerics%geobj_coord_tfm%hmin=1
        gnumerics%geobj_coord_tfm%rhmin=1
        gnumerics%geobj_coord_tfm%offvec=0
        gnumerics%geobj_coord_tfm%nqtfm=1
        ! and position transform
        gnumerics%position_coord_tfm%scale=1
        gnumerics%position_coord_tfm%offset=0
        gnumerics%position_coord_tfm%matrix(1,:)=(/1,0,0/)
        gnumerics%position_coord_tfm%matrix(2,:)=(/0,1,0/)
        gnumerics%position_coord_tfm%matrix(3,:)=(/0,0,1/)
        gnumerics%position_coord_tfm%ntfm=1
     end if
     ! quantise relevant quantities
     call powcal_quant(self,gnumerics)
  end if

end subroutine powcal_init
!---------------------------------------------------------------------
!> quantise relevant quantities
subroutine powcal_quant(self,gnumerics)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< object data structure
  type(numerics_t), intent(in) :: gnumerics   !< object control data structure

  !! local
  character(*), parameter :: s_name='powcal_quant' !< subroutine name
  type(posveclis_t) :: zposl !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zposq !< local variable
  integer(ki4) :: inqtfm   !< save transform type
  integer(ki4) :: inkn1   !< number of knots in 1-direction
  integer(ki4) :: inkn2   !< number of knots in 2-direction
  integer(ki4) :: inkn   !< maximum number of knots in either direction
  real(kr8), dimension(3)  :: zpar !< scaling parameters
  real(kr8) :: zh1 !< scaling parameter
  real(kr8) :: zh2 !< scaling parameter
  real(kr8) :: zh3 !< scaling parameter
  real(kr8) :: zo1 !< offset parameter
  real(kr8) :: zo3 !< offset parameter
  real(kr8) :: zepsdif !<  \f$ \epsilon \f$ for \f$ \zeta \f$ boundary crossing
  real(kr8) :: zrnns !< \f$ \zeta \f$ scaling
  real(kr8) :: zarnns  !< \f$ \zeta \f$ component scaling
  integer(ki4) :: idir !< coordinate direction \f$ (R,Z,\xi) \f$ in order \f$ 1,2,3 \f$

  zpar=0
  self%powres%beq%domlen=0

  self%powres%geobjl%tfmdata=gnumerics%position_coord_tfm
  self%powres%geobjl%quantfm=gnumerics%geobj_coord_tfm
  ! rotate and quantise results geometry positions
  call position_tfmlis(self%powres%geobjl%posl,self%powres%geobjl%tfmdata)
  call position_qtfmlis(self%powres%geobjl%posl,self%powres%geobjl%quantfm)

  calcn_type: select case (self%n%caltype)
  case('afws')
     ! psi-theta-zeta case
     ! multiplies RHS of ode
     self%powres%qfac=self%powres%geobjl%quantfm%hmin(2)/self%powres%geobjl%quantfm%hmin(3)

     ! limits on psi-theta interpolation, also timestep control
     allocate(zposl%pos(6),stat=status)
     call log_alloc_check(m_name,s_name,10,status)
     zposl%pos(1)%posvec=0 ! scaled below
     zposl%pos(2)%posvec=(/self%powres%beq%rjac%org1, self%powres%beq%rjac%org2, 0._kr8/)
     zposl%pos(3)%posvec=(/self%powres%psimin, self%powres%thetamin, 0._kr8/)
     zposl%pos(4)%posvec=(/self%powres%psimax, self%powres%thetamax, self%odes%n%tmax/)
     zposl%pos(5)%posvec=(/self%powres%beq%psiaxis, 0._kr8, 0._kr8/)
     zposl%pos(6)%posvec=(/self%powres%beq%psiqbdry, 0._kr8, 0._kr8/)
     zposl%np=6
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     self%powres%beq%rjac%org1=zposl%pos(2)%posvec(1)
     self%powres%beq%rjac%org2=zposl%pos(2)%posvec(2)
     self%powres%psimin=zposl%pos(3)%posvec(1)
     self%powres%thetamin=zposl%pos(3)%posvec(2)
     self%powres%psimax=zposl%pos(4)%posvec(1)
     self%powres%thetamax=zposl%pos(4)%posvec(2)
     self%powres%beq%domlen(1)=self%powres%psimax-self%powres%psimin
     self%powres%beq%domlen(2)=self%powres%thetamax-self%powres%thetamin
     self%odes%n%tmax=zposl%pos(4)%posvec(3)
     self%powres%beq%psiaxis=zposl%pos(5)%posvec(1)
     self%powres%beq%psiqbdry=zposl%pos(6)%posvec(1)

     zpos%posvec=(/self%powres%beq%rjac%h1, self%powres%beq%rjac%h2 , self%odes%n%dt0/)
     inqtfm=self%powres%geobjl%quantfm%nqtfm
     ! case 1 just scales
     self%powres%geobjl%quantfm%nqtfm=1
     zposq=position_qtfm(zpos,self%powres%geobjl%quantfm)
     self%powres%beq%rjac%h1=zposq%posvec(1)
     self%powres%beq%rjac%rh1=1/zposq%posvec(1)
     self%powres%beq%rjac%rh2=1/zposq%posvec(2)
     self%powres%beq%rjac%h2=zposq%posvec(2)
     self%odes%n%dt0=zposq%posvec(3)
     self%powres%geobjl%quantfm%nqtfm=inqtfm
     deallocate(zposl%pos)

     ! now knots
     inkn1=self%powres%beq%rjac%n1p+self%powres%beq%rjac%nord
     inkn2=self%powres%beq%rjac%n2p+self%powres%beq%rjac%nord

     inkn=max(inkn1,inkn2)
     allocate(zposl%pos(inkn),stat=status)
     call log_alloc_check(m_name,s_name,20,status)
     do j=1,inkn
        zposl%pos(j)%posvec=0
     end do
     do j=1,inkn1
        zposl%pos(j)%posvec(1)=self%powres%beq%rjac%knot1(j)
     end do
     do j=1,inkn2
        zposl%pos(j)%posvec(2)=self%powres%beq%rjac%knot2(j)
     end do
     zposl%np=inkn
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     do j=1,inkn1
        self%powres%beq%rjac%knot1(j)=zposl%pos(j)%posvec(1)
     end do
     do j=1,inkn2
        self%powres%beq%rjac%knot2(j)=zposl%pos(j)%posvec(2)
     end do
     deallocate(zposl%pos)

  case default
     ! R-Z-xi case

     ! limits on R-Z interpolation
     call spl2d_ptscale(self%powres%beq%psi,self%powres%geobjl%quantfm)
     call spl2d_ptscale(self%powres%beq%dpsidr,self%powres%geobjl%quantfm)
     call spl2d_ptscale(self%powres%beq%dpsidz,self%powres%geobjl%quantfm)
     fldspec_type1: select case (self%powres%beq%n%fldspec)
     case (3)
        call spl2d_ptscale(self%powres%beq%rispldr,self%powres%geobjl%quantfm)
        call spl2d_ptscale(self%powres%beq%rispldz,self%powres%geobjl%quantfm)
     end select fldspec_type1

     ! domain control
     allocate(zposl%pos(5),stat=status)
     call log_alloc_check(m_name,s_name,30,status)
     zposl%pos(1)%posvec=(/0._kr8, 0._kr8, self%powres%beq%n%ximin/)
     zposl%pos(2)%posvec=(/0._kr8, 0._kr8, self%powres%beq%n%ximax/)
     zepsdif=const_epsbdry*(self%powres%beq%n%ximax-self%powres%beq%n%ximin)
     zposl%pos(3)%posvec=(/0._kr8, 0._kr8, self%powres%beq%n%ximin+zepsdif/)
     zposl%pos(4)%posvec=(/0._kr8, 0._kr8, self%powres%beq%n%ximax-zepsdif/)
     zposl%pos(5)%posvec=(/0._kr8, self%odes%n%termcon(3), 0._kr8/)
     zposl%np=5
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     self%powres%beq%n%ximin=zposl%pos(1)%posvec(3)
     self%powres%beq%n%ximax=zposl%pos(2)%posvec(3)
     self%powres%beq%domlen(3)=self%powres%beq%n%ximax-self%powres%beq%n%ximin
     self%powres%beq%ximinp=zposl%pos(3)%posvec(3)
     self%powres%beq%ximaxm=zposl%pos(4)%posvec(3)
     self%odes%n%termcon(3)=zposl%pos(5)%posvec(2)
     !
     inqtfm=self%powres%geobjl%quantfm%nqtfm
     ! case 1 only scales, as required here
     self%powres%geobjl%quantfm%nqtfm=1

     zpar(3)=self%powres%beq%nzets
     zposl%pos(1)%posvec=(/0._kr8, 0._kr8, zpar(3)*self%odes%n%tmax/)
     zposl%pos(2)%posvec=(/0._kr8, 0._kr8, zpar(3)*self%odes%n%dt0/)
     zposl%np=2
     call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
     self%odes%n%tmax=zposl%pos(1)%posvec(3)
     self%odes%n%dt0=zposl%pos(2)%posvec(3)
     ! restore
     self%powres%geobjl%quantfm%nqtfm=inqtfm
     deallocate(zposl%pos)

     if (self%n%ltermplane) then
        ! scale termination planes
        allocate(zposl%pos(self%n%termp%ntermplane),stat=status)
        call log_alloc_check(m_name,s_name,31,status)
        do j=1,self%n%termp%ntermplane
           zposl%pos(j)%posvec=0
           idir=self%n%termp%termplanedir(j,1)
           zposl%pos(j)%posvec(idir)=self%n%termp%termplane(j)
        end do
        zposl%np=self%n%termp%ntermplane
        call position_qtfmlis(zposl,self%powres%geobjl%quantfm)
        do j=1,self%n%termp%ntermplane
           idir=self%n%termp%termplanedir(j,1)
           self%n%termp%termplane(j)=zposl%pos(j)%posvec(idir)
        end do
        deallocate(zposl%pos)
     end if

     ! multiplies RHS of ode
     self%powres%qfaca=self%powres%geobjl%quantfm%hmin(3)/self%powres%geobjl%quantfm%hmin

     if (self%powres%beq%n%vacfile=='null') then
        if (self%powres%beq%n%mrip/=0) then
           ! set ripple parameters (before scaling)
           zdum=beq_ripple_h1(1._kr8,0._kr8,0._kr8,-1,&
 &         self%powres%beq%n%mrip,self%powres%beq%ivac,self%powres%beq%n%arip)
           ! set scale function arguments
           zh1=self%powres%geobjl%quantfm%hmin(1)
           zh2=self%powres%geobjl%quantfm%hmin(2)
           zh3=self%powres%geobjl%quantfm%hmin(3)
           zdum=beq_ripple_h1(1._kr8,0._kr8,0._kr8,-1,-1,zh1,zh3)
           zo1=self%powres%geobjl%quantfm%offvec(1)
           zo3=self%powres%geobjl%quantfm%offvec(3)
           zdum=beq_ripple_h1(1._kr8,0._kr8,0._kr8,-1,-2,zo1,zo3)
           zrnns=real(self%powres%beq%n%mrip,kr8)/self%powres%beq%nzets
           zarnns=self%powres%beq%n%arip*zrnns
           zdum=beq_ripple_h1(1._kr8,0._kr8,0._kr8,-1,-3,zrnns,zarnns)
           call nrsolve_spl2dfn(self%powres%beq%psi,self%powres%beq%dpsidr,&
 &         self%powres%beq%dpsidz,zdum,&
 &         zh1,zh2,idum,0)
        else
           fldspec_type: select case (self%powres%beq%n%fldspec)
           case (3)
              ! no 3-cpt so in fact value of qfaca(3) does not matter
           case default
              ! working with f=I in 3-component, so need 1/h_1 factor for extra R
              self%powres%qfaca(3)=1._kr8/self%powres%geobjl%quantfm%hmin(1)
           end select fldspec_type
        end if
     else
        ! ripple defined by spl3d vacuum field
        call spl3d_ptscale(self%powres%beq%vacfld,self%powres%geobjl%quantfm)
        ! here scale consistent with xi update
        self%powres%qfaca(3)=self%powres%beq%nzets
     end if

     ! normalise termination controls
     self%odes%n%termcon(2)=self%odes%n%termcon(2)*self%powres%beq%domlen(3)

  end select calcn_type

end subroutine powcal_quant
!---------------------------------------------------------------------
!> namelist input of controls
subroutine powcal_readcon(selfn,kin)

  !! arguments
  type(pnumerics_t), intent(out) :: selfn   !< object control data structure
  integer(ki4), intent(in) :: kin   !< input channel for object data structure


  !! local
  character(*), parameter :: s_name='powcal_readcon' !< subroutine name
  real(kr8) :: power_split !< power split ion to electron direction
  real(kr8) :: decay_length !< power decay length at outer midplane (m)
  real(kr8) :: power_loss !< power crossing the separatrix (W)
  real(kr8) :: diffusion_length !< (m)
  real(kr8) :: q_parallel0 !< parallel reference flux (W/sqm)
  integer(ki4) :: refine_level   !< level of refinement to be used
  integer(ki4) :: shadow_control   !< control shadowing
  character(len=80) :: calculation_type !< calculation type
  logical :: termination_planes !<  termination planes present
  integer(ki4) :: analytic_launch_type !< analytic launch type
  integer(ki4) :: power_on_launch_geo !< power on launch geometry
  integer(ki4) :: power_on_shadow_geo !< power on shadow geometry
  !! powcal parameters
  namelist /powcalparameters/ &
 &power_split, decay_length, power_loss, refine_level, shadow_control, &
 &q_parallel0, &
 &diffusion_length, &
 &calculation_type, &
 &termination_planes, &
 &analytic_launch_type, power_on_launch_geo, power_on_shadow_geo

  !! set default powcal parameters
  power_split=0.5_kr8
  decay_length=.012_kr8
  diffusion_length=0
  power_loss=7.5e+06_kr8
  refine_level=1
  shadow_control=0
  !! ignore q_parallel0 if negative
  q_parallel0=-1_kr8

  calculation_type='unset'
  termination_planes=.FALSE.
  analytic_launch_type=0
  power_on_launch_geo=0
  power_on_shadow_geo=0

  !!read powcal parameters
  read(kin,nml=powcalparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading powcal parameters')
     print '("Fatal error reading powcal parameters")'
  end if

  !! check for valid data
  if(power_split<0.OR.power_split>1) &
 &call log_error(m_name,s_name,11,error_fatal,'power_split must be >=0 and <=1')
  if(decay_length<=0) &
 &call log_error(m_name,s_name,12,error_fatal,'decay_length (m) must be >0')
  if(diffusion_length<0) &
 &call log_error(m_name,s_name,13,error_fatal,'diffusion_length (m) must be >=0')
  if(q_parallel0<=0.AND.power_loss<=0) &
 &call log_error(m_name,s_name,14,error_fatal,'power_loss (W) must be >0')
  ! need to check against value in input file
  if(refine_level<1.OR.refine_level>3) &
 &call log_error(m_name,s_name,15,error_fatal,'refine_level must be 1, 2 or 3')

  if(analytic_launch_type<0.OR.analytic_launch_type>10) &
 &call log_error(m_name,s_name,23,error_fatal,'analytic_launch_type must be >=0 and <=10')
  if(power_on_launch_geo<0.OR.power_on_launch_geo>10) &
 &call log_error(m_name,s_name,24,error_fatal,'power_on_launch_geo must be >=0 and <=10')
  if(power_on_shadow_geo<0.OR.power_on_shadow_geo>10) &
 &call log_error(m_name,s_name,25,error_fatal,'power_on_shadow_geo must be >=0 and <=10')

  !! store values
  selfn%f=power_split
  selfn%lmid=decay_length
  selfn%ploss=power_loss
  selfn%sigma=diffusion_length
  selfn%nlevel=refine_level
  selfn%shadow=shadow_control
  selfn%qpara0=q_parallel0

  selfn%caltype=calculation_type
  selfn%ltermplane=termination_planes
  selfn%nanalau=analytic_launch_type
  selfn%nlaupow=power_on_launch_geo
  selfn%nshapow=power_on_shadow_geo

  if (selfn%ltermplane) then
     call termplane_readcon(selfn%termp,kin)
  end if

end subroutine powcal_readcon
!---------------------------------------------------------------------
!> loop over entries in results
subroutine powcal_move(self,gshadl,btree)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(geobjlist_t), intent(inout) :: gshadl   !< shadow data structure
  type(btree_t), intent(inout) :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powcal_move' !< subroutine name
  type(powelt_t) :: zelt   !< power element

  ! check for shadowing
  if (self%powres%beq%n%vacfile=='null') then
     if (self%powres%beq%n%mrip/=0) then
        ! analytic ripple model
        do i=1,self%powres%npowe
           zelt%ie=i
           do j=imlevel,inlevel
              zelt%je=j
              call powelt_move2(zelt,self,gshadl,btree)
           end do
        end do
     else
        ! no ripple
        calcn_type: select case (self%powres%beq%n%fldspec)
        case (3)
           do i=1,self%powres%npowe
              zelt%ie=i
              do j=imlevel,inlevel
                 zelt%je=j
                 call powelt_move1(zelt,self,gshadl,btree)
              end do
           end do
        case default
           do i=1,self%powres%npowe
              zelt%ie=i
              do j=imlevel,inlevel
                 zelt%je=j
                 call powelt_move(zelt,self,gshadl,btree)
              end do
           end do
        end select calcn_type
     end if
  else if (self%n%ltermplane) then
     ! ripple defined by vacuum file
     do i=1,self%powres%npowe
        zelt%ie=i
        do j=imlevel,inlevel
           zelt%je=j
           call powelt_move4(zelt,self,gshadl,btree)
        end do
     end do
  else
     ! ripple defined by vacuum file
     do i=1,self%powres%npowe
        zelt%ie=i
        do j=imlevel,inlevel
           zelt%je=j
           call powelt_move3(zelt,self,gshadl,btree)
        end do
     end do
  end if

end subroutine powcal_move
!---------------------------------------------------------------------
!> coordinate calculation of power deposition (=powres_XX)
subroutine powcal_calc(self,gshadl)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure

  !! local
  character(*), parameter :: s_name='powcal_calc' !< subroutine name
  type(powelt_t) :: zelt   !< powelt data instance

  ! calculate power deposition
  do i=1,self%powres%npowe
     zelt%ie=i
     do j=imlevel,inlevel
        zelt%je=j
        call powelt_dep(zelt,self,gshadl)
     end do
  end do

end subroutine powcal_calc
!---------------------------------------------------------------------
!> coordinate calculation of power statistics (=powres_XX)
subroutine powcal_stat(self,gshadl)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure

  !! local
  character(*), parameter :: s_name='powcal_stat' !< subroutine name
  type(powelt_t) :: zelt   !< powelt data instance

  ! calculate power fluctuation statistic
  zelt%je=inlevel
  do i=1,self%powres%npowe
     zelt%ie=i
     call powelt_stat(zelt,self)
  end do
  ! calculate power average
  zelt%je=inlevel
  do i=1,self%powres%npowe
     zelt%ie=i
     call powelt_avg(zelt,self)
  end do

end subroutine powcal_stat
!---------------------------------------------------------------------
!> write (vtk) powcal data
subroutine powcal_writev(self,kchar,kplot)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  character(*), intent(in) :: kchar !< case
  integer(ki4), intent(inout) :: kplot   !< output channel for vis. data

  !! local
  character(*), parameter :: s_name='powcal_writev' !< subroutine name
  type(posvecl_t), dimension(:), allocatable :: workpos !< save objlist variable
  integer(ki4), dimension(:), allocatable :: work !< save objlist variable
  integer(ki4) :: ing   !< save objlist variable

  plot_type: select case (kchar)
  case('psi-theta-zeta')
     ! rearrange nodl, keep copy of structure in work and ing
     call log_alloc_check(m_name,s_name,11,status)
     ing=self%powres%geobjl%ng
     ! ??reset number of objects to be output according to refine_level statistics input
     self%powres%geobjl%ng=ing*powelt_table(self%powres%n%nlevel,2)/powelt_table(infilelevel,2)
     self%powres%geobjl%ng=self%powres%npowe
     allocate(work(3*self%powres%geobjl%ng), stat=status)
     work=self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)
     call geobjlist_nodlmv(self%powres%geobjl,infilelevel,1,self%powres%npowe)
     call geobjlist_writev(self%powres%geobjl,'geometry',kplot)
     ! restore nodl
     self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)=work
     deallocate(work)
     call vfile_rscalarwrite(self%powres%powa,self%powres%geobjl%ng,'Q-avg','CELL',kplot,0)
     call vfile_rscalarwrite(self%powres%pows,self%powres%geobjl%ng,'Q-dev','CELL',kplot,0)
     self%powres%geobjl%ng=ing

  case('allpsi-theta-zeta')
     ! reset number of objects to be output according to refine_level input
     ing=self%powres%geobjl%ng
     self%powres%geobjl%ng=ing*powelt_table(self%powres%n%nlevel,2)/powelt_table(infilelevel,2)
     if (infilelevel>self%powres%n%nlevel) then
        allocate(work(3*self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name,s_name,31,status)
        work=self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)
        call geobjlist_nodlmv(self%powres%geobjl,infilelevel,self%powres%n%nlevel,self%powres%npowe)
     end if
     call geobjlist_writev(self%powres%geobjl,'geometry',kplot)
     call vfile_rscalarwrite(abs(self%powres%pow),self%powres%geobjl%ng,'Q','CELL',kplot,0)
     if (infilelevel>self%powres%n%nlevel) then
        ! restore nodl
        self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)=work
     end if
     ! restore number of objects
     self%powres%geobjl%ng=ing

  case('cartesian')
     ! replace R-Z-xi positions with Cartesian
     allocate(workpos(self%powres%geobjl%np), stat=status)
     call log_alloc_check(m_name,s_name,21,status)
     workpos=self%powres%geobjl%posl%pos
     self%powres%geobjl%posl%pos=self%powres%vecx%pos
     ! rearrange nodl, keep copy of structure in work and ing
     ing=self%powres%geobjl%ng
     ! reset number of objects to be output according to refine_level statistics input
     self%powres%geobjl%ng=ing*powelt_table(self%n%nlevel,2)/powelt_table(infilelevel,2)
     self%powres%geobjl%ng=self%powres%npowe
     allocate(work(3*self%powres%geobjl%ng), stat=status)
     call log_alloc_check(m_name,s_name,22,status)
     work=self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)
     call geobjlist_nodlmv(self%powres%geobjl,infilelevel,1,self%powres%npowe)
     call geobjlist_writev(self%powres%geobjl,'geometry',kplot)
     self%powres%geobjl%posl%pos=workpos
     deallocate(workpos)
     ! restore nodl
     self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)=work
     deallocate(work)
     call vfile_rscalarwrite(self%powres%powa,self%powres%geobjl%ng,'Q-avg','CELL',kplot,0)
     call vfile_rscalarwrite(self%powres%pows,self%powres%geobjl%ng,'Q-dev','CELL',kplot,0)
     self%powres%geobjl%ng=ing

  case('allcartesian')
     ! replace R-Z-xi positions with Cartesian
     allocate(workpos(self%powres%geobjl%np), stat=status)
     call log_alloc_check(m_name,s_name,40,status)
     workpos=self%powres%geobjl%posl%pos
     self%powres%geobjl%posl%pos=self%powres%vecx%pos
     ! reset number of objects to be output according to refine_level input
     ing=self%powres%geobjl%ng
     self%powres%geobjl%ng=ing*powelt_table(self%n%nlevel,2)/powelt_table(infilelevel,2)
     if (infilelevel>self%n%nlevel) then
        allocate(work(3*self%powres%geobjl%ng), stat=status)
        call log_alloc_check(m_name,s_name,41,status)
        work=self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)
        call geobjlist_nodlmv(self%powres%geobjl,infilelevel,self%n%nlevel,self%powres%npowe)
     end if
     call geobjlist_writev(self%powres%geobjl,'geometry',kplot)
     call vfile_rscalarwrite(self%powres%pow,self%powres%geobjl%ng,'Qs','CELL',kplot,0)
     if (allocated(self%powres%psista)) then
        call vfile_rscalarwrite(self%powres%psista,self%powres%geobjl%ng,'psista','CELL',kplot,0)
     end if
     if (infilelevel>self%n%nlevel) then
        ! restore nodl
        self%powres%geobjl%nodl(1:3*self%powres%geobjl%ng)=work
     end if
     self%powres%geobjl%posl%pos=workpos
     deallocate(workpos)
     ! restore number of objects
     self%powres%geobjl%ng=ing

  end select plot_type

end subroutine powcal_writev
!---------------------------------------------------------------------
!> set refinement level
subroutine powcal_refine(self,gshadl,btree)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(geobjlist_t), intent(inout), optional :: gshadl   !< shadow data structure
  type(btree_t), intent(inout), optional :: btree !< btree data

  !! local
  character(*), parameter :: s_name='powcal_refine' !< subroutine name

  calcn_type: select case (self%n%caltype)
  case('afws')
  case('msum')
     ! need launch and shadow obj arrays  (=weight)s
     allocate(self%powres%geobjl%obj(self%powres%geobjl%ng), stat=status)
     self%powres%geobjl%nwset=2
     allocate(gshadl%obj(gshadl%ng), stat=status)
     gshadl%nwset=2
     call log_alloc_check(m_name,s_name,1,status)
  case('msus')
     ! need launch obj array (=weight)
     allocate(self%powres%geobjl%obj(self%powres%geobjl%ng), stat=status)
     self%powres%geobjl%nwset=2
     call log_alloc_check(m_name,s_name,2,status)
  end select calcn_type

  imlevel=1
  do level=1,self%n%nlevel
     inlevel=powelt_table(level,2)
     if (self%n%shadow/=0) then
        call powcal_move(self,gshadl,btree)
     end if
     call powcal_calc(self,gshadl)
     imlevel=inlevel+1
  end do

  call powcal_stat(self,gshadl)

end subroutine powcal_refine
!---------------------------------------------------------------------
!> write powcal data
subroutine powcal_write(self,kout)

  !! arguments
  type(powcal_t), intent(in) :: self   !< powcal data structure
  integer(ki4), intent(in) :: kout   !< output channel for powcal data structure

  !! local
  character(*), parameter :: s_name='powcal_write' !< subroutine name

  write(kout,*,iostat=status) 'npowe'
  call log_write_check(m_name,s_name,1,status)
  write(kout,*,iostat=status) self%powres%npowe
  call log_write_check(m_name,s_name,2,status)
  write(kout,*,iostat=status) 'npow'
  call log_write_check(m_name,s_name,3,status)
  write(kout,*,iostat=status) self%powres%geobjl%ng
  call log_write_check(m_name,s_name,4,status)
  !X      write(kout,*,iostat=status) 'rblfac'
  !X      call log_write_check(m_name,s_name,5,status)
  !X      write(kout,*,iostat=status) self%powres%rblfac
  !X      call log_write_check(m_name,s_name,5,status)
  !X      write(kout,*,iostat=status) 'fpfac'
  !X      call log_write_check(m_name,s_name,6,status)
  !X      write(kout,*,iostat=status) self%powres%fpfac
  !X      call log_write_check(m_name,s_name,7,status)
  !
  write(kout,*,iostat=status) 'pows'
  call log_write_check(m_name,s_name,10,status)
  write(kout,*,iostat=status) self%powres%pows
  call log_write_check(m_name,s_name,11,status)
  write(kout,*,iostat=status) 'pow'
  call log_write_check(m_name,s_name,12,status)
  write(kout,*,iostat=status) self%powres%pow
  call log_write_check(m_name,s_name,13,status)

  !X      write(kout,*,iostat=status) 'f'
  !X      call log_write_check(m_name,s_name,20,status)
  !X      write(kout,*,iostat=status) self%n%f
  !X      call log_write_check(m_name,s_name,21,status)
  !X      write(kout,*,iostat=status) 'lmid'
  !X      call log_write_check(m_name,s_name,22,status)
  !X      write(kout,*,iostat=status) self%n%lmid
  !X      call log_write_check(m_name,s_name,23,status)
  !X      write(kout,*,iostat=status) 'ploss'
  !X      call log_write_check(m_name,s_name,24,status)
  !X      write(kout,*,iostat=status) self%n%ploss
  !X      call log_write_check(m_name,s_name,25,status)
  write(kout,*,iostat=status) 'nlevel'
  call log_write_check(m_name,s_name,26,status)
  write(kout,*,iostat=status) self%n%nlevel
  call log_write_check(m_name,s_name,27,status)
  write(kout,*,iostat=status) 'shadow'
  call log_write_check(m_name,s_name,28,status)
  write(kout,*,iostat=status) self%n%shadow
  call log_write_check(m_name,s_name,29,status)
  write(kout,*,iostat=status) 'start_control'
  call log_write_check(m_name,s_name,30,status)
  write(kout,*,iostat=status) self%odes%n%nstartcon
  call log_write_check(m_name,s_name,31,status)
  write(kout,*,iostat=status) 'termination_control'
  call log_write_check(m_name,s_name,32,status)
  write(kout,*,iostat=status) self%odes%n%ntermcon
  call log_write_check(m_name,s_name,33,status)
  write(kout,*,iostat=status) 'termination_parameters'
  call log_write_check(m_name,s_name,34,status)
  write(kout,*,iostat=status) self%odes%n%termcon
  call log_write_check(m_name,s_name,35,status)
  write(kout,*,iostat=status) 'analytic_launch_type'
  call log_write_check(m_name,s_name,36,status)
  write(kout,*,iostat=status) self%n%nanalau
  call log_write_check(m_name,s_name,37,status)
  write(kout,*,iostat=status) 'power_on_launch_geo'
  call log_write_check(m_name,s_name,38,status)
  write(kout,*,iostat=status) self%n%nlaupow
  call log_write_check(m_name,s_name,39,status)
  write(kout,*,iostat=status) 'power_on_shadow_geo'
  call log_write_check(m_name,s_name,40,status)
  write(kout,*,iostat=status) self%n%nshapow
  call log_write_check(m_name,s_name,41,status)
  !X      write(kout,*,iostat=status) 'sigma'
  !X      call log_write_check(m_name,s_name,42,status)
  !X      write(kout,*,iostat=status) self%n%sigma
  !X      call log_write_check(m_name,s_name,43,status)

  call edgprof_write(self%edgprof,kout)

end subroutine powcal_write
!---------------------------------------------------------------------
!> transfer power to shadow wall geometry
subroutine powcal_xferpow(self,gshadl)

  !! arguments
  type(powcal_t), intent(inout) :: self   !< powcal data structure
  type(geobjlist_t), intent(inout) :: gshadl   !< shadow data structure

  !! local
  character(*), parameter :: s_name='powcal_xferpow' !< subroutine name
  integer(ki4) :: iobj   !< object number

  if (gshadl%nwset==0) then
     allocate(gshadl%obj(gshadl%ng), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     gshadl%nwset=2
  end if
  ! zero assigned power
  do j=1,gshadl%ng
     gshadl%obj(j)%weight=0._kr4
  end do
  ! assign power to saved object number
  do j=1,self%powres%geobjl%ng
     iobj=self%powres%geobjl%obj(j)%weight+0.5
     if (iobj>0) then
        gshadl%obj(iobj)%weight=self%powres%pow(j)
     end if
  end do

end subroutine powcal_xferpow
!---------------------------------------------------------------------
!> delete powcal_t
subroutine powcal_delete(self)

  !! arguments
  type(powcal_t), intent(inout) :: self !<  powcal data
  !! local

  call position_deleteveclis(self%powres%vecb)
  call position_deleteveclis(self%powres%vecx)
  call geobjlist_delete(self%powres%geobjl)
  deallocate(self%powres%pows)
  deallocate(self%powres%pow)
  calcn_type: select case (self%n%caltype)
  case('msus')
     deallocate(self%powres%psista)
  case('afws')
     call beq_deletepart(self%powres%beq)
  case default
     call beq_deleteplus(self%powres%beq)
  end select calcn_type
  call edgprof_delete(self%edgprof)

end subroutine powcal_delete

end module powcal_m
