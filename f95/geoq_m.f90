module geoq_m

  use const_kind_m
  use log_m
  use const_numphys_h
  use position_h
  use geobjlist_h
  use position_m
  use control_h
  use ls_m
  use btree_m
  use geobj_m
  use query_m
  use beq_h
  use posang_h
  use posang_m
  use geobjlist_m
  use beq_m

  implicit none
  private

! public subroutines
  public :: &
  geoq_delete,   & !< delete  geoq data structure
  geoq_beqscale,   & !< scale geobjlist according to beq data
  geoq_read,   & !< read (vtk)  geoq data structure
  geoq_init,   & !< initialise geometry+field quantities
  geoq_psilimiter,   & !< calculate limits of limiter object
  geoq_psisilh,   & !< calculate \f$ \psi \f$ of silhouette object
  geoq_dsilhcont,   & !< distance between silhouette and flux contour
  geoq_writev, &    !< write (vtk)  geoq data structure
  geoq_writeg    !< write (gnuplot)  geoq data structure

! public types
!> data structure describing geometrical objects and equilibrium field
  type, public :: geoq_t
     type(geobjlist_t) :: objl !< geobj coord data and useful bits
     type(beq_t) :: beq !< equilibrium field
  end type geoq_t

!public variables

! private types

! private variables
  character(*), parameter :: m_name='geoq_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer   :: status   !< error status
  integer(ki4) :: nin   !< input channel for geobj data
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  type(posveclis_t) :: rposl   !< list of position data
  integer(ki4) :: inpos !< number of position vectors to be read
  integer(ki4) :: inobj !< number of geobj records to be read
  integer(ki4) :: inls !< number of entries in list
  integer(ki4) :: idum !< dummy integer
  real(kr8), parameter :: deltal=1 !< silhoutte sampling length (mm)
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> delete geoq_t
subroutine geoq_delete(self)

  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  call geobjlist_delete(self%objl)
  call beq_delete(self%beq)

end subroutine geoq_delete
!---------------------------------------------------------------------
!> scale geobjlist according to beq data
subroutine geoq_beqscale(self,kt)

  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  integer(ki4), intent(in) :: kt    !<  type of transform

  !! local
  character(*), parameter :: s_name='geoq_beqscale' !< subroutine name
  real(kr8), dimension(3)   :: zpar !< scaling parameters
  type(geobjlist_t) :: zobjl !< geobj coord data and useful bits
  type(posang_t) :: zposang !< position and vector involving angles

  zpar=0
  if (self%beq%n%xiopt==2) then
     ! copy to workspace
     zobjl%np=self%objl%np
     !! allocate storage
     if(zobjl%np>0) then
        allocate(zobjl%posl%pos(zobjl%np),stat=status)
        call log_alloc_check(m_name,s_name,1,status)
     else
        call log_error(m_name,s_name,2,error_fatal,'No data')
     end if
     zobjl%posl%np=self%objl%posl%np
     !! copy positions
     do j=1,zobjl%np
        zobjl%posl%pos(j)%posvec=self%objl%posl%pos(j)%posvec
     end do
     print '("number of geobj coordinates copied = ",i10)',zobjl%np
     call log_value("number of geobj coordinates copied ",zobjl%np)

     do j=1,zobjl%np
        ! transform positions to R-Z-zeta space
        zposang%pos=self%objl%posl%pos(j)%posvec
        zposang%opt=0 ; zposang%units=-3
        call posang_invtfm(zposang,0)
        self%objl%posl%pos(j)%posvec=zposang%pos
     end do

     ! scale in zeta
     zpar(2)=0
     zpar(3)=real(self%beq%nzets)
     call geobjlist_spectfm(zobjl,kt,zpar,3)

     ! put back into geobjlist
     do j=1,zobjl%np
        self%objl%posl%pos(j)%posvec=zobjl%posl%pos(j)%posvec
     end do

     deallocate(zobjl%posl%pos)
  end if

end subroutine geoq_beqscale
!---------------------------------------------------------------------
!> read geobj coordinates
subroutine geoq_read(self,infileo,infileb)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  character(*),intent(in) :: infileo !< name of geobjlist input file
  character(*),intent(in) :: infileb !< name of beq input file

  !! local
  character(*), parameter :: s_name='geoq_read' !< subroutine name
  integer(ki4) :: ifldspec !< field specification
  integer(ki4):: ipsibig !< flag whether psi overlarge by 2pi

  call geobjlist_read(self%objl,infileo)
  ifldspec=self%beq%n%fldspec
  ipsibig=abs(self%beq%n%psibig)
  call log_value("psi 2pi too big if unity ",ipsibig)
  call beq_readequil(self%beq,infileb,ifldspec,ipsibig)
  call log_error(m_name,s_name,70,log_info,'geoq input data files read')

end subroutine geoq_read
!---------------------------------------------------------------------
!> initialise geometry+field quantities
!> read geobj coordinates
subroutine geoq_init(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_init' !< subroutine name

  if(self%beq%n%bdryopt/=3.AND.self%beq%n%bdryopt/=6) then
     ! calculate limits based on geometry, unless are silhouette options (3,6)
     call geoq_psilimiter(self)
  end if
  ! psi boundary definition
  if(self%beq%n%bdryopt==1.OR.self%beq%n%bdryopt==5.OR.self%beq%n%bdryopt==9) then
     ! user defined
     self%beq%psibdry=self%beq%n%psiref
     self%beq%psiltr=self%beq%n%psiref
  else if(self%beq%n%bdryopt==2) then
     ! boundary based on nodes of geometry only
     if (beq_rsig()>0) then
        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psiltr)
     else
        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psiltr)
     end if
  else if(self%beq%n%bdryopt==3) then
     ! boundary based on geometry sampled at rate given by parameter deltal
     call geoq_psisilh(self)
!B     if (beq_rsig()>0) then
!B        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psiltr)
!B     else
!B        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psiltr)
!B     end if
!B     ! geoq_override for specific project
!B     if (BEQ_OVERRIDE_ITER) then
!B        call log_error(m_name,s_name,49,log_info,'override for ITER')
        self%beq%psibdry=self%beq%psiltr
!B     end if
  else if(self%beq%n%bdryopt==6) then
     ! boundary based on geometry sampled at rate given by parameter deltal
     call geoq_psisilh(self)
     if (beq_rsig()>0) then
        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psiltr)
     else
        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psiltr)
     end if
  else if(self%beq%n%bdryopt==4.OR.self%beq%n%bdryopt==8) then
       ! boundary from X-point
       call beq_psix(self%beq)
!B     if (beq_rsig()>0) then
!B        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psixpt)
!B     else
!B        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psixpt)
!B     end if
!B     ! geoq_override for specific project
!B     if (BEQ_OVERRIDE_ITER) then
!B     call log_error(m_name,s_name,50,log_info,'override for ITER')
        self%beq%psibdry=self%beq%psixpt
        self%beq%psiltr=self%beq%psixpt
!B     end if
  else if(self%beq%n%bdryopt==7.OR.self%beq%n%bdryopt==10) then
     ! boundary from X-point
     call beq_psix(self%beq)
     if (beq_rsig()>0) then
        self%beq%psibdry=min(self%beq%psiqbdry,self%beq%psixpt)
     else
        self%beq%psibdry=max(self%beq%psiqbdry,self%beq%psixpt)
     end if
  else if(self%beq%n%bdryopt==11.OR.self%beq%n%bdryopt==12) then
     ! boundary from EQDSK regardless
     self%beq%psibdry=self%beq%psiqbdry
  end if

  if (beq_rsig()*(self%beq%psiqbdry-self%beq%psibdry)<0) then
     call log_error(m_name,s_name,1,error_warning,'boundary psi outside separatrix')
  end if

  call log_error(m_name,s_name,80,log_info,'boundary psi value set')
  call log_value("psi boundary reference",self%beq%psibdry)
  call log_value("psi limiter or X-point",self%beq%psiltr)
  call log_value("psi plasma boundary",self%beq%psiqbdry)

end subroutine geoq_init
!---------------------------------------------------------------------
!> calculate limits of limiter object
subroutine geoq_psilimiter(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_psilimiter' !< subroutine name
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zrmin    !<   \f$ R \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zrmax    !<   \f$ R \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zzmin    !<   \f$ Z \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zzmax    !<   \f$ Z \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zpsi    !<  \f$ \psi(R,Z) \f$
  real(kr8) :: zf    !<  \f$ \f(\psi) \f$
  type(posang_t) :: posang !< position and vector involving anglea
  real(kr8) :: zpsimin    !<  least value of \f$ \psi \f$ on geometry
  real(kr8) :: zpsimax    !<  largest value of \f$ \psi \f$ on geometry
  real(kr8) :: ztheta    !<  \f$ \theta(R,Z) \f$
  real(kr8) :: zthetamin    !<  least value of \f$ \theta \f$ on geometry
  real(kr8) :: zthetamax    !<  largest value of \f$ \theta \f$ on geometry
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zsr    !< \f$  r_i \f$
  real(kr8) :: cylj    !<  current estimated in cylindrical approx.

  zpsimin=1.E+8
  zpsimax=-1.E+8
  zthetamin=1.E+8
  zthetamax=-1.E+8
  do j=1,self%objl%np
     ! transform positions to R-Z-zeta space
     posang%pos=self%objl%posl%pos(j)%posvec
     posang%opt=0 ; posang%units=-3
     call posang_invtfm(posang,0)
     zr=posang%pos(1)
     zz=posang%pos(2)
     ! then to flux coordinates
     !! evaluate psi
     !      call spl2d_eval(self%beq%psi,zr,zz,zpsi)
     !   zpsimin=min(zpsi,zpsimin)
     call posang_psitfm(posang,self%beq)
     zpsi=posang%pos(1)
     ztheta=posang%pos(2)
     if (zpsi<zpsimin) then
        zpsimin=zpsi
        zrmin=zr ; zzmin=zz
     end if
     !   zpsimax=max(zpsi,zpsimax)
     if (zpsi>zpsimax) then
        zpsimax=zpsi
        zrmax=zr ; zzmax=zz
     end if
     zthetamin=min(zthetamin,ztheta)
     zthetamax=max(zthetamax,ztheta)
  end do

  if (beq_rsig()>0) then
     ! equiv. (self%psiaxis<self%psiqbdry)
     ! psi increasing outward
     self%beq%psiltr=zpsimin
     zr=zrmin ; zz=zzmin
     self%beq%psiotr=zpsimax
  else
     self%beq%psiltr=zpsimax
     zr=zrmax ; zz=zzmax
     self%beq%psiotr=zpsimin
  end if

  self%beq%thetagmax=zthetamax
  self%beq%thetagmin=zthetamin

  call spl2d_eval(self%beq%dpsidr,zr,zz,zdpdr)
  call spl2d_eval(self%beq%dpsidz,zr,zz,zdpdz)
  call spl2d_eval(self%beq%psi,zr,zz,zpsi)
  self%beq%rbdry=zr
  self%beq%bpbdry=(1/zr)*sqrt( max(0.,(zdpdr**2+zdpdz**2)) )

  ! evaluate I aka f at psi
  call spleval(self%beq%f,self%beq%mr,self%beq%psiaxis,self%beq%psiqbdry,zpsi,zf,1)
  self%beq%btotbdry=sqrt( max(0.,(self%beq%bpbdry**2+(zf/zr)**2)) )

  call log_error(m_name,s_name,1,log_info,'Reference boundary values')
  call log_value("SMITER-GEOQ psiltr ",self%beq%psiltr)
  call log_value("SMITER-GEOQ rbdry ",self%beq%rbdry)
  call log_value("SMITER-GEOQ zbdry ",zz)
  call log_value("SMITER-GEOQ bpbdry ",self%beq%bpbdry)
  call log_value("SMITER-GEOQ btotbdry ",self%beq%btotbdry)
  zsr=sqrt( max(0.,(zr-self%beq%n%rcen)**2+(zz-self%beq%n%zcen)**2) )
  cylj=abs(zsr*self%beq%bpbdry/2.e-7)
  call log_value("Estimated cylindrical current ",cylj)

  call log_error(m_name,s_name,2,log_info,'Object limits ')
  call log_value("OBJECT LIMIT psimax ",zpsimax)
  call log_value("OBJECT LIMIT psimin ",zpsimin)
  call log_value("OBJECT LIMIT thetamax ",zthetamax)
  call log_value("OBJECT LIMIT thetamin ",zthetamin)

end subroutine geoq_psilimiter
!---------------------------------------------------------------------
!> calculate psi of silhouette object
subroutine geoq_psisilh(self)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data

  !! local
  character(*), parameter :: s_name='geoq_psisilh' !< subroutine name
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zrmin    !<   \f$ R \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zrmax    !<   \f$ R \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  real(kr8) :: zzmin    !<   \f$ Z \f$ at minimum \f$ \psi \f$ on geometry
  real(kr8) :: zzmax    !<   \f$ Z \f$ at maximum \f$ \psi \f$ on geometry
  real(kr8) :: zpsi    !<  \f$ \psi(R,Z) \f$
  real(kr8) :: zf    !<  \f$ \f(\psi) \f$
  type(posang_t) :: posang !< position and vector involving anglea
  real(kr8) :: zpsimin    !<  least value of \f$ \psi \f$ on geometry
  real(kr8) :: zpsimax    !<  largest value of \f$ \psi \f$ on geometry
  real(kr8), dimension(3)  :: zv1 !<   3-D Cartesian vector
  real(kr8), dimension(3)  :: zvd !<   3-D Cartesian vector difference
  real(kr8) :: zline    !<  length of line
  integer(ki4) :: inline    !<  length of line in segments
  integer(ki4) :: iobj   !< object pointer
  integer(ki4) :: inode   !< point (node) pointer
  integer(ki4) :: ityp   !< object type
  integer(ki4) :: inumpts   !< number of points in object
  real(kr8) :: zdpdr    !<  \f$ \frac{\partial\psi}{\partial R} \f$
  real(kr8) :: zdpdz    !<  \f$ \frac{\partial\psi}{\partial Z} \f$
  real(kr8) :: zsr    !< \f$  r_i \f$
  real(kr8) :: cylj    !<  current estimated in cylindrical approx.

  zpsimin=1.E+8
  zpsimax=-1.E+8
  loop_object: do j=1,self%objl%ng
     iobj=self%objl%obj2(j)%ptr
     ityp=self%objl%obj2(j)%typ
     inumpts=geobj_entry_table(ityp)
     do k=1,inumpts-1
        ! length of line
        inode=self%objl%nodl(iobj)
        zv1=self%objl%posl%pos(inode)%posvec
        inode=self%objl%nodl(iobj+1)
        zvd=self%objl%posl%pos(inode)%posvec-zv1
        zline=sqrt( max(0.,zvd(1)**2+zvd(2)**2+zvd(3)**2) )
        inline=zline/deltal
        !dbg   if (j<=10) write(*,*) inode,zvd,inline   !dbg
        loop_segment: do l=1,max(inline,1)
           ! transform positions to R-Z-zeta space
           posang%pos=zv1+(l-1)*(zvd/inline)
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           zr=posang%pos(1)
           zz=posang%pos(2)
           ! evaluate psi
           call spl2d_eval(self%beq%psi,zr,zz,zpsi)
           !dbg   if (j<=10) write(*,*) zr,zz,zpsi  !dbg
           !   zpsimin=min(zpsi,zpsimin)
           if (zpsi<zpsimin) then
              zpsimin=zpsi
              zrmin=zr ; zzmin=zz
           end if
           !   zpsimax=max(zpsi,zpsimax)
           if (zpsi>zpsimax) then
              zpsimax=zpsi
              zrmax=zr ; zzmax=zz
           end if
        end do loop_segment
     end do
  end do loop_object

  !dbg   write(*,*) zrmin,zrmax,zpsimin,zpsimax !dbg
  if (beq_rsig()>0) then
     ! equiv. (self%psiaxis<self%psiqbdry)
     ! psi increasing outward
     self%beq%psiltr=zpsimin
     zr=zrmin ; zz=zzmin
     self%beq%psiotr=zpsimax
  else
     self%beq%psiltr=zpsimax
     zr=zrmax ; zz=zzmax
     self%beq%psiotr=zpsimin
  end if
  write(*,*) 'psiltr,psiotr', self%beq%psiltr,self%beq%psiotr !dbg
  call spl2d_eval(self%beq%dpsidr,zr,zz,zdpdr)
  call spl2d_eval(self%beq%dpsidz,zr,zz,zdpdz)
  self%beq%rbdry=zr
  self%beq%bpbdry=(1/zr)*sqrt( max(0.,(zdpdr**2+zdpdz**2)) )

  ! evaluate I aka f at psi
  zpsi=self%beq%psiotr
  call spleval(self%beq%f,self%beq%mr,self%beq%psiaxis,self%beq%psiqbdry,zpsi,zf,1)
  self%beq%btotbdry=sqrt( max(0.,(self%beq%bpbdry**2+(zf/zr)**2)) )
  !dbg   write(*,*) self%beq%psiltr,self%beq%psiotr !dbg
  call log_error(m_name,s_name,1,log_info,'Reference boundary values')
  call log_value("SMITER-GEOQ psiltr ",self%beq%psiltr)
  call log_value("SMITER-GEOQ rbdry ",self%beq%rbdry)
  call log_value("SMITER-GEOQ zbdry ",zz)
  call log_value("SMITER-GEOQ bpbdry ",self%beq%bpbdry)
  call log_value("SMITER-GEOQ btotbdry ",self%beq%btotbdry)
  zsr=sqrt( max(0.,(zr-self%beq%n%rcen)**2+(zz-self%beq%n%zcen)**2) )
  cylj=abs(zsr*self%beq%bpbdry/2.e-7)
  call log_value("Estimated cylindrical current ",cylj)

end subroutine geoq_psisilh
!---------------------------------------------------------------------
!> calculate distance between silhouette and flux contour as function of silhouette position
subroutine geoq_dsilhcont(self,prc,pzc,prs,pzs,knear,pdist)
  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  real(kr8), dimension(:), intent(in) :: prc !< \f$ R \f$ values of contour
  real(kr8), dimension(:), intent(in) :: pzc !< \f$ Z \f$ values of contour
  real(kr8), dimension(:), allocatable, intent(out) :: prs !< \f$ R \f$ values of silhouette
  real(kr8), dimension(:), allocatable, intent(out) :: pzs !< \f$ Z \f$ values of silhouette
  integer(ki4), dimension(:), allocatable, intent(out) :: knear !< nearest segment of silhouette
  real(kr8), dimension(:), allocatable, intent(out) :: pdist !< distances of contour from silhouette

  !! local
  character(*), parameter :: s_name='geoq_dsilhcont' !< subroutine name
  real(kr8) :: zr    !<   \f$ R(x,y,z) \f$
  real(kr8) :: zz    !<   \f$ Z(x,y,z) \f$
  type(posang_t) :: posang !< position and vector involving anglea
  real(kr8) :: zdist    !<  initial least value of distance
  real(kr8) :: zdistsq    !<  initial least value of distance squared
  real(kr8) :: zdsq    !<  distance squared
  real(kr8), dimension(3)  :: zv1 !<   3-D Cartesian vector
  real(kr8), dimension(3)  :: zvd !<   3-D Cartesian vector difference
  real(kr8) :: zline    !<  length of line
  integer(ki4) :: inc    !<  length of contour in segments
  integer(ki4) :: inseg    !<  length of silhouette in segments
  integer(ki4) :: inline    !<  length of line in segments
  integer(ki4) :: il    !<  loope counter for line in segments
  integer(ki4) :: iobj   !< object pointer
  integer(ki4) :: inode   !< point (node) pointer
  integer(ki4) :: ityp   !< object type
  integer(ki4) :: inumpts   !< number of points in object

  inc=size(prc)
  ! work out size of distance array
  inseg=0
  loop_size: do j=1,self%objl%ng
     iobj=self%objl%obj2(j)%ptr
     ityp=self%objl%obj2(j)%typ
     inumpts=geobj_entry_table(ityp)
     do k=1,inumpts-1
        ! length of line
        inode=self%objl%nodl(iobj)
        zv1=self%objl%posl%pos(inode)%posvec
        inode=self%objl%nodl(iobj+1)
        zvd=self%objl%posl%pos(inode)%posvec-zv1
        zline=sqrt( max(0.,zvd(1)**2+zvd(2)**2+zvd(3)**2) )
        inline=zline/deltal
        inseg=inseg+inline
     end do
  end do loop_size

  !! allocate distance/silhouette position storage
  if(inseg>0) then
     allocate(prs(inseg),pzs(inseg),knear(inseg),pdist(inseg),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No silhouette position data')
  end if

  ! calculate positions for distance array
  il=0
  loop_object: do j=1,self%objl%ng
     iobj=self%objl%obj2(j)%ptr
     ityp=self%objl%obj2(j)%typ
     inumpts=geobj_entry_table(ityp)
     do k=1,inumpts-1
        ! length of line
        inode=self%objl%nodl(iobj)
        zv1=self%objl%posl%pos(inode)%posvec
        inode=self%objl%nodl(iobj+1)
        zvd=self%objl%posl%pos(inode)%posvec-zv1
        zline=sqrt( max(0.,zvd(1)**2+zvd(2)**2+zvd(3)**2) )
        inline=zline/deltal
        loop_segment: do l=1,inline
           ! transform positions to R-Z-zeta space
           posang%pos=zv1+(l-1)*(zvd/inline)
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           il=il+1
           ! save
           prs(il)=posang%pos(1)
           pzs(il)=posang%pos(2)
        end do loop_segment
     end do
  end do loop_object

  ! now calculate distance
  loop_segments: do j=1,inseg
     zr=prs(j)
     zz=pzs(j)
     ! evaluate distance to contour
     zdist=1.E+8
     zdistsq=zdist**2
     do k=1,inc
        ! squared distance between points
        zdsq=(zr-prc(k))**2+(zz-pzc(k))**2
        if (zdsq<zdistsq) then
           zdistsq=zdsq
           knear(j)=k
        end if
     end do
     pdist(j)=sqrt( max(0._kr8,zdistsq) )
  end do loop_segments

  call log_error(m_name,s_name,1,log_info,'Nearest contour segments calculated')

end subroutine geoq_dsilhcont
!---------------------------------------------------------------------
!> write (vtk)  geoq data structure
subroutine geoq_writev(self,kchar,kplot)

  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  character(*), intent(in) :: kchar  !< case
  integer(ki4) :: kplot   !< output channel for vis. data


  !! local
  character(*), parameter :: s_name='geoq_writev' !< subroutine name
  type(posang_t) :: posang !< position and vector involving angles
  integer(ki4) :: isum !< sum workspace
  integer(ki2) :: inn !< number of nodes defining geobj
  integer(ki4) :: iplall !< flag plot 'all' option
  integer(ki4) :: irzzeta !< flag plot 'all' option
  real(kr4), dimension(3) :: znormal !< unit normal vector
  real(kr4), dimension(3) :: zbary !< vector of barycentre of geobj
  real(kr4) :: zmag !< magnitude of normal vector
  type(geobj_t) :: iobj !< object
  logical, parameter :: output_cartv=.TRUE. !< local variable
  type(posvecl_t) :: zpos !< local variable
  type(posvecl_t) :: zpostfm !< local variable
  type(tfmdata_t) :: ztfmdata !< local variable

  ! preamble
  irzzeta=0
  iplall=0
  plot_init_type: select case (kchar)
  case('geometry')
  case('all')
     iplall=1
  case('frzzeta')
     irzzeta=1
     iplall=1
  case('frzxi')
     irzzeta=2
     iplall=1
  case('allcart')
  end select plot_init_type

  plot_type: select case (kchar)
  case('geometry','all','frzzeta','frzxi')
     write(kplot,'(''DATASET UNSTRUCTURED_GRID'')')
     write(kplot,'(''POINTS '',I8, '' float'')') self%objl%np
     ! transform positions
     if (irzzeta==0) then
        do j=1,self%objl%np
           ! transform positions to psi-theta-zeta space and write
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           call posang_psitfm(posang,self%beq)
           call posang_writev(posang,kplot)
        end do
     else if (irzzeta==1) then
        do j=1,self%objl%np
           ! transform positions to R-Z-zeta space and write
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           call posang_writev(posang,kplot)
        end do
     else if (irzzeta==2) then
        ! default identity
        ztfmdata%ntfm=2
        ztfmdata%matrix(1,:)=(/1.,  0.,  0. /)
        ztfmdata%matrix(2,:)=(/0.,  1.,  0. /)
        ztfmdata%matrix(3,:)=(/0.,  0.,  1. /)
        ztfmdata%scale=(/1.,  1.,  1. /)
        ztfmdata%offset=(/0.,  0.,  0. /)
        ztfmdata%scale(3)=real(self%beq%nzets)

        do j=1,self%objl%np
           ! transform positions to R-Z-xi space and write
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           zpos%posvec=posang%pos
           zpostfm=position_tfm(zpos,ztfmdata)
           posang%pos=zpostfm%posvec
           call posang_writev(posang,kplot)
        end do
     end if
     write(kplot, '('' '')')

     ! count entries in CELLS
     isum=0
     do j=1,self%objl%ng
        inn=geobj_entry_table(self%objl%obj2(j)%typ)
        isum=isum+inn
     end do

     ! output CELL data
     write(kplot,'(''CELLS '',I8,I8)') self%objl%ng,self%objl%ng+isum
     i=1
     do j=1,self%objl%ng
        inn=geobj_entry_table(self%objl%obj2(j)%typ)
        write(kplot,'(8(1X,I8))') inn,(self%objl%nodl(i+ij-1)-1,ij=1,inn)
        i=i+inn
     end do
     write(kplot, '('' '')')

     ! output CELL types
     write(kplot,'(''CELL_TYPES '',I8)') self%objl%ng
     do j=1,self%objl%ng
        write(kplot,'(1X,I8)') self%objl%obj2(j)%typ
     end do
     write(kplot, '('' '')')

     write(kplot,'(''POINT_DATA '',I8)') self%objl%np

     ! output cartesian vectors (mm) usually

     if (output_cartv) then
        write(kplot,'(''VECTORS Xcart float'')')
        do j=1,self%objl%np
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0
           call posang_writev(posang,kplot)
        end do
     end if

     if(iplall==1)then

        write(kplot, '('' '')')
        ! output cartesian B vectors
        write(kplot,'(''VECTORS Bcart float'')')
        do j=1,self%objl%np
           ! position
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           !BB! first position to cyl polars to get psi derivs
           !BB         call posang_invtfm(posang,0)
           !BB         zr=posang%pos(1)
           !BB         zz=posang%pos(2)
           !BB         call spl2d_eval(self%beq%dpsidr,zr,zz,zdpdr)
           !BB         call spl2d_eval(self%beq%dpsidz,zr,zz,zdpdz)
           !BB         call spl2d_eval(self%beq%psi,zr,zz,zpsi)
           !BB! evaluate I aka f at psi
           !BB         call spleval(self%beq%f,self%beq%mr,self%beq%psiaxis,self%beq%psiqbdry,zpsi,zf,1)
           !BB! B in toroidal-cyl polars
           !BB         posang%vec(1)=-zdpdz/zr ; posang%vec(2)=zdpdr/zr ; posang%vec(3)=zf/zr
           !BB         posang%opt=17
           !BB! vec to cartesians and output
           !BB         call posang_tfm(posang,-3)
           call beq_b(self%beq,posang,0)
           call posang_writev(posang,kplot,2)
        end do

     end if

  case('allcart')
     write(kplot,'(''DATASET UNSTRUCTURED_GRID'')')
     write(kplot,'(''POINTS '',I8, '' float'')') self%objl%np
     do j=1,self%objl%np
        ! write positions in Cartesians (mm)
        posang%pos=self%objl%posl%pos(j)%posvec
        posang%opt=0 ; posang%units=-3
        call posang_writev(posang,kplot)
     end do
     write(kplot, '('' '')')

     ! count entries in CELLS
     isum=0
     do j=1,self%objl%ng
        inn=geobj_entry_table(self%objl%obj2(j)%typ)
        isum=isum+inn
     end do

     ! output CELL data
     write(kplot,'(''CELLS '',I8,I8)') self%objl%ng,self%objl%ng+isum
     i=1
     do j=1,self%objl%ng
        inn=geobj_entry_table(self%objl%obj2(j)%typ)
        write(kplot,'(8(1X,I8))') inn,(self%objl%nodl(i+ij-1)-1,ij=1,inn)
        i=i+inn
     end do
     write(kplot, '('' '')')

     ! output CELL types
     write(kplot,'(''CELL_TYPES '',I8)') self%objl%ng
     do j=1,self%objl%ng
        write(kplot,'(1X,I8)') self%objl%obj2(j)%typ
     end do
     write(kplot, '('' '')')

     ! output psi scalar
     write(kplot,'(''POINT_DATA '',I8)') self%objl%np
     write(kplot,'(''SCALARS Psi float'')')
     write(kplot,'(''LOOKUP_TABLE default'')')
     do j=1,self%objl%np
        ! transform positions to psi-theta-zeta space and write psi
        posang%pos=self%objl%posl%pos(j)%posvec
        posang%opt=0 ; posang%units=-3
        call posang_invtfm(posang,0)
        call posang_psitfm(posang,self%beq)
        write(kplot, cfmtbs ) posang%pos(1)
     end do

     write(kplot, '('' '')')
     ! output cartesian B vectors
     write(kplot,'(''VECTORS Bcart float'')')
     do j=1,self%objl%np
        ! field at nodes
        posang%pos=self%objl%posl%pos(j)%posvec
        posang%opt=0 ; posang%units=-3
        call beq_b(self%beq,posang,0)
        call posang_writev(posang,kplot,2)
     end do

     write(kplot, '('' '')')
     write(kplot,'(''CELL_DATA '',I8)') self%objl%ng
     ! output normal vectors
     write(kplot,'(''VECTORS Normal float'')')
     do j=1,self%objl%ng
        iobj%geobj=self%objl%obj2(j)%ptr
        iobj%objtyp=self%objl%obj2(j)%typ
        call geobj_normal(iobj,self%objl%posl,self%objl%nodl,znormal,zmag)
        write(kplot, cfmtb1v ) znormal
     end do
     write(kplot, '('' '')')
     ! output cartesian B vectors as cell data
     write(kplot,'(''VECTORS Bcell float'')')
     do j=1,self%objl%ng
        ! calculate barycentre
        iobj%geobj=self%objl%obj2(j)%ptr
        iobj%objtyp=self%objl%obj2(j)%typ
        call geobj_centre(iobj,self%objl%posl,self%objl%nodl,zbary)
        ! evaluate field at barycentre
        posang%pos=zbary
        posang%opt=0 ; posang%units=-3
        call beq_b(self%beq,posang,0)
        call posang_writev(posang,kplot,2)
     end do

  end select plot_type

end subroutine geoq_writev

!---------------------------------------------------------------------
!> write (gnu)  geoq data structure
subroutine geoq_writeg(self,kchar,kout)

  !! arguments
  type(geoq_t), intent(inout) :: self !< geometrical objects and equilibrium data
  character(*), intent(in) :: kchar  !< case
  integer(ki4) :: kout   !< output channel for gnuplot data

  !! local
  character(*), parameter :: s_name='geoq_writeg' !< subroutine name
  type(posang_t) :: posang !< position and vector involving angles

  plot_type: select case (kchar)
  case('gnusil')
        ! positions in R-Z-zeta space
        do j=1,self%objl%np
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           write(kout,'(1x,i9,'//cfmt2v,iostat=status) &
&          j,posang%pos(1),posang%pos(2),posang%pos(3)
           call log_write_check(m_name,s_name,1,status)
        end do
  case('gnusilm')
        ! positions in psi-theta-zeta space
        do j=1,self%objl%np
           posang%pos=self%objl%posl%pos(j)%posvec
           posang%opt=0 ; posang%units=-3
           call posang_invtfm(posang,0)
           call posang_psitfm(posang,self%beq)
           write(kout,'(1x,i9,'//cfmt2v,iostat=status) &
&          j,posang%pos(1),posang%pos(2),posang%pos(3)
           call log_write_check(m_name,s_name,2,status)
        end do
  end select plot_type

end subroutine geoq_writeg

end module geoq_m
