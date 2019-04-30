module mtest_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use control_h
  use position_h
  use position_m
  use vfile_m
  use geobjlist_h

  implicit none
  private

! public subroutines
  public ::  &
  mtest_init, & !< create  mtest data structure
  mtest_readv, & ! read mtest position data
  mtest_set, & ! allocate rest of mtest data
  mtest_writev, & ! write in visualisation format
  mtest_writeptsv, & ! write intersection points in visualisation format
  mtest_delete ! deallocate mtest array

  type, public :: mtest_t
     type(mtnumerics_t) :: m !< data for mtest
     type(posveclis_t) :: posl !< list of positions
     type(posveclis_t) :: posln !< list of new positions
     integer(ki4), dimension(:), allocatable :: nodl !< list of nodes
     integer(ki4), dimension(:), allocatable :: nodln !< list of new nodes
     integer(ki4), dimension(:), allocatable :: objl !< list of objects hit
     integer(ki4) :: np !< number of entries in posl
     integer(ki4) :: npn !< number of entries in posln
     integer(ki4) :: ntrk   !< numbers of tracks
     integer(ki4) :: nptrk   !< number of positions + tracks
  end type mtest_t

! private types

! private variables
  character(*), parameter :: m_name='mtest_m' !< module name
  integer   :: status   !< error status
  integer(ki4) :: i!< loop counters
  integer(ki4) :: j!< loop counters
  integer(ki4) :: k!< loop counters
  integer(ki4) :: l !< loop counters
  character(len=80) :: icbuf !< local variable
  logical :: unitused !< flag to test unit is available
  logical :: iltest !< logical flag
  integer, save :: nmt   !< output channel for mtest position data
  integer(ki4) :: iopt   !< optionf for vfile_rvectorread
  integer(ki4) :: ip   !< number of positions
  integer(ki4) :: iph   !< numbers of tracks
  integer(ki4) :: ipph   !< number of positions + tracks
  real(kr4), dimension(:), allocatable :: work !< work array for nominal velocity
  type(posveclis_t) :: pwork !< list of working positions

  contains
!---------------------------------------------------------------------
!> create  mtest data structure
subroutine mtest_init(self,numerics)

  !! arguments
  type(mtest_t), intent(inout) :: self   !< object data structure
  type(mtnumerics_t), intent(in) :: numerics   !< object control data structure

  !! local
  character(*), parameter :: s_name='mtest_init' !< subroutine name

  ! initialise control structure
  self%m=numerics

end subroutine mtest_init
!---------------------------------------------------------------------
!> read in (vtk) mtest position data
subroutine mtest_readv(self,geobjl,infile,kinmt)

  !! arguments
  type(mtest_t), intent(inout) :: self   !< mtest position data
  type(geobjlist_t), intent(inout)  :: geobjl   !< geobjl defining mtest
  character(*),intent(in) :: infile !< name of input file
  integer :: kinmt   !< input channel for mtest position data

  logical isnumb
  external isnumb

  !! local
  character(*), parameter :: s_name='mtest_readv' !< subroutine name
  integer(ki4), dimension(3) :: iadim !<  dimensions of vector list data
  integer(ki4) :: inres !< number of result vectors read
  type(geobj_t) :: iobj !< geo object
  real(kr4), dimension(3) :: zcentr !< centre
  real(kr4), dimension(3) :: znormal !< normal
  real(kr4) :: zmag !< magnitude of normal vector

  file_content: select case (self%m%vcontent)
  case (0,1)
     call position_readlis(self%posl,infile,nmt)
     ip=self%posl%np

  case (2)
     iph=geobjl%np
     ip=2*iph
     allocate(self%posl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,1,status)
     icbuf=' '
     iopt=1
     call vfile_rvectorread(work,iph,iadim,icbuf,self%m%fldnam,kinmt,iopt)
     do j=1,iph
        self%posl%pos(2*j-1)=geobjl%posl%pos(j)
        self%posl%pos(2*j)%posvec(1)=work(3*j-2)
        self%posl%pos(2*j)%posvec(2)=work(3*j-1)
        self%posl%pos(2*j)%posvec(3)=work(3*j)
     end do
  case (3)
     iph=geobjl%ng
     ip=2*iph
     allocate(self%posl%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,5,status)
     do j=1,iph
        iobj%geobj=geobjl%obj2(j)%ptr
        iobj%objtyp=geobjl%obj2(j)%typ
        call geobj_centre(iobj,geobjl%posl,geobjl%nodl,zcentr)
        call geobj_normal(iobj,geobjl%posl,geobjl%nodl,znormal,zmag)
        self%posl%pos(2*j-1)%posvec=zcentr
        self%posl%pos(2*j)%posvec=znormal
     end do
  end select file_content

  self%np=ip
  call log_value("number of data read ",self%np)

end subroutine mtest_readv
!---------------------------------------------------------------------
!> set up in mtest data
subroutine mtest_set(self)

  !! arguments
  type(mtest_t), intent(inout) :: self   !< mtest data

  !! local
  character(*), parameter :: s_name='mtest_set' !< subroutine name

  ip=self%posl%np
  iph=ip/2
  ipph=ip+iph
  !! allocate storage
  if(iph>0) then
     allocate(self%posln%pos(iph),self%nodl(iph),self%nodln(iph),self%objl(iph), &
 &   stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No position data')
  end if
  self%nptrk=ipph
  self%ntrk=iph

  file_content: select case (self%m%vcontent)
  case default

  case (1,3)
     ! alternate position is in fact velocity or normal
     do j=1,ip,2
        self%posl%pos(j+1)%posvec=self%posl%pos(j)%posvec+self%posl%pos(j+1)%posvec*self%m%mdt
     end do

  case (2)
     ! alternate block position is in fact velocity
     !! allocate storage
     allocate(pwork%pos(ip), stat=status)
     call log_alloc_check(m_name,s_name,3,status)
     do j=1,iph
        pwork%pos(2*j-1)=self%posl%pos(j)
        pwork%pos(2*j)%posvec=self%posl%pos(j)%posvec+self%posl%pos(iph+j)%posvec*self%m%mdt
     end do
     self%posl%pos=pwork%pos
     deallocate(pwork%pos)
  end select file_content

end subroutine mtest_set
!!---------------------------------------------------------------------
!> output mtest position vector
subroutine mtest_writev(self,numerics,kchar,kout)

  !! arguments
  type(mtest_t), intent(in) :: self   !< mtest position data
  type(mtnumerics_t), intent(in) :: numerics   !< numeric controls
  character(*),intent(in):: kchar !< control output (dummy)
  integer, intent(in) :: kout   !< output channel for mtest position data


  !! local
  character(*), parameter :: s_name='mtest_writev' !< subroutine name
  type(posvecl_t) :: zpos !< local variable

  ip=self%np
  iph=self%ntrk
  ipph=self%nptrk
  !! plot list of all positions
  write(kout,'(''DATASET UNSTRUCTURED_GRID'')')
  write(kout,'(''POINTS '',I8, '' float'')') ipph
  ! start and end-points
  do j=1,ip
     pos_type: select case (kchar)
     case('real space')
        zpos=position_invqtfm(self%posl%pos(j),numerics%n%geobj_coord_tfm)
     case default
        !     case('quantised space')
        zpos=self%posl%pos(j)
     end select pos_type
     call position_writev(zpos,kout)
  end do
  ! intersection points
  do j=1,iph
     pos_type2: select case (kchar)
     case('real space')
        zpos=position_invqtfm(self%posln%pos(j),numerics%n%geobj_coord_tfm)
     case default
        !     case('quantised space')
        zpos=self%posln%pos(j)
     end select pos_type2
     call position_writev(zpos,kout)
  end do
  write(kout, '('' '')')
  ! cells as polylines
  write(kout,'(''CELLS '',I8,1X,I8)') iph,ipph+iph
  do j=0,iph-1
     write(kout,'(''3'',3(1X,I8))') 2*j,ip+j,(2*j+1)
  end do

  write(kout,'(''CELL_TYPES '',I8)') iph
  do j=1,iph
     write(kout,'(''4'')')
  end do
  write(kout, '('' '')')

  !! intersecting objects
  write(kout,'(''POINT_DATA'',I8)') ipph
  write(kout,'(''SCALARS objects int 1'')')
  write(kout,'(''LOOKUP_TABLE default'')')
  plot_type: select case (kchar)
  case default
     !     write(kout,'((I8,1x))') ( self%nodl(i),  i=1,iph)
     write(kout,'((I8,1x))') ( 0,  i=1,ip)
     write(kout,'((I8,1x))') ( self%objl(i),  i=1,iph)
  end select plot_type


end subroutine mtest_writev
!!---------------------------------------------------------------------
!> output mtest intersection points
subroutine mtest_writeptsv(self,numerics,kchar,kout)

  !! arguments
  type(mtest_t), intent(in) :: self   !< mtest position data
  type(mtnumerics_t), intent(in) :: numerics   !< numeric controls
  character(*),intent(in):: kchar !< control output (dummy)
  integer, intent(in) :: kout   !< output channel for mtest position data

  !! local
  character(*), parameter :: s_name='mtest_writeptsv' !< subroutine name
  type(posvecl_t) :: zpos !< local variable

  ip=self%np
  iph=self%ntrk
  ipph=self%nptrk
  !! plot list of all positions
  write(kout,'(''DATASET UNSTRUCTURED_GRID'')')
  write(kout,'(''POINTS '',I8, '' float'')') iph
  ! intersection points
  do j=1,iph
     pos_type2: select case (kchar)
     case('real space')
        zpos=position_invqtfm(self%posln%pos(j),numerics%n%geobj_coord_tfm)
     case default
        !     case('quantised space')
        zpos=self%posln%pos(j)
     end select pos_type2
     call position_writev(zpos,kout)
  end do
  write(kout, '('' '')')

  ! points as nodal vertices
  write(kout,'(''CELLS 1 '',I8)') iph+1
  write(kout,'(I8)') iph,  (i-1,i=1,iph)
  write(kout,'(''CELL_TYPES 1'')')
  write(kout,'(''2'')')
  write(kout, '('' '')')

  !! intersecting objects
  write(kout,'(''POINT_DATA'',I8)') iph
  write(kout,'(''SCALARS objects int 1'')')
  write(kout,'(''LOOKUP_TABLE default'')')
  plot_type: select case (kchar)
  case default
     write(kout,'((I8,1x))') ( self%objl(i),  i=1,iph)
  end select plot_type


end subroutine mtest_writeptsv
!---------------------------------------------------------------------
!> delete mtest_t
subroutine mtest_delete(self)

  !! arguments
  type(mtest_t), intent(inout) :: self !< mtest data

  deallocate(self%posl%pos)
  deallocate(self%posln%pos)
  deallocate(self%nodl)
  deallocate(self%nodln)
  deallocate(self%objl)

end subroutine mtest_delete


end module mtest_m
