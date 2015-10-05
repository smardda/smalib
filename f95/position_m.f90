module position_m

  use const_kind_m
  use log_m
  use position_h

  implicit none
  private

! public subroutines
  public :: &
 &position_inbox, & !< test whether within real box
 &position_innbox, & !< test whether within normalised box
 &position_asvecb, & !< encode as binary vector
 &position_tfm, & !< transform position vector
 &position_invtfm, & !< inverse transform position vector
 &position_qtfm, & !< quantise position vector
 &position_invqtfm, & !< inverse quantise position vector
 &position_readv, & !< read in visualisation format
 &position_writev, & !< write in visualisation format
 &position_readcon, & !< read position control data
 &position_readlis, & !< read position list data
 &position_readveclis, & !< read position list data as vector
 &position_deleteveclis, & !< delete list of vectors
 &position_writelis, & !< write position list data
 &position_lenlis, & !< physical length of position list data
 &position_tfmlis, & !< transform position list
 &position_invtfmlis, & !< inverse transform position list
 &position_qtfmlis, & !< quantise position list
 &position_invqtfmlis ! inverse quantise position list


!public variables

! private types

! private variables
  character(*), parameter :: m_name='position_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer   :: status   !< error status
  integer(ki4) :: nin   !< input channel for position list data
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  logical :: iltest !< logical flag


  contains
!---------------------------------------------------------------------
!> is position in real box
  pure logical function position_inbox(self,box)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  real(kr4), dimension(3,2), intent(in) :: box !< bounding box corners


!! local
  character(*), parameter :: s_name='position_inbox' !< subroutine name
  logical :: linint   !< position in interval
  logical :: linbox   !< position in interval
  integer(ki4) :: jj   !< loop counter

  linbox=.true.
  do jj=1,3
     linint=( box(jj,1)<=self%posvec(jj) )
     linint=linint.and.( self%posvec(jj)<box(jj,2) )
     linbox=linbox.and.linint
  end do

  position_inbox=linbox

end function position_inbox
!---------------------------------------------------------------------
!> is position in logical box
  pure logical function position_innbox(self,kcorn,kext)

!! arguments
  type(posvecl_t), intent(in) :: self   !< quantised position data
  integer(ki2), dimension(3), intent(in) :: kcorn !< bounding box corner
  integer(ki2), dimension(3), intent(in) :: kext !< bounding box etents


!! local
  character(*), parameter :: s_name='position_innbox' !< subroutine name
  logical :: linbox   !< position in interval
  integer(ki4) :: jj   !< loop counter
  integer(ki2), dimension(3) :: ivec !< integer position data
  integer(ki2) :: ie !< integer element
  integer(ki2) :: iesum !< integer sum

  if (any(ivec<0)) then
     linbox=.false.
  else
     ivec=int(self%posvec)
     iesum=0
     do jj=1,3
        ie=ishft( ieor(kcorn(jj),ivec(jj)), -kext(jj) )
        iesum=iesum+ie
     end do

     linbox=(iesum==0)
  end if
  position_innbox=linbox

end function position_innbox
!---------------------------------------------------------------------
!> convert transformed position to binary representation (3 X integer*2)
  pure type(posvecb_t) function position_asvecb(self)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data


!! local
  character(*), parameter :: s_name='position_asvecb' !< subroutine name
  integer(ki2), dimension(3) :: ivecb   !< binary rep of position data
  ivecb=int(self%posvec)

  position_asvecb%posb=ivecb

end function position_asvecb
!---------------------------------------------------------------------
!> transform in position space
  pure type(posvecl_t) function position_tfm(self,tfmdata)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform


!! local
  character(*), parameter :: s_name='position_tfm' !< subroutine name
  real(kr4), dimension(3) :: zvec   !< transformed position data
  integer(ki4) :: intfm   !< number of transform
  integer(ki4) :: jj   !< loop

  intfm=tfmdata%ntfm
  transform_number :  select case (intfm)
  case(1)
     zvec=self%posvec*tfmdata%scale
  case(2)
     zvec=self%posvec*tfmdata%scale+tfmdata%offset
  case(3)
     do jj=1,3
        zvec(jj)=dot_product(self%posvec,tfmdata%matrix(jj,:))+tfmdata%offset(jj)
     end do
  case(4)
     do jj=1,3
        zvec(jj)=tfmdata%offset(jj)+dot_product(self%posvec-tfmdata%offset,tfmdata%matrix(jj,:))
     end do
  end select transform_number

  position_tfm%posvec=zvec

end function position_tfm
!---------------------------------------------------------------------
!> inverse transform in position space
  pure type(posvecl_t) function position_invtfm(self,tfmdata)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform


!! local
  character(*), parameter :: s_name='position_invtfm' !< subroutine name
  real(kr4), dimension(3) :: zvec   !< transformed position data
  integer(ki4) :: intfm   !< number of transform

  intfm=tfmdata%ntfm
  transform_number :  select case (intfm)
  case(1)
     zvec=self%posvec/tfmdata%scale
  case(2)
     zvec=(self%posvec-tfmdata%offset)/tfmdata%scale
  case(3)
!   call log_error(m_name,s_name,3,error_fatal,'inverse position transform not defined')
  end select transform_number

  position_invtfm%posvec=zvec

end function position_invtfm
!---------------------------------------------------------------------
!> quantise in position space
  pure type(posvecl_t) function position_qtfm(self,qtfmdata)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  type(quantfm_t), intent(in) :: qtfmdata   !< data defining transform


!! local
  character(*), parameter :: s_name='position_qtfm' !< subroutine name
  real(kr4), dimension(3) :: zvec   !< quantised position data
  integer(ki4) :: inqtfm   !< number of transform

  inqtfm=qtfmdata%nqtfm
  transform_number :  select case (inqtfm)
  case(1)
     zvec=self%posvec*qtfmdata%rhmin
  case(2)
     zvec=self%posvec*qtfmdata%rhmin+qtfmdata%offvec
  end select transform_number

  position_qtfm%posvec=zvec

end function position_qtfm
!---------------------------------------------------------------------
!> inverse quantise in position space
  pure type(posvecl_t) function position_invqtfm(self,qtfmdata)

!! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  type(quantfm_t), intent(in) :: qtfmdata   !< data defining transform


!! local
  character(*), parameter :: s_name='position_invqtfm' !< subroutine name
  real(kr4), dimension(3) :: zvec   !< quantised position data
  integer(ki4) :: inqtfm   !< number of transform

  inqtfm=qtfmdata%nqtfm
  transform_number :  select case (inqtfm)
  case(1)
     zvec=self%posvec*qtfmdata%hmin
  case(2)
     zvec=(self%posvec-qtfmdata%offvec)*qtfmdata%hmin
  end select transform_number

  position_invqtfm%posvec=zvec

end function position_invqtfm
!---------------------------------------------------------------------
!> read in position data
subroutine position_readv(self,kin)

  !! arguments
  type(posvecl_t), intent(out) :: self   !< position data
  integer(ki4), intent(in) :: kin   !< input channel for position data


  !! local
  character(*), parameter :: s_name='position_readv' !< subroutine name

  read(kin,*,iostat=status) self%posvec
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading Particle position')
  end if

end subroutine position_readv
!---------------------------------------------------------------------
!> output position vector
subroutine position_writev(self,kplot)

  !! arguments
  type(posvecl_t), intent(in) :: self   !< position data
  integer(ki4), intent(in) :: kplot   !< output channel for position data


  !! local
  character(*), parameter :: s_name='position_writev' !< subroutine name

  write(kplot,cfmtbv1,iostat=status) self%posvec
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error writing Particle position')
  end if


end subroutine position_writev
!---------------------------------------------------------------------
!> read position control data
subroutine position_readcon(ztfmdata,kin)

  !! arguments
  type(tfmdata_t), intent(out) :: ztfmdata   !< position transform numeric controls
  integer(ki4) :: kin  !< local variable

  !! local
  character(*), parameter :: s_name='position_readcon' !< subroutine name
  real(kr4), dimension(3,3):: position_matrix  !< local variable
  real(kr4), dimension(3):: position_scale  !< local variable
  real(kr4), dimension(3):: position_offset  !< local variable
  integer(ki4):: position_transform  !< local variable

  !! position parameters
  namelist /positionparameters/ &
 &position_matrix , &
 &position_scale , &
 &position_offset , &
 &position_transform

  !! set default position parameters
  ! default identity
  position_transform=1
  position_matrix(1,:)=(/1.,  0.,  0. /)
  position_matrix(2,:)=(/0.,  1.,  0. /)
  position_matrix(3,:)=(/0.,  0.,  1. /)
  position_scale=(/1.,  1.,  1. /)
  position_offset=(/0.,  0.,  0. /)

  !!read position parameters
  read(kin,nml=positionparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading position parameters')
     print '("Fatal error reading position parameters")'
  end if


  !! check for valid data
  if(all(position_scale==0)) &
 &call log_error(m_name,s_name,2,error_fatal,'position_scale must be > 0')
  if(position_transform<1) &
 &call log_error(m_name,s_name,3,error_fatal,'position_transform must be > 0')

  if (position_transform==3) then
     ! FZK's transform
     position_matrix(1,:)=(/0.,  0.,  1. /)
     position_matrix(2,:)=(/-1.,  0.,  0. /)
     position_matrix(3,:)=(/0.,  -1.,  0. /)
     position_scale=(/1.,  1.,  1. /)
     ! working in velocity space, position offset can be ignored
     !     position_offset=(/50.,  0.,  -390. /)
     position_offset=(/0.,  0.,  0. /)
  end if
  !! store values
  ztfmdata%matrix=position_matrix
  ztfmdata%scale=position_scale
  ztfmdata%offset=position_offset
  ztfmdata%ntfm=position_transform

end subroutine position_readcon
!---------------------------------------------------------------------
!> read (vtk) list of position coordinates
subroutine position_readlis(self,infile,kin,kopt)
  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  character(*),intent(in) :: infile !< name of input file
  integer(ki4), intent(out) :: kin   !< input channel for position list data
  integer(ki4), intent(in), optional :: kopt   !< options

  !! local
  character(*), parameter :: s_name='position_readlis' !< subroutine name
  logical :: unitused !< flag to test unit is available
  integer(ki4) :: inpos    !< local variable

  logical :: isnumb !< local variable
  external isnumb

  if(present(kopt)) then
     !! assume unit already open and reading infile
     nin=kin
  else

     !! get file unit
     do i=99,1,-1
        inquire(i,opened=unitused)
        if(.not.unitused)then
           kin=i
           exit
        end if
     end do
     nin=kin

     !! open file
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,1,error_fatal,'Error opening position list data file')
     else
        call log_error(m_name,s_name,1,log_info,'Position list data file opened')
     end if

  end if
  !! File unit sorted out
  !! read header information
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,1,error_fatal,'Error reading header data')
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:6)=='POINTS') then
           iltest=isnumb(ibuf2,inpos,7)
           self%np=inpos
           exit
        end if
     end if
  end do

  !! allocate position storage
  if(self%np>0) then
     allocate(self%pos(self%np), &
 &   stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,2,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,3,error_fatal,'No position data')
  end if

  !! read coordinates
  do j=1,self%np
     call position_readv(self%pos(j),nin)
  end do
  print '("number of position coordinates read = ",i10)',self%np
  call log_value("number of position coordinates read ",self%np)

end subroutine position_readlis
!---------------------------------------------------------------------
!> read (vtk) list of vectors
subroutine position_readveclis(self,infile,kcname,kin,kopt)
  !! arguments
  type(posveclis_t), intent(inout) :: self !< vector list data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of field required
  integer(ki4), intent(inout) :: kin   !< input channel for vector list data
  integer(ki4), intent(in), optional :: kopt   !< options

  !! local
  character(*), parameter :: s_name='position_readveclis' !< subroutine name
  logical :: unitused !< flag to test unit is available
  integer(ki4) :: invec   !< number of vectors
  character(len=80) :: vname !< name of vector
  integer(ki4) :: islen   !< length of vector field name
  integer(ki4) :: islen2   !< length of required vector field name

  logical :: isnumb !< local variable
  external isnumb

  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(present(kopt)) then
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening vector list data file')
        else
           call log_error(m_name,s_name,1,log_info,'Vector list data file opened')
           nin=kin
        end if
     else
        nin=kin
     end if
  else

     !! get file unit
     do i=99,1,-1
        inquire(i,opened=unitused)
        if(.not.unitused)then
           kin=i
           exit
        end if
     end do
     nin=kin

     !! open file
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        call log_error(m_name,s_name,2,error_fatal,'Error opening vector list data file')
     else
        call log_error(m_name,s_name,2,log_info,'Vector list data file opened')
     end if

  end if

  !! File unit now sorted, get to where point data begin
  if (self%np==0) then
     !! read local header information to get number of vectors
     do
        read(nin,fmt='(a)',iostat=status) ibuf1
        !!eof
        if(status<0) then
           exit
           !! error
        else if (status>0) then
           call log_error(m_name,s_name,3,error_fatal,'Error reading header data')
        else
           ibuf2=adjustl(ibuf1)
           if(ibuf2(1:10)=='POINT_DATA') then
              iltest=isnumb(ibuf2,invec,11)
              self%np=invec
              exit
           end if
        end if
     end do
  end if

  !! find vector header
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,error_fatal,'Error reading header data')
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:7)=='VECTORS') then
           ibuf1=adjustl(ibuf2(8:))
           islen=max(2,scan(ibuf1,' '))-1
           vname=ibuf1(:islen)
           if (vname(:islen)==kcname(:islen2)) then
              call log_value("Vector field found ",vname)
              exit
           else
              call log_value("Skipped vector field ",vname)
           end if
           !         else if(ibuf2(1:10)=='POINT_DATA') then
           !            iltest=isnumb(ibuf2,invec,11)
           !            self%np=invec
        end if
     end if
  end do

  !! allocate vector storage
  if(self%np>0) then
     allocate(self%pos(self%np), &
 &   stat=status)
     !! check successful allocation
     if(status/=0) then
        call log_error(m_name,s_name,5,error_fatal,'Unable to allocate memory')
     end if
  else
     call log_error(m_name,s_name,6,error_fatal,'No vector data')
  end if

  !! read coordinates
  do j=1,self%np
     call position_readv(self%pos(j),nin)
  end do
  print '("number of vectors read = ",i10)',self%np
  call log_value("number of vectors read ",self%np)

end subroutine position_readveclis
!---------------------------------------------------------------------
!> delete list of vectors
subroutine position_deleteveclis(self)
  !! arguments
  type(posveclis_t), intent(inout) :: self !< vector list data

  !! local
  character(*), parameter :: s_name='position_deleteveclis' !< subroutine name

  deallocate(self%pos)

end subroutine position_deleteveclis
!---------------------------------------------------------------------
!> write (vtk)  list of position coordinates
subroutine position_writelis(self,kchar,kplot)

  !! arguments
  type(posveclis_t), intent(in) :: self !< position list data
  character(*), intent(in) :: kchar  !< case, specifies objects
  integer(ki4), intent(in) :: kplot   !< output channel for position list data
  !     type(beq_t), intent(inout), optional :: pbeq   !> beq data structure

  !! local
  character(*), parameter :: s_name='position_writelis' !< subroutine name
  integer(ki4) :: iseg   !< number of segments in track

  !      if(present(pbeq)) then
  !! convert positions to cartesians
  !      end if

  !! plot list of all positions
  write(kplot,'(''DATASET UNSTRUCTURED_GRID'')')
  write(kplot,'(''POINTS '',I8, '' float'')') self%np
  do j=1,self%np
     call position_writev(self%pos(j),kplot)
  end do

  write(kplot, '('' '')')
  !      if(present(pbeq)) then
  !! convert positions back to psi-theta-zeta
  !      end if

  plot_type: select case (kchar)
  case('track')
     ! output positions as a track
     iseg=self%np-1
     write(kplot,'(''CELLS '',I8,1X,I8)') iseg, 3*iseg
     do i=1,iseg
        write(kplot,'(''2 '',I8,1X,I8)')  i-1,i
     end do
     write(kplot, '('' '')')
     write(kplot,'(''CELL_TYPES '',I8)') iseg
     do i=1,iseg
        write(kplot,'(''3'')')
     end do
     write(kplot, '('' '')')

  case default
     ! output positions as points
     write(kplot,'(''CELLS 1 '',I8)') self%np+1
     write(kplot,'(I8)') self%np,  (i-1,i=1,self%np)
     write(kplot,'(''CELL_TYPES 1'')')
     write(kplot,'(''2'')')
     write(kplot, '('' '')')

  end select plot_type

end subroutine position_writelis
!---------------------------------------------------------------------
!> physical length of list of position coordinates
subroutine position_lenlis(self,length)

  !! arguments
  type(posveclis_t), intent(in) :: self !< position list data
  real(kr4), intent(out) :: length   !< physical length

  !! local
  character(*), parameter :: s_name='position_lenlis' !< subroutine name
  real(kr8) :: zlensum   !< sum of physical length (working)
  real(kr8) :: zlen   !< physical length (working)
  real(kr8), dimension(3)  :: zlenvec !< physical length vector (working)


  !! sum
  zlensum=0
  do j=1,self%np-1
     zlenvec=self%pos(j+1)%posvec-self%pos(j)%posvec
     zlen=sqrt(max(0._kr8,zlenvec(1)**2+zlenvec(2)**2+zlenvec(3)**2))
     zlensum=zlensum+zlen
  end do

  length=zlensum

end subroutine position_lenlis
!---------------------------------------------------------------------
!> transform list of position coordinates
subroutine position_tfmlis(self,tfmdata)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='position_tfmlis' !< subroutine name
  type(posvecl_t) :: zpos !< local variable

  !! transform list of all positions
  do j=1,self%np
     zpos=position_tfm(self%pos(j),tfmdata)
     self%pos(j)=zpos
  end do

end subroutine position_tfmlis
!---------------------------------------------------------------------
!> inverse transform list of position coordinates
subroutine position_invtfmlis(self,tfmdata)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  type(tfmdata_t), intent(in) :: tfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='position_invtfmlis' !< subroutine name
  type(posvecl_t) :: zpos !< local variable


  !! transform list of all positions
  do j=1,self%np
     zpos=position_invtfm(self%pos(j),tfmdata)
     self%pos(j)=zpos
  end do

end subroutine position_invtfmlis
!---------------------------------------------------------------------
!> transform quantised list of position coordinates
subroutine position_qtfmlis(self,qtfmdata)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  type(quantfm_t), intent(in) :: qtfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='position_qtfmlis' !< subroutine name
  !      type(posvecl_t) position_qtfm
  type(posvecl_t) :: zpos !< local variable


  !! transform list of all positions
  do j=1,self%np
     zpos=position_qtfm(self%pos(j),qtfmdata)
     self%pos(j)=zpos
  end do

end subroutine position_qtfmlis
!---------------------------------------------------------------------
!> inverse transform quantised list of position coordinates
subroutine position_invqtfmlis(self,qtfmdata)

  !! arguments
  type(posveclis_t), intent(inout) :: self !< position list data
  type(quantfm_t), intent(in) :: qtfmdata   !< data defining transform

  !! local
  character(*), parameter :: s_name='position_invqtfmlis' !< subroutine name
  !      type(posvecl_t) position_invqtfm
  type(posvecl_t) :: zpos !< local variable


  !! transform list of all positions
  do j=1,self%np
     zpos=position_invqtfm(self%pos(j),qtfmdata)
     self%pos(j)=zpos
  end do

end subroutine position_invqtfmlis

end module position_m

