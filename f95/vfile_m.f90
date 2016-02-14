module vfile_m

  use const_kind_m
  use log_m

  implicit none
  private

! public subroutines
  public :: &
 &vfile_init, & !< initialise vtk file
 &vfile_rscalarread, & !< read vtk real scalars
 &vfile_iscalarread, & !< read vtk integer scalars
 &vfile_rscalarwrite, & !< write vtk real scalars
 &vfile_iscalarwrite, & !< write vtk integer scalars
 &vfile_close !< close vtk file

! public types

  integer(ki4), parameter, public :: vfile_made_up_data = 1 !< local variable

! private variables
  character(*), parameter :: m_name='vfile_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: nin  !< file unit for input
  integer(ki4) :: nplot  !< file unit for output
  integer(ki4) :: status  !< status flag
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> initialise vtk file
subroutine vfile_init(fplot,descriptor,kplot)

  !! arguments
  character(len=*), intent(in) :: fplot !< file name root
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer(ki4), intent(inout) :: kplot   !< unit number

  !! local
  character(*), parameter :: s_name='vfile_init' !< subroutine name
  logical :: unitused !< flag to test unit is available

  !! open file

  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        kplot=i
        exit
     end if
  end do

  open(unit=kplot,file=trim(fplot)//'.vtk')

  !! write vtk header
  write(kplot,'(''# vtk DataFile Version 2.0'')')
  write(kplot,'(A)') descriptor
  write(kplot,'(''ASCII'')')
  write(kplot,'(''           '')')

  nplot=kplot

end subroutine vfile_init
!---------------------------------------------------------------------
!> read vtk real scalars
subroutine vfile_rscalarread(self,kp,infile,kcname,kin,kopt)
  !! arguments
  real(kr4), dimension(:), allocatable, intent(inout) :: self !< real scalar list data
  integer(ki4), intent(inout) :: kp   !< size of scalar list data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of field required
  integer(ki4), intent(inout) :: kin   !< input channel for scalar list data
  !> if positive on input, assume unit open and do not terminate execution if
  !! data is missing, return positive code if trouble, else zero
  !! if zero on input, missing data is a fatal error
  integer(ki4), intent(inout) :: kopt   !< .

  !! local
  character(*), parameter :: s_name='vfile_rscalarread' !< subroutine name
  logical :: unitused !< flag to test unit is available
  integer :: ierror   !<  whether error is to be fatal or not
  integer(ki4) :: insca   !< number of scalars
  character(80) :: sname !< name of scalar
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name

  logical :: isnumb !< local variable
  external isnumb

  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(kopt>0) then
     ierror=error_warning
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening scalar list data file')
        else
           call log_error(m_name,s_name,1,log_info,'Scalar list data file opened')
        end if
     end if
     nin=kin
  else
     ierror=error_fatal
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
        call log_error(m_name,s_name,2,error_fatal,'Error opening scalar list data file')
     else
        call log_error(m_name,s_name,2,log_info,'Scalar list data file opened')
     end if

  end if

  !! File unit now sorted, get to where point data begin
  !! read local header information
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
           iltest=isnumb(ibuf2,insca,11)
           kp=insca
           exit
        else if(ibuf2(1:8)=='POLYGONS') then
           iltest=isnumb(ibuf2,insca,9)
           kp=insca
           exit
        else if(ibuf2(1:9)=='CELL_DATA') then
           iltest=isnumb(ibuf2,insca,10)
           kp=insca
           exit
        end if
     end if
  end do

  !! find scalar header
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,ierror,'Error reading header data')
        kopt=4
        return
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:7)=='SCALARS') then
           ibuf1=adjustl(ibuf2(8:))
           islen=max(2,scan(ibuf1,' '))-1
           sname=ibuf1(:islen)
           if (sname(:islen)==kcname(:islen2)) then
              call log_value("Scalar field found ",sname)
              exit
           else
              call log_value("Skipped scalar field ",sname)
           end if
        end if
     end if
  end do

  !! allocate scalar storage
  if(kp>0) then
     allocate(self(kp),stat=status)
     call log_alloc_check(m_name,s_name,5,status)
  else
     call log_error(m_name,s_name,6,ierror,'No scalar data')
     kopt=6
     return
  end if
  !! skip LOOKUP table default
  read(nin,fmt='(a)',iostat=status) ibuf1
  !! read coordinates
  read(nin,*,iostat=status) (self(j),j=1,kp)
  call log_read_check(m_name,s_name,7,status)
  print '("number of scalars read = ",i10)',kp
  call log_value("number of scalars read ",kp)
  kopt=0

end subroutine vfile_rscalarread
!---------------------------------------------------------------------
!> read vtk integer scalars
subroutine vfile_iscalarread(kself,kp,infile,kcname,kin,kopt)
  !! arguments
  integer(ki4), dimension(:), allocatable, intent(inout) :: kself !< integer scalar list data
  integer(ki4), intent(inout) :: kp   !< size of scalar list data
  character(*),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kcname !< name of field required
  integer(ki4), intent(inout) :: kin   !< input channel for scalar list data
  !> if positive on input, assume unit open and do not terminate execution if
  !! data is missing, return positive code if trouble, else zero
  !! if zero on input, missing data is a fatal error
  integer(ki4), intent(inout) :: kopt   !< .

  !! local
  character(*), parameter :: s_name='vfile_iscalarread' !< subroutine name
  logical :: unitused !< flag to test unit is available
  integer :: ierror   !<  whether error is to be fatal or not
  integer(ki4) :: insca   !< number of scalars
  character(80) :: sname !< name of scalar
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name
  logical :: ilfound !< flag whether field name found

  logical :: isnumb !< local variable
  external isnumb

  ibuf1=adjustl(kcname)
  islen2=max(2,scan(ibuf1,' '))-1
  if(kopt>0) then
     ierror=error_warning
     !! assume unit already open and reading infile
     if (kin==0) then
        inquire(file=infile,number=kin,iostat=status)
        if(status/=0.OR.kin==-1)then
           !! error opening file
           call log_error(m_name,s_name,1,error_fatal,'Error opening scalar list data file')
        else
           call log_error(m_name,s_name,1,log_info,'Scalar list data file opened')
        end if
     end if
     nin=kin
  else
     ierror=error_fatal
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
        call log_error(m_name,s_name,2,error_fatal,'Error opening scalar list data file')
     else
        call log_error(m_name,s_name,2,log_info,'Scalar list data file opened')
     end if

  end if

  !! File unit now sorted, get to where point data begin
  !! read local header information
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
           iltest=isnumb(ibuf2,insca,11)
           kp=insca
           exit
        else if(ibuf2(1:8)=='POLYGONS') then
           iltest=isnumb(ibuf2,insca,9)
           kp=insca
           exit
        else if(ibuf2(1:9)=='CELL_DATA') then
           iltest=isnumb(ibuf2,insca,10)
           kp=insca
           exit
        end if
     end if
  end do

  !! find scalar header
  ilfound=.FALSE.
  do
     read(nin,fmt='(a)',iostat=status) ibuf1
     !!eof
     if(status<0) then
        exit
        !! error
     else if (status>0) then
        call log_error(m_name,s_name,4,ierror,'Error reading header data')
        kopt=4
        return
     else
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:7)=='SCALARS') then
           ibuf1=adjustl(ibuf2(8:))
           islen=max(2,scan(ibuf1,' '))-1
           sname=ibuf1(:islen)
           if (sname(:islen)==kcname(:islen2)) then
              call log_value("Scalar field found ",sname)
              ilfound=.TRUE.
              exit
           else
              call log_value("Skipped scalar field ",sname)
           end if
        end if
     end if
  end do

  !! allocate scalar storage
  if(ilfound) then
     if(kp>0) then
        allocate(kself(kp),stat=status)
        call log_alloc_check(m_name,s_name,5,status)
     else
        call log_error(m_name,s_name,6,ierror,'No scalar data')
        kopt=6
        return
     end if
  else if(vfile_made_up_data/=0) then
     if(kp>0) then
        allocate(kself(kp),stat=status)
        call log_alloc_check(m_name,s_name,7,status)
     else
        call log_error(m_name,s_name,8,ierror,'No scalar data')
        kopt=8
        return
     end if
  else
     call log_error(m_name,s_name,9,ierror,'No scalar data')
     kopt=9
     return
  end if
  !! data
  if(ilfound) then
     !! skip LOOKUP table default
     read(nin,fmt='(a)',iostat=status) ibuf1
     !! read coordinates
     read(nin,*,iostat=status) (kself(j),j=1,kp)
     call log_read_check(m_name,s_name,10,status)
     print '("number of scalars read = ",i10)',kp
     call log_value("number of scalars read ",kp)
  else if(vfile_made_up_data/=0) then
     kself=vfile_made_up_data
     call log_error(m_name,s_name,10,error_warning,'Made up data')
     call log_value("number of scalars made up ",kp)
  end if
  kopt=0

end subroutine vfile_iscalarread
!---------------------------------------------------------------------
!> write vtk real scalars
subroutine vfile_rscalarwrite(self,kp,kcname,kctyp,kplot,kheader)
  !! arguments
  real(kr4), dimension(kp), intent(in) :: self !< real scalar list data
  integer(ki4), intent(in) :: kp   !< size of scalar list data
  character(*),intent(in) :: kcname !< name of field
  character(*),intent(in) :: kctyp !< type of data
  integer(ki4), intent(in) :: kplot   !< output channel for scalar list data
  integer(ki4), intent(in) :: kheader   !< header options

  !! local
  character(*), parameter :: s_name='vfile_rscalarwrite' !< subroutine name
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name

  ibuf1=adjustl(kctyp)
  islen=len_trim(ibuf1)
  ibuf1=adjustl(kcname)
  islen2=len_trim(ibuf1)

  !! output scalar header if kheader is unity
  if(kheader==1) then
     write(kplot,'(A,''_DATA'',I8)',iostat=status) kctyp(1:islen),kp
  end if
  call log_write_check(m_name,s_name,1,status)
  write(kplot,'(''SCALARS '',A,'' float 1'')',iostat=status), kcname(1:islen2)
  call log_write_check(m_name,s_name,2,status)
  write(kplot,'(''LOOKUP_TABLE default'')',iostat=status)
  call log_write_check(m_name,s_name,3,status)

  !! write scalars
  do j=1,kp
     write(kplot,fmt=cfmtbs,iostat=status) self(j)
     call log_write_check(m_name,s_name,5,status)
  end do

  call log_value("number of scalars written ",kp)

end subroutine vfile_rscalarwrite
!---------------------------------------------------------------------
!> write vtk real scalars
subroutine vfile_iscalarwrite(self,kp,kcname,kctyp,kplot,kheader)
  !! arguments
  integer(ki4), dimension(kp), intent(in) :: self !< integer scalar list data
  integer(ki4), intent(in) :: kp   !< size of scalar list data
  character(*),intent(in) :: kcname !< name of field
  character(*),intent(in) :: kctyp !< type of data
  integer(ki4), intent(in) :: kplot   !< output channel for scalar list data
  integer(ki4), intent(in) :: kheader   !< header options

  !! local
  character(*), parameter :: s_name='vfile_iscalarwrite' !< subroutine name
  integer(ki4) :: islen   !< length of scalar field name
  integer(ki4) :: islen2   !< length of required scalar field name

  ibuf1=adjustl(kctyp)
  islen=len_trim(ibuf1)
  ibuf1=adjustl(kcname)
  islen2=len_trim(ibuf1)

  !! output scalar header
  if(kheader==1) then
     write(kplot,'(A,''_DATA'',I8)',iostat=status) kctyp(1:islen),kp
  end if
  call log_write_check(m_name,s_name,1,status)
  write(kplot,'(''SCALARS '',A,'' int '')',iostat=status), kcname(1:islen2)
  call log_write_check(m_name,s_name,2,status)
  write(kplot,'(''LOOKUP_TABLE default'')',iostat=status)
  call log_write_check(m_name,s_name,3,status)

  !! write scalars
  do j=1,kp
     write(kplot,fmt='(I9)',iostat=status) self(j)
     call log_write_check(m_name,s_name,5,status)
  end do

  call log_value("number of scalars written ",kp)

end subroutine vfile_iscalarwrite
!---------------------------------------------------------------------
!> close vis plot file on unit nplot
subroutine vfile_close

  close(nplot)

end subroutine vfile_close

end module vfile_m
