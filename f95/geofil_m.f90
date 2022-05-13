!> @addtogroup groupname4
!> @{
module geofil_m
!> @}
  use log_m
  use misc_m
  use const_numphys_h
  use const_kind_m
  use dcontrol_h
  use dcontrol_m
  use gfcontrol_h
  use gfcontrol_m
  use geofil_h
  use geobj_m
  use geobjlist_h
  use geobjlist_m
  use vfile_m
  use indict_m

  implicit none
  private

! public subroutines
  public :: &
  geofil_read,  & !< read data from file
  geofil_objaddcon,  & !< control addition of geobjlists
  geofil_initwrite, & !< open new file, making up name
  geofil_writev, &  !< write out object as vtk
  geofil_augment, & !< augment data associated with geobjlist
  geofil_closewrite, & !< close write file
  geofil_delete, & !< delete object
  geofil_close !< close file

! private variables
  character(*), parameter :: m_name='geofil_m' !< module name
  integer  :: status   !< error status
  integer, save  :: ningf=0     !< geometry file unit number
  integer, save  :: ninscal=0     !< scalx file unit number
  integer, save  :: noutgf=0      !< output file unit number
  character(len=80), save :: controlfile !< control file name
  character(len=80), save :: outputfile !< output file name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: iopt !< option
  integer(ki4), dimension(:), allocatable :: iwork !< work array
  integer(ki4) :: imaxstat   !< number of outputs for geofil data structure
  real(kr8), dimension(:), allocatable :: work !< 1D work array
  character(len=80) :: iched  !< local variable

  contains
!---------------------------------------------------------------------
!> initialise object with geobjl data
subroutine geofil_read(self,file,numerics)

  !! arguments
  type(geofil_t), intent(out) :: self !< type which data will be assigned to
  type(gffiles_t), intent(in)      :: file      !< names of files
  type(gfnumerics_t), intent(in)   :: numerics  !< numerical parameters
  !! local
  character(*), parameter :: s_name='geofil_read' !< subroutine name

  self%geobjl%ngtype=2
  ! just one file with Body data
  call geobjlist_read(self%geobjl,file%vtk,leave_open=.true.)
  ningf=0
  iopt=1
  ! set nkey here because vtk file may not contain key or body data
  self%nscag=self%geobjl%ng
  call vfile_iscalarread(self%scag,self%nscag,file%vtk,numerics%namekey,ningf,iopt)

  if (4 <=iopt.AND.iopt<=9) then
    ! no scalar data in file, fix up
    call log_value("Fixing up file for missing body/cell labels replacement ",1)
    allocate(self%scag(self%nscag),stat=status)
    call log_alloc_check(m_name,s_name,1,status)
    self%scag=1
  else  if (iopt/=0) then
    call log_error(m_name,m_name,iopt,error_fatal,'Corrupt vtk file')
  end if

  self%n=numerics

end subroutine geofil_read
!---------------------------------------------------------------------
subroutine geofil_objaddcon(self)
  !! arguments
  type(geofil_t), intent(inout) :: self !< geometrical objects and file data

  !! local
  character(*), parameter :: s_name='geofil_objaddcon' !< subroutine name
  integer :: iin      !< local control file unit number
  type(dfiles_t) :: file !< file names specifying \f$ (R,Z) \f$ geometry
  type(dnumerics_t) :: numerics !< control numerics
  integer(ki4) :: jpla !<  number of object planes to add to geobjlist
  integer(ki2par) :: igcode !< integer scalar geometry code
  real(kr4) :: zetamin   !<  minimum \f$ \zeta \f$ of any point
  real(kr4) :: zetamax   !<  maximum \f$ \zeta \f$ of any point
  logical :: filedata !< data file specifies geometry

  filedata=.FALSE.

  if (self%n%calcangle) then
     ! find angular extent of geometry
     call geobjlist_angext(self%geobjl,zetamin,zetamax)
     ! get unit for input
     call gfcontrol_getunit(iin)
  else
     zetamin=0
     zetamax=0
  end if


  ! objects based on datvtkparameters input
  do j=1,size(self%n%objadd)
     igcode=j-1
     desc_type: select case (igcode)
     case(GEOBJ_ABSORB, GEOBJ_INVISI, GEOBJ_ERRLOS, GEOBJ_CUTOUT, GEOBJ_SKYLIT)
        do jpla=1,self%n%objadd(igcode)
           call dcontrol_readnum(numerics,iin,filedata)
           ! test
           if (numerics%descode-igcode/=0) then
              call log_error(m_name,s_name,1,error_warning,'Object description does not match')
              call log_value("Requested description code",igcode)
              call log_value("Found description code",numerics%descode)
           end if
           if (filedata) then
              call dcontrol_readprogfiles(file,iin)
              call dcontrol_readatfile(file,numerics)
              call dcontrol_lines2d(numerics)
           end if
           numerics%minang=zetamin
           numerics%maxang=zetamax
           call geobjlist_objadd(self%geobjl,numerics,igcode)
           call dcontrol_delete(numerics)
        end do
     case(GEOBJ_PERMAX, GEOBJ_PERMIN, GEOBJ_ESCAPE)
        ! TO DO to write
     case default
     end select desc_type
  end do

end subroutine geofil_objaddcon
!---------------------------------------------------------------------
!> open new file, making up name
subroutine geofil_initwrite(fileroot,channel)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer, intent(out),optional :: channel   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='geofil_initwrite' !< subroutine name
  character(len=80) :: outputfile !< output file name
  !! logical :: unitused !< flag to test unit is available
  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then
  !! if (present(channel)) channel=i exit end if end do noutgf=i

  !! open file
  outputfile=trim(fileroot)//"_geofil.out"
  call log_value("Control data file",trim(outputfile))
  call misc_getfileunit(noutgf)
  open(unit=noutgf,file=outputfile,status='NEW',iostat=status)
!  uncomment these to disable lock file
!  if(status/=0)then
!     open(unit=noutgf,file=outputfile,status='REPLACE',iostat=status)
!  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open new output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open new output data file')
     stop
  end if
  if (present(channel)) channel=noutgf

end subroutine geofil_initwrite
!---------------------------------------------------------------------
!> write object data as vtk
subroutine geofil_writev(self,select,channel)

  !! arguments
  type(geofil_t), intent(in) :: self   !< object data structure
  character(*), intent(in) :: select  !< case
  integer, intent(in) :: channel   !< output channel for geofil data structure

  !! local
  character(*), parameter :: s_name='geofil_writev' !< subroutine name
  integer :: iout   !< output channel for geofil data structure
  integer(ki4) :: imaxstat   !< number of outputs for geofil data structure
  integer(ki4) :: ikey !< local variable
  integer(ki4) :: indx !< local variable

  iout=channel

  plot_type: select case(select)
  case('full')
     call geobjlist_writev(self%geobjl,'full',iout)
     call vfile_iscalarwrite(self%scag,self%nscag,self%n%namekey,'CELL',iout,0)
     call vfile_close

  case('geometry')
     ! just write geometry
     call geobjlist_writev(self%geobjl,'geometry',iout)
     call vfile_close

  case default

  end select plot_type

  call log_error(m_name,s_name,10,log_info,'vtk file produced')

end subroutine geofil_writev
!---------------------------------------------------------------------
!> augment data associated with geobjlist
subroutine geofil_augment(self,select)

  !! arguments
  type(geofil_t), intent(inout) :: self   !< object data structure
  character(*), intent(in) :: select  !< case

  !! local
  character(*), parameter :: s_name='geofil_augment' !< subroutine name
  integer(ki4) :: iaugment !< local variable
  integer(ki4) :: inkey !< local variable
  integer(ki4) :: imaxblk !< local variable
  integer(ki4), dimension(:), allocatable :: iscag !< new scalar array

  iaugment=0
  var_type: select case(select)
  case('scag')
     ! simple case of just one associated scalar associated with geometry
     if (self%geobjl%ng>self%nscag) then

        iaugment=1

        inkey=self%geobjl%ng
        allocate(iscag(inkey),stat=status)
        call log_alloc_check(m_name,s_name,1,status)
        imaxblk=maxval(self%scag)+1
        iscag(1:self%nscag)=self%scag
        iscag(self%nscag+1:inkey)=imaxblk
        ! now move  back to subroutine arguments
        deallocate(self%scag)
        allocate(self%scag(inkey),stat=status)
        call log_alloc_check(m_name,s_name,2,status)
        self%scag=iscag
        self%nscag=inkey
        deallocate(iscag)

     end if

  end select var_type

  if (iaugment>0) call log_error(m_name,s_name,1,log_info,'vtk file augmented')

end subroutine geofil_augment
!---------------------------------------------------------------------
!> close write file
subroutine geofil_closewrite

  !! local
  character(*), parameter :: s_name='geofil_closewrite' !< subroutine name

  !! close file
  close(unit=noutgf,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine geofil_closewrite
!---------------------------------------------------------------------
!> delete object
subroutine geofil_delete(self)

  !! arguments
  type(geofil_t), intent(inout) :: self !< module object
  !! local
  character(*), parameter :: s_name='geofil_delete' !< subroutine name

  call geobjlist_delete(self%geobjl)
  deallocate(self%scag)

end subroutine geofil_delete
!---------------------------------------------------------------------
!> close file
subroutine geofil_close

  !! local
  character(*), parameter :: s_name='geofil_close' !< subroutine name

  !! close file
  close(unit=ningf,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close control data file')
     stop
  end if

end subroutine geofil_close

end module geofil_m
