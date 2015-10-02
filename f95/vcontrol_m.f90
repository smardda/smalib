module vcontrol_m

  use const_kind_m
  use log_m
  use position_h
  use vcontrol_h
  use position_m

  implicit none
  private

! public subroutines
  public :: &
 &vcontrol_init, &
 &vcontrol_read


! private variables
  character(*), parameter :: m_name='vcontrol_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4)  :: nin      !< control file unit number
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine vcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='vcontrol_init' !< subroutine name
  logical :: unitused !< flag to test unit is available
  character(len=80) :: controlfile !< control file name


  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        nin=i
        exit
     end if
  end do

  !! open file
  controlfile=trim(fileroot)//".ctl"
  root=fileroot
  call log_value("Control data file",trim(controlfile))
  open(unit=nin,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if


end  subroutine vcontrol_init
!---------------------------------------------------------------------
!> read data for this run
subroutine vcontrol_read(file,numerics)

  !! arguments
  type(vfiles_t), intent(out) :: file !< file names
  type(vnumerics_t), intent(out) :: numerics !< input numerical parameters

  !!local
  character(*), parameter :: s_name='vcontrol_read' !< subroutine name
  character(len=80) :: vtk_input_file  !< vtk format geometry data input file
  character(len=80) :: vtk_output_file  !< vtk format geometry data output file
  logical :: filefound !< true of file exists
  integer(ki4) :: number_of_panels !< number of panels for which transform defined
  integer(ki4) :: number_of_transforms !< number of panels for which transform defined
  integer(ki4), parameter :: max_number_of_panels=100 !< maximum number of panels allowed
  integer(ki4), dimension(max_number_of_panels) :: panel_transform !< number of transform to apply
  integer(ki4), dimension(2,max_number_of_panels) :: panel_bodies !< bodies defining a panel
  type(tfmdata_t) :: ztfmdata   !< position transform numeric controls

  !! file names
  namelist /vtkfiles/ &
 &vtk_input_file, &
 &vtk_output_file

  !! misc parameters
  namelist /miscparameters/ &
 &number_of_panels,number_of_transforms

  !! panelarray parameters
  namelist /panelarrayparameters/ &
 &panel_bodies,panel_transform

  !! read input file names
  vtk_input_file='null'
  read(nin,nml=vtkfiles,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
     print '("Fatal error reading input filenames")'
  end if

  file%vtkdata   = vtk_input_file
  file%vtkout   = vtk_output_file

  !!check file exists

  call log_value("geometry data file, vtk_input_file",trim(file%vtkdata))
  if(file%vtkdata/='null') then
     inquire(file=vtk_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find Beq field data file, ",a)',vtk_input_file
        call log_error(m_name,s_name,3,error_fatal,'Beq field data file not found')
     end if
  end if

  !! create output file names from root

  !! output file
  file%vout = trim(root)//"_v.out"
  filefound=.false.
  !     inquire(file=file%vout,exist=filefound)
  !     call log_value("vtktfm output file",trim(file%vout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%vout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,4,error_fatal,'Output file already exists')
  end if

  !! set default misc parameters
  number_of_panels = 1
  number_of_transforms = 1

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error reading misc parameters')
     print '("Fatal error reading misc parameters")'
  end if

  !! check for valid data
  if(number_of_panels<=0) &
 &call log_error(m_name,s_name,11,error_fatal,'number of panels must be >0')
  if(number_of_panels>max_number_of_panels) then
     call log_value("max_number_of_panels parameter",max_number_of_panels)
     call log_error(m_name,s_name,12,error_fatal,'too many panels: increase parameter')
  end if
  if(number_of_transforms<=0) &
 &call log_error(m_name,s_name,13,error_fatal,'number of transforms must be >0')

  numerics%npans=number_of_panels
  numerics%vptfm%ntfmtyp=number_of_transforms

  !      allocate(panel_transform(number_of_panels), stat=status)
  !      call log_alloc_check(m_name,s_name,13,status)
  panel_transform=0
  panel_bodies=0

  !!read panelarray parameters
  read(nin,nml=panelarrayparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,20,error_fatal,'Error reading array parameters')
     print '("Fatal error reading array parameters")'
  end if

  !! check for valid panel data
  if(minval(panel_transform(:number_of_panels))<0) &
 &call log_error(m_name,s_name,21,error_fatal,'negative transforms not allowed')

  !! check for valid panel bodies
  if(minval(panel_bodies(:,:number_of_panels))<=0) &
 &call log_error(m_name,s_name,22,error_fatal,'bodies must have positive identifiers')
  !! allocate and assign
  !! panel transform parameters
  allocate(numerics%pantfm(number_of_panels), stat=status)
  call log_alloc_check(m_name,s_name,23,status)
  numerics%pantfm=panel_transform(:number_of_panels)
  !! panel bodies
  allocate(numerics%panbod(2,number_of_panels), stat=status)
  call log_alloc_check(m_name,s_name,29,status)
  numerics%panbod=panel_bodies(:,:number_of_panels)
  !DPRT     write(*,*) numerics%panbod !DPRT
  !      do j=1,2
  !      idum(j)=geobjl%nodl(j)
  !      end do

  allocate(numerics%vptfm%scale(3,number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,25,status)
  allocate(numerics%vptfm%offset(3,number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,26,status)
  allocate(numerics%vptfm%matrix(3,3,number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,27,status)
  allocate(numerics%vptfm%ntfm(number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,28,status)
  !! transform array parameters
  do j=1,number_of_transforms
     call log_value('read namelist',j)
     call position_readcon(ztfmdata,nin)
     numerics%vptfm%scale(:,j)=ztfmdata%scale
     numerics%vptfm%offset(:,j)=ztfmdata%offset
     numerics%vptfm%matrix(:,:,j)=ztfmdata%matrix
     numerics%vptfm%ntfm(j)=ztfmdata%ntfm
  end do

end  subroutine vcontrol_read

end module vcontrol_m
