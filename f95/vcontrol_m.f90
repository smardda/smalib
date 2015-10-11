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
  integer(ki4) :: islen !< string length
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
  integer(ki4), parameter  :: maximum_number_of_files=100 !< maximum number of files allowed
  character(len=80), dimension(maximum_number_of_files) :: vtk_input_file !< vtk format geometry data input file
  integer(ki4), dimension(maximum_number_of_files) :: number_of_copies !< number of copies to make of input file
  character(len=80), dimension(maximum_number_of_files) :: vtk_label  !< labels vtk format geometry data input file
  integer(ki4), dimension(maximum_number_of_files) :: vtk_label_number !< numeric label for input file (not used)
  character(len=20) :: angle_units  !< units, either radian(s) or degree(s)
  character(len=80) :: vtk_output_file  !< vtk format geometry data output file
  integer(ki4) :: number_of_panels !< number of panels for which bodies defined
  integer(ki4) :: number_of_transforms !< number of panels for which transform defined
  integer(ki4) :: max_number_of_files !< number >= no. of files for which transform defined
  integer(ki4) :: max_number_of_panels !< number >= no. of panels for which transform defined
  integer(ki4) :: max_number_of_transforms !< number >= no. of panels for which transform defined
  integer(ki4), parameter  :: maximum_number_of_panels=100 !< maximum number of panels allowed
  integer(ki4), dimension(maximum_number_of_panels) :: panel_transform !< number of transform to apply
  character(len=20), dimension(maximum_number_of_panels) :: transform_id !< id of transform to apply
  integer(ki4), dimension(2,maximum_number_of_panels) :: panel_bodies !< bodies defining a panel
  logical :: filefound !< true of file exists
  type(tfmdata_t) :: ztfmdata   !< position transform numeric controls
  integer(ki4) :: infil=0 !< number of files detected
  integer(ki4) :: inpan=0 !< number of panels detected
  integer(ki4) :: intfm=0 !< number of transforms detected
  integer(ki4) :: iflag=0 !< allow position_readcon to return on error

  !> misc parameters, unusually comes first
  namelist /miscparameters/ &
 &max_number_of_files, angle_units, &
 &max_number_of_panels,max_number_of_transforms,&
 &number_of_panels,number_of_transforms

  !! file names
  namelist /vtkfiles/ &
 &vtk_input_file, number_of_copies, &
 &vtk_label, vtk_label_number, &
 &vtk_output_file

  !! panelarray parameters
  namelist /panelarrayparameters/ &
 &panel_bodies,panel_transform, &
 &transform_id

  !---------------------------------------------------------------------
  !! set default misc parameters
  max_number_of_panels = 1
  max_number_of_files = 1
  max_number_of_transforms = 1
  number_of_panels = 0
  number_of_transforms = 0
  angle_units = 'radians'

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading misc parameters')
     print '("Fatal error reading misc parameters")'
  end if

  !! check for valid data
  if(max_number_of_files<=0) &
 &call log_error(m_name,s_name,2,error_fatal,'max number of files must be >0')
  if(max_number_of_files>maximum_number_of_files) then
     call log_value("maximum_number_of_files parameter",maximum_number_of_files)
     call log_error(m_name,s_name,3,error_fatal,'too many files: increase parameter')
  end if
  !!added for upwards compatibility
  if (number_of_panels>0) max_number_of_panels=number_of_panels
  if(max_number_of_panels<=0) &
 &call log_error(m_name,s_name,11,error_fatal,'max number of panels must be >0')
  if(max_number_of_panels>maximum_number_of_panels) then
     call log_value("maximum_number_of_panels parameter",maximum_number_of_panels)
     call log_error(m_name,s_name,12,error_fatal,'too many panels: increase parameter')
  end if
  if (number_of_transforms>0) max_number_of_transforms=number_of_transforms
  if(max_number_of_transforms<=0) &
 &call log_error(m_name,s_name,13,error_fatal,'max number of transforms must be >0')
  if(angle_units(1:6)/='degree'.AND.angle_units(1:6)/='radian') then
     call log_value("angle_units not recognised",angle_units)
     call log_error(m_name,s_name,14,error_warning,'radians assumed')
  end if
  numerics%angles=angle_units(1:6)

  !---------------------------------------------------------------------
  !! read input file names and associated data
  vtk_input_file='null'
  number_of_copies=1
  vtk_label='null'
  vtk_label_number=0
  read(nin,nml=vtkfiles,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,21,error_fatal,'Error reading input filenames')
     print '("Fatal error reading input filenames")'
  end if

  ! check number of files against actuals
  do j=1,maximum_number_of_files
     if (trim(vtk_input_file(j))=='null') exit
     infil=j
  end do
  file%nvtkdata=infil
  call log_value("actual number of files",infil)

  !!check files exist
  do j=1,infil
     call log_value("geometry data file number",j)
     call log_value("geometry data file, vtk_input_file",trim(vtk_input_file(j)))
     if(vtk_input_file(j)/='null') then
        inquire(file=vtk_input_file(j),exist=filefound)
        if(.not.filefound) then
           !! error opening file
           print '("Fatal error: Unable to find field data file, ",a)',vtk_input_file(j)
           call log_error(m_name,s_name,22,error_fatal,'field data file not found')
        end if
        if (number_of_copies(j)<0) then
           call log_error(m_name,s_name,24,error_fatal,'Number of copies is negative')
        end if
     end if
  end do

  ! allocate and assign files
  allocate(file%vtkdata(infil), file%vtkcopies(infil),&
 &file%vtklabel(infil), stat=status)
  call log_alloc_check(m_name,s_name,23,status)
  do j=1,infil
     file%vtkdata(j)   = vtk_input_file(j)
     file%vtklabel(j)   = vtk_label(j)
     file%vtkcopies(j)   = number_of_copies(j)
  end do

!! strip any final '.vtk' string
     islen=len_trim(vtk_output_file)
     if (islen>4) then
        if (vtk_output_file(islen-3:islen)=='.vtk') then
           file%vtkout = vtk_output_file(:islen-4)
        end if
     end if

  !> create output file names from root
  !! output file
  file%vout = trim(root)//"_v.out"
  filefound=.false.
  call log_value("vtktfm output file",trim(file%vout))
  inquire(file=file%vout,exist=filefound)
  if(filefound) then
     !! error opening file
     print '("Warning: output file ",a, " already exists")',trim(file%vout)
     !    print '("Remove it and restart run")'
     call log_error(m_name,s_name,31,error_warning,'Output file already exists')
  end if

  !---------------------------------------------------------------------
  panel_transform=-1
  panel_bodies=0
  transform_id='null'

  !!read panelarray parameters
  read(nin,nml=panelarrayparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading array parameters")'
     call log_error(m_name,s_name,32,error_fatal,'Error reading array parameters')
  end if

  !! check for valid panel data
  ! check number of panels against actuals
  do j=1,maximum_number_of_panels
     if (panel_bodies(1,j)<=0) exit
     inpan=j
  end do
  call log_value("actual number of panels",inpan)
  if (inpan==0) &
 &call log_error(m_name,s_name,33,error_fatal,'No panels present')
  numerics%npans=inpan
  do j=1,inpan
     if (panel_transform(j)<0.AND.transform_id(j)=='null') then
        call log_value("panel number",j)
        call log_error(m_name,s_name,34,error_fatal,'Transform not defined')
     end if
     if (panel_bodies(2,j)<=0) panel_bodies(2,j)=panel_bodies(1,j)
  end do

  !! allocate and assign
  !! panel transform parameters
  allocate(numerics%pantfm(inpan), stat=status)
  call log_alloc_check(m_name,s_name,43,status)
  numerics%pantfm=panel_transform(:inpan)
  !! panel bodies
  allocate(numerics%panbod(2,inpan), stat=status)
  call log_alloc_check(m_name,s_name,49,status)
  numerics%panbod=panel_bodies(:,:inpan)

  !---------------------------------------------------------------------
  ! define transformations
  allocate(numerics%vptfm%scale(3,max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,45,status)
  allocate(numerics%vptfm%offset(3,max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,46,status)
  allocate(numerics%vptfm%matrix(3,3,max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,47,status)
  allocate(numerics%vptfm%ntfm(max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,48,status)
!ID  allocate(numerics%vptfm%id(max_number_of_transforms), stat=status)
  call log_alloc_check(m_name,s_name,49,status)
  !! transform array parameters
  do j=1,max_number_of_transforms
     call log_value('read namelist',j)
     call position_readcon(ztfmdata,nin,iflag)
     if (iflag/=0) exit
     intfm=j
     numerics%vptfm%scale(:,j)=ztfmdata%scale
     numerics%vptfm%offset(:,j)=ztfmdata%offset
!F11 write(*,*) 'ztfmdata%offset',ztfmdata%offset
     numerics%vptfm%matrix(:,:,j)=ztfmdata%matrix
     numerics%vptfm%ntfm(j)=ztfmdata%ntfm
!ID     numerics%vptfm%id(j)=ztfmdata%id
  end do
  call log_value('No transforms present',intfm)

  ! assign transform numbers
  do j=1,inpan
     if (transform_id(j)(1:4)=='unit') then
        numerics%pantfm(j)=0
        cycle
     end if
!ID     do k=1,intfm
!ID        if (transform_id(j)==numerics%vptfm%id(k)) then
!ID           numerics%pantfm(j)=numerics%vptfm%ntfm(k)
!ID           exit
!ID        end if
!ID     end do
  end do
  !! final check
  do j=1,inpan
     if (numerics%pantfm(j)<0.OR.&
 &   numerics%pantfm(j)>intfm) then
        call log_value("panel number",j)
        call log_error(m_name,s_name,50,error_fatal,'Transform not defined')
     end if
  end do
  numerics%vptfm%ntfmtyp=intfm
!F11 write(*,*) 'numerics%pantfm(j)',(numerics%pantfm(j),j=1,inpan)

end subroutine vcontrol_read

end module vcontrol_m
