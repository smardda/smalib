!> @addtogroup groupname4
!> @{
module gfcontrol_m
!> @}
  use const_kind_m
  use log_m
  use gfcontrol_h
  use geobj_m
  use misc_m

  implicit none
  private

! public subroutines
  public :: &
 &gfcontrol_init, & !< open input control data file
 &gfcontrol_getunit,  & !< get unit number
 &gfcontrol_close, & !< close input control data file
 &gfcontrol_read !< read data for this object

! public variables

! private variables
  character(*), parameter :: m_name='gfcontrol_m' !< module name
  integer  :: status   !< error status
  integer, save  :: nin      !< control file unit number
  integer  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: islen !< string length
  character(len=80) :: root !< file root


  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine gfcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='gfcontrol_init' !< subroutine name
  !! logical :: unitused !< flag to test unit is available
  character(len=80) :: controlfile !< control file name

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then nin=i exit end if end do

  !! open file
  controlfile=trim(fileroot)//".ctl"
  root=fileroot
  call log_value("Control data file",trim(controlfile))
  call misc_getfileunit(nin)
  open(unit=nin,file=controlfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open control file, ",a)',controlfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open control data file')
     stop
  end if

end  subroutine gfcontrol_init
!---------------------------------------------------------------------
!> get unit number of input
subroutine gfcontrol_getunit(kunit)

  !! arguments
  integer, intent(out) :: kunit    !< log unit number

  kunit=nin

end subroutine gfcontrol_getunit
!---------------------------------------------------------------------
!> close input control data file
subroutine gfcontrol_close

  !! local
  character(*), parameter :: s_name='gfcontrol_close' !< subroutine name

  !! close file unit
  close(unit=nin)

end  subroutine gfcontrol_close
!---------------------------------------------------------------------
!> read data for this run
subroutine gfcontrol_read(file,numerics,plot,channel)

  !! arguments
  type(gffiles_t), intent(out) :: file !< file names
  type(gfnumerics_t), intent(inout) :: numerics !< object control data structure
  type(gfplots_t), intent(out) :: plot !< plot controls
  integer, intent(inout), optional :: channel   !< input channel for object control data structures

  !! local
  character(*), parameter :: s_name='gfcontrol_read' !< subroutine name
  character(len=80) :: vtk_input_file  !< vtk format geometry data input file
  character(len=80) :: vtk_output_file  !< vtk format geometry data output file
  logical :: filefound !< true of file exists
  integer(ki4) :: ndummy !< dummy namelist variable
  logical :: plot_vtk !< vtk plot selector
  logical :: centre_line_cutout !< flag whether centre line is to be a cutout
  integer(ki4) :: absorber_objects !< number of absorbers defined geometrically
  integer(ki4) :: invisible_objects !< number of invisible objects defined geometrically
  integer(ki4) :: skylight_objects !< number of skylights defined geometrically
  integer(ki4) :: beancan_objects !< number of beancan surfaces defined geometrically
  integer(ki4) :: cutout_objects !< number of cutouts defined geometrically
  logical :: angular_segment !< angular_segment
  character(len=80) :: scalar_key !< key to be used to sort scalar variable
  integer(ki4) :: ierr  !< for checking statistics

  !! file names
  namelist /inputfiles/ &
 &vtk_input_file,&
 &vtk_output_file

  !! misc parameters
  namelist /miscparameters/ &
 &ndummy

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_vtk

  !! geofil parameters
  namelist /geofilparameters/ &
 &skylight_objects,&
 &absorber_objects,&
 &invisible_objects,&
 &beancan_objects,&
 &cutout_objects,&
 &centre_line_cutout,&
 &scalar_key,&
 &angular_segment

  if(present(channel).AND.channel/=0) then
     !! assume unit already open and reading infile
     nin=channel
  end if

  !! read input file names
  vtk_input_file='null'
  vtk_output_file='null'
  read(nin,nml=inputfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_getunit(ilog)
     write(ilog,nml=inputfiles)
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%vtk   = vtk_input_file
  file%vtkout   = vtk_output_file

  !!check files exist
  call log_value("geometry data file, vtk_input_file",trim(file%vtk))
  if(file%vtk/='null') then
     inquire(file=vtk_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find geometry data file, ",a)',vtk_input_file
        call log_error(m_name,s_name,2,error_fatal,'vtk data file not found')
     end if
  end if

  if(vtk_output_file/='null') then
     !! strip any final '.vtk' string
     islen=len_trim(vtk_output_file)
     if (islen>4) then
        if (vtk_output_file(islen-3:islen)=='.vtk') then
           file%vtkout = vtk_output_file(:islen-4)
        end if
     end if
  else
     file%vtkout = trim(root)
  end if
  call log_value("geometry data file, vtk_output_file",trim(file%vtkout))

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading misc parameters")'
     call log_getunit(ilog)
     write(ilog,nml=miscparameters)
     call log_error(m_name,s_name,3,error_fatal,'Error reading misc parameters')
  end if

  plot_vtk = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     print '("Fatal error reading plot selections")'
     call log_getunit(ilog)
     write(ilog,nml=plotselections)
     call log_error(m_name,s_name,4,error_fatal,'Error reading plot selections')
  end if

  !! store values
  plot%vtk = plot_vtk

  !! set default geofil parameters
  absorber_objects=0
  invisible_objects=0
  skylight_objects=0
  beancan_objects=0
  cutout_objects=0
  centre_line_cutout=.FALSE.
  angular_segment=.TRUE.
  scalar_key='Body'

  !!read geofil parameters
  read(nin,nml=geofilparameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading geofil parameters")'
     call log_getunit(ilog)
     write(ilog,nml=geofilparameters)
     call log_error(m_name,s_name,10,error_fatal,'Error reading geofil parameters')
  end if

  !! check for valid data
  if(absorber_objects<0) &
 &call log_error(m_name,s_name,12,error_fatal,'absorber_objects must be non-negative integer')
  if(invisible_objects<0) &
 &call log_error(m_name,s_name,13,error_fatal,'invisible_objects must be non-negative integer')
  if(skylight_objects<0) &
 &call log_error(m_name,s_name,14,error_fatal,'skylight_objects must be non-negative integer')
  if(beancan_objects<0) &
 &call log_error(m_name,s_name,15,error_fatal,'beancan_objects must be non-negative integer')
  if(cutout_objects<0) &
 &call log_error(m_name,s_name,16,error_fatal,'cutout_objects must be non-negative integer')

  !! store values
  numerics%objadd=0
  numerics%objadd(GEOBJ_ABSORB)=absorber_objects
  numerics%objadd(GEOBJ_INVISI)=invisible_objects
  numerics%objadd(GEOBJ_SKYLIT)=skylight_objects
  numerics%objadd(GEOBJ_ERRLOS)=beancan_objects
  numerics%objadd(GEOBJ_CUTOUT)=cutout_objects
  numerics%calcangle=angular_segment
  numerics%namekey=scalar_key


end  subroutine gfcontrol_read

end module gfcontrol_m
