module dcontrol_m

  use const_kind_m
  use const_numphys_h
  use log_m
  use dcontrol_h

  implicit none
  private

! public subroutines
  public :: &
 &dcontrol_init, & !< open input control data file
 &dcontrol_close, & !< close input control data file
 &dcontrol_read !< read data for this run

! private variables
  character(*), parameter :: m_name='dcontrol_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4), save  :: nin      !< control file unit number
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  real(kr8) :: zdum !< dummy real
  character(len=80) :: root !< file root

  contains
!---------------------------------------------------------------------
!> open input control data file
subroutine dcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='dcontrol_init' !< subroutine name
  character(len=80) :: controlfile !< control file name

  call  misc_getfileunit(nin)

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

end  subroutine dcontrol_init
!---------------------------------------------------------------------
!> close input control data file
subroutine dcontrol_close

  !! local
  character(*), parameter :: s_name='dcontrol_close' !< subroutine name

  !! close file unit
  close(unit=nin,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close file unit, ",i5)',nin
     call log_error(m_name,s_name,1,error_fatal,'Cannot close data file')
     stop
  end if

end  subroutine dcontrol_close
!---------------------------------------------------------------------
!> read data for this run
subroutine dcontrol_read(file,numerics,plot)

  !! arguments
  type(dfiles_t), intent(out) :: file !< file names
  type(dnumerics_t), intent(out) :: numerics !< control numerics
  type(dplots_t), intent(out), optional :: plot !< plot selectors

  !!local
  character(*), parameter :: s_name='dcontrol_read' !< subroutine name
  logical, parameter :: noplotnml=.TRUE. !< shut off plot namelist read
  logical :: filefound !< true if file exists

  character(len=80) :: r_input_file  !< position coordinate 1 data input file
  character(len=80) :: z_input_file !< position coordinate 2 data input file
  character(len=80) :: rz_input_file  !< position coordinate 1 & 2 data input file
  character(len=80) :: transform_type  !< transform type
  real(kr8) :: start_angle !< start angle
  real(kr8) :: finish_angle !< finish angle
  real(kr8), dimension(3) :: start_position !< start position
  real(kr8), dimension(3) :: finish_position !< finish position
  integer(ki4) :: number_of_divisions !< number of divisions
  character(len=20) :: angle_units  !< units, either radian(s) or degree(s)
  real(kr4) :: angfac=1. !< angles default is radians

  logical :: plot_vtk !< vtk plot selector
  logical :: plot_gnu !< gnuplot plot selector

  !! file names
  namelist /progfiles/ &
 &r_input_file, &
 &z_input_file, &
 &rz_input_file

  !! transformation parameters
  namelist /datvtk_parameters/ &
 &transform_type,&
 &start_angle, finish_angle, &
 &start_position, finish_position, &
 &number_of_divisions

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_vtk, &
 &plot_gnu

  !---------------------------------------------------------------------
  !! read input file names
  r_input_file='null'
  z_input_file='null'
  rz_input_file='null'
  read(nin,nml=progfiles,iostat=status)
  if(status/=0) then
     print '("Fatal error reading input filenames")'
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
  end if

  file%rfile = r_input_file
  file%zfile = z_input_file
  file%rzfile = rz_input_file

  !!check files exist

  call log_value("datvtk data file, r_input_file",trim(file%rfile))
  if(file%rfile/='null') then
     inquire(file=r_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find datvtk data file, ",a)',r_input_file
        call log_error(m_name,s_name,2,error_fatal,'datvtk data file not found')
     end if
  end if
  call log_value("datvtk data file, r_input_file",trim(file%zfile))
  if(file%zfile/='null') then
     inquire(file=z_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find datvtk data file, ",a)',z_input_file
        call log_error(m_name,s_name,3,error_fatal,'datvtk data file not found')
     end if
  end if
  if(file%rzfile/='null') then
     inquire(file=rz_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find datvtk data file, ",a)',rz_input_file
        call log_error(m_name,s_name,4,error_fatal,'datvtk data file not found')
     end if
  end if

  !! create output file names from root

  !!vtk file
  file%vtk = trim(root)//"_vtk"
  !!gnu file roots
  file%gnu = trim(root)//"_gnu"

  !---------------------------------------------------------------------
  !! set default datvtk parameters
  angle_units='degree'
  transform_type = 'rotate'
  start_angle = 0
  finish_angle = 90
  start_position = 0
  finish_position = (/0,0,1/)
  number_of_divisions = 10

  !!read datvtk parameters
  read(nin,nml=datvtk_parameters,iostat=status)
  if(status/=0) then
     print '("Fatal error reading datvtk parameters")'
     call log_error(m_name,s_name,10,error_fatal,'Error reading datvtk parameters')
  end if

  !! check for valid data
  if(transform_type/='rotate'.AND.transform_type/='translate') then
     call log_error(m_name,s_name,11,error_fatal,'transform type not recognised')
  end if
  ! Check angles
  !! set up factor for angles
  if (angle_units(1:6)=='degree') angfac=const_degrad
  if(start_angle*angfac<-2.000001_kr8*const_pid.OR.start_angle*angfac>2.000001_kr8*const_pid) then
     call log_error(m_name,s_name,20,error_fatal,'start angle not in range -360 to +360')
  end if
  if(finish_angle*angfac<-2.000001_kr8*const_pid.OR.finish_angle*angfac>2.000001_kr8*const_pid) then
     call log_error(m_name,s_name,21,error_fatal,'finish angle not in range -360 to +360')
  end if
  ! positive integer parameter
  if(number_of_divisions<=0) then
     call log_error(m_name,s_name,30,error_fatal,'number of divisions must be positive')
  end if

  numerics%tfm=transform_type
  numerics%stang=start_angle*angfac
  numerics%finang=finish_angle*angfac
  numerics%stpos=start_position
  numerics%finpos=finish_position
  numerics%div=2*((number_of_divisions+1)/2)

  !---------------------------------------------------------------------
  !! set default plot selections
  plot_vtk = .false.
  plot_gnu = .false.

  if (noplotnml) then
     !!read plot selections
     read(nin,nml=plotselections,iostat=status)
     if(status/=0) then
        print '("Fatal error reading plot selections")'
        call log_error(m_name,s_name,50,error_fatal,'Error reading plot selections')
     end if
  end if

  if (present(plot)) then
     !! store values
     plot%vtk = plot_vtk
     plot%gnu = plot_gnu
  end if

  call dcontrol_close

  !---------------------------------------------------------------------
  !! read silhouette files
  if (rz_input_file=='null') then
     call log_value("r input file",trim(r_input_file))
     open(unit=nin,file=r_input_file,status='OLD',iostat=status)
     call log_open_check(m_name,s_name,40,status)

     i=0
     do
        read(nin,*,iostat=status) zdum
        if(status/=0) exit
        i=i+1
     end do
     allocate(numerics%r(i),numerics%z(i),stat=status)
     call log_alloc_check(m_name,s_name,41,status)
     numerics%npos=i

     !! read into store values
     rewind(nin,iostat=status)
     call log_read_check(m_name,s_name,42,status)
     do j=1,i
        read(nin,*,iostat=status) numerics%r(j)
        call log_read_check(m_name,s_name,43,status)
     end do
     call log_value("number of values read from r input file",i)

     call dcontrol_close

     call log_value("z input file",trim(z_input_file))
     open(unit=nin,file=z_input_file,status='OLD',iostat=status)
     call log_open_check(m_name,s_name,44,status)

     i=0
     !! read into store values
     do j=1,numerics%npos
        i=i+1
        read(nin,*,iostat=status) numerics%z(j)
        call log_read_check(m_name,s_name,46,status)
     end do
     call log_value("number of values read from z input file",i)

  else
     call log_value("rz input file",trim(rz_input_file))
     open(unit=nin,file=rz_input_file,status='OLD',iostat=status)
     call log_open_check(m_name,s_name,50,status)

     i=0
     do
        read(nin,*,iostat=status) zdum, zdum
        if(status/=0) exit
        i=i+1
     end do
     allocate(numerics%r(i),numerics%z(i),stat=status)
     call log_alloc_check(m_name,s_name,51,status)
     numerics%npos=i

     !! read into store values
     rewind(nin,iostat=status)
     call log_read_check(m_name,s_name,52,status)
     do j=1,i
        read(nin,*,iostat=status) numerics%r(j),numerics%z(j)
        call log_read_check(m_name,s_name,53,status)
     end do
     call log_value("number of values read from rz input file",i)
  end if

  call dcontrol_close

end  subroutine dcontrol_read

subroutine misc_getfileunit(kunit)
  integer(ki4), intent(out) :: kunit !< file unit

  integer(ki4) :: i !< loop counter
  logical :: unitused !< flag to test unit is available
  !! get file unit
  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        kunit=i
        exit
     end if
  end do
end subroutine misc_getfileunit

end module dcontrol_m
