module bcontrol_m

  use const_kind_m
  use log_m
  use beq_h
  use beq_m

  implicit none
  private

! public subroutines
  public :: &
 &bcontrol_init, &
 &bcontrol_read


! public types
  type, public :: bfiles_t
   !! file names
     character(len=80)  :: geoqout       !< output data
     character(len=80)  :: log            !< log file
     character(len=80)  :: eqdsk         !< equilibrium type eqtype file name
     character(len=80)  :: equil         !< generic equilibrium type file name
     character(len=80)  :: vtkdata         !< vtk input data file
     character(len=80)  :: laudata         !< launch data
     character(len=80)  :: cartv   !< \f$ \bf{B} \f$ in Cartesians
     character(len=80)  :: allcartv   !< all \f$ \bf{B} \f$ in Cartesians
     character(len=80)  :: fptv    !< vtk plot file of functions of \f$ (\psi,\theta) \f$
     character(len=80)  :: allgeoq  !< plot file of all vtk info needed for SMITER-GEOQ
     character(len=80)  :: frzv    !< vtk plot file of functions of \f$ (R,Z) \f$
     character(len=80)  :: frzzeta    !< vtk plot file of functions of \f$ (R,Z,\zeta) \f$
     character(len=80)  :: frzxi    !< vtk plot file of functions of \f$ (R,Z,\xi) \f$
     character(len=80)  :: geofld  !< vtk plot file of geometry and field
     character(len=80)  :: geofldq  !< vtk plot file of quantised geometry and field
     character(len=80)  :: gnu !< gnuplot file of \f$ \psi(R,Z) \f$
     character(len=80)  :: gnusil !< gnuplot file of silhouette as function of \f$ (R,Z) \f$
     character(len=80)  :: gnusilm !< gnuplot file of silhouette as function of \f$ (\psi,\theta) \f$
  end type bfiles_t


  type, public :: bplots_t
   !! vtk plot output selectors
     logical  :: cartv   !< \f$ \bf{B} \f$ in Cartesians
     logical  :: allcartv   !< all \f$ \bf{B} \f$ in Cartesians
     logical  :: fptv    !< vtk plot file of functions of \f$ (\psi,\theta) \f$
     logical  :: allgeoq    !< plot file of all vtk info needed for SMITER-GEOQ
     logical  :: frzv    !< vtk plot file of functions of \f$ (R,Z) \f$
     logical  :: frzzeta    !< vtk plot file of functions of \f$ (R,Z,\zeta) \f$
     logical  :: frzxi    !< vtk plot file of functions of \f$ (R,Z,\xi) \f$
     logical  :: geofld  !< vtk plot file of geometry and field
     logical  :: geofldq  !< vtk plot file of quantised geometry and field
     logical  :: gnu !< gnuplot file of \f$ \psi(R,Z) \f$
     logical  :: gnusil !< gnuplot file of silhouette as function of \f$ (R,Z) \f$
     logical  :: gnusilm !< gnuplot file of silhouette as function of \f$ (\psi,\theta) \f$
     logical  :: gnuptz !< gnuplot file of \f$ fns(\psi,\theta) \f$
  end type bplots_t


! private variables
  character(*), parameter :: m_name='bcontrol_m' !< module name
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
subroutine bcontrol_init(fileroot)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  !! local
  character(*), parameter :: s_name='bcontrol_init' !< subroutine name
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


end  subroutine bcontrol_init
!---------------------------------------------------------------------
!> read data for this run
subroutine bcontrol_read(file,numerics,plot)

  !! arguments
  type(bfiles_t), intent(out) :: file !< file names
  type(bnumerics_t), intent(out) :: numerics !< input numerical parameters
  type(bplots_t), intent(out) :: plot !< vtk plot selectors

  !!local
  character(*), parameter :: s_name='bcontrol_read' !< subroutine name
  character(len=80) :: eqdsk_input_file  !< EQDSK G format data input file
  character(len=80) :: equil_input_file  !< generic equilibrium input file
  character(len=80) :: vtk_input_file  !< vtk format geometry data input file
  character(len=80) :: lau_input_file  !< launch data input file
  character(len=80) :: icsuf  !< eqdsk file suffix
  logical :: filefound !< true of file exists
  integer(ki4) :: dummy_number !< dummy
  integer(ki4) :: indot !< position of dot in filename
  integer(ki4) :: ilent !< local variable
  logical :: plot_cartv !< vtk plot selector
  logical :: plot_allcartv !< vtk plot selector
  logical :: plot_fptv !< vtk plot selector
  logical :: plot_allgeoq !< vtk plot selector
  logical :: plot_frzv !< vtk plot selector
  logical :: plot_frzzeta !< vtk plot selector
  logical :: plot_frzxi !< vtk plot selector
  logical :: plot_geofld !< vtk plot selector
  logical :: plot_geofld_quantised !< vtk plot selector
  logical :: plot_gnu !< gnu plot selector
  logical :: plot_gnuptz !< gnuptz plot selector
  logical :: plot_gnusil !< gnusil plot selector
  logical :: plot_gnusilm !< gnusilm plot selector

  !! file names
  namelist /inputfiles/ &
 &eqdsk_input_file,vtk_input_file,lau_input_file, &
 &equil_input_file

  !! misc parameters
  namelist /miscparameters/ &
 &dummy_number

  !! plot selection parameters
  namelist /plotselections/ &
 &plot_cartv, &
 &plot_allcartv, &
 &plot_fptv, &
 &plot_allgeoq, &
 &plot_frzv, &
 &plot_frzzeta, &
 &plot_frzxi, &
 &plot_geofld, &
 &plot_geofld_quantised, &
 &plot_gnu, &
 &plot_gnuptz, &
 &plot_gnusil, &
 &plot_gnusilm

  !! read input file names
  vtk_input_file='null'
  lau_input_file='null'
  equil_input_file='null'
  eqdsk_input_file='null'
  read(nin,nml=inputfiles,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,1,error_fatal,'Error reading input filenames')
     print '("Fatal error reading input filenames")'
  end if

  if ( equil_input_file=='null') then
     equil_input_file   = eqdsk_input_file
  end if
  file%eqdsk   = equil_input_file
  file%equil   = equil_input_file
  file%vtkdata   = vtk_input_file
  file%laudata   = lau_input_file

  !!check files exist

  call log_value("beq field data file, equil_input_file",trim(file%equil))
  inquire(file=equil_input_file,exist=filefound)
  if(.not.filefound) then
     !! error opening file
     print '("Fatal error: Unable to find Beq field data file, ",a)',equil_input_file
     call log_error(m_name,s_name,2,error_fatal,'Beq field data file not found')
  end if

  ! set equilibrium type depending on suffix
  indot=index(equil_input_file,'.',.TRUE.)
  if (indot>80) then
     call log_error(m_name,s_name,2,error_fatal,'Beq field data file has no suffix')
  end if
  icsuf=equil_input_file(indot+1:)
  ilent=len_trim(icsuf)
  call lowor(icsuf,1,ilent)
  numerics%eqtype=icsuf

  call log_value("geometry data file, vtk_input_file",trim(file%vtkdata))
  if(file%vtkdata/='null') then
     inquire(file=vtk_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find Beq geometry data file, ",a)',vtk_input_file
        call log_error(m_name,s_name,3,error_fatal,'Beq geometry data file not found')
     end if
  end if

  if(file%laudata/='null') then
     inquire(file=lau_input_file,exist=filefound)
     if(.not.filefound) then
        !! error opening file
        print '("Fatal error: Unable to find launch data file, ",a)',lau_input_file
        call log_error(m_name,s_name,4,error_fatal,'launch data file not found')
     end if
  end if

  !! create output file names from root

  !! output file
  file%geoqout = trim(root)//"_geoq.out"
  inquire(file=file%geoqout,exist=filefound)
  call log_value("geoq output file",trim(file%geoqout))
  if(filefound) then
     !! error opening file
     print '("Fatal error: Output file ",a, " already exists")',trim(file%geoqout)
     print '("Remove it and restart run")'
     call log_error(m_name,s_name,5,error_fatal,'Output file already exists')
  end if

  !!vtk file roots
  file%cartv     =trim(root)//"_cartv"
  file%allcartv     =trim(root)//"_allcartv"
  file%fptv    =trim(root)//"_fptv"
  file%allgeoq    =trim(root)//"_allgeoq"
  file%frzv    =trim(root)//"_frzv"
  file%frzzeta    =trim(root)//"_frzzeta"
  file%frzxi    =trim(root)//"_frzxi"
  file%geofld =trim(root)//"_geofld"
  file%geofldq    =trim(root)//"_geofldq"
  !!gnu file roots
  file%gnu     =trim(root)//"_gnu"
  file%gnusil     =trim(root)//"_gnusil"
  file%gnusilm     =trim(root)//"_gnusilm"

  !! set default misc parameters
  dummy_number = 10

  !!read misc parameters
  read(nin,nml=miscparameters,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,10,error_fatal,'Error reading misc parameters')
     print '("Fatal error reading misc parameters")'
  end if

  !! check for valid data
  if(dummy_number<0) &
 &call log_error(m_name,s_name,11,error_fatal,'dummy_number must be >=0')

  !W     numerics%ndummy = dummy_number

  !! set default plot selections
  plot_cartv = .false.
  plot_allcartv = .false.
  plot_fptv = .false.
  plot_allgeoq = .false.
  plot_frzv = .false.
  plot_frzzeta = .false.
  plot_frzxi = .false.
  plot_geofld = .false.
  plot_geofld_quantised = .false.
  plot_gnu = .false.
  plot_gnuptz = .false.
  plot_gnusil = .false.
  plot_gnusilm = .false.

  !!read plot selections
  read(nin,nml=plotselections,iostat=status)
  if(status/=0) then
     call log_error(m_name,s_name,12,error_fatal,'Error reading plot selections')
     print '("Fatal error reading plot selections")'
  end if

  !! store values
  plot%cartv   = plot_cartv
  plot%allcartv   = plot_allcartv
  plot%fptv    = plot_fptv
  plot%allgeoq    = plot_allgeoq
  plot%frzv     = plot_frzv
  plot%frzzeta     = plot_frzzeta
  plot%frzxi     = plot_frzxi
  plot%geofld    = plot_geofld
  plot%geofldq     = plot_geofld_quantised
  plot%gnu     = plot_gnu
  plot%gnuptz     = plot_gnuptz
  plot%gnusil     = plot_gnusil
  plot%gnusilm     = plot_gnusilm

  call beq_readcon(numerics,nin)

end  subroutine bcontrol_read

end module bcontrol_m
