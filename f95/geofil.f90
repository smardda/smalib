program geofil_p

  use date_time_m
  use const_kind_m
  use log_m
  use misc_m
  use clock_m
  use const_numphys_h
  use spl2d_m
  use spl3d_m
  use position_h
  use fmesh_h
  use bods_h
  use geobjlist_h
  use geobjlist_m
  use geofil_h
  use geofil_m
  use position_m
  use control_h
  use vcontrol_h
  use dcontrol_h
  use gfcontrol_h
  use dcontrol_m
  use gfcontrol_m
  use beq_h
  use posang_h
  use ls_m
  use btree_m
  use li_m
  use ld_m
  use dbtree_h
  use dbtree_m
  use geobj_m
  use query_m
  use posang_m
  use datline_h
  use datline_m
  use stack_m
  use indict_m
  use vfile_m
  use pcle_h
  use pcle_m
  use termplane_h
  use mtest_m

  implicit none

! Local variables
  character(*), parameter :: m_name='geofil' !< module name
  type(gffiles_t)     :: file      !< names of files
  type(gfplots_t)     :: plot      !< plot controls
  type(gfnumerics_t)  :: numerics  !< numerical control parameters for run
  type(dnumerics_t)  :: dnumerics  !< numerical control parameters
  type(geofil_t)  :: geofil      !< geometrical objects plus from vtk file
  type(geobjlist_t)  :: geobjl      !< geometrical objects
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=20) :: buf1 !< buffer input
  character(len=20) :: buf2 !< buffer input
  character(len=256) :: vtkdesc !< descriptor line for vtk files
  character(len=80),save :: iched !< vtk field file descriptor

  type(bods_t) :: bods !< type of bods for geometrical objects
  integer:: nplot !< unit for vtk files
  integer:: nread !< unit for dat files
  integer:: nin !< unit for other data
  integer(ki4):: i !< loop variable
  integer(ki4):: j !< loop variable
  integer(ki4) :: islen   !< length of input field filename
  integer(ki4), dimension(2):: idum !< dummy array
  integer(ki4) :: inarg   !< number of command line arguments
  integer(ki4) :: inda   !< index of dash in command line argument
  integer(ki4) :: iopt   !< option for reading bods
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'geofil run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'geofil: add surfaces to vtk geometry file'
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run geofil type at the command line:'
     print *, '   geofil fileroot '
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)

!--------------------------------------------------------------------------
!! lock file open check

  call clock_start(31,'lockfile time')
  call geofil_initwrite(fileroot)
  call clock_stop(31)

!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'gfcontrol_init time')
  call gfcontrol_init(fileroot)
  call gfcontrol_read(file,numerics,plot)
  call clock_stop(2)

!--------------------------------------------------------------------------
!! read vtk data

  call clock_start(4,'geofil_read time')
  call geofil_read(geofil,file,numerics)
  call clock_stop(4)

!--------------------------------------------------------------------------
!! do the main work

  call clock_start(5,'geofil_readcon time')
!call geofil_readcon(geofil)
  call geofil_objaddcon(geofil)
  call clock_stop(5)

!--------------------------------------------------------------------------
!! output file

  call clock_start(30,'outfile_init time')
  iched=geofil%geobjl%hed
  call geobjlist_makehedline(geofil%geobjl,iched,vtkdesc)
  call vfile_init(file%vtkout,vtkdesc,nplot)
  call geofil_augment(geofil,'scag')
  call geofil_writev(geofil,'full',nplot)
  call clock_stop(30)

!--------------------------------------------------------------------------
!! lock file close

  call clock_start(32,'lockfile time')
  call geofil_closewrite
  call clock_stop(32)

!--------------------------------------------------------------------------
!! cleanup and closedown
  call geofil_delete(geofil)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete
!--------------------------------------------------------------------------

end program geofil_p
