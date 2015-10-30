program datvtkprog

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use control_h
  use vcontrol_h
  use log_m
  use clock_m
  use spl2d_m
  use spl3d_m
  use beq_h
  use posang_h
  use posang_m

  use position_h
  use geobjlist_h
  use position_m
  use ls_m
  use btree_m
  use geobj_m
  use query_m
  use geobjlist_m

  use vcontrol_m
  use vfile_m
  use gfile_m
  use datline_h
  use datline_m
  use dfile_m
  use bods_m
  use stack_m

  implicit none


! Local variables
  character(*), parameter :: m_name='datvtkprog' !< module name
  type(vfiles_t)     :: file      !< names of files
  type(vnumerics_t)  :: numerics  !< numerical control parameters
  type(geobjlist_t)  :: geobjl      !< geometrical objects
  type(geobjlist_t)  :: geobj2      !< geometrical objects
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=3) :: optarg='nxx' !< optional argument
  character(len=256) :: vtkdesc !< descriptor line for vtk files
  character(len=80),save :: iched !< vtk field file descriptor

  integer(ki4), dimension(:), allocatable :: bods !< array of bodies for geometrical objects
  integer(ki4):: nplot !< unit for vtk files
  integer(ki4):: nread !< unit for dat files
  integer(ki4):: nin !< unit for other data
  integer(ki4):: nscal !< number of scalars (body identifiers)
  integer(ki4):: i !< loop variable
  integer(ki4):: j !< loop variable
  integer(ki4) :: islen   !< length of input field filename
  integer(ki4), dimension(2):: idum !< dummy array
  logical :: iltest !< logical flag
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'datvtkprog run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'datvtkprog: convert dat to vtk geometry file'
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run datvtkprog type at the command line:'
     print *, '   datvtk fileroot [opt]'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
     iltest=(command_argument_count()==2)
     if (iltest) then
        call get_command_argument(2,value=optarg)
     end if
!! strip any final '.dat' string
     islen=len_trim(fileroot)
     if (islen>4) then
        if (fileroot(islen-3:islen)=='.dat') then
           fileroot(islen-3:islen)='    '
        end if
     end if
  end if

! strip leading dash

!! start log
  call log_init(fileroot,timestamp)
!  write(*,*) 'this output helps gfortran sometimes', fileroot
!--------------------------------------------------------------------------
!sk  !! read control file (nominally)
!sk
!sk      call clock_start(2,'vcontrol_init time')
!sk      call vcontrol_init(fileroot)
!  open(file='rpiu.dat',unit=nread)
!sk      call vcontrol_read(file,numerics)
!sk      call clock_stop(2)
!--------------------------------------------------------------------------
!! do the main work
!! read  geobjl data

  call clock_start(4,'geobjlist_(d)read time')
  input_type: select case(optarg(1:1))
  case('v')
     call geobjlist_read(geobjl,trim(fileroot)//'.vtk',iched)
     inquire(file=trim(fileroot)//'.vtk',number=nread)
  case default
     iched='converted dat file '//fileroot//'.dat'
     call dfile_init(fileroot,nread)
     call geobjlist_dread(geobjl,nread)
  end select input_type
  call clock_stop(4)
!--------------------------------------------------------------------------
!! shell/reorient  geobjl data as needed

  call clock_start(5,'geobjlist_orientri time')
  transform_type: select case(optarg(2:2))
  case('o')
     call geobjlist_orientri(geobjl)
  case('s')
     call geobjlist_shelltets(geobjl,geobj2)
     call geobjlist_delete(geobjl)
     call geobjlist_copy(geobj2,geobjl,1)
     call geobjlist_delete(geobj2)
     call geobjlist_orientri(geobjl)
  end select transform_type
  call clock_stop(5)
!--------------------------------------------------------------------------
!! output file

  call clock_start(30,'outfile_init time')
!     write(*,*) 'iched=',iched
  call bods_init(bods,geobjl,1)
  divide_type: select case(optarg(3:3))
  case('d','s')
!! write out as separate files
     call bods_write(bods,size(bods),geobjl,fileroot,'none',1)
  case default
     call vfile_init(trim(fileroot)//'_out',iched,nplot)
     call geobjlist_writev(geobjl,'geometry',nplot)
     call vfile_iscalarwrite(bods,size(bods),'Body','CELL',nplot,1)
  end select divide_type
  call clock_stop(30)
!--------------------------------------------------------------------------
!! cleanup and closedown
  call geobjlist_delete(geobjl)

  call clock_stop(1)
  call clock_summary

  call vfile_close
  call log_close
  call clock_delete
!--------------------------------------------------------------------------

end program datvtkprog
