program vtktfmprog

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use control_h
  use vcontrol_h
  use log_m
  use clock_m
  use posang_h
  use posang_m
  use spl2d_m

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
  use stack_m

  use bods_m

  implicit none


! Local variables
  character(*), parameter :: m_name='vtktfmprog' !< module name
  type(vfiles_t)     :: file      !< names of files
  type(vnumerics_t)  :: numerics  !< numerical control parameters
  type(geobjlist_t)  :: geobjl      !< geometrical objects
  type(geobjlist_t)  :: igeobjl      !< more geometrical objects
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=80),save :: iched !< vtk field file descriptor

  integer(ki4), dimension(:), allocatable :: bods !< array of bodies for geometrical objects
  integer(ki4), dimension(:), allocatable :: ibods !< array of more bodies for geometrical objects
  integer(ki4):: nplot !< unit for vtk files
  integer(ki4):: nin !< unit for other data
  integer(ki4):: cpstart !< start number of copies
  integer(ki4):: nscal !< number of scalars (body identifiers)
  integer(ki4):: iopt=0 !< option for bodies
  integer(ki4):: j !< loop variable
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'vtktfmprog run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'vtktfmprog: transform vtk geometry file'
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run vtktfmprog type at the command line:'
     print *, '   vtktfm fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)
!--------------------------------------------------------------------------
!! read control file (nominally)

  call clock_start(2,'vcontrol_init time')
  call vcontrol_init(fileroot)
  call vcontrol_read(file,numerics)
  call clock_stop(2)
!--------------------------------------------------------------------------
!! read  geobjl data

  call clock_start(4,'geobjlist_read time')
  if (file%nvtkdata==1) then
! just one file with Body data
     call geobjlist_read(geobjl,file%vtkdata(1),iched)
!     write(*,*) 'first',(geobjl%nodl(j),j=1,20)
     nin=0
     call vfile_iscalarread(bods,nscal,file%vtkdata(1),'Body',nin,1)
  else
     call geobjlist_read(geobjl,file%vtkdata(1),iched)
     call bods_init(bods,geobjl,101)
     call geobjlist_close()
     cpstart=2
     do j=1,file%nvtkdata
        call geobjlist_read(igeobjl,file%vtkdata(j),iched)
!     write(*,*) 'first',(geobjl%nodl(j),j=1,20)
        call geobjlist_cumulate(geobjl,igeobjl,cpstart,file%vtkcopies(j),iopt)
        call bods_cumulate(bods,igeobjl,j,cpstart,file%vtkcopies(j),iopt)
        cpstart=1
        call geobjlist_delete(igeobjl)
        call geobjlist_close()
     end do
  end if
!! check nscal vs objl
  call clock_stop(4)
!--------------------------------------------------------------------------
!! do the main work

  call clock_start(6,'position transform time')
!     line below needed for debugging version to run
  write(*,*) 'mysterious output',(geobjl%nodl(j),j=1,2)
!     write(*,*) 'second',((numerics%panbod(i,j),i=1,2),j=1,9)
  call geobjlist_paneltfm(geobjl,bods,numerics)
  call clock_stop(6)
!--------------------------------------------------------------------------
!! output file

  call clock_start(30,'outfile_init time')
  call vfile_init(file%vtkout,iched,nplot)
!     write(*,*) 'iched=',iched
  call geobjlist_writev(geobjl,'geometry',nplot)
  nscal=ubound(bods,1)
  call vfile_iscalarwrite(bods,nscal,'Body','CELL',nplot,1)
  call clock_stop(30)
!--------------------------------------------------------------------------
!! cleanup and closedown
  call geobjlist_delete(geobjl)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete
!--------------------------------------------------------------------------

end program vtktfmprog
