!> @addtogroup groupname0
!> @{
program move_p
!> @}
  use const_kind_m
  use const_numphys_h
  use control_h
  use date_time_m
  use control_m
  use log_m
  use clock_m
  use position_h
  use position_m
  use geobj_m
  use ls_m
  use btree_m
  use mtest_m
  use pcle_h
  use pcle_m
  use geobjlist_h
  use geobjlist_m
  use vcontrol_h
  use dcontrol_h
  use posang_h
  use li_m
  use ld_m
  use dbtree_h
  use dbtree_m
  use bods_h
  use termplane_h
  use vfile_m
  use outfile_m
  use hdsfile_m
  use misc_m
  use query_m
  use posang_m
  use datline_h
  use datline_m
  use stack_m
  use indict_m
  use fmesh_h
  use beq_h
  use spl2d_m
  use spl3d_m

  implicit none

! Local variables
  character(*), parameter :: m_name='move' !< module name
  type(files_t)     :: file      !< names of files
  type(mtnumerics_t)  :: numerics  !< numerical parameters for mtest, includes hdsgen parameters as %n
  type(plots_t)     :: plot      !< diagnostic plot selectors
  type(btree_t)  :: btree      !< binary tree
  type(mtest_t)  :: mtest      !< mtest, set of rays to trace
  type(geobjlist_t) :: geobjl !< geometrical objects
  type(geobjlist_t) :: geobjm !< geometrical objects for mtest
  type(date_time_t) :: timestamp !< timestamp of run
  type(date_time_t) :: timeold !< timestamp of preceding HDSGEN run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=80) :: icbuf !< local variable
  integer(ki4):: nplot !< unit for vfiles
  integer(ki4):: nouth !< unit for hdsfiles
  integer(ki4):: nmt !< unit for mtest

!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'move run time')

!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'move: test particle move using hybrid data structure'
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long

!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run move type at the command line:'
     print *, '   move fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)

!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'control_init time')
  call control_init(fileroot)
  call control_mread(file,numerics,plot)
  call clock_stop(2)

!--------------------------------------------------------------------------
!! read geobjm file to describe rays

  if (numerics%qfilesuf=='vtk') then
!geobjm%ngtype=0
     icbuf=' '
     call geobjlist_read(geobjm,file%qrydata,icbuf,nmt)
  end if

!--------------------------------------------------------------------------
!! mtest data

  call clock_start(3,'mtest_init time')
  call mtest_init(mtest,numerics)
  call mtest_readv(mtest,geobjm,file%qrydata,nmt)
  call mtest_set(mtest)
  call clock_stop(3)

!--------------------------------------------------------------------------
!! read  HDS  data

  call clock_start(4,'hds_init time')
  call hdsfile_dinit(file%hdsdata,timeold)
  call hdsfile_read(geobjl,numerics%n,btree)
  call clock_stop(4)

!--------------------------------------------------------------------------
!! read  geobjl data

  call clock_start(5,'geobjlist_init time')
  call geobjlist_init(geobjl,file%vtkdata,numerics%n)
  call clock_stop(5)

!--------------------------------------------------------------------------
!!evaluate particle moves on mtest set

  call clock_start(9,'moves evaluation time')
  call geobjlist_stepmove(geobjl,btree,mtest)
  call clock_stop(9)

!--------------------------------------------------------------------------
!! data diagnostics

!!plot HDS
  if(plot%hdsm) then
     call clock_start(10,'vfile_hds time')
     call vfile_init(file%hdsm,'hds in real space',nplot)
     call btree_writev(btree,numerics%n,'hds',nplot)
     call vfile_close
     call clock_stop(10)
  end if

!!plot moves
  if(plot%movq) then
     call clock_start(17,'vfile_mov time')
     call vfile_init(file%movq,'moves in quantised space',nplot)
! in quantised space
     call mtest_writev(mtest,numerics,'quantised space',nplot)
     call vfile_close
     call clock_stop(17)
  end if

!! intersection points output file

  if(plot%ptsq) then
     call clock_start(18,'ptsfile_init time')
     call vfile_init(file%ptsq,'intersection points in quantised space',nplot)
! in quantised space
     call mtest_writeptsv(mtest,numerics,'quantised space',nplot)
     call vfile_close
     call clock_stop(18)
  end if

!!plot all geobjs
  if(plot%geobjq) then
     call clock_start(16,'vfile_allgeobj time')
     call vfile_init(file%geobjq,'all geobj',nplot)
! in position space
     call position_invqtfmlis(geobjl%posl,geobjl%quantfm)
     call geobjlist_writev(geobjl,'all transformed',nplot)
     call position_qtfmlis(geobjl%posl,geobjl%quantfm)
     call vfile_close
     call clock_stop(16)
  end if

!--------------------------------------------------------------------------
!! output file

  call clock_start(20,'outfile_init time')
  call outfile_minit(file,timestamp)
  call outfile_mwrite(timeold)
  call outfile_close
  call clock_stop(20)

!--------------------------------------------------------------------------
!! mtest/move output file

  if(plot%mov) then
     call clock_start(21,'movfile_init time')
     call vfile_init(file%mov,'moves',nplot)
     call mtest_writev(mtest,numerics,'real space',nplot)
     call vfile_close
     call clock_stop(21)
  end if

!--------------------------------------------------------------------------
!! intersection points output file

  if(plot%pts) then
     call clock_start(22,'ptsfile_init time')
     call vfile_init(file%pts,'intersection points',nplot)
     call mtest_writeptsv(mtest,numerics,'real space',nplot)
     call vfile_close
     call clock_stop(22)
  end if

!--------------------------------------------------------------------------
!! cleanup and closedown
  call geobjlist_delete(geobjl)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete

!--------------------------------------------------------------------------

end program move_p_lib
