program hdsgen

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
  use posang_h
  use ls_m
  use btree_m
  use mtest_m
  use pcle_h
  use pcle_m
  use geobjlist_h
  use geobjlist_m
  use vfile_m
  use outfile_m
  use hdsfile_m
  use query_m
  use posang_m
  use datline_h
  use datline_m
  use stack_m
  use spl2d_m

  implicit none


! Local variables
  character(*), parameter :: m_name='hdsgen' !< module name
  type(files_t)     :: file      !< names of files
  type(numerics_t)  :: numerics  !< numerical parameters
  type(plots_t)     :: plot      !< diagnostic plot selectors
  type(btree_t)  :: btree      !< binary tree
  type(geobjlist_t) :: geobjl !< geometrical objects
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run

  integer(ki4):: nplot !< unit for vfiles
  integer(ki4):: nouth !< unit for hdsfiles

!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(20)
  call clock_start(1,'hdsgen run time')

!--------------------------------------------------------------------------
!! print header

  print *, '--------------------------------------------------'
  print *, 'hdsgen: map data to hybrid data structure'
  print *, '--------------------------------------------------'
  print '(a)', timestamp%long

!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run hdsgen type at the command line:'
     print *, '   hdsgen fileroot'
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
  call control_read(file,numerics,plot)
  call clock_stop(2)

!--------------------------------------------------------------------------
!! initialise  geobjl data

  call clock_start(3,'geobjlist_init time')
  call geobjlist_init(geobjl,file%vtkdata,numerics)
  call clock_stop(3)

!--------------------------------------------------------------------------
!! initialise btree

  call clock_start(4,'btreelist_init time')
  call btree_init(btree,numerics)
  call clock_stop(4)

!--------------------------------------------------------------------------
!! initial btree data diagnostics

! NOT NEEDED

!--------------------------------------------------------------------------
!!sort geobjl onto bins

  call clock_start(9,'geobjlist_sort time')
  call geobjlist_step(geobjl,btree)
  numerics%geobj_coord_tfm=geobjl%quantfm
  call clock_stop(9)

!--------------------------------------------------------------------------
!! sorted btree and geobjl data diagnostics


!!plot HDS
  if(plot%hds) then
     call clock_start(10,'vfile_hds time')
     call vfile_init(file%hdsv,'hds lowest',nplot)
     call btree_writev(btree,numerics,'hds lowest',nplot)
     call vfile_close
     call clock_stop(10)
  end if

!!plot HDS quantised (default)
  if(plot%hdsbin) then
     call clock_start(11,'vfile_hdsq time')
     call vfile_init(file%hdsq,'hds quantised',nplot)
     call btree_writev(btree,numerics,'hds quantised',nplot)
     call vfile_close
     call clock_stop(11)
  end if

!!plot assigned geobj
  if(plot%geobj) then
     call clock_start(14,'vfile_assignedgeobj time')
     call vfile_init(file%geobj,'assigned geobj',nplot)
! in quantised space
     call geobjlist_writev(geobjl,'assigned quantised',nplot)
     call vfile_close
     call clock_stop(14)
  end if

!!plot unassigned geobjs
  if(plot%lostgeobj) then
     call clock_start(15,'vfile_unassignedgeobj time')
     call vfile_init(file%lostgeobj,'unassigned geobj',nplot)
! in quantised space
     call geobjlist_writev(geobjl,'unassigned quantised',nplot)
     call vfile_close
     call clock_stop(15)
  end if

!!plot all geobjs
  if(plot%allgeobjq) then
     call clock_start(16,'vfile_allgeobjq time')
     call vfile_init(file%allgeobjq,'all geobjq',nplot)
! in quantised space
     call geobjlist_writev(geobjl,'all quantised',nplot)
     call vfile_close
     call clock_stop(16)
  end if

!!plot density geobjs
  if(plot%densitygeobj) then
     call clock_start(17,'vfile_densitygeobj time')
     call vfile_init(file%densitygeobj,'density geobj',nplot)
! in quantised space
     call geobjlist_writev(geobjl,'density quantised',nplot)
     call vfile_close
     call clock_stop(17)
  end if

!--------------------------------------------------------------------------
!! output file

  call clock_start(20,'outfile_init time')
  call outfile_init(file,timestamp)
  call outfile_write(geobjl,btree)
  call outfile_close
  call clock_stop(20)

!--------------------------------------------------------------------------
!! hds output file

  call clock_start(21,'hdsfile_init time')
  call hdsfile_init(file%hdslist,timestamp,nouth)
  call hdsfile_write(geobjl,numerics,btree)
  call hdsfile_close
  call clock_stop(21)

!--------------------------------------------------------------------------
!! cleanup and closedown
  call geobjlist_delete(geobjl)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete

!--------------------------------------------------------------------------

end program hdsgen
