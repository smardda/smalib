program vtktfm_p

#ifdef WITH_MPI
  use mpi
#endif
  use const_kind_m
  use const_numphys_h
  use date_time_m
  use control_h
  use vcontrol_h
  use dcontrol_h
  use log_m
  use clock_m
  use position_h
  use fmesh_h
  use beq_h
  use posang_h
  use posang_m
  use spl2d_m
  use spl3d_m

  use geobjlist_h
  use position_m
  use ls_m
  use li_m
  use ld_m
  use dbtree_h
  use dbtree_m
  use btree_m
  use geobj_m
  use query_m
  use geobjlist_m
  use indict_m

  use vcontrol_m
  use vfile_m
  use gfile_m

  use datline_h
  use datline_m
  use stack_m

  use bods_h
  use bods_m
  use scontrol_h
  use scontrol_m

  implicit none


! Local variables
  character(*), parameter :: m_name='vtktfm' !< module name
  type(vfiles_t)     :: file      !< names of files
  type(vnumerics_t)  :: numerics  !< numerical control parameters
  type(geobjlist_t)  :: geobjl      !< geometrical objects
  type(geobjlist_t)  :: igeobjl      !< more geometrical objects
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=80),save :: iched !< vtk field file descriptor

  type(bods_t) :: bods !< type bods for geometrical objects
  integer(ki4), dimension(:), allocatable :: ibods !< array of more bodies for geometrical objects
  integer(ki4):: nplot !< unit for vtk files
  integer(ki4):: nin !< unit for other data
  integer(ki4):: cpstart !< start number of copies
  integer(ki4):: nscal !< number of scalars (body identifiers)
  integer(ki4):: iopt=1 !< option for bodies
  integer(ki4):: icall !< call number for bods_cumulate2
  integer(ki4):: j !< loop variable
#ifdef WITH_MPI
  integer error, rank
  call MPI_Init ( error )
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, error)
  if (rank .eq. 0) then
#endif
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'vtktfm run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'vtktfm: transform vtk geometry file'
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run vtktfm type at the command line:'
     print *, '   vtktfm fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)
!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'vcontrol_init time')
  call vcontrol_init(fileroot)
  call vcontrol_read(file,numerics)
  call clock_stop(2)
!--------------------------------------------------------------------------
!! read  geobjl data

  call clock_start(4,'geobjlist_read time')
  geobjl%ngtype=2
! write(*,*) "fn,iopt=",file%nvtkdata, iopt
  if (file%nvtkdata==1) then
! just one file with Body data
     call geobjlist_read(geobjl,file%vtkdata(1),iched)
!     write(*,*) 'first',(geobjl%nodl(j),j=1,20)
     nin=0
     call vfile_iscalarread(bods%list,nscal,file%vtkdata(1),numerics%name,nin,iopt) !W
     !write(*,*) "fn,iopt=",file%nvtkdata, iopt
     if (iopt==0) then
        numerics%npans=maxval(bods%list)
! read data OK, but may still suppress
        if (numerics%same) then
           call bods_initlist(bods,geobjl,numerics%nvalue) !W
        end if
     else if (4<=iopt.AND.iopt<=9) then
! no body data in file, fix up
        call log_value("Fixing up file for missing body/cell labels replacement ",numerics%nvalue)
        call bods_initlist(bods,geobjl,numerics%nvalue)
     else
        call log_error(m_name,m_name,iopt,error_fatal,'Corrupt vtk file')
     end if
     call bods_init(bods,geobjl,numerics) !W
  else if (numerics%preserve) then
     call geobjlist_read(geobjl,file%vtkdata(1),iched)
     nin=0
     call vfile_iscalarread(bods%list,bods%nbod,file%vtkdata(1),numerics%name,nin,iopt) !W
     call bods_init(bods,geobjl,numerics) !W
     call geobjlist_close()
     cpstart=2
     do j=1,file%nvtkdata
        icall=j-1
        igeobjl%ngtype=2
        call geobjlist_read(igeobjl,file%vtkdata(j),iched)
        nin=0
        call vfile_iscalarread(ibods,nscal,file%vtkdata(j),numerics%name,nin,iopt) !W
!     write(*,*) 'first',(geobjl%nodl(j),j=1,20)
        call geobjlist_cumulate(geobjl,igeobjl,cpstart,file%vtkcopies(j),iopt)
        call bods_cumulate2(bods,ibods,igeobjl,j,cpstart,file%vtkcopies(j),icall)
        cpstart=1
        call geobjlist_delete(igeobjl)
        call geobjlist_close()
     end do
  else
     call geobjlist_read(geobjl,file%vtkdata(1),iched)
     call bods_initlist(bods,geobjl,101)
     call bods_init(bods,geobjl,numerics) !W
     call geobjlist_close()
     cpstart=2
     do j=1,file%nvtkdata
        igeobjl%ngtype=2
        call geobjlist_read(igeobjl,file%vtkdata(j),iched)
!     write(*,*) 'first',(geobjl%nodl(j),j=1,20)
        call geobjlist_cumulate(geobjl,igeobjl,cpstart,file%vtkcopies(j),iopt)
        call bods_cumulate(bods,igeobjl,j,cpstart,file%vtkcopies(j),0)
        cpstart=1
        call geobjlist_delete(igeobjl)
        call geobjlist_close()
     end do
  end if
  call clock_stop(4)
!--------------------------------------------------------------------------
!! do the main work

  call clock_start(6,'position transform time')
!     line below needed for debugging version to run
!  write(*,*) 'this output helps gfortran sometimes',(geobjl%nodl(j),j=1,2)
!     write(*,*) 'second',((numerics%panbod(i,j),i=1,2),j=1,9)
  if (numerics%paneltfm) then
     call geobjlist_paneltfm(geobjl,bods,numerics)
  end if
  if (numerics%extract) then
     call geobjlist_extract(geobjl,bods%list,numerics)
  end if
  call clock_stop(6)
!--------------------------------------------------------------------------
!! output file(s)

  call clock_start(30,'outfile_init time')
  if (numerics%split) then
! write out as separate files
     call bods_write(bods,geobjl,fileroot,'none',numerics%name,1)
  else
     call vfile_init(file%vtkout,iched,nplot)
!     write(*,*) 'iched=',iched
     call geobjlist_writev(geobjl,'geometry',nplot)
     call vfile_iscalarwrite(bods%list,bods%nbod,numerics%name,'CELL',nplot,1)
     call vfile_close
  end if
  call clock_stop(30)
!--------------------------------------------------------------------------
!! cleanup and closedown
  call geobjlist_delete(geobjl)
  call bods_delete(bods)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete
!--------------------------------------------------------------------------
#ifdef WITH_MPI
  end if
  call MPI_Finalize ( error )
#endif
end program vtktfm_p
