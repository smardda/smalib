program geoq_p

  use const_kind_m
  use const_numphys_h
  use date_time_m
  use control_h
  use dcontrol_h
  use log_m
  use clock_m
  use beq_h
  use beqan_h
  use posang_h
  use bcontrol_m
  use vcontrol_h
  use mcontrol_h
  use mcontrol_m
  use beq_m
  use beqan_m
  use posang_m

  use position_h
  use geobjlist_h
  use position_m
  use ls_m
  use btree_m
  use geobj_m
  use stack_m
  use query_m
  use datline_h
  use datline_m
  use geobjlist_m
  use geoq_m

  use spl2d_m
  use spl3d_m
  use boutfile_m
  use moutfile_m
  use vfile_m
  use gfile_m

  implicit none


! Local variables
  character(*), parameter :: m_name='geoq' !< module name
  type(bfiles_t)     :: file      !< names of files
  type(bnumerics_t)  :: numerics  !< numerical control parameters
  type(bplots_t)     :: plot      !< diagnostic plot selectors
  type(geoq_t)  :: geoq      !< equilibrium B + geometrical objects
  type(beqan_t)  :: beqan      !< analytic equilibrium B
  type(date_time_t) :: timestamp !< timestamp of run
  character(len=80) :: fileroot !< reference name for all files output by run
  character(len=80) :: fileq !< equil file name
  character(len=256) :: vtkdesc !< descriptor line for vtk files
  character(len=80) :: mapfld !< field output in geoqm file

  integer(ki4):: nplot !< unit for vtk files
  integer(ki4):: nprint !< unit for gnuplot files
  integer(ki4):: nin=0 !< unit for other data
  integer(ki4):: nana=0 !< unit for analytic field data
  integer(ki4):: ifldspec !< field specification
  integer(ki4):: ipsibig !< flag whether psi overlarge by 2pi
  real(kr8):: zivac !< value of I in field file
  character(len=80) :: ibuf !< character workspace
!--------------------------------------------------------------------------
!! initialise timing

  call date_time_init(timestamp)
  call clock_init(30)
  call clock_start(1,'geoq run time')
!--------------------------------------------------------------------------
!! print header

  print *, '----------------------------------------------------'
  print *, 'geoq: combined geometry and equilibrium processing '
  print *, '----------------------------------------------------'
  print '(a)', timestamp%long
!--------------------------------------------------------------------------
!! get file root from arg
  if(command_argument_count()<1) then
!! no file root specified
     print *, 'Fatal error: no file root name specified.'
     print *, 'To run geoq type at the command line:'
     print *, '   geoq fileroot'
     stop
  else
!!get fileroot
     call get_command_argument(1,value=fileroot)
  end if

!! start log
  call log_init(fileroot,timestamp)
!--------------------------------------------------------------------------
!! read control file

  call clock_start(2,'bcontrol_init time')
  call bcontrol_init(fileroot)
  call bcontrol_read(file,numerics,plot)
  call clock_stop(2)
!--------------------------------------------------------------------------
!! beq data

  call clock_start(3,'beq_init time')
  fileq=file%equil
  ifldspec=numerics%fldspec
  ipsibig=numerics%psibig
  select case (numerics%eqtype)
  case ('ana')
! need additional data to specify field analytically
     call beqan_init(file%equil,nana)
     call beqan_readcon(beqan)
     call beq_readana(geoq%beq,beqan,ifldspec)
     call beqan_initwrite(fileroot)
     call beqan_write(beqan)
     call beqan_closewrite()
     call beqan_delete(beqan)
  case ('equ')
     call beq_readequ(geoq%beq,fileq,ifldspec)
  case default
     call beq_readequil(geoq%beq,fileq,ifldspec,ipsibig)
  end select
  call beq_move(geoq%beq,numerics)
  call beq_init(geoq%beq,numerics)
  call clock_stop(3)
!--------------------------------------------------------------------------
!! read  geobjl data

  call clock_start(4,'geobjlist_init time')
  call geobjlist_read(geoq%objl,file%vtkdata)
  call clock_stop(4)
!--------------------------------------------------------------------------
!! do the main work

  call clock_start(6,'beqs evaluation time')

  select case (ifldspec)
  case (2)
     call beq_centre(geoq%beq)
     call geoq_init(geoq)
     call beq_bdryrb(geoq%beq)
     call beq_ipsi(geoq%beq)
     call beq_xilt(geoq%beq)
  case (3)
     call beq_centre(geoq%beq)
     call geoq_init(geoq)
     call beq_bdryrb(geoq%beq)
     call beq_ipsi(geoq%beq)
     call beq_rispld(geoq%beq)
     call beq_xilt(geoq%beq)
  case default
     call beq_centre(geoq%beq)
     call geoq_init(geoq)
     call beq_psilt(geoq%beq)
     call beq_thetalt(geoq%beq)
     call beq_rextrema(geoq%beq)
     call beq_rpsi(geoq%beq)
     call beq_rjpsi(geoq%beq)
     call beq_bdryrb(geoq%beq)
     call beq_ipsi(geoq%beq)
  end select
  call clock_stop(6)
!--------------------------------------------------------------------------
!! data diagnostics

  select case (ifldspec)
  case (1)
     mapfld='all'
     call beq_dia(geoq%beq)
  case default
     mapfld='frzxi'
! these options are not valid so are suppressed
     plot%geoqvolm=.FALSE.
     plot%gnuptz=.FALSE.
     plot%gnum=.FALSE.
     plot%gnusilm=.FALSE.
     plot%geoqx=.FALSE.
! this one not implemented (or needed?)
     plot%frzv=.FALSE.
  end select
!! optionally read in ripple data
  if (geoq%beq%n%vacfile/='null') then
     if(plot%geoqx.OR.plot%geoqvolx.OR. &
 &   plot%geoqm.OR.plot%geofldx.OR.plot%frzzeta.OR.plot%geoqvolm) then
        call moutfile_read(geoq%beq%n%vacfile,zivac,nin)
        call spl3d_read(geoq%beq%vacfld,geoq%beq%n%vacfile,nin)
     end if
  end if
!--------------------------------------------------------------------------
!!plots
!!plot cartesian B
  if(plot%geoqx) then
     call clock_start(10,'vfile_geoqx time')
     call vfile_init(file%geoqx,'B in cartesian space',nplot)
     call beq_writev(geoq%beq,'cartesian',nplot)
     call vfile_close
     call clock_stop(10)
  end if

!!plot cartesian all B
  if(plot%geoqvolx) then
     call clock_start(11,'vfile_geoqvolx time')
     call vfile_init(file%geoqvolx,'B in all cartesian space',nplot)
     call beq_writev(geoq%beq,'allcartesian',nplot)
     call vfile_close
     call clock_stop(11)
  end if

!!plots for psi-theta-zeta space
  if(plot%geoqvolm) then
     call clock_start(16,'vfile_geoqvolm time')
     call vfile_init(file%geoqvolm,'fields in psi-theta space',nplot)
     call beq_writev(geoq%beq,'psi-theta',nplot)
     call vfile_close
     call clock_stop(16)
  end if

!!plots for R-Z space
  if(plot%frzv) then
     call clock_start(17,'vfile_frzv time')
     call vfile_init(file%frzv,'fields in R-Z space',nplot)
     call beq_writev(geoq%beq,'R-Z',nplot)
     call vfile_close
     call clock_stop(17)
  end if

!! gnu plots
  if(plot%gnu) then
     call clock_start(18,'gfile_gnu time')
     call gfile_init(trim(file%gnu),'psi sample in R-Z space',nprint)
     call spl2d_writeg(geoq%beq%psi,'sampl',nprint)
     call gfile_close
     call gfile_init(trim(file%gnu)//'_dR','dpsi/dR sample in R-Z space',nprint)
     call spl2d_writeg(geoq%beq%dpsidr,'sampl',nprint)
     call gfile_close
     call gfile_init(trim(file%gnu)//'_dZ','dpsi/dZ sample in R-Z space',nprint)
     call spl2d_writeg(geoq%beq%dpsidz,'sampl',nprint)
     call gfile_close
     call clock_stop(18)
  end if

  if(plot%gnum) then
     call clock_start(19,'gfile_gnum time')
     call gfile_init(trim(file%gnu)//'_R','R sample in psi-theta space',nprint)
     call spl2d_writeg(geoq%beq%r,'sampl',nprint)
     call gfile_close
     call gfile_init(trim(file%gnu)//'_Z','Z sample in psi-theta space',nprint)
     call spl2d_writeg(geoq%beq%z,'sampl',nprint)
     call gfile_close
     call gfile_init(trim(file%gnu)//'_RJ','R/J sample in psi-theta space',nprint)
     call spl2d_writeg(geoq%beq%rjac,'sampl',nprint)
     call gfile_close
     call gfile_init(trim(file%gnu)//'_RZ','RZ contours',nprint)
     call beq_writeg(geoq%beq,'R-Z',nprint)
     call gfile_close
     call clock_stop(19)
  end if

!!plot all geobjs
  if(plot%geoqm) then
     call clock_start(20,'vfile_geoqm time')
     write(vtkdesc,'(A3,1X,A18,I3)') mapfld(1:3),'Integer_Parameter=',geoq%objl%nparam
     call vfile_init(file%geoqm,vtkdesc,nplot)
! in mapped position space
     call geoq_writev(geoq,mapfld,nplot)
     call vfile_close
     call clock_stop(20)
  end if

  if(plot%geofldx) then
     call clock_start(22,'vfile_geofldx time')
     call vfile_init(file%geofldx,'allcart',nplot)
     ibuf=geoq%beq%n%vacfile
!BB! in Cartesians, omit ripple contribution
!BB         geoq%beq%n%vacfile='null'
     call geoq_writev(geoq,'allcart',nplot)
     geoq%beq%n%vacfile=ibuf
     call vfile_close
     call clock_stop(22)
  end if

  if(plot%frzzeta) then
     call clock_start(23,'vfile_frzzeta time')
     call vfile_init(file%frzzeta,'frzzeta',nplot)
     call geoq_writev(geoq,'frzzeta',nplot)
     call vfile_close
     call clock_stop(23)
  end if

  !! silhouette of geometry in (R,Z)
  if(plot%gnusil) then
     call clock_start(25,'vfile_gnusil time')
     call gfile_init(file%gnusil,'gnusil',nprint)
     call geoq_writeg(geoq,'gnusil',nprint)
     call gfile_close
     call clock_stop(25)
  end if

  !! silhouette of geometry in (psi,theta)
  if(plot%gnusilm) then
     call clock_start(26,'vfile_gnusilm time')
     call gfile_init(file%gnusilm,'gnusilm',nprint)
     call geoq_writeg(geoq,'gnusilm',nprint)
     call gfile_close
     call clock_stop(26)
  end if
!--------------------------------------------------------------------------
!! output file

  call clock_start(30,'outfile_init time')
  call boutfile_init(file,timestamp)
  call boutfile_write(geoq%beq,timestamp)
  call boutfile_close
  call clock_stop(30)
!--------------------------------------------------------------------------
!! cleanup and closedown
  call geoq_delete(geoq)

  call clock_stop(1)
  call clock_summary

  call log_close
  call clock_delete
!--------------------------------------------------------------------------

end program geoq_p
