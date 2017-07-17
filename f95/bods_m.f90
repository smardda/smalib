module bods_m

  use const_kind_m
  use log_m
  use const_numphys_h
  use position_h
  use geobjlist_h
  use position_m
  use control_h
  use vcontrol_h
  use posang_h
  use ls_m
  use btree_m
  use geobj_m
  use query_m
  use posang_m
  use datline_h
  use datline_m
  use vfile_m
  use geobjlist_m

  implicit none
  private

! public subroutines
  public :: &
 &bods_init, & !< initialise dat vtk file bodies
 &bods_getlist, & !< get list of different bodies in bods array
 &bods_write, & !< write out bodies in separate files
 &bods_cumulate !< cumulate body information

! public types
!integer(ki4), parameter, public :: bods_made_up_data=1

! private variables
  character(*), parameter :: m_name='bods_m' !< module name
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: status  !< status flag
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> initialise dat vtk file bodies
subroutine bods_init(self,geobjl,blockno)

  !! arguments
  integer(ki4), dimension(:), allocatable, intent(out) :: self !< bodies
  type(geobjlist_t), intent(in) :: geobjl   !< geobj list data
  integer(ki4), intent(in) :: blockno !< first block number

  !! local
  character(*), parameter :: s_name='bods_init' !< subroutine name
  integer(ki4)  :: iobod !< number of objects in a body

  !! allocate array
  iobod=geobjl%ng
  if(iobod>0) then
     allocate(self(iobod),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No bodies data')
  end if

  !! copy from geobjlist and integer convert
  if (allocated(geobjl%obj)) then
     do j=1,iobod
        self(j)=geobjl%obj(j)%weight+0.5
     end do
  else
     do j=1,iobod
        self(j)=blockno
     end do
  end if

end subroutine bods_init
!---------------------------------------------------------------------
!> get list of different bodies in bods array
subroutine bods_getlist(self,selfo,kno)

  !! arguments
  integer(ki4), dimension(:), intent(in) :: self !< bodies
  integer(ki4), dimension(:), allocatable, intent(out) :: selfo !< bodies short list
  integer(ki4), intent(out) :: kno !< bodies short list dimension

  !! local
  character(*), parameter :: s_name='bods_getlist' !< subroutine name
  integer(ki4)  :: inbod !< number of entries in self
  integer(ki4), dimension(:), allocatable :: iwork !< scratch bodies list
  integer(ki4)  :: inew !< new entry in scratch list

  !! allocate array
  inbod=ubound(self,1)
  if(inbod>0) then
     allocate(iwork(inbod),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No bodies data')
  end if

  !! find different entries
  kno=1
  iwork(kno)=self(1)
  do j=2,inbod
     inew=1
     do k=1,kno
        if (self(j)==iwork(k)) then
           inew=0
           exit
        end if
     end do
     if (inew==1) then
        kno=kno+1
        iwork(kno)=self(j)
     end if
  end do
  !! copy to smaller output array
  allocate(selfo(kno),stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  selfo=iwork(1:kno)
  
  deallocate(iwork)

end subroutine bods_getlist
!---------------------------------------------------------------------
!> write out bodies in separate files
subroutine bods_write(self,knbod,geobjl,fileroot,kctyp,kcname,kheader)
  !! arguments
  integer(ki4), dimension(knbod), intent(in) :: self !< bodies list
  integer(ki4), intent(in) :: knbod   !< size of bodies list data
  type(geobjlist_t), intent(in) :: geobjl   !< geobj list data
  character(*),intent(in) :: fileroot !< file root
  character(*),intent(in) :: kctyp !< type of points compression (none)
  character(*),intent(in) :: kcname !< name of cell attribute
  integer(ki4), intent(in) :: kheader   !< header options

  !! local
  character(*), parameter :: s_name='bods_write' !< subroutine name
  type(geobjlist_t) :: glsmall   !< geobj list data small
  integer(ki4) :: infile   !< number of different bodies files
  integer(ki4)  :: inpt !< number of points in generic geobjlist
  integer(ki4)  :: inobj !< number of objects in generic geobjlist
  integer(ki4)  :: innod !< number of nodal data in generic geobjlist
  integer(ki4)  :: inptsmall !< number of points in small geobjlist
  integer(ki4)  :: inobjsmall !< number of objects in small geobjlist
  integer(ki4)  :: innodsmall !< number of nodal data in small geobjlist
  integer(ki4)  :: iostart!< end points of object list if bodies' objects contiguous
  integer(ki4)  :: ioend !< end points of object list if bodies' objects contiguous
  integer(ki4) :: iplot   !< output channel for bodies list data
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: iobj !< local variable
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: ityp !< local variable
  character(len=5) :: ic5  !< local variable
  character(len=80) :: icplot  !< local variable
  character(len=80) :: iched  !< local variable

  infile=max(1,maxval(self))

  ! write out all points
  inptsmall=geobjl%np

  do i=1,infile

     !! make up icplot from fileroot, iched from fileroot
     write(ic5,'(I5.5)') i
     icplot=trim(fileroot)//'_'//ic5
     iched='original filename root '//fileroot
     call vfile_init(icplot,iched,iplot)

     !! construct small geobjl, first allocate storage

     inobjsmall=0
     innodsmall=0
     do j=1,geobjl%ng

        if (self(j)==i) then
           ioend=j
           inobjsmall=inobjsmall+1
           !     iobj=geobjl%obj2(j)%ptr
           ityp=geobjl%obj2(j)%typ
           inumpts=geobj_entry_table(ityp)
           innodsmall=innodsmall+inumpts
        end if
     end do
     iostart=ioend+1-inobjsmall

     call geobjlist_iinit(glsmall,inptsmall,inobjsmall,innodsmall,2,1)

     !! fill small geobjl with data

     !! copy all points
     glsmall%posl=geobjl%posl

     inobj=0
     innd=1
     do j=1,geobjl%ng

        if (self(j)==i) then
           inobj=inobj+1
           iobj=geobjl%obj2(j)%ptr
           ityp=geobjl%obj2(j)%typ
           inumpts=geobj_entry_table(ityp)
           do k=1,inumpts
              glsmall%nodl(innd-1+k)=geobjl%nodl(iobj-1+k)
           end do
           glsmall%obj2(inobj)%ptr=innd
           glsmall%obj2(inobj)%typ=ityp
           innd=innd+inumpts
        end if
     end do

     call geobjlist_ptcompress(glsmall)

     call geobjlist_writev(glsmall,'geometry',iplot)

     call vfile_iscalarwrite(self(iostart:ioend),inobjsmall,kcname,'CELL',iplot,kheader)

     call vfile_close

     call geobjlist_delete(glsmall)

  end do

end subroutine bods_write

!---------------------------------------------------------------------
!> cumulate body information
subroutine bods_cumulate(self,geobjl,nfile,start,copy,kopt)

  !! arguments
  integer(ki4), dimension(:), allocatable, intent(inout) :: self !< local variable
  type(geobjlist_t), intent(in) :: geobjl   !< geobj list data
  integer(ki4), intent(in) :: nfile   !< number of file
  integer(ki4), intent(in) :: start !< start number of copies required
  integer(ki4), intent(in) :: copy !< stop number of copies required
  integer(ki4), intent(in) :: kopt   !< if nonzero, try to copy weights

  !! local
  character(*), parameter :: s_name='bods_cumulate' !< subroutine name
  integer(ki4)  :: icopy !< actual number of copies required
  integer(ki4)  :: ibod !< number of objects in self
  integer(ki4)  :: inbod !< number of objects in geobjl
  integer(ki4), dimension(:), allocatable :: iwork !< integer work array
  integer(ki4) :: ing !< number of objects in geobjlist

  if (copy<=0) return
  icopy=copy+1-start
  !! allocate array
  ibod=ubound(self,1)
  ing=geobjl%ng
  if(ing<0) &
 &call log_error(m_name,s_name,2,error_fatal,'No bodies data')

  inbod=ibod+ing*icopy
  !set up replacement bodies array and destroy old bodies
  allocate(iwork(inbod), stat=status)
  call log_alloc_check(m_name,s_name,6,status)
  do j=1,ibod
     iwork(j)=self(j)
  end do
  i=ibod+1
  !! copy from geobjlist and integer convert
  do k=start,copy
     do j=1,ing
        if (allocated(geobjl%obj)) then
           iwork(i)=100*(100*nfile+k)+geobjl%obj(j)%weight+0.5
        else
           iwork(i)=100*nfile+k
        end if
        i=i+1
     end do
  end do
  deallocate(self)
  allocate(self(inbod), stat=status)
  call log_alloc_check(m_name,s_name,7,status)
  self=iwork
  deallocate(iwork)

end subroutine bods_cumulate

end module bods_m
