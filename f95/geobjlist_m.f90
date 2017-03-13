module geobjlist_m

  use const_kind_m
  use log_m
  use const_numphys_h
  use position_h
  use geobjlist_h
  use position_m
  use control_h
  use vcontrol_h
  use dcontrol_h
  use posang_h
  use ls_m
  use btree_m
  use geobj_m
  use query_m
  use posang_m
  use datline_h
  use datline_m
  use stack_m

  implicit none
  private

! public subroutines
  public :: &
  geobjlist_init,   & !< read in geobjlist data structure
  geobjlist_delete,   & !< delete  geobjlist data structure
  geobjlist_close,   & !< close  geobjlist data file
  geobjlist_read,   & !< read (vtk)  geobjlist data structure
  geobjlist_addcube,   & !< add cube(s) to  geobjlist data structure
  geobjlist_writev,   & !< write (vtk)  geobjlist data structure
  geobjlist_nodlmv,   & !< rearrange (vtk)  geobjlist data structure
  geobjlist_step, &   !< control processing of  geobjlist data structure
  geobjlist_stepquery, &   !< density processing of geobjlist data structure
  geobjlist_getbb,   & !< bb of  geobjlist data structure
  geobjlist_tfm,   & !< transform geobjlist data structure
  geobjlist_spectfm,   & !< special transform of geobjlist data structure
  geobjlist_qtfm,   & !< quantise geobjlist data structure
  geobjlist_paneltfm,   & !< transform positions on a panel-by-panel basis
  geobjlist_bin,    & !< bin geobjlist data structure
  geobjlist_mbin,   & !< bin many-one geobjlist data structure
  geobjlist_dread,  & !< read dat file
  geobjlist_iinit,  & !< initialise geobjlist from internal data
  geobjlist_ptcompress,  & !< compress points
  geobjlist_orientri,  & !< orient triangles consistently
  geobjlist_shelltets, & !< shell set of tetrahedra
  geobjlist_copy, &  !< copy geobjlist to another
  geobjlist_cumulate, &  !< append geobjlist to first
  geobjlist_create, & !< create some 2-D geobjlists
  geobjlist_create3d !< create geobjlists by translation/rotation

! public types

!public variables

  integer(ki4), public :: geobjlist_maxchildn !< local variable
! private types

! private variables
  character(*), parameter :: m_name='geobjlist_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  character(len=132) :: bigbuf !<big buffer for input/output
  integer   :: status   !< error status
  integer(ki4), save :: nin   !< input channel for geobj data
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  integer(ki4) :: jj !< loop counter
  type(posveclis_t) :: rposl   !< list of position data
  real(kr4), dimension(3)  :: xlbb !< lower \f$ x \f$ of geobj bounding box
  real(kr4), dimension(3)  :: xubb !< upper \f$ x \f$ of geobj bounding box
  real(kr4), dimension(3)  :: origin !< origin of current hoc block
  real(kr4), dimension(3) :: zwork !< work vector
  integer(ki4) :: inpos !< number of position vectors to be read
  integer(ki4) :: inobj !< number of geobj records to be read
  integer(ki4) :: iempt !< first empty entry lin list
  integer(ki4) :: inls !< number of entries in list
  integer(ki4) :: idum !< dummy integer
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> read in geobjlist data structure
subroutine geobjlist_init(self,vtkfile,numerics)

  !! arguments
  type(geobjlist_t), intent(out) :: self   !< geobj list data
  character(*), intent(in)       :: vtkfile !< name of input file
  type(numerics_t), intent(inout) :: numerics !< numerical control data


  !! local
  character(*), parameter :: s_name='geobjlist_init' !< subroutine name

  self%ngtype=numerics%geomtype
  self%tolerance=max(numerics%maxtolerance,0.0)
  self%nbdcub=numerics%nbdcub
  self%cornflag=numerics%cornflag
  self%lowcorner=numerics%lowcorner
  self%upcorner=numerics%upcorner
  self%dilen=numerics%dilen
  self%dolen=numerics%dolen
  self%minmaxtolerance=max(numerics%mintolerance,0.0)
  self%minobj=numerics%mingeobjinbin
  self%nquant=numerics%nquante
  self%tfmdata=numerics%position_coord_tfm
  self%quantfm=numerics%geobj_coord_tfm
  self%ngunassigned=0
  self%nwset=0

  !! read coords
  call geobjlist_read(self,vtkfile)
  numerics%ngeobj=self%ng
  call log_error(m_name,s_name,1,log_info,'geobj data read')

  ! allocate array
  status=0
  if (.not.allocated(rposl%pos)) allocate(rposl%pos(1), stat=status)
  call log_alloc_check(m_name,s_name,2,status)

end subroutine geobjlist_init
!---------------------------------------------------------------------
!> delete geobjlist_t
subroutine geobjlist_delete(self)

  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list coord data
  !! local

  deallocate(self%posl%pos)
  if (self%ngtype==1.OR.self%nwset/=0) then
     !  write(*,*) self%ngtype,self%nwset
     !    if (allocated(self%obj)) deallocate(self%obj)
     deallocate(self%obj)
  end if
  if (self%ngtype/=1) then
     deallocate(self%nodl)
     deallocate(self%obj2)
  end if

end subroutine geobjlist_delete
!---------------------------------------------------------------------
!> close geobjlist file unit
subroutine geobjlist_close(self)

  !! arguments
  type(geobjlist_t), intent(inout), optional :: self !< geobj list coord data
  !! local
  character(*), parameter :: s_name='geobjlist_close' !< subroutine name

  !! close file
  close(unit=nin,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close file unit, ",i5)',nin
     call log_error(m_name,s_name,1,error_fatal,'Cannot close data file')
     stop
  end if

end subroutine geobjlist_close
!---------------------------------------------------------------------
!> read geobj coordinates
subroutine geobjlist_read(self,infile,kched,kin)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  character(*),intent(in) :: infile !< name of input file
  character(len=80),intent(out), optional :: kched !< field file header
  integer(ki4), intent(inout), optional :: kin   !< input channel for object data structure

  !! local
  character(*), parameter :: s_name='geobjlist_read' !< subroutine name
  logical :: unitused !< flag to test unit is available
  logical :: ilnumb !< local variable
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: inobj !< local variable
  integer(ki4) :: insto !< local variable
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: ityp !< local variable
  integer(ki4) :: ivtktyp !< local variable
  integer(ki4) :: ir   !< counter for file header read
  integer(ki4) :: k !< loop counter
  integer(ki4) :: ieq !< position of equals in string
  integer(ki4) :: iieq !< position of another or same equals in string
  integer(ki4) :: isubstr !< start of substring in string
  integer(ki4), dimension(geobj_max_entry_table) :: inodl  !< local variable
  integer(ki4) :: ifmt   !< formatting of position vectors

  logical :: isnumb !< local variable
  external isnumb

  if(present(kin).AND.kin>0) then
     !! assume unit already open and reading infile
     nin=kin
     rewind(nin)
  else
     !! get file unit
     do i=99,1,-1
        inquire(i,opened=unitused)
        if(.not.unitused)then
           nin=i
           exit
        end if
     end do
     if (present(kin)) kin=nin

     !! open file
     open(unit=nin,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0)then
        !! error opening file
        write(*,*) 'infile=',infile
        call log_error(m_name,s_name,1,error_fatal,'Error opening geobj data file')
     else
        call log_error(m_name,s_name,2,log_info,'geometrical object  data file opened')
     end if
  end if

  !! first set of reads to determine format of position vectors
  ir=0
  do
     ir=ir+1
     read(nin,fmt='(a)',iostat=status) ibuf1
     call log_read_check(m_name,s_name,70,status)
     !! look for key at start of line
     ibuf2=adjustl(ibuf1)
     if(ibuf2(1:6)=='POINTS') exit
  end do
  if (ir>20) then
     call log_error(m_name,s_name,71,error_fatal,'Error reading header data')
  end if
  read(nin,fmt='(a)',iostat=status) bigbuf
  call log_read_check(m_name,s_name,72,status)
  call misc_countnos(bigbuf,ifmt)
  if (ifmt<0.OR.ifmt>4) then
     call log_error(m_name,s_name,73,error_fatal,'Error cannot count points on input line')
  end if
  rewind(nin)

  !! read header information
  ivtktyp=1
  ir=0
  do
     ir=ir+1
     read(nin,fmt='(a)',iostat=status) ibuf1
     if (status/=0) then
        !!eof or error
        call log_error(m_name,s_name,3,error_fatal,'Error reading header data')
     end if
     if (ir==2) then
        if (present(kched)) kched=ibuf1
     end if
     !! look for keys at start of line
     ibuf2=adjustl(ibuf1)
     if(ibuf2(1:7)=='DATASET') then
        ibuf1=adjustl(ibuf2(9:))
        if(trim(ibuf1)=='POLYDATA') then
           !
           ivtktyp=2
        else if(trim(ibuf1)=='UNSTRUCTURED_GRID') then
           !
           ivtktyp=3
        else
           ! error
           call log_error(m_name,s_name,4,error_fatal,'Unrecognised data set type')
        end if
     else if(ibuf2(1:6)=='POINTS') then
        iltest=isnumb(ibuf1,inpos,7)
        self%np=inpos
        exit
     end if
     !! look for keys embedded in line if "=" present
     ieq=index(ibuf2,'=')
     if (ieq/=0) then
        isubstr=index(ibuf2,'Integer_Parameter=')
        if (isubstr/=0) then
           ibuf1=ibuf2(isubstr:)
           iieq=index(ibuf1,'=')
           read(ibuf1(iieq+2:),'(I3)') self%nparam
        end if
     end if
  end do

  if (self%ngtype==1) then
     ! overwrite ivtktyp so only points processed
     ivtktyp=1
  end if

  !! allocate position storage
  if(self%np>0) then
     allocate(self%posl%pos(self%np), stat=status)
     call log_alloc_check(m_name,s_name,5,status)
  else
     call log_error(m_name,s_name,6,error_fatal,'No position data')
  end if

  self%posl%np=self%np
  !! read coordinates
  call position_readonlylis(self%posl,nin,ifmt)
  !  do j=1,self%np
  !     call position_readv(self%posl%pos(j),nin)
  !  end do
  ! end position_readlis
  print '("number of geobj coordinates read = ",i10)',self%np
  call log_value("number of geobj coordinates read ",self%np)

  self%nwset=0
  geobject_typer : select case (ivtktyp)
  case(1)
     ! skip to POINT_DATA, ignore CELL and CELL_TYPES etc.
     do
        read(nin,fmt='(a)',iostat=status) ibuf1
        if (status/=0) then
           !!eof or error
           call log_error(m_name,s_name,10,error_fatal,'Error reading header data')
        end if
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:10)=='POINT_DATA') then
           iltest=isnumb(ibuf1,inobj,11)
           if (inobj/=inpos) then
              call log_error(m_name,s_name,11,error_fatal,'Numbers of data disagree')
           end if
           exit
        end if
     end do

     self%ng=inobj
     ! check geobj storage
     if(self%ng>0) then
        allocate(self%obj(self%ng), stat=status)
        call log_alloc_check(m_name,s_name,12,status)
        self%nwset=1
     else
        call log_error(m_name,s_name,13,error_fatal,'No geobj data')
     end if
     ! skip 1 or 2 records
     read(nin,fmt='(a)',iostat=status) ibuf1
     if (status/=0) then
        !!eof or error
        call log_error(m_name,s_name,14,error_fatal,'Error reading skip 1 data')
     end if
     ibuf2=adjustl(ibuf1)
     if(ibuf2(1:7)=='SCALARS') then
        ! expect a look-up table line
        read(nin,fmt='(a)',iostat=status) ibuf1
        if (status/=0) then
           !!eof or error
           call log_error(m_name,s_name,15,error_fatal,'Error reading skip 2 data')
        end if
     end if

     do j=1,self%ng
        read(nin,*,iostat=status) self%obj(j)%weight
        if(status/=0) then
           if(status<0) then
              call log_error(m_name,s_name,16,error_warning,'EOF encountered, ng reset')
              self%ng=j
              call log_value("number of points ",self%ng)
              exit
           else
              call log_error(m_name,s_name,17,error_fatal,'Error in geobj read')
           end if
        end if
     end do


     ! General surface grid of polygons
     !
  case(2)
     innd=1
     inobj=0
     insto=0
     ! read POLYGONS and CELL_TYPES
     do
        read(nin,fmt='(a)',iostat=status) ibuf1
        !!eof
        if (status/=0) then
           !!eof or error
           call log_error(m_name,s_name,20,error_fatal,'Error reading keyword data')
        end if
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:8)=='POLYGONS') then
           ibuf1=ibuf2(10:)
           read(ibuf1,*,iostat=status) inobj,insto
           if(status/=0) then
              call log_error(m_name,s_name,21,error_warning,'Error in keyword data')
           end if
           exit
        end if
     end do

     ! check geobj storage
     if(inobj>0) then
        self%ng=inobj
        allocate(self%obj2(self%ng), stat=status)
        call log_alloc_check(m_name,s_name,22,status)
     else
        call log_error(m_name,s_name,23,error_fatal,'No geobj data')
     end if
     if(insto>0) then
        self%nnod=insto-inobj
        allocate(self%nodl(self%nnod), stat=status)
        call log_alloc_check(m_name,s_name,24,status)
     else
        call log_error(m_name,s_name,25,error_fatal,'No geobj data')
     end if

     do j=1,self%ng
        read(nin,fmt='(a)',iostat=status) ibuf1
        if(status/=0) then
           if(status<0) then
              call log_error(m_name,s_name,30,error_warning,'EOF encountered, ng reset')
              self%ng=j
              exit
           else
              call log_error(m_name,s_name,31,error_fatal,'Error in geobj read')
           end if
        else
           ilnumb=isnumb(ibuf1,inumpts,1)
           if(.not.ilnumb.or.inumpts<=0) then
              call log_error(m_name,s_name,32,error_warning,'Corrupt data encountered, ng reset')
              self%ng=j
              exit
           else
              read(ibuf1,*,iostat=status) idum, (inodl(k),k=1,inumpts)
              do k=1,inumpts
                 self%nodl(innd-1+k)=inodl(k)+1
              end do
              self%obj2(j)%ptr=innd
              self%obj2(j)%typ=inumpts
              innd=innd+inumpts
           end if
        end if
     end do
     innd=innd-1
     print '("number of node pointers read = ",i10)',innd
     call log_value("number of node pointers read ",innd)
     if (innd.ne.self%nnod) then
        call log_error(m_name,s_name,33,error_warning,'Not all pointers read in, nnod reset')
        self%nnod=innd
     end if
     !
     ! now set cell types

     do j=1,self%ng
        ityp=self%obj2(j)%typ
        if (ityp==geobj_entry_table(5)) then
           ! triangles
           self%obj2(j)%typ=VTK_TRIANGLE
        else if (ityp==geobj_entry_table(9)) then
           ! quadrilaterals
           self%obj2(j)%typ=VTK_QUAD
        else
           ! error
           call log_error(m_name,s_name,34,error_fatal,'Not all polydata recognised reading data')
        end if
     end do
     ! General case of unstructured grid
     !
  case(3)
     innd=1
     inobj=0
     insto=0
     ! read CELLS and CELL_TYPES
     do
        read(nin,fmt='(a)',iostat=status) ibuf1
        !!eof
        if (status/=0) then
           !!eof or error
           call log_error(m_name,s_name,40,error_fatal,'Error reading keyword data')
        end if
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:5)=='CELLS') then
           ibuf1=ibuf2(6:)
           read(ibuf1,*,iostat=status) inobj,insto
           if(status/=0) then
              call log_error(m_name,s_name,41,error_warning,'Error in keyword data')
           end if
           exit
        end if
     end do

     ! check geobj storage
     if(inobj>0) then
        self%ng=inobj
        allocate(self%obj2(self%ng), stat=status)
        call log_alloc_check(m_name,s_name,42,status)
     else
        call log_error(m_name,s_name,43,error_fatal,'No geobj data')
     end if
     if(insto>0) then
        self%nnod=insto-inobj
        allocate(self%nodl(self%nnod), stat=status)
        call log_alloc_check(m_name,s_name,44,status)
     else
        call log_error(m_name,s_name,45,error_fatal,'No geobj data')
     end if

     do j=1,self%ng
        read(nin,fmt='(a)',iostat=status) ibuf1
        if(status/=0) then
           if(status<0) then
              call log_error(m_name,s_name,50,error_warning,'EOF encountered, ng reset')
              self%ng=j
              exit
           else
              call log_error(m_name,s_name,51,error_fatal,'Error in geobj read')
           end if
        else
           ilnumb=isnumb(ibuf1,inumpts,1)
           if(.not.ilnumb.or.inumpts<=0) then
              call log_error(m_name,s_name,52,error_warning,'Corrupt data encountered, ng reset')
              self%ng=j
              exit
           else
              read(ibuf1,*,iostat=status) idum, (inodl(jj),jj=1,inumpts)
              do jj=1,inumpts
                 self%nodl(innd-1+jj)=inodl(jj)+1
              end do
              self%obj2(j)%ptr=innd
              innd=innd+inumpts
           end if
        end if
     end do
     innd=innd-1
     print '("number of node pointers read = ",i10)',innd
     call log_value("number of node pointers read ",innd)
     if (innd.ne.self%nnod) then
        call log_error(m_name,s_name,53,error_warning,'Not all pointers read in, nnod reset')
        self%nnod=innd
     end if
     !
     ! Now read CELL_TYPES
     do
        read(nin,fmt='(a)',iostat=status) ibuf1
        if (status/=0) then
           !!eof or error
           call log_error(m_name,s_name,60,error_fatal,'Error reading keyword data')
        end if
        ibuf2=adjustl(ibuf1)
        if(ibuf2(1:10)=='CELL_TYPES') then
           ! ignore rest of line
           exit
        end if
     end do

     do j=1,self%ng
        read(nin,*,iostat=status) ityp
        if(status/=0) then
           if(status<0) then
              call log_error(m_name,s_name,61,error_warning,'EOF encountered, ng reset')
              self%ng=j
              exit
           else
              call log_error(m_name,s_name,62,error_fatal,'Error in cell types read')
           end if
        else
           self%obj2(j)%typ=ityp
        end if
     end do
  end select geobject_typer


  print '("number of geobj read = ",i10)',self%ng
  call log_value("number of geobj read ",self%ng)
  call log_error(m_name,s_name,70,log_info,'geobjlist read in from data file')

  !w write(*,*) (self%nodl(k),k=1,100)
  call geobjlist_addcube(self)
  !w write(*,*) 'after'
  !w write(*,*) (self%nodl(k),k=1,100)

end subroutine geobjlist_read
!---------------------------------------------------------------------
!> add cube(s) to  geobjlist data structure
subroutine geobjlist_addcube(self)

  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data

  !! local
  character(*), parameter :: s_name='geobjlist_addcube' !< subroutine name
  type(posveclis_t) :: zposl !< augmented list of positions
  type(geobj2_t), dimension(:), allocatable :: iobj2 !< augmented list of object variables
  integer(ki4) :: inp !< number of entries in augmented list of positions
  integer(ki4) :: ing !< number of entries in augmented list of objects
  integer(ki4) :: ioffp !< offset for addition to zposl
  integer(ki4) :: ioffg !< offset for addition to iobj2
  integer(ki4), dimension(:), allocatable :: inodl !< augmented list of nodes
  integer(ki4) :: innod !< number of entries in augmented list of nodes
  integer(ki4) :: ioffn !< offset for addition to inodl
  integer(ki4) :: idelta !< offset of positions for addition to inodl
  real(kr4), dimension(3) :: xlbb !< min corner of geometry
  real(kr4), dimension(3) :: xubb !< max corner of geometry
  real(kr4), dimension(3)  :: zlow !< temporary lower corner of bounding cube
  real(kr4), dimension(3)  :: zup !< temporary upper corner of bounding cube
  real(kr4) :: zdelta   !< temporary cube separation
  integer(ki4), dimension(36), parameter :: icube= &!< vtk description of cube surface triangles
  (/2,      0,        3,        3,        0,        1,        1, &
5     ,        0,        5,        4,        0,        5,        1, &
6     ,        6,        1,        3,        6,        3,        7, &
7     ,        3,        2,        0,        4,        2,        4, &
7     ,        2,        7,        4,        6,        6,        4, &
5     /)
  !     check external cubes present
  !write(*,*) 'self%nbdcub',self%nbdcub
  if (self%nbdcub==0) return
  ! allocate new geobjlist position storage
  inp=self%np+8*self%nbdcub
  allocate(zposl%pos(inp),  stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  ! allocate new geobjlist object description storage
  innod=self%nnod+12*3*self%nbdcub
  allocate(inodl(innod), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  ! allocate new geobjlist object type storage
  ing=self%ng+12*self%nbdcub
  allocate( iobj2(ing), stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  ! find bounds on model
  xlbb=self%posl%pos(1)%posvec
  xubb=xlbb
  do j=2,self%np
     do i=1,3
        xlbb(i)=min(xlbb(i),self%posl%pos(j)%posvec(i))
        xubb(i)=max(xubb(i),self%posl%pos(j)%posvec(i))
     end do
  end do
  ! generate and add positions to new
  ioffp=0
  do j=1,self%nbdcub
     if (self%cornflag==1) then
        do i=1,3
           zlow(i)=min(xlbb(i),self%lowcorner(i))
           zup(i)=max(xubb(i),self%upcorner(i))
        end do
     else
        if (j==1) zdelta=self%dilen
        if (j==2) zdelta=self%dilen+self%dolen
        do i=1,3
           zlow(i)=xlbb(i)-zdelta
           zup(i)=xubb(i)+zdelta
        end do
     end if
     zposl%pos(ioffp+1)%posvec=zlow
     zposl%pos(ioffp+2)%posvec=(/zup(1),zlow(2),zlow(3)/)
     zposl%pos(ioffp+3)%posvec=(/zlow(1),zup(2),zlow(3)/)
     zposl%pos(ioffp+4)%posvec=(/zup(1),zup(2),zlow(3)/)
     zposl%pos(ioffp+5)%posvec=(/zlow(1),zlow(2),zup(3)/)
     zposl%pos(ioffp+6)%posvec=(/zup(1),zlow(2),zup(3)/)
     zposl%pos(ioffp+7)%posvec=zup
     zposl%pos(ioffp+8)%posvec=(/zlow(1),zup(2),zup(3)/)
     ioffp=ioffp+8
  end do
  ! add objects to new
  ioffn=0
  ioffg=0
  do j=1,self%nbdcub
     idelta=(j-1)*8
     do i=1,36
        inodl(ioffn+i)=icube(i)+idelta+1
     end do
     do i=1,12
        iobj2(ioffg+i)%ptr=ioffn+1
        iobj2(ioffg+i)%typ=VTK_TRIANGLE
        ioffn=ioffn+3
     end do
     ioffg=ioffg+12
  end do
  ! add rest to new
  do j=1,self%np
     zposl%pos(ioffp+j)%posvec=self%posl%pos(j)%posvec
  end do
  do j=1,self%ng
     iobj2(ioffg+j)%ptr=self%obj2(j)%ptr+ioffn
     iobj2(ioffg+j)%typ=self%obj2(j)%typ
  end do
  do j=1,self%nnod
  inodl(ioffn+j)=self%nodl(j)+8*self%nbdcub
  end do
  ! replace self with new
  deallocate(self%posl%pos)
  deallocate(self%nodl)
  deallocate(self%obj2)
  ! reallocate geobjlist position storage
  allocate(self%posl%pos(inp), stat=status)
  call log_alloc_check(m_name,s_name,11,status)
  ! reallocate geobjlist object description storage
  allocate(self%nodl(innod), stat=status)
  call log_alloc_check(m_name,s_name,12,status)
  ! reallocate geobjlist object type storage
  allocate(self%obj2(ing), stat=status)
  call log_alloc_check(m_name,s_name,13,status)
  ! copy new arrays and indices
  do j=1,inp
     self%posl%pos(j)%posvec=zposl%pos(j)%posvec
  end do
  do j=1,ing
     self%obj2(j)%ptr=iobj2(j)%ptr
     self%obj2(j)%typ=iobj2(j)%typ
  end do
  self%nodl=inodl

  self%nnod=innod
  self%np=inp
  self%posl%np=inp
  self%ng=ing
  ! tidy up
  deallocate(zposl%pos)
  deallocate(inodl)
  deallocate(iobj2)

end subroutine geobjlist_addcube
!---------------------------------------------------------------------
!> write (vtk)  geobjlist data structure
subroutine geobjlist_writev(self,kchar,kplot)

  !! arguments
  type(geobjlist_t), intent(in) :: self !< geobj list data
  character(*), intent(in) :: kchar  !< case
  integer(ki4), intent(in) :: kplot   !< output channel for vis. data


  !! local
  character(*), parameter :: s_name='geobjlist_writev' !< subroutine name
  integer(ki4) :: isum !< sum workspace
  integer(ki2) :: inn !< number of nodes

  plot_type: select case (kchar)
  case('geometry')

     write(kplot,'(''DATASET UNSTRUCTURED_GRID'')')
     write(kplot,'(''POINTS '',I8, '' float'')') self%np
     do j=1,self%np
        call position_writev(self%posl%pos(j),kplot)
     end do
     write(kplot, '('' '')')

     ! count entries in CELLS
     isum=0
     do j=1,self%ng
        inn=geobj_entry_table(self%obj2(j)%typ)
        isum=isum+inn
     end do

     ! output CELL data
     write(kplot,'(''CELLS '',I8,I8)') self%ng,self%ng+isum
     i=1
     do j=1,self%ng
        inn=geobj_entry_table(self%obj2(j)%typ)
        write(kplot,'(8(1X,I8))') inn,(self%nodl(i+ij-1)-1,ij=1,inn)
        i=i+inn
     end do
     write(kplot, '('' '')')

     ! output CELL types
     write(kplot,'(''CELL_TYPES '',I8)') self%ng
     do j=1,self%ng
        write(kplot,'(1X,I8)') self%obj2(j)%typ
     end do

     if (self%nwset==2) then
        write(kplot, '('' '')')
        !! weights
        write(kplot,'(''CELL_DATA'',I8)') self%ng
        write(kplot,'(''SCALARS weights float 1'')')
        write(kplot,'(''LOOKUP_TABLE default'')')
        write(kplot,cfmtbs1) (self%obj(i)%weight, i=1,self%ng)
     end if

  case default
     !! plot list of all quantised positions
     ! output points only
     write(kplot,'(''DATASET UNSTRUCTURED_GRID'')')
     write(kplot,'(''POINTS '',I8, '' float'')') self%np
     do j=1,self%np
        call position_writev(self%posl%pos(j),kplot)
     end do
     write(kplot, '('' '')')

     write(kplot,'(''CELLS 1 '',I8)') self%np+1
     write(kplot,'(I8)') self%np,  (i-1,i=1,self%np)
     write(kplot,'(''CELL_TYPES 1'')')
     write(kplot,'(''2'')')
     write(kplot, '('' '')')

     if (self%ngtype==1) then
        !! weights
        write(kplot,'(''POINT_DATA'',I8)') self%np
        write(kplot,'(''SCALARS weights float 1'')')
        write(kplot,'(''LOOKUP_TABLE default'')')
        write(kplot,cfmtbs1) (self%obj(i)%weight, &
 &      i=1,self%np)
     end if

  end select plot_type

end subroutine geobjlist_writev
!---------------------------------------------------------------------
!> rearrange (vtk)  geobjlist data structure
subroutine geobjlist_nodlmv(self,infilelevel,knlevel,kpowe)

  !! arguments
  type(geobjlist_t), intent(inout) :: self   !< geobj list data
  integer(ki4), intent(in)       :: infilelevel !< input refinement level
  integer(ki4), intent(in)       :: knlevel !< output refinement level
  integer(ki4), intent(in)       :: kpowe !< number of level 1 elements

  !! local
  character(*), parameter :: s_name='geobjlist_nodlmv' !< subroutine name
  integer(ki4) :: ie   !< element number
  integer(ki4) :: ide   !< entry in element array
  integer(ki4) :: iee   !< entry in element array

  if (knlevel==1) then
     self%ng=kpowe
     if (infilelevel==2) then
        do i=1,kpowe
           ie=3*(i-1)
           iee=3*(kpowe+ie)
           self%nodl(ie+1)=self%nodl(iee+1)
           self%nodl(ie+2)=self%nodl(iee+4)
           self%nodl(ie+3)=self%nodl(iee+7)
        end do
     else if (infilelevel==3) then
        do i=1,kpowe
           ie=3*(i-1)
           iee=3*4*(kpowe+ie)
           self%nodl(ie+1)=self%nodl(iee+10)
           self%nodl(ie+2)=self%nodl(iee+19)
           self%nodl(ie+3)=self%nodl(iee+28)
        end do
     end if
  else if (knlevel==2) then
     self%ng=4*kpowe
     if (infilelevel==3) then
        do i=1,kpowe
           ie=3*(i-1)
           iee=3*4*(kpowe+ie)
           self%nodl(ie+1)=self%nodl(iee+1)
           self%nodl(ie+2)=self%nodl(iee+4)
           self%nodl(ie+3)=self%nodl(iee+7)
           ide=3*(kpowe+ie)
           self%nodl(ide+1)=self%nodl(iee+10)
           self%nodl(ide+2)=self%nodl(iee+7)
           self%nodl(ide+3)=self%nodl(iee+4)
           self%nodl(ide+4)=self%nodl(iee+19)
           self%nodl(ide+5)=self%nodl(iee+1)
           self%nodl(ide+6)=self%nodl(iee+7)
           self%nodl(ide+7)=self%nodl(iee+28)
           self%nodl(ide+8)=self%nodl(iee+4)
           self%nodl(ide+9)=self%nodl(iee+1)
        end do
     end if
  end if

end subroutine geobjlist_nodlmv
!---------------------------------------------------------------------
!>< control sorting of  coordinates into bins
subroutine geobjlist_step(self,btree)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  type(btree_t), intent(inout) :: btree !< btree data
  !! local
  character(*), parameter :: s_name='geobjlist_step' !< subroutine name

  ! rotate and translate positions as needed

  !     call geobjlist_tfm(self,1)
  call position_tfmlis(self%posl,self%tfmdata)

  ! get bounding box and quantise positions
  call geobjlist_getbb(self,btree)

  !     call geobjlist_qtfm(self,1)
  call position_qtfmlis(self%posl,self%quantfm)

  ! create btree according to quantised positions
  if (btree%nttype==1) then
     ! special for BSP
     call geobjlist_bin(self,btree)
  else
     call geobjlist_mbin(self,btree)
  end if


end subroutine geobjlist_step
!---------------------------------------------------------------------
!> control calculation of density
subroutine geobjlist_stepquery(self,btree,query)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  type(btree_t), intent(inout) :: btree !< btree data
  type(queryset_t), intent(inout) :: query   !< query position data
  !! local
  character(*), parameter :: s_name='geobjlist_stepquery' !< subroutine name

  ! rotate and quantise geobj positions
  !     call geobjlist_tfm(self,1)
  call position_tfmlis(self%posl,self%tfmdata)
  call position_qtfmlis(self%posl,self%quantfm)
  ! only quantise query positions
  call position_qtfmlis(query%posl,self%quantfm)

  ! calculate density at quantised positions
  call geobjlist_query(self,btree,query)

end subroutine geobjlist_stepquery
!---------------------------------------------------------------------
!> transform positions
subroutine geobjlist_tfm(self,kt)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  integer(ki4), intent(in) :: kt   !< transform or inverse

  !! local
  character(*), parameter :: s_name='geobjlist_tfm' !< subroutine name
  type(posvecl_t) :: zpos !< local variable
  type(tfmdata_t) :: ztfmdata !< local variable

  ztfmdata%ntfm=self%tfmdata%ntfm
  ztfmdata%scale=self%tfmdata%scale
  ztfmdata%offset=self%tfmdata%offset
  ztfmdata%matrix=self%tfmdata%matrix

  if (kt==1) then
     do j=1,self%np
        zpos=position_tfm(self%posl%pos(j),ztfmdata)
        self%posl%pos(j)=zpos
     end do
  else
     do j=1,self%np
        zpos=position_invtfm(self%posl%pos(j),ztfmdata)
        self%posl%pos(j)=zpos
     end do
  end if

end subroutine geobjlist_tfm
!---------------------------------------------------------------------
!> special transform of positions
subroutine geobjlist_spectfm(self,kt,pars,ktyp)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  integer(ki4), intent(in) :: kt   !< transform or inverse
  real(kr8), intent(in), dimension(:) :: pars  !< parameters for special transform
  integer(ki4), intent(in) :: ktyp   !< type of special transform

  !! local
  character(*), parameter :: s_name='geobjlist_spectfm' !< subroutine name
  type(posvecl_t) :: zpos !< local variable
  type(tfmdata_t) :: ztfmdata !< local variable

  ! default identity
  ztfmdata%ntfm=2
  ztfmdata%matrix(1,:)=(/1.,  0.,  0. /)
  ztfmdata%matrix(2,:)=(/0.,  1.,  0. /)
  ztfmdata%matrix(3,:)=(/0.,  0.,  1. /)
  ztfmdata%scale=(/1.,  1.,  1. /)
  ztfmdata%offset=(/0.,  0.,  0. /)
  select case (ktyp)
  case (1,2,3)
     ztfmdata%offset(2)=pars(2)*const_pid
     ztfmdata%scale(3)=pars(3)
  end select

  if (kt==1) then
     do j=1,self%np
        zpos=position_tfm(self%posl%pos(j),ztfmdata)
        self%posl%pos(j)=zpos
     end do
  else
     do j=1,self%np
        zpos=position_invtfm(self%posl%pos(j),ztfmdata)
        self%posl%pos(j)=zpos
     end do
  end if

end subroutine geobjlist_spectfm
!---------------------------------------------------------------------
!> set up bounding box
subroutine geobjlist_getbb(self,btree)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  type(btree_t), intent(inout) :: btree !< btree data

  !! local
  character(*), parameter :: s_name='geobjlist_getbb' !< subroutine name
  real(kr4) :: margin !< margin added to bb
  real(kr4) :: maxside !< max side of bb
  real(kr4) :: minside !< min side of bb
  real(kr4) :: zhmin !< min side of bb for quantising
  real(kr4) :: zmind  !< local variable
  real(kr4) :: zpow2  !< local variable
  real(kr4), dimension(3) :: zlbb !< min corner of bb for quantising
  real(kr4), dimension(3) :: zubb !< max corner of bb for quantising
  integer(ki2), dimension(3) :: ixyz  !< local variable
  integer(ki2) :: lmin  !< local variable
  integer(ki2) :: imxyz  !< local variable
  integer(ki4) :: ipow2  !< local variable
  integer(ki4):: iquant !< max scaled extent of btree
  integer(ki4):: imtype  !< local variable

  !! compute xbb
  xlbb=self%posl%pos(1)%posvec
  xubb=xlbb
  do j=2,self%np
     do i=1,3
        xlbb(i)=min(xlbb(i),self%posl%pos(j)%posvec(i))
        xubb(i)=max(xubb(i),self%posl%pos(j)%posvec(i))
     end do
  end do
  self%coordbb(:,1)=xlbb
  self%coordbb(:,2)=xubb
  !!$    print '("coordbb=", 6(1x,g12.5))', ((self%coordbb(i,j),i=1,3),j=1,2)


  ! set binning bb
  maxside=0.0
  do i=1,3
     if(maxside<xubb(i)-xlbb(i))then
        maxside=xubb(i)-xlbb(i)
        l=i
     end if
  end do

  margin=maxside*epsilon(maxside)
  self%tolerance=max(self%tolerance,margin)+margin
  self%minmaxtolerance=max(self%minmaxtolerance,margin)+margin

  !! set bb with margin
  do i=1,3
     zlbb(i)=xlbb(i)-self%tolerance
     zubb(i)=xubb(i)+self%tolerance
  end do

  ! adjust depending on type of tree
  if (btree%nttype==1) then
     ! BSP
     ixyz=(/2,1,1/)
     geobjlist_maxchildn=2
  else if (btree%nttype==2) then
     ! standard octree
     ixyz=(/2,2,2/)
     geobjlist_maxchildn=8
  else if (btree%nttype==3) then
     !special top octree
     if (btree%nttalg==0) then
        ! default, assumes hx=hy=hz unspecified
        !! define top level box sizes
        !           lmin=minloc(zubb-zlbb)
        !           minside=zubb(lmin)-zlbb(lmin)
        lmin=1
        minside=zubb(1)-zlbb(1)
        do i=2,3
           if(minside>zubb(i)-zlbb(i))then
              minside=zubb(i)-zlbb(i)
              lmin=i
           end if
        end do

        do i=1,3
           ixyz(i)=int((margin+zubb(i)-zlbb(i))/minside,ki2)
        end do
        ! redefine upper bounds
        do i=1,3
           zubb(i)=zlbb(i)+float(ixyz(i))*minside
        end do
     else if (btree%nttalg==1) then
        ! assumes hx,hy,hz specified
        zwork=(zubb-zlbb)/btree%hxyz
        zmind=min(zwork(1),zwork(2),zwork(3))
        zpow2=2.
        do j=2,31
           ij=j
           if (zpow2.gt.zmind) exit
           zpow2=2.*zpow2
        end do
        btree%ndepth=ij
        zwork=(zubb-zlbb+2.*margin)/(btree%hxyz*zpow2)
        do i=1,3
           ixyz(i)=max(1,ceiling(zwork(i)))
           zubb(i)=zlbb(i)+float(ixyz(i))*btree%hxyz(i)*zpow2
        end do
     else if (btree%nttalg==2) then
        ! assumes nx,ny,nz specified
        ixyz=btree%nxyz
        !           zwork=(zubb-zlbb)*(1+epsilon(1.))/ixyz
        !           minside=min(zwork(1),zwork(2),zwork(3))
        !           zubb=zlbb+float(ixyz)*minside
     end if
     geobjlist_maxchildn=max(ixyz(1)*ixyz(2)*ixyz(3),8)
     ! test quantising not excessive
     imxyz=max(ixyz(1),ixyz(2),ixyz(3))
     ipow2=2
     do j=1,31
        ij=j
        if (ipow2.ge.imxyz) exit
        ipow2=2*ipow2
     end do
     if (self%nquant+ij.gt.ki2bits-2) then
        ! too large for type ki2 integers
        call log_error(m_name,s_name,3,error_warning,'quantising number too large, reset')
        iquant=ki2bits-2-ij
        self%nquant=iquant
        btree%exten(:,1)=(/iquant,iquant,iquant/)
     end if
     ! test depth not too great
     if (btree%ndepth>self%nquant) then
        call log_error(m_name,s_name,4,error_warning,'depth too great for quantising_number, reset')
        btree%ndepth=self%nquant
     end if
  end if


  !! set bb (with margin)
  self%binbb(:,1)=zlbb
  self%binbb(:,2)=zubb

  !!set quantising transform
  self%quantfm%nqtfm=2
  self%quantfm%hmin=(/0.,0.,0./)
  zhmin=0.
  imtype=btree%mtype

  !! loop always executes at least once, even if imtype=margin_type=0
  do j=0,min(imtype,1)
     ! effect is to expand  box by 0, h or approx 3h+2h^2, depending on imtype
     if ( (btree%nttype==3).AND.(btree%nttalg==0) ) then
        ! multi-octree, automatic quantisation
        ! on second pass zhmin is non-zero,
        self%binbb(:,1)=self%binbb(:,1)-imtype*zhmin
        self%binbb(:,2)=self%binbb(:,2)+imtype*zhmin
     else
        self%binbb(:,1)=self%binbb(:,1)-imtype*self%quantfm%hmin
        self%binbb(:,2)=self%binbb(:,2)+imtype*self%quantfm%hmin
     end if
     !
     if (btree%nttype==1) then
        ! BSP
        self%quantfm%hmin=(self%binbb(:,2)-self%binbb(:,1))/2**self%nquant
     else if (btree%nttype==2) then
        ! octree
        self%quantfm%hmin=(self%binbb(:,2)-self%binbb(:,1))/2**self%nquant
     else if (btree%nttype==3) then
        ! multi-octree
        if (btree%nttalg==0) then
           ! multi-octree, automatic quantisation
           zhmin=(self%binbb(lmin,2)-self%binbb(lmin,1))/2**self%nquant
           self%quantfm%hmin=(/zhmin,zhmin,zhmin/)
        else if (btree%nttalg==1) then
           ! multi-octree, using hxyz
           self%quantfm%hmin=(self%binbb(:,2)-self%binbb(:,1))/(ixyz*2**self%nquant)
        else if (btree%nttalg==2) then
           ! multi-octree, using nxyz
           self%quantfm%hmin=(self%binbb(:,2)-self%binbb(:,1))/(ixyz*2**self%nquant)
        end if
     end if
  end do
  self%quantfm%rhmin=1./self%quantfm%hmin
  self%quantfm%offvec=-self%binbb(:,1)*self%quantfm%rhmin

  btree%nxyz=ixyz

end subroutine geobjlist_getbb
!---------------------------------------------------------------------
!> quantise positions
subroutine geobjlist_qtfm(self,kt)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  integer(ki4), intent(in) :: kt   !< transform or inverse

  !! local
  character(*), parameter :: s_name='geobjlist_qtfm' !< subroutine name
  type(posvecl_t) :: zpos !< local variable
  type(quantfm_t) :: zquantfm !< local variable

  zquantfm%nqtfm=self%quantfm%nqtfm
  zquantfm%hmin=self%quantfm%hmin
  zquantfm%rhmin=self%quantfm%rhmin
  zquantfm%offvec=self%quantfm%offvec

  if (kt==1) then
     do j=1,self%np
        zpos=position_qtfm(self%posl%pos(j),zquantfm)
        self%posl%pos(j)=zpos
     end do

  else
     do j=1,self%np
        zpos=position_invqtfm(self%posl%pos(j),zquantfm)
        self%posl%pos(j)=zpos
     end do
  end if

end subroutine geobjlist_qtfm
!---------------------------------------------------------------------
!> transform positions on a panel-by-panel basis
subroutine geobjlist_paneltfm(self,kbods,numerics)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  integer(ki4), dimension(*), intent(in) :: kbods !< integer scalar list data of bodies for each point
  type(vnumerics_t), intent(inout) :: numerics !< input numerical parameters

  !! local
  character(*), parameter :: s_name='geobjlist_paneltfm' !< subroutine name
  type(posvecl_t) :: zpos !< scratch position vector
  type(posvecl_t) :: zpos1 !< scratch position vector
  type(tfmdata_t) :: ztfmdata !< local variable
  type(posang_t) :: zposang !< local variable
  integer(ki4) :: ibod   !< body identifier
  integer(ki4) :: ibod1  !< body identifier
  integer(ki4) :: ipan   !< panel identifier
  integer(ki4) :: inbod   !< number of bodies
  integer(ki4) :: inpan   !< number of panels
  integer(ki4) :: ipantfm   !< transform or inverse
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: ityp !< type of object (5 for triangle)
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: ipt !< number of point (starting at 1)
  integer(ki4) :: itfm !< type of transform
  integer(ki4), dimension(:), allocatable :: ipansum !< number of points in panel
  real(kr4), dimension(:,:), allocatable :: zpansum !< sum over points in panel
  !BP   integer(ki4), dimension(:), allocatable :: ibodpan !< body-to-panel array
  real(kr4) :: zthetar !< rotation angle \f$ \theta_r \f$ for poloidal tilt
  real(kr4) :: zctr !< \f$ \cos(\theta_r) \f$
  real(kr4) :: zmctr !< \f$ 1-\cos(\theta_r) \f$
  real(kr4) :: zstr !< \f$ \sin(\theta_r) \f$
  real(kr4) :: zetab !< toroidal angle \f$ \zeta_b \f$ for poloidal tilt
  real(kr4) :: zczb !< \f$ \cos(\zeta_b) \f$
  real(kr4) :: zszb !< \f$ \sin(\zeta_b) \f$
  real(kr4) :: zphir !< rotation angle \f$ \phi_r \f$ for toroidal tilt
  real(kr4) :: zcpr !< \f$ \cos(\phi_r) \f$
  real(kr4) :: zspr !< \f$ \sin(\phi_r) \f$
  real(kr4) :: zroffset !< offset in minor radius
  real(kr4) :: ztheta    !<  \f$ \theta \f$ defined relative to geometry
  real(kr4) :: zrmajor !< centre of geometry
  real(kr4) :: zzmajor !< centre of geometry
  real(kr4) :: angfac=1. !< angles in radians

  !! set up marker (weight) array for each point, so only processed once
  allocate(self%obj(self%np), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  self%nwset=1
  self%obj%weight=-1.
  !! set up posang for mm to mm transform of position to R,Z
  zposang%vec=0.
  zposang%units=3
  !! set up factor for angles
  if (numerics%angles=='degree') angfac=const_degrad

  inpan=numerics%npans
  !! set up panel summing arrays
  allocate(zpansum(3,inpan), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  zpansum=0.
  allocate(ipansum(inpan), stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  ipansum=0
  !! set up body-to-panel array
  inbod=2*inpan
  !BP   allocate(ibodpan(inbod), stat=status)
  call log_alloc_check(m_name,s_name,4,status)
  !! set up panel transformation arrays
  allocate(numerics%vpantfm%scale(3,inpan), stat=status)
  call log_alloc_check(m_name,s_name,5,status)
  allocate(numerics%vpantfm%offset(3,inpan), stat=status)
  call log_alloc_check(m_name,s_name,6,status)
  allocate(numerics%vpantfm%matrix(3,3,inpan), stat=status)
  call log_alloc_check(m_name,s_name,7,status)
  allocate(numerics%vpantfm%ntfm(inpan), stat=status)
  call log_alloc_check(m_name,s_name,8,status)

  !BP   do i=1,inpan
  !BP      do k=1,2
  !BP         ibod=numerics%panbod(k,i)
  !BP         ibodpan(ibod)=i
  !BP      end do
  !BP   end do
  !BP   !      write(*,*) 'ibodpan',ibodpan

  !! loop over objects to get panel centroids (if needed)
  do j=1,self%ng
     !w write(*,*)'j=',j
     ibod=kbods(j)
     !BP      !w    ipan=ibodpan(ibod)
     ipan=indict(inpan,numerics%panbod,ibod)
     if (ipan==0) then
        print '("Fatal error in list of panel_bodies ")'
        print '("Check number of copies - and try again")'
        call log_error(m_name,s_name,9,error_fatal,'No matching object')
     end if
     ! transform number for this panel
     ipantfm=numerics%pantfm(ipan)
     if (ipantfm>0) then
        ! type of transform
        itfm=numerics%vptfm%ntfm(ipantfm)
        if (itfm==6.OR.itfm==7.OR.itfm==42) then
           innd=self%obj2(j)%ptr
           ityp=self%obj2(j)%typ
           inumpts=geobj_entry_table(ityp)
           !! loop over points defining object
           do jj=1,inumpts
              ipt=self%nodl(innd-1+jj)
              !w if (ipt<=0) then
              !w write(*,*) 'pt<=0', j,ibod,ipan,ipantfm,itfm,innd,ityp,inumpts
              !w write(*,*) (self%nodl(k),k=1,100)
              !w end if
              if (self%obj(ipt)%weight<0.) then
                 ! add point to sum (once)
                 zpansum(:,ipan)=zpansum(:,ipan)+self%posl%pos(ipt)%posvec
                 ipansum(ipan)=ipansum(ipan)+1
                 ! mark point as summed
                 self%obj(ipt)%weight=1.
              end if

           end do
        end if
     end if

  end do


  !! process panel information
  do i=1,inpan
     ! panel centroid (if present)
     if (ipansum(i)>0) then
        zpansum(:,i)=zpansum(:,i)/max(1,ipansum(i))
     end if
     ! transform for this panel
     ipantfm=numerics%pantfm(i)
     if (ipantfm>0) then
        ! type of transform
        itfm=numerics%vptfm%ntfm(ipantfm)
        ! default is to copy
        numerics%vpantfm%ntfm(i)=itfm
        numerics%vpantfm%scale(:,i)=numerics%vptfm%scale(:,ipantfm)
        numerics%vpantfm%offset(:,i)=numerics%vptfm%offset(:,ipantfm)
        numerics%vpantfm%matrix(:,:,i)=numerics%vptfm%matrix(:,:,ipantfm)
        ! whether special processing required
        if (itfm==6) then
           !special poloidal tilting transformation
           ! theta displacement from first entry in offset vector
           zthetar=numerics%vpantfm%offset(1,i)*angfac
           zstr=sin(zthetar)
           zctr=cos(zthetar)
           zmctr=1-zctr
           ! convert centroid to (R,Z,zeta)
           zposang%opt=0
           zposang%pos=zpansum(:,i)
           call posang_invtfm(zposang,3)
           zetab=zposang%pos(3)
           zszb=sin(zetab)
           zczb=cos(zetab)
           ! set up transform
           numerics%vpantfm%ntfm(i)=4
           ! offset now centroid
           numerics%vpantfm%offset(:,i)=zpansum(:,i)
           ! rotation matrix from zthetar and zetab
           numerics%vpantfm%matrix(1,1,i)=zctr+zczb**2*zmctr
           numerics%vpantfm%matrix(2,1,i)=zczb*zszb*zmctr
           numerics%vpantfm%matrix(3,1,i)=-zszb*zstr
           numerics%vpantfm%matrix(1,2,i)=zczb*zszb*zmctr
           numerics%vpantfm%matrix(2,2,i)=zctr+zszb**2*zmctr
           numerics%vpantfm%matrix(3,2,i)=zczb*zstr
           numerics%vpantfm%matrix(1,3,i)=zszb*zstr
           numerics%vpantfm%matrix(2,3,i)=-zczb*zstr
           numerics%vpantfm%matrix(3,3,i)=zctr
        else if (itfm==7) then
           !special toroidal tilting transformation
           ! phi displacement from first entry in offset vector
           zphir=numerics%vpantfm%offset(1,i)*angfac
           zspr=sin(zphir)
           zcpr=cos(zphir)
           ! set up transform
           numerics%vpantfm%ntfm(i)=4
           ! offset now centroid
           numerics%vpantfm%offset(:,i)=zpansum(:,i)
           ! rotation matrix from zphir
           numerics%vpantfm%matrix(:,:,i)=0
           numerics%vpantfm%matrix(1,1,i)=zcpr
           numerics%vpantfm%matrix(2,1,i)=zspr
           numerics%vpantfm%matrix(1,2,i)=-zspr
           numerics%vpantfm%matrix(2,2,i)=zcpr
           numerics%vpantfm%matrix(3,3,i)=1
        else if (itfm==12) then
           ! this offset is an angle, so scale
           zphir=numerics%vpantfm%offset(3,i)*angfac
           numerics%vpantfm%offset(3,i)=zphir
        else if (itfm==22) then
           ! displace in Z only
           numerics%vpantfm%offset(1,i)=0
           numerics%vpantfm%offset(2,i)=0
           numerics%vpantfm%ntfm(i)=2
        else if (itfm==42) then
           ! move in minor radius
           ! recalculate transformation as Cartesian vector
           ! offset(2) should be 6210mm for ITER, second =0
           zrmajor=numerics%vptfm%offset(2,ipantfm)
           zzmajor=numerics%vptfm%offset(3,ipantfm)
           ! 1. convert centroid to (R,Z,zeta)
           zposang%opt=0
           zposang%pos=zpansum(:,i)
           call posang_invtfm(zposang,3)
           zpos%posvec=zposang%pos
           ! 2. calculate transform offset in (R,Z) from r displacement
           ! displacement in minor radius
           zroffset=numerics%vptfm%offset(1,ipantfm)
           ! polar angle theta of centroid
           ztheta=atan2(zpos%posvec(2)-zzmajor,zpos%posvec(1)-zrmajor)
           ! 3. set up transform
           ztfmdata%ntfm=2
           ztfmdata%scale=numerics%vptfm%scale(:,ipantfm)
           ztfmdata%offset(1)=zroffset*cos(ztheta)
           ztfmdata%offset(2)=zroffset*sin(ztheta)
           ztfmdata%offset(3)=0.
           ztfmdata%matrix=numerics%vptfm%matrix(:,:,ipantfm)
           ! and apply
           zpos1=position_tfm(zpos,ztfmdata)
           ! 4. back to cartesians
           zposang%opt=1
           zposang%pos=zpos1%posvec
           call posang_tfm(zposang,3)
           ! 5. difference determines new displacement vector for panel
           numerics%vpantfm%ntfm(i)=2
           numerics%vpantfm%offset(:,i)=zposang%pos-zpansum(:,i)
        end if
     end if

  end do


  self%obj%weight=-1.
  !! loop over objects to transform
  ! generate output whenever body number changes
  ibod1=-999
  do j=1,self%ng

     ibod=kbods(j)
     !BP      !w    ipan=ibodpan(ibod)
     ipan=indict(inpan,numerics%panbod,ibod)
     if (ipan==0) then
        print '("Fatal error in list of panel_bodies ")'
        print '("Check number of copies - and try again")'
        call log_error(m_name,s_name,19,error_fatal,'No matching object')
     end if
     ! transform for this panel
     ipantfm=numerics%pantfm(ipan)
     if (ipantfm>0) then
        ! type of transform
        itfm=numerics%vpantfm%ntfm(ipan)
        ztfmdata%ntfm=itfm
        ztfmdata%scale=numerics%vpantfm%scale(:,ipan)
        ztfmdata%offset=numerics%vpantfm%offset(:,ipan)
        ztfmdata%matrix=numerics%vpantfm%matrix(:,:,ipan)
        innd=self%obj2(j)%ptr
        ityp=self%obj2(j)%typ
        inumpts=geobj_entry_table(ityp)
        if (ibod1/=ibod) then
           ibod1=ibod
           !F11            write(*,*) 'ibod1',ibod1 !F11
           !F11            write(*,*) 'ipantfm', ipantfm !F11
           !F11            write(*,*) 'ipan', ipan !F11
           !F11            write(*,*) 'itfm', itfm !F11
           !F11            write(*,*) 'ztfmdata%scale', ztfmdata%scale !F11
           !F11            write(*,*) 'ztfmdata%offset', ztfmdata%offset !F11
           !F11            write(*,*) 'ztfmdata%matrix', ztfmdata%matrix !F11
           !F11            write(*,*) 'innd', innd !F11
           !F11            write(*,*) 'ityp', ityp !F11
           call log_value("body number",ibod1)
           call log_value("transform number",itfm)
           call log_value("----- offset-1",ztfmdata%offset(1))
           call log_value("----- offset-2",ztfmdata%offset(2))
           call log_value("----- offset-3",ztfmdata%offset(3))
        end if
        !! loop over points defining object
        do jj=1,inumpts
           ipt=self%nodl(innd-1+jj)
           if (self%obj(ipt)%weight<0.) then
              ! transform point (once)
              if (itfm<=4) then
                 ! standard transformation
                 zpos=self%posl%pos(ipt)
                 self%posl%pos(ipt)=position_tfm(zpos,ztfmdata)
              else if (itfm==12) then
                 ! offset in (R,Z,zeta)
                 ! first express in (R,Z,zeta) coordinates
                 zposang%opt=0
                 zposang%pos=self%posl%pos(ipt)%posvec
                 call posang_invtfm(zposang,3)
                 zpos%posvec=zposang%pos
                 ! transform
                 ztfmdata%ntfm=2
                 zpos1=position_tfm(zpos,ztfmdata)
                 ! back to cartesians
                 zposang%opt=1
                 zposang%pos=zpos1%posvec
                 call posang_tfm(zposang,3)
                 zpos%posvec=zposang%pos
                 self%posl%pos(ipt)=zpos
              end if
              ! mark point as transformed
              self%obj(ipt)%weight=1.
           end if

        end do
     end if

  end do

  deallocate(zpansum)
  deallocate(ipansum)
  !BP   deallocate(ibodpan)
  deallocate(self%obj)
  self%nwset=0

end subroutine geobjlist_paneltfm
!---------------------------------------------------------------------
!> sort coordinates into bins
subroutine geobjlist_bin(self,btree)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  type(btree_t), intent(inout) :: btree !< btree data
  !! local
  character(*), parameter :: s_name='geobjlist_bin' !< subroutine name
  integer(ki4) :: first  !< local variable
  integer(ki4) :: last  !< local variable
  type(geobj_t) :: iobj !< geo object
  integer(ki2) :: id  !< local variable
  integer(ki2) :: idepth  !< local variable
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: iextn  !< local variable
  integer(ki2), dimension(3,2) :: inbox !< new box
  integer(ki2), dimension(3) :: icorn !< new corner
  integer(ki2), dimension(3) :: icornn !< new corner
  integer(ki4) :: inext  !< local variable
  integer(ki4) :: inadr  !< local variable
  integer(ki4) :: ileaf  !< local variable
  logical :: ilsplit !< logical flag
  integer(ki4), dimension(8) :: iadra  !< local variable
  integer(ki4) :: iadr0  !< local variable
  integer(ki4) :: iadr1  !< local variable
  integer(ki4) :: in0  !< local variable
  integer(ki4) :: in1  !< local variable
  integer(ki4) :: iadr00  !< local variable
  integer(ki4) :: iadr01  !< local variable
  integer(ki4) :: iadd  !< local variable

  !      logical geobj_innbox
  !      external geobj_innbox

  ! Note that btree%desc(1,k),btree%pter(2,k),btree%pter(3,k)
  ! have different meaning depending whether node is terminal:
  ! Entry           |   Not terminal        |           Terminal         |
  ! btree%desc(1,k) | Splitting direction   | 0                          |
  ! btree%pter(2,k) | Pointer to child node | Address of objects in cell |
  ! btree%pter(3,k) | Pointer to child node | -1 (skip)  or -2 (empty)   |


  first=1
  last=1
  ilsplit=.false.

  !! loop over depth
  loop_99 : do j=1,btree%ndepth
     idepth=j
     id=1+mod(j-1,3)
     ! first entry in list is 0 (empty node)
     iempt=2
     iltest=.true.
     ! loop over volumes down to level j
     loop_9 : do k=1,last

        !! test for empty node
        ileaf=btree%pter(3,k)
        if (ileaf==-2.OR. (k<first.AND.ileaf/=-1) ) then
           ! no object address to worry about
           cycle
        else
           inadr=btree%pter(2,k)
           inls=btree%objectls%list(inadr,2)
           if (ileaf==-1) then
              ! no need to process further, except for object address
              ilsplit=.false.
           else
              if (inls<=self%minobj) then
                 ! no need to process further, except for object address
                 ilsplit=.false.
                 ileaf=-1
              else
                 !! test for criterion P
                 ! new box array
                 iext=btree%exten(:,btree%desc(2,k))
                 iextn=iext
                 iextn(id)=iext(id)-1
                 ! update corner
                 inbox(:,1)=btree%corner(:,k)
                 inbox(:,2)=inbox(:,1)+2**iextn
                 icorn=btree%corner(:,k)
                 icornn=icorn
                 icornn(id)=icorn(id)+2**iextn(id)

                 ! loop over list, assign to new ls
                 iadr00=iempt
                 iadr0=iempt+1
                 iadr1=1
                 in0=0
                 in1=0
                 loop_list: do l=1,inls

                    inext=btree%objectls%list(inadr+l,2)
                    iobj%geobj=inext
                    iobj%objtyp=self%ngtype
                    if ( geobj_innbox( iobj,self%posl,icorn,iextn )  ) then
                       ! add to 0 ls (iadr0)
                       call ls_add(btree%objectls,iadr0,0,inext)
                       iadr0=iadr0+1
                       in0=in0+1
                    else if ( geobj_innbox( iobj,self%posl,icornn,iextn )  ) then
                       ! add to 1 ls (iadr1)
                       call ls_add(btree%objectls,iadr1,1,inext)
                       iadr1=iadr1+1
                       in1=in1+1
                    else
                       call log_error(m_name,s_name,1,error_warning,'object unassigned')
                       self%ngunassigned=self%ngunassigned+1
                    end if

                 end do loop_list

                 ! test whether object numbers changed
                 if (in0==inls .AND. in1==inls ) then
                    ileaf=-1
                    ilsplit=.false.
                 else
                    ilsplit=.true.
                 end if
              end if
           end if
        end if

        if (ilsplit) then
           ! criterion satisfied
           !  complete 0 list (with number of entries)
           call ls_add(btree%objectls,iadr00,0,in0)
           ! add 1 list to end of 0
           iadd=iadr0
           iadr01=iadd
           !  complete 1 list
           call ls_add(btree%objectls,iadd,0,in1)
           iadd=iadd+1
           ! then add to 0 list
           call ls_copy(btree%objectls,1,1,iadd,0,in1)
           iadd=iadd+in1

           iempt=iadd
           ! not able to quit on this pass
           iltest=.false.

           ! update btree
           if (in0==0) then
              iadra(1)=-2
           else
              iadra(1)=iadr00
           end if
           if (in1==0) then
              iadra(2)=-2
           else
              iadra(2)=iadr01
           end if
           call btree_add(btree,id,k,iadra)
        else
           ! put list of entries including no. entries at end of 0 list
           iadd=iempt
           iadr00=iadd
           call ls_copy(btree%objectls,inadr,2,iadd,0,inls+1)
           iempt=iadd+inls+1
           ! ensure node gets no more processing
           btree%pter(2,k)=iadr00
           btree%pter(3,k)=ileaf
        end if


     end do loop_9
     first=last+1
     last=btree%nt
     ! copy the 0 list back to 2 list, overwriting
     call ls_copy(btree%objectls,1,0,1,2,iempt-1)
     btree%objectls%nlist=iempt-1

     if (iltest) then
        exit
     end if
  end do loop_99

  if (.not.iltest) then
     call log_error(m_name,s_name,2,error_warning, &
 &   'binary tree not completed-increase quantising and/or depth and/or minobj')
  end if

  btree%ndepth=idepth

end subroutine geobjlist_bin
!---------------------------------------------------------------------
!> sort objects into bins (one object may go in many bins)
subroutine geobjlist_mbin(self,btree)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  !> btree data
  !! Note that btree%desc(1,k),btree%pter(2,k),btree%pter(3,k)
  !! have different meaning depending whether node is terminal:
  !! Entry           |   Not terminal        |           Terminal         |
  !! btree%desc(1,k) | Splitting direction   | 0                          |
  !! btree%pter(2,k) | Pointer to child node | Address of objects in cell |
  !! btree%pter(3,k) | Pointer to child node | -1 (skip)  or -2 (empty)   |
  type(btree_t), intent(inout) :: btree  !< local variable

  !! local
  character(*), parameter :: s_name='geobjlist_mbin' !< subroutine name
  integer(ki4) :: first  !< local variable
  integer(ki4) :: last  !< local variable
  integer(ki2) :: ittype  !< local variable
  type(geobj_t) :: iobj !< geo object
  integer(ki2) :: id  !< local variable
  integer(ki2), dimension(3) :: ixyz  !< local variable
  integer(ki2) :: idepth  !< local variable
  integer(ki2), dimension(3) :: iextn  !< local variable
  real(kr4), dimension(3,2) :: box !< bounding box corners
  integer(ki2), dimension(3,geobjlist_maxchildn) :: icorna !< new corner array
  integer(ki4) :: inext  !< local variable
  integer(ki4) :: inadr  !< local variable
  integer(ki4) :: ileaf  !< local variable
  logical :: ilsplit !< logical flag
  logical :: ilall !< logical flag
  integer(ki4), dimension(geobjlist_maxchildn) :: iadra  !< local variable
  integer(ki4), dimension(geobjlist_maxchildn) :: ina  !< local variable
  integer(ki4) :: iadr0 !< loop variable
  integer(ki4) :: in0  !< local variable
  integer(ki4) :: in1  !< local variable
  integer(ki4) :: inl1  !< local variable
  integer(ki4) :: iadr00  !< local variable
  integer(ki4) :: iadr01  !< local variable
  integer(ki4) :: iadd  !< local variable
  integer(ki4), dimension(3) :: ibsiz !< actual size of box
  integer(ki2) :: ichildn  !< local variable
  integer(ki2) :: id1  !< local variable

  first=1
  last=1
  ilsplit=.false.
  ittype=btree%nttype

  !! loop over depth
  loop_99 : do j=1,btree%ndepth
     idepth=j
     id=1+mod(j-1,3)
     if (ittype==1 ) then
        ixyz=(/1,1,1/)
        ixyz(id)=2
     else
        if ( j==1 ) then
           ixyz=btree%nxyz
        else
           ixyz=(/2,2,2/)
        end if
     end if

     ! first entry in list is 0 (empty node)
     iempt=2
     iltest=.true.
     ! loop over volumes down to level j
     loop_9 : do k=1,last

        !! test for empty node
        ileaf=btree%pter(3,k)
        if (ileaf==-2.OR. (k<first.AND.ileaf/=-1) ) then
           ! no object address to worry about
           cycle
        else
           inadr=btree%pter(2,k)
           inls=btree%objectls%list(inadr,2)
           if (ileaf==-1) then
              ! no need to process further, except for object address
              ilsplit=.false.
           else
              if (inls<=self%minobj) then
                 ! no need to process further, except for object address
                 ilsplit=.false.
                 ileaf=-1
              else
                 !! test for criterion P
                 ! define child boxes
                 call btree_mdefn(btree,ixyz,k,id1,iextn,icorna,ichildn,ibsiz)
                 ! loop over child boxes
                 loop_child: do i=1,ichildn
                    box(:,1)=float(icorna(:,i))
                    box(:,2)=float(icorna(:,i)+ibsiz)
                    if (i==1) then
                       ! set 1-ls to 1 =list of unassigned
                       do l=1,inls
                          call ls_add(btree%objectls,l+2,1,1)
                       end do
                       iadr0=iempt
                    end if
                    ! loop over list, assign to new ls
                    in0=0
                    ! address of number of entries
                    iadr01=iadr0

                    iadr0=iadr0+1

                    loop_list: do l=1,inls

                       inext=btree%objectls%list(inadr+l,2)
                       iobj%geobj=self%obj2(inext)%ptr
                       iobj%objtyp=self%obj2(inext)%typ
                       if ( geobj_inbox( iobj,self%posl,self%nodl,box )  ) then
                          ! add to 0-ls (iadr0)
                          call ls_add(btree%objectls,iadr0,0,inext)
                          iadr0=iadr0+1
                          in0=in0+1
                          ! mark as assigned in 1-ls
                          call ls_add(btree%objectls,l+2,1,0)
                       end if

                    end do loop_list
                    iadra(i)=iadr01
                    ina(i)=in0
                    !  complete 0 list (with number of entries)
                    call ls_add(btree%objectls,iadr01,0,in0)
                 end do loop_child

                 ! check all assigned
                 inl1=0
                 do l=1,inls
                    in1=btree%objectls%list(inadr+l,1)
                    if (in1>0) then
                       call log_error(m_name,s_name,1,error_warning,'object unassigned')
                       inl1=inl1+1
                    end if
                 end do
                 self%ngunassigned=self%ngunassigned+inl1

                 ! test whether object numbers changed
                 ilall=.true.
                 do i=1,ichildn
                    ilall=(ilall.and.(ina(i)==inls))
                 end do
                 if (ilall) then
                    ileaf=-1
                    ilsplit=.false.
                 else
                    ilsplit=.true.
                 end if
              end if
           end if
        end if

        if (ilsplit) then
           ! criterion satisfied
           iadd=iadr0-1
           iempt=iadr0
           ! not able to quit on this pass
           iltest=.false.
           ! update btree
           do i=1,ichildn
              if (ina(i)==0) iadra(i)=-2
           end do
           call btree_madd(btree,ixyz,k,iadra,id1,iextn,icorna,ichildn,ibsiz)
        else
           ! put list of entries including no. entries at end of 0 list
           iadr00=iempt
           call ls_copy(btree%objectls,inadr,2,iadr00,0,inls+1)
           iadd=iadr00+inls
           iempt=iadd+1
           ! ensure node gets no more processing
           btree%pter(2,k)=iadr00
           btree%pter(3,k)=ileaf
        end if


     end do loop_9
     first=last+1
     last=btree%nt
     ! copy the 0 list back to 2 list, overwriting
     call ls_copy(btree%objectls,1,0,1,2,iadd)
     btree%objectls%nlist=iempt-1

     if (iltest) then
        exit
     end if
  end do loop_99

  if (.not.iltest) then
     call log_error(m_name,s_name,2,error_fatal, &
 &   'binary tree not completed-try increasing limit_geobj_in_bin')
  end if

  btree%ndepth=idepth

end subroutine geobjlist_mbin
!---------------------------------------------------------------------
!> read dat file
subroutine geobjlist_dread(self,kread)
  !! arguments
  type(geobjlist_t), intent(out) :: self   !< geobj list data
  integer(ki4), intent(inout) :: kread   !< input channel for dat file
  !integer(ki4), intent(in), optional :: kopt   !< options

  !! local
  character(*), parameter :: s_name='geobjlist_dread' !< subroutine name
  integer(ki4) :: itri !< local variable
  integer(ki4) :: iquad !< local variable
  integer(ki4) :: ipt !< local variable
  integer(ki4) :: iskip !< local variable
  integer(ki4) :: ig !< local variable
  integer(ki4) :: inod !< local variable
  integer(ki4) :: inn !< local variable
  integer(ki4) :: inpt !< number of points in geobjlist
  integer(ki4) :: inobj !< number of objects in geobjlist
  integer(ki4) :: innod !< number of nodal data in geobjlist
  integer(ki4), dimension(8) :: iarray !< local variable
  real(kr4), dimension(8) :: zarray !< local variable
  character(len=8) :: field!< first field
  character(len=8) :: fieldin  !< first field
  integer(ki4) :: islen   !< length of first field entry
  type(datline_t) :: datl !< local variable

  ! two main loops

  !! first count objects
  itri=0
  iquad=0
  ipt=0

  iskip=0
  start_first_loop_over_file : do
     read(kread,fmt='(a80)',iostat=status) ibuf2
     call log_read_check(m_name,s_name,1,status)
     ibuf1=adjustl(ibuf2)
     !D   write(*,*) iskip, ibuf1(1:10) !D
     if (ibuf1(1:10)=='BEGIN BULK') goto 1
     iskip=iskip+1
  end do start_first_loop_over_file
  !error exit
  call log_error(m_name,s_name,1,error_fatal,'BEGIN BULK not found')

1     continue
  first_loop_over_file : do
     read(kread,fmt='(a80)',iostat=status) ibuf2
     call log_read_check(m_name,s_name,2,status)
     ibuf1=adjustl(ibuf2)
     if (ibuf1(1:1)=='$') cycle
     if (ibuf1(1:7)=='ENDDATA') goto 2

     ! call datline_entry(ibuf2,field)
     fieldin=ibuf2(1:8)
     islen=len_trim(fieldin)
     if (fieldin(islen:islen)=='*') fieldin(islen:islen)=' '
     field=adjustl(fieldin)

     call datline_read(datl,ibuf2,kread)

     entry_type: select case ( field )
     case ('CTRIA3  ')
        itri=itri+1
     case ('CQUAD4  ')
        iquad=iquad+1
     case ('GRID    ')
        ipt=ipt+1
     end select entry_type

  end do first_loop_over_file

  ! warning
  call log_error(m_name,s_name,1,error_warning,'ENDDATA not found')

2     continue
  rewind(kread,iostat=status)
  call log_read_check(m_name,s_name,3,status)

  !! allocation of geobjlist storage
  ig=itri+iquad
  inod=3*itri+4*iquad
  if (ig*ipt==0) then
     ! error exit
     call log_error(m_name,s_name,1,error_fatal,'no objects found')
  end if
  call geobjlist_iinit(self,ipt,ig,inod,2,1)
  !      self%ngtype=2
  !      self%np=ipt
  !      self%ng=ig
  !      self%nnod=inod
  !      self%ngunassigned=0
  !      allocate(self%posl%pos(self%np), stat=status)
  !      call log_alloc_check(m_name,s_name,3,status)
  !      allocate(self%obj(self%ng), stat=status)
  !      allocate(self%obj2(self%ng), stat=status)
  !      call log_alloc_check(m_name,s_name,4,status)
  !      allocate(self%nodl(self%nnod), stat=status)
  !      call log_alloc_check(m_name,s_name,5,status)

  !! second loop  place objects in geobjlist
  innod=0
  inpt=0
  inobj=0

  start_second_loop_over_file : do j=1,iskip+1
     read(kread,fmt='(a80)',iostat=status) ibuf2
     call log_read_check(m_name,s_name,6,status)
  end do start_second_loop_over_file

  second_loop_over_file : do
     read(kread,fmt='(a80)',iostat=status) ibuf2
     call log_read_check(m_name,s_name,7,status)
     ibuf1=adjustl(ibuf2)
     if (ibuf1(1:1)=='$') cycle
     if (ibuf1(1:7)=='ENDDATA') exit

     call datline_read(datl,ibuf2,kread)
     call datline_fieldtype(datl)
     call datline_split(datl)

     !     write(*,*) 'datl%flds8(1)', datl%flds8(1)
     entry_type2: select case ( datl%flds8(1) )
     case ('CTRIA3  ', 'CQUAD4  ', 'CTRIA6  ')
        read(datl%flds8(1)(6:6),'(I1)') inn
        do i=1,inn+2
           if (datl%type=='large') then
              read(datl%flds16(i),*) iarray(i)
           else
              read(datl%flds8(i+1),*) iarray(i)
           end if
        end do
        inobj=inobj+1
        self%obj2(inobj)%ptr=innod+1
        !C 1st object numbers are ignored, assumed consecutive
        ! these are body numbers
        self%obj(inobj)%weight=iarray(2)
        self%nodl(innod+1:innod+inn)=iarray(3:2+inn)
        innod=innod+inn
        select case (inn)
        case(3)
           self%obj2(inobj)%typ=VTK_TRIANGLE
        case(4,8)
           self%obj2(inobj)%typ=VTK_QUAD
        case(6)
           self%obj2(inobj)%typ=VTK_TRIANGLE
        end select

     case ('GRID    ')
        do i=1,2
           if (datl%type=='large') then
              read(datl%flds16(i),*) iarray(i)
           else
              read(datl%flds8(i+1),*) iarray(i)
           end if
        end do
        do i=1,3
           if (datl%type=='large') then
              read(datl%flds16(i+2),*) zarray(i)
           else
              read(datl%flds8(i+3),*) zarray(i)
           end if
        end do
        inpt=inpt+1
        !C implicit assumption that iarray(1)=inpt
        !C coordinate system set by iarray(2)
        self%posl%pos(inpt)%posvec=zarray(1:3)
     case ('CORD2R  ')
        !C TO DO
     end select entry_type2

  end do second_loop_over_file

  print '("number of objects read = ",i10)',inobj
  print '("number of nodes read = ",i10)',innod
  print '("number of points read = ",i10)',inpt
  call log_value("number of objects read ",inobj)
  call log_value("number of nodes read ",innod)
  call log_value("number of points read ",inpt)

end subroutine geobjlist_dread
!---------------------------------------------------------------------
!> initialise geobjlist from dat file data
subroutine geobjlist_iinit(self,knpt,knobj,knnod,kgtype,kopt)

  !! arguments
  type(geobjlist_t), intent(out) :: self   !< geobj list data
  integer(ki4), intent(in)  :: knpt !< number of points in geobjlist
  integer(ki4), intent(in)  :: knobj !< number of objects in geobjlist
  integer(ki4), intent(in)  :: knnod !< number of nodal data in geobjlist
  integer(ki4), intent(in) :: kgtype   !< default type of object
  integer(ki4), intent(in) :: kopt   !< if nonzero, allocate weight as well

  !! local
  character(*), parameter :: s_name='geobjlist_iinit' !< subroutine name

  self%np=knpt
  self%ng=knobj
  self%nnod=knnod
  self%ngtype=kgtype
  self%ngunassigned=0
  allocate(self%posl%pos(self%np), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  if (kgtype==1.OR.kopt/=0) then
     allocate(self%obj(self%ng), stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     self%nwset=min(1,kopt)
  end if
  if (kgtype/=1) then
     allocate(self%nodl(self%nnod), stat=status)
     call log_alloc_check(m_name,s_name,3,status)
     allocate(self%obj2(self%ng), stat=status)
     call log_alloc_check(m_name,s_name,4,status)
     self%nwset=kopt
  end if

end subroutine geobjlist_iinit
!---------------------------------------------------------------------
!> copy geobjlist to another
subroutine geobjlist_copy(selfin,selfout,kopt)

  !! arguments
  type(geobjlist_t), intent(in) :: selfin  !< geobj list data
  type(geobjlist_t), intent(out) :: selfout   !< geobj list data
  integer(ki4), intent(in) :: kopt   !< if nonzero, try to copy weights

  !! local
  character(*), parameter :: s_name='geobjlist_copy' !< subroutine name

  selfout%np=selfin%np
  selfout%ng=selfin%ng
  selfout%nnod=selfin%nnod
  selfout%ngtype=selfin%ngtype
  selfout%ngunassigned=selfin%ngunassigned
  allocate(selfout%posl%pos(selfout%np), stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  !     selfout%posl=selfin%posl
  do j=1,selfout%np
     selfout%posl%pos(j)%posvec=selfin%posl%pos(j)%posvec
  end do

  if (selfin%ngtype==1.OR.kopt/=0) then
     allocate(selfout%obj(selfout%ng), stat=status)
     call log_alloc_check(m_name,s_name,2,status)
     selfout%nwset=min(1,kopt)
     if (selfin%nwset/=0.OR.selfin%ngtype==1) then
        selfout%obj(1:selfout%ng)%weight=selfin%obj(1:selfout%ng)%weight
     end if
  end if
  if (selfin%ngtype/=1) then
     allocate(selfout%nodl(selfout%nnod), stat=status)
     call log_alloc_check(m_name,s_name,3,status)
     selfout%nodl=selfin%nodl
     allocate(selfout%obj2(selfout%ng), stat=status)
     call log_alloc_check(m_name,s_name,4,status)
     selfout%obj2(1:selfout%ng)%ptr=selfin%obj2(1:selfout%ng)%ptr
     selfout%obj2(1:selfout%ng)%typ=selfin%obj2(1:selfout%ng)%typ
     selfout%nwset=kopt
  end if

end subroutine geobjlist_copy
!---------------------------------------------------------------------
!> append geobjlist to first
subroutine geobjlist_cumulate(self,selfin,start,copy,kopt)

  !! arguments
  type(geobjlist_t), intent(inout) :: self !< module object
  type(geobjlist_t), intent(in) :: selfin !< module object to be added
  integer(ki4), intent(in) :: start !< start number of copies required
  integer(ki4), intent(in) :: copy !< stop number of copies required
  integer(ki4), intent(in) :: kopt   !< if nonzero, try to copy weights
  !! local
  character(*), parameter :: s_name='geobjlist_cumulate' !< subroutine name
  real(kr4), dimension(:), allocatable :: rwork !< real work array
  real(kr4), dimension(:,:), allocatable :: rwork2 !< real 2D work array
  integer(ki4), dimension(:), allocatable :: iwork !< integer work array
  type(geobj1_t), dimension(:), allocatable :: zobj !< type variable
  type(geobj2_t), dimension(:), allocatable :: zobj2 !< type variable
  integer(ki4)  :: icopy !< actual number of copies required
  integer(ki4) :: ip !< number of points in initial geobjlist
  integer(ki4) :: inn !< number of nodeal data in initial geobjlist
  integer(ki4) :: inpt !< number of points in geobjlist
  integer(ki4) :: inobj !< number of objects in geobjlist
  integer(ki4) :: innod !< number of nodal data in geobjlist

  if (copy<=0) return
  icopy=copy+1-start
  !set up replacement posl array and destroy old posl
  inpt=self%np+selfin%np*icopy
  allocate(rwork2(3,inpt),stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  do j=1,self%np
     rwork2(:,j)=self%posl%pos(j)%posvec
  end do
  i=self%np+1
  do k=start,copy
     do j=1,selfin%np
        rwork2(:,i)=selfin%posl%pos(j)%posvec
        i=i+1
     end do
  end do
  deallocate(self%posl%pos)
  allocate(self%posl%pos(inpt), stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  do j=1,inpt
     self%posl%pos(j)%posvec=rwork2(:,j)
  end do
  deallocate(rwork2)

  inobj=self%ng+selfin%ng*icopy
  if (self%ngtype==1.OR.kopt/=0) then
     !set up replacement obj array and destroy old obj
     self%nwset=min(1,kopt)
     if (selfin%nwset/=0.OR.selfin%ngtype==1) then
        allocate(zobj(inobj),stat=status)
        call log_alloc_check(m_name,s_name,3,status)
        !set up replacement weight array and destroy old weight
        do j=1,self%ng
           zobj(j)%weight=self%obj(j)%weight
        end do
        i=self%ng+1
        do k=start,copy
           do j=1,selfin%ng
              zobj(i)%weight=selfin%obj(j)%weight
              i=i+1
           end do
        end do
        deallocate(self%obj)
        allocate(self%obj(inobj), stat=status)
        call log_alloc_check(m_name,s_name,4,status)
        do j=1,inobj
           self%obj(j)=zobj(j)
        end do
        deallocate(zobj)
        !     else
        !        deallocate(self%obj)
        !        allocate(self%obj(inobj), stat=status)
        !        call log_alloc_check(m_name,s_name,5,status)
     end if
  end if

  innod=self%nnod+selfin%nnod*icopy
  !F11   write(*,*)'inn,innod', self%nnod,innod !F11
  if (self%ngtype/=1) then
     !set up replacement nodl array and destroy old nodl
     allocate(iwork(innod), stat=status)
     call log_alloc_check(m_name,s_name,6,status)
     do j=1,self%nnod
        iwork(j)=self%nodl(j)
     end do
     i=self%nnod+1
     do k=start,copy
        ip=self%np+(k-start)*selfin%np
        do j=1,selfin%nnod
           iwork(i)=selfin%nodl(j)+ip
           i=i+1
        end do
     end do
     deallocate(self%nodl)
     allocate(self%nodl(innod), stat=status)
     call log_alloc_check(m_name,s_name,7,status)
     self%nodl=iwork
     deallocate(iwork)
     !set up replacement obj2 array and destroy old obj2
     allocate(zobj2(inobj),stat=status)
     call log_alloc_check(m_name,s_name,8,status)
     do j=1,self%ng
        zobj2(j)%ptr=self%obj2(j)%ptr
        zobj2(j)%typ=self%obj2(j)%typ
     end do
     i=self%ng+1
     do k=start,copy
        inn=self%nnod+(k-start)*selfin%nnod
        do j=1,selfin%ng
           zobj2(i)%ptr=selfin%obj2(j)%ptr+inn
           zobj2(i)%typ=selfin%obj2(j)%typ
           i=i+1
        end do
     end do
     deallocate(self%obj2)
     allocate(self%obj2(inobj), stat=status)
     call log_alloc_check(m_name,s_name,9,status)
     do j=1,inobj
        self%obj2(j)%ptr=zobj2(j)%ptr
        self%obj2(j)%typ=zobj2(j)%typ
     end do
     deallocate(zobj2)
     self%nwset=kopt
  end if

  self%np=inpt
  self%ng=inobj
  self%nnod=innod

end subroutine geobjlist_cumulate
!---------------------------------------------------------------------
!> create some 2-D geobjlists
subroutine geobjlist_create(self,kchar,prc,pzc,prs,pzs,knear,pdist)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  character(*),intent(in) :: kchar !< type of geobj to create
  real(kr8), dimension(:), intent(in) :: prc !< \f$ R \f$ values of contour
  real(kr8), dimension(:), intent(in) :: pzc !< \f$ Z \f$ values of contour
  real(kr8), dimension(:), intent(in) :: prs !< \f$ R \f$ values of silhouette
  real(kr8), dimension(:), intent(in) :: pzs !< \f$ Z \f$ values of silhouette
  integer(ki4), dimension(:), intent(in) :: knear !< nearest segment of silhouette
  real(kr8), dimension(:), intent(in) :: pdist !< distances of contour from silhouette

  !! local
  character(*), parameter :: s_name='geobjlist_create' !< subroutine name
  integer(ki4) :: inc    !<  length of contour in segments
  integer(ki4) :: inseg    !<  length of silhouette in segments
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: inobj !< number of geometrical objects
  integer(ki4) :: insto !< total storage
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: k !< loop counter
  integer(ki4), dimension(geobj_max_entry_table) :: inodl !< one object as nodes


  inc=size(prc)
  inseg=size(prs)
  plot_type: select case (kchar)
  case ('points')
     ! points and weights
     self%ngtype=1
     self%nwset=1
     self%np=inseg
     !! allocate storage
     if(self%np>0) then
        !        if (allocated(self%posl%pos)) deallocate(self%posl%pos)
        allocate(self%posl%pos(self%np),self%obj(self%np),stat=status)
        call log_alloc_check(m_name,s_name,1,status)
     else
        call log_error(m_name,s_name,2,error_fatal,'No silhouette data')
     end if
     self%posl%np=self%np
     !! define positions
     do j=1,inseg
        self%posl%pos(j)%posvec(1)=prs(j)
        self%posl%pos(j)%posvec(2)=0
        self%posl%pos(j)%posvec(3)=pzs(j)
     end do
     print '("number of geobj coordinates created = ",i10)',self%np
     call log_value("number of geobj coordinates created ",self%np)
     !! define weights
     self%ng=self%np
     do j=1,inseg
        self%obj(j)%weight=pdist(j)
     end do

  case ('triangles')
     ! surface grid of triangles, joining silhouette to contour
     self%ngtype=2
     self%nwset=1
     self%np=inseg+inc

     !! allocate position storage
     if(self%np>0) then
        !        if (allocated(self%posl%pos)) deallocate(self%posl%pos)
        allocate(self%posl%pos(self%np),self%obj(self%np),stat=status)
        call log_alloc_check(m_name,s_name,11,status)
     else
        call log_error(m_name,s_name,12,error_fatal,'No position data')
     end if
     self%posl%np=self%np
     !! define positions
     do j=1,inseg
        self%posl%pos(j)%posvec(1)=prs(j)
        self%posl%pos(j)%posvec(2)=0
        self%posl%pos(j)%posvec(3)=pzs(j)
     end do
     do j=inseg+1,self%np
        self%posl%pos(j)%posvec(1)=prc(j-inseg)
        self%posl%pos(j)%posvec(2)=0
        self%posl%pos(j)%posvec(3)=pzc(j-inseg)
     end do
     ! end of position creation
     print '("number of geobj coordinates created = ",i10)',self%np
     call log_value("number of geobj coordinates created ",self%np)

     innd=1
     ! CELLS and CELL_TYPES
     inobj=inseg-1
     inumpts=3
     insto=inobj*(1+inumpts)
     ! create geobj storage
     self%ng=inobj
     allocate(self%obj2(self%ng), stat=status)
     call log_alloc_check(m_name,s_name,14,status)
     self%nnod=insto-inobj
     allocate(self%nodl(self%nnod), stat=status)
     call log_alloc_check(m_name,s_name,15,status)
     do j=1,inobj
        ! triangle nodes (array positions, starting at unity)
        inodl(1)=j
        inodl(2)=j+1
        inodl(3)=inseg+knear(j)
        do k=1,inumpts
           self%nodl(innd-1+k)=inodl(k)
        end do
        self%obj2(j)%ptr=innd
        self%obj2(j)%typ=inumpts
        innd=innd+inumpts
     end do
     innd=innd-1
     ! now set cell type as triangle
     do j=1,self%ng
        self%obj2(j)%typ=VTK_TRIANGLE
     end do

     print '("number of geobj created = ",i10)',self%ng
     call log_value("number of geobj created ",self%ng)

  end select plot_type

  call log_error(m_name,s_name,70,log_info,'geobjlist created')

end subroutine geobjlist_create
!---------------------------------------------------------------------
!> create geobjlists by translation/rotation
subroutine geobjlist_create3d(self,numerics)
  !! arguments
  type(geobjlist_t), intent(out) :: self !< geobj list data
  type(dnumerics_t), intent(in) :: numerics !< input numerical parameters

  !! local
  character(*), parameter :: s_name='geobjlist_create3d' !< subroutine name
  integer(ki4) :: nseg    !<  length of curve in segments
  integer(ki4) :: n    !<  length of curve in points
  integer(ki4) :: m    !<  number of repeats of curve
  integer(ki4) :: mp1    !<  number of repeat segments
  real(kr8) :: zeta !< value of \f$ \zeta \f$
  real(kr8) :: delzeta !< increment of \f$ \zeta \f$
  type(posang_t) :: zposang   !< posang data structure
  type(posvecl_t) :: zpos1   !< one position data
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: inobj !< number of geometrical objects
  integer(ki4) :: insto !< total storage
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: io !< loop counter
  integer(ki4) :: ip !< loop counter
  ! integer(ki4), dimension(geobj_max_entry_table) :: inodl !< one object as nodes

  n=numerics%npos
  nseg=n-1
  m=numerics%div
  mp1=m+1
  self%ngtype=2
  self%nwset=0
  self%np=n*mp1

  !! allocate position storage
  if(self%np>0) then
     !        if (allocated(self%posl%pos)) deallocate(self%posl%pos)
     allocate(self%posl%pos(self%np),stat=status)
     call log_alloc_check(m_name,s_name,1,status)
  else
     call log_error(m_name,s_name,2,error_fatal,'No position data')
  end if
  self%posl%np=self%np
  ! surface grid of triangles
  !! define positions
  tfm_type: select case (numerics%tfm)
  case ('rotate')
     zeta=numerics%stang
     delzeta=(numerics%finang-numerics%stang)/m
     ip=0
     do j=1,mp1
        do i=1,n
           ip=ip+1
           zpos1%posvec(1)=numerics%r(i)
           zpos1%posvec(2)=numerics%z(i)
           zpos1%posvec(3)=zeta
           zposang%pos=zpos1%posvec
           zposang%opt=1 ; zposang%units=0
           call posang_tfm(zposang,0)
           self%posl%pos(ip)%posvec=zposang%pos
        end do
        zeta=zeta+delzeta
     end do

  case ('translate')
  end select tfm_type
  ! end of position creation
  print '("number of geobj coordinates created = ",i10)',self%np
  call log_value("number of geobj coordinates created ",self%np)

  innd=1
  ! CELLS and CELL_TYPES
  inobj=2*nseg*m
  inumpts=3
  insto=inobj*(1+inumpts)
  ! create geobj storage
  self%ng=inobj
  allocate(self%obj2(self%ng), stat=status)
  call log_alloc_check(m_name,s_name,14,status)
  self%nnod=insto-inobj
  allocate(self%nodl(self%nnod), stat=status)
  call log_alloc_check(m_name,s_name,15,status)
  io=0
  ip=0
  do j=0,m-1,2
     do i=1,nseg
        ! triangle nodes (array positions, starting at unity)
        self%nodl(ip+1)=j*n+i
        self%nodl(ip+2)=j*n+i+1
        self%nodl(ip+3)=(j+1)*n+i
        io=io+1
        innd=ip+1
        self%obj2(io)%ptr=innd
        self%obj2(io)%typ=inumpts
        ip=ip+inumpts
        self%nodl(ip+1)=j*n+i+1
        self%nodl(ip+2)=(j+1)*n+i+1
        self%nodl(ip+3)=(j+1)*n+i
        io=io+1
        innd=ip+1
        self%obj2(io)%ptr=innd
        self%obj2(io)%typ=inumpts
        ip=ip+inumpts
        self%nodl(ip+1)=(j+1)*n+i
        self%nodl(ip+2)=(j+1)*n+i+1
        self%nodl(ip+3)=(j+2)*n+i+1
        io=io+1
        innd=ip+1
        self%obj2(io)%ptr=innd
        self%obj2(io)%typ=inumpts
        ip=ip+inumpts
        self%nodl(ip+1)=(j+1)*n+i
        self%nodl(ip+2)=(j+2)*n+i+1
        self%nodl(ip+3)=(j+2)*n+i
        io=io+1
        innd=ip+1
        self%obj2(io)%ptr=innd
        self%obj2(io)%typ=inumpts
        ip=ip+inumpts
     end do
  end do
  ! now set cell type as triangle
  do j=1,self%ng
     self%obj2(j)%typ=VTK_TRIANGLE
  end do

  print '("number of geobj created = ",i10)',self%ng
  call log_value("number of geobj created ",self%ng)

  call log_error(m_name,s_name,70,log_info,'geobjlist created')

end subroutine geobjlist_create3d
!---------------------------------------------------------------------
!> compress points
subroutine geobjlist_ptcompress(self)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data

  !! local
  character(*), parameter :: s_name='geobjlist_ptcompress' !< subroutine name
  real(kr4), dimension(:,:), allocatable :: rwork2 !< real 2D work array
  integer(ki4), dimension(:), allocatable :: iwork !< integer work array
  integer(ki4),dimension(1) :: idums !< dummy
  integer(ki4) :: inpt !< index of points
  integer(ki4) :: inptc !< index of compressed points
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: ityp !< type of object (5 for triangle)
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: iw !< entry in work array
  integer(ki4) :: inode !< node number
  integer(ki4) :: ista !< start of loop

  allocate(iwork(2*self%nnod), stat=status)
  call log_alloc_check(m_name,s_name,1,status)

  inpt=0
  ! set up list of points needed for objects
  do j=1,self%ng
     innd=self%obj2(j)%ptr
     ityp=self%obj2(j)%typ
     inumpts=geobj_entry_table(ityp)
     !! loop over points defining object
     do k=1,inumpts
        inpt=inpt+1
        iwork(inpt)=self%nodl(innd-1+k)
     end do
  end do

  !! order points (using SLATEC)
  call isort(iwork, idums, inpt, 1)

  !! compress list of points, removing duplicates
  inptc=1
  iw=iwork(1)
  do j=2,inpt
     if (iwork(j)==iw) cycle
     iw=iwork(j)
     inptc=inptc+1
     iwork(inptc)=iw
  end do

  !set up replacement posl array and destroy old posl
  allocate(rwork2(3,self%np),stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  do j=1,self%np
     rwork2(:,j)=self%posl%pos(j)%posvec
  end do
  deallocate(self%posl%pos)
  ! copy back sorted and compressed version
  allocate(self%posl%pos(inptc),stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  do j=1,inptc
     self%posl%pos(j)%posvec=rwork2(:,iwork(j))
  end do
  self%posl%np=inptc
  self%np=inptc
  ! sort out nodes array
  ista=1
  do k=1,self%nnod
     inode=self%nodl(k)
     do j=ista,inptc
        if (inode==iwork(j)) then
           self%nodl(k)=j
           ista=j
           goto 1
        end if
     end do
     do j=1,ista-1
        if (inode==iwork(j)) then
           self%nodl(k)=j
           ista=j
           goto 1
        end if
     end do
     call log_error(m_name,s_name,1,error_fatal,'No matching node')
1        continue
  end do

end subroutine geobjlist_ptcompress
!---------------------------------------------------------------------
!> orient triangles consistently
subroutine geobjlist_orientri(self)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data

  !! local
  character(*), parameter :: s_name='geobjlist_orientri' !< subroutine name
  integer(ki4), dimension(:), allocatable :: imark !< mark objects as processed
  integer(ki4), dimension(:), allocatable :: iptob !< points to objects
  integer(ki4), dimension(:), allocatable :: iobjs !< objects corresponding to edges
  integer(ki4), dimension(:), allocatable :: indx !< index sorted array
  integer(ki4), dimension(:), allocatable :: indxfp !< index sorted array first points
  integer(ki4),dimension(2) :: iedge !< edge as two points
  integer(ki4),dimension(2) :: iobjm !< matching objects
  integer(ki4) :: ip !< point number
  integer(ki4) :: iob !< object number
  integer(ki4) :: inptx !< index of first point array
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: ityp !< type of object (5 for triangle)
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: ikp !< k+1
  integer(ki4) :: ista !< start of loop over objects
  integer(ki4) :: ifp !< first point working
  !W      integer(ki4) :: ij !< diagnostic j !W
  integer(ki4) :: ipt1 !< index of first point of face
  integer(ki4) :: ip1 !< index of first point of face
  integer(ki4) :: ipp1 !< index of first point of face
  integer(ki4) :: ipt2 !< index of second point of face
  integer(ki4) :: ip2 !< index of second point of face
  integer(ki4) :: ipp2 !< index of second point of face
  integer(ki4) :: i1 !< object index
  integer(ki4) :: i2 !< object index
  integer(ki4) :: imatch !< no. of matches
  integer(ki4) :: iki4=1 !< type for stack
  !     logical stack_empty ! not needed

  !! produce index of points to objects
  allocate(iptob(self%nnod),iobjs(self%nnod),stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  allocate(indx(self%nnod),stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  allocate(indxfp(self%np+1),stat=status)
  call log_alloc_check(m_name,s_name,3,status)
  ip=0
  do j=1,self%ng
     ityp=self%obj2(j)%typ
     if (ityp==VTK_TRIANGLE) then
        ! triangles only
        innd=self%obj2(j)%ptr
        inumpts=geobj_entry_table(ityp)
        !! loop over points defining object, add corresponding edges
        do k=1,inumpts
           ip=ip+1
           indx(ip)=ip
           iobjs(ip)=j
           iptob(ip)=self%nodl(innd-1+k)
        end do
     end if
  end do

  !W      write(*,*) 'iptob',iptob !W

  call isort(iptob,indx,self%nnod,2)

  !W      write(*,*) 'iptob',iptob !W
  !W      write(*,*) 'indx',indx !W
  !! index of first points
  inptx=1
  ifp=iptob(1)
  indxfp(1)=ifp
  do j=2,self%nnod
     if (iptob(j)==ifp) cycle
     ifp=iptob(j)
     inptx=inptx+1
     indxfp(inptx)=j
  end do
  ! use distance to next entry to define number of entries
  indxfp(self%np+1)=self%nnod+1
  !W     write(*,*) 'indxfp',indxfp !W

  !! to business
  allocate(imark(self%ng),stat=status)
  call log_alloc_check(m_name,s_name,5,status)
  imark=0
  call stack_init(iki4,2)

  ista=0

  ! loop over edges on stack
  stack_loop: do

     if (stack_empty(iki4,2)) then
        ! try to add edge to stack
        do j=ista+1,self%ng
           ij=j
           if (imark(j)==0) then
              ityp=self%obj2(j)%typ
              if (ityp==VTK_TRIANGLE) then
                 ! triangles only
                 innd=self%obj2(j)%ptr
                 inumpts=geobj_entry_table(ityp)
                 !! loop over points defining object, add corresponding edges
                 do k=1,inumpts
                    ikp=1+mod(k,inumpts)
                    iedge(1)=self%nodl(innd-1+k)
                    iedge(2)=self%nodl(innd-1+ikp)
                    call stack_add(iki4,2,iedge)
                    imark(j)=1
                    ista=j
                 end do
                 goto 1
              end if
           end if
        end do
        ! finished, all objects marked
        goto 9
     end if

     ! there is edge on stack
1        continue
     !W      write(*,*) 'ij',ij !W
     !W      write(*,*) 'imark',imark !W

     !! get an edge
     call stack_get(iki4,2,iedge)

     !W      if (maxval(iedge)>self%nnod) then !W
     !W      write(*,*) 'iki4',iki4 !W
     !W      call stack_get(iki4,2,iedge) !W
     !W      end if !W
     !! find adjacent triangles
     ! objects associated with edge start point
     ipt1=iedge(1)
     ip1=indxfp(ipt1)
     ipp1=indxfp(ipt1+1)-1
     ! objects associated with edge end point
     ipt2=iedge(2)
     ip2=indxfp(ipt2)
     ipp2=indxfp(ipt2+1)-1
     ! seek matching objects in list of each point
     imatch=0
     iobjm=0
     do i=ip1,ipp1
        i1=indx(i)
        do j=ip2,ipp2
           i2=indx(j)
           !W      write(*,*) 'objs', iobjs(i1),iobjs(i2) !W
           if (iobjs(i1)==iobjs(i2)) then
              imatch=imatch+1
              iobjm(imatch)=iobjs(i1)
              if (imatch==2) goto 2
           end if
        end do
     end do
     ! if no matching edge, something is corrupt
     if (imatch==0) call log_error(m_name,s_name,1,error_fatal,'No matching object')
     ! one matching edge, lies on boundary, will be OK
     if (imatch==1) cycle ! stack_loop
2        continue

     ! two matches found
     adjacent_objects: select case ( imark(iobjm(1))+imark(iobjm(2)) )
     case (0)
        ! self-explanatory error, something is corrupt
        call log_error(m_name,s_name,2,error_fatal,'One object should be marked')
     case (1)
        ! find unmarked (new) triangle
        do i=1,2
           ip=i
           if (imark(iobjm(ip))==0) exit
        end do

        iob=iobjm(ip)
        innd=self%obj2(iob)%ptr
        ityp=self%obj2(iob)%typ
        inumpts=geobj_entry_table(ityp)
        !! loop over points defining object, find matching
        do k=1,inumpts
           ikp=1+mod(k,inumpts)
           if (self%nodl(innd-1+k)==ipt1) then
              if (self%nodl(innd-1+ikp)==ipt2) then
                 ! swop edge direction and hence triangle orientation
                 self%nodl(innd-1+k)=ipt2
                 self%nodl(innd-1+ikp)=ipt1
                 exit
              end if
           end if
        end do
        ! add other edges of object to stack
        do k=1,inumpts
           ikp=1+mod(k,inumpts)
           if (self%nodl(innd-1+k)==ipt2) cycle
           iedge(1)=self%nodl(innd-1+k)
           iedge(2)=self%nodl(innd-1+ikp)
           call stack_add(iki4,2,iedge)
        end do
        imark(iob)=1
     case (2)
        ! both objects already marked, skip
     end select adjacent_objects

  end do stack_loop

9     continue
  call stack_delete(iki4)

end subroutine geobjlist_orientri
!---------------------------------------------------------------------
!> shell set of tetrahedra
subroutine geobjlist_shelltets(self,geobjtri)
  !! arguments
  type(geobjlist_t), intent(in) :: self !< geobj list data
  type(geobjlist_t), intent(out) :: geobjtri !< geobj list data for triangles

  !! local
  character(*), parameter :: s_name='geobjlist_shelltets' !< subroutine name
  integer(ki4), dimension(:), allocatable :: imark !< mark objects as processed
  integer(ki4), dimension(:), allocatable :: iptob !< points to objects
  integer(ki4), dimension(:), allocatable :: iobjs !< objects corresponding to faces
  integer(ki4), dimension(:), allocatable :: indx !< index sorted array
  integer(ki4), dimension(:), allocatable :: indxfp !< index sorted array first points
  integer(ki4), dimension(:), allocatable :: iwork !< array of surface objects
  integer(ki4),dimension(3) :: iface !< face as three points
  integer(ki4),dimension(3) :: ifacen !< new face
  integer(ki4),dimension(2) :: iobjm !< matching objects
  integer(ki4) :: ip !< point number
  integer(ki4) :: iob !< object number
  integer(ki4) :: inptx !< index of first point array
  integer(ki4) :: innd !< position of first entry for object in nodl
  integer(ki4) :: ityp !< type of object (10 for tetrahedron)
  integer(ki4) :: inumpts !< length of object in nodl array
  integer(ki4) :: ikp !< k+1 mod inumpts
  integer(ki4) :: ikpp !< k+2 mod inumpts
  integer(ki4) :: ista !< start of loop over objects
  integer(ki4) :: ifp !< first point working
  !W      integer(ki4) :: ij !< diagnostic j !W
  integer(ki4) :: ipt1 !< index of first point of face
  integer(ki4) :: ip1 !< index of first point of face
  integer(ki4) :: ipp1 !< index of first point of face
  integer(ki4) :: ipt2 !< index of second point of face
  integer(ki4) :: ip2 !< index of second point of face
  integer(ki4) :: ipp2 !< index of second point of face
  integer(ki4) :: ipt3 !< index of third point of face
  integer(ki4) :: ip3 !< index of third point of face
  integer(ki4) :: ipp3 !< index of third point of face
  integer(ki4) :: i1 !< object index
  integer(ki4) :: i2 !< object index
  integer(ki4) :: i3 !< object index
  integer(ki4) :: iptmin !< min index of face
  integer(ki4) :: iptmax !< max index of face
  integer(ki4) :: ipnmin !< min index of new face
  integer(ki4) :: ipnmax !< max index of new face
  integer(ki4) :: imatch !< no. of matches
  integer(ki4) :: inw=0 !< number of entries in work
  integer(ki4) :: inobj !< number of surface objects
  integer(ki4) :: iki4=1 !< type for stack

  !! produce index of points to objects
  allocate(iptob(self%nnod),iobjs(self%nnod),stat=status)
  call log_alloc_check(m_name,s_name,1,status)
  allocate(indx(self%nnod),stat=status)
  call log_alloc_check(m_name,s_name,2,status)
  allocate(indxfp(self%np+1),stat=status)
  call log_alloc_check(m_name,s_name,3,status)

  allocate(iwork(2*self%nnod),stat=status)
  call log_alloc_check(m_name,s_name,4,status)
  ip=0
  do j=1,self%ng
     ityp=self%obj2(j)%typ
     if (ityp==10) then
        ! tetrahedra only
        innd=self%obj2(j)%ptr
        inumpts=geobj_entry_table(ityp)
        !! loop over points defining object, index by object
        do k=1,inumpts
           ip=ip+1
           indx(ip)=ip
           iobjs(ip)=j
           iptob(ip)=self%nodl(innd-1+k)
        end do
     end if
  end do

  !W      write(*,*) 'iptob',iptob !W

  call isort(iptob,indx,self%nnod,2)

  !W      write(*,*) 'iptob',iptob !W
  !W      write(*,*) 'indx',indx !W
  !! index of first points
  inptx=1
  ifp=iptob(1)
  indxfp(1)=ifp
  do j=2,self%nnod
     if (iptob(j)==ifp) cycle
     ifp=iptob(j)
     inptx=inptx+1
     indxfp(inptx)=j
  end do
  ! use distance to next entry to define number of entries
  indxfp(self%np+1)=self%nnod+1
  !W      write(*,*) 'indxfp',indxfp !W

  !! to business
  allocate(imark(self%ng),stat=status)
  call log_alloc_check(m_name,s_name,5,status)
  imark=0
  call stack_init(iki4,3)

  ista=0

  ! loop over faces on stack
  stack_loop: do

     if (stack_empty(iki4,3)) then
        ! try to add face to stack
        do j=ista+1,self%ng
           ij=j
           if (imark(j)==0) then
              ityp=self%obj2(j)%typ
              if (ityp==10) then
                 ! tetrahedra only
                 innd=self%obj2(j)%ptr
                 inumpts=geobj_entry_table(ityp)
                 !! loop over points defining object, add corresponding faces
                 do k=1,inumpts
                    ikp=1+mod(k,inumpts)
                    ikpp=1+mod(ikp,inumpts)
                    iface(1)=self%nodl(innd-1+k)
                    iface(2)=self%nodl(innd-1+ikp)
                    iface(3)=self%nodl(innd-1+ikpp)
                    call stack_add(iki4,3,iface)
                    imark(j)=1
                    ista=j
                 end do
                 goto 1
              end if
           end if
        end do
        ! finished, all objects marked
        goto 9
     end if

     ! there is face on stack
1        continue
     !W      write(*,*) 'ij',ij !W
     !W      write(*,*) 'imark',imark !W

     !! get an face
     call stack_get(iki4,3,iface)

     !W      if (maxval(iface)>self%nnod) then !W
     !W      write(*,*) 'iki4',iki4 !W
     !W      call stack_get(iki4,3,iface) !W
     !W      end if !W
     !! find adjacent triangles
     ! objects associated with face 1 point
     ipt1=iface(1)
     ip1=indxfp(ipt1)
     ipp1=indxfp(ipt1+1)-1
     ! objects associated with face 2 point
     ipt2=iface(2)
     ip2=indxfp(ipt2)
     ipp2=indxfp(ipt2+1)-1
     ! objects associated with face 3 point
     ipt3=iface(3)
     ip3=indxfp(ipt3)
     ipp3=indxfp(ipt3+1)-1
     ! seek matching objects in list of each point
     imatch=0
     iobjm=0
     do i=ip1,ipp1
        i1=indx(i)
        do j=ip2,ipp2
           i2=indx(j)
           !W      write(*,*) 'objs', iobjs(i1),iobjs(i2) !W
           if (iobjs(i1)==iobjs(i2)) then
              do k=ip3,ipp3
                 i3=indx(k)
                 if (iobjs(i1)==iobjs(i3)) then
                    imatch=imatch+1
                    iobjm(imatch)=iobjs(i1)
                    if (imatch==2) goto 2
                 end if
              end do
           end if
        end do
     end do
     ! if no matching face, something is corrupt
     if (imatch==0) call log_error(m_name,s_name,1,error_fatal,'No matching object')
     ! one matching face, add to boundary geobjl
     if (imatch==1)  then
        do k=1,3
           inw=inw+1
           iwork(inw)=iface(k)
        end do
        cycle ! stack_loop
     end if
2        continue

     ! two matches found
     adjacent_objects: select case ( imark(iobjm(1))+imark(iobjm(2)) )
     case (0)
        ! self-explanatory error, something is corrupt
        call log_error(m_name,s_name,2,error_fatal,'One object should be marked')
     case (1)
        ! find unmarked (new) tet
        do i=1,2
           ip=i
           if (imark(iobjm(ip))==0) exit
        end do

        iob=iobjm(ip)
        innd=self%obj2(iob)%ptr
        ityp=self%obj2(iob)%typ
        inumpts=geobj_entry_table(ityp)
        !! ignore matching face, add other faces to stack
        !! loop over points defining object
        do k=1,inumpts
           ikp=1+mod(k,inumpts)
           ikpp=1+mod(ikp,inumpts)
           ifacen(1)=self%nodl(innd-1+k)
           ifacen(2)=self%nodl(innd-1+ikp)
           ifacen(3)=self%nodl(innd-1+ikpp)
           ipnmax=maxval(ifacen)
           ipnmin=minval(ifacen)
           iptmax=maxval(iface)
           iptmin=minval(iface)
           if ( ipnmax==iptmax.AND.ipnmin==iptmin.AND.sum(ifacen)==sum(iface) ) then
              ! do not add
           else
              ! add other faces of object to stack
              call stack_add(iki4,3,ifacen)
           end if
        end do
        imark(iob)=1
     case (2)
        ! both objects already marked, skip
     end select adjacent_objects

  end do stack_loop

9     continue
  call stack_delete(iki4)

  ! create new geobjl
  ityp=VTK_TRIANGLE
  inumpts=geobj_entry_table(ityp)
  inobj=inw/inumpts
  call geobjlist_iinit(geobjtri,self%np,inobj,inw,2,1)
  ! copy node list and points list
  geobjtri%nodl=iwork
  geobjtri%posl%pos=self%posl%pos
  innd=1
  do j=1,inobj
     ! triangle definitions
     geobjtri%obj2(j)%ptr=innd
     geobjtri%obj2(j)%typ=ityp
     innd=innd+inumpts
  end do
  ! remove unnecessary points
  call geobjlist_ptcompress(geobjtri)

end subroutine geobjlist_shelltets
!---------------------------------------------------------------------
!> perform calculation of density
subroutine geobjlist_query(self,btree,query)
  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  type(btree_t), intent(inout) :: btree !< btree data
  type(queryset_t), intent(inout) :: query   !< query position data

  !! local
  character(*), parameter :: s_name='geobjlist_query' !< subroutine name
  type(geobj_t) :: iobj !< geo object
  integer(ki4) :: inode  !< local variable

  ! loop over query positions
  do j=1,query%np
     iobj%geobj=j
     iobj%objtyp=self%ngtype
     call btree_find(btree,iobj,query,inode)
     call geobjlist_querynode(self,btree,inode,query%which,query%res(j))
  end do
  ! any additional smoothing

end subroutine geobjlist_query

subroutine geobjlist_querynode(self,btree,knode,kchar,pres)

  !! arguments
  type(geobjlist_t), intent(inout) :: self !< geobj list data
  type(btree_t), intent(in) :: btree   !< binary tree data
  integer(ki4), intent(in) :: knode   !< node
  character(*),intent(in):: kchar !< control output
  real(kr4), intent(out) :: pres   !< query result

  !! local
  character(*), parameter :: s_name='geobjlist_querynode' !< subroutine name
  integer(ki4) :: in  !< local variable
  integer(ki4) :: iadr  !< local variable
  integer(ki4) :: iaobj  !< local variable
  real(kr4) :: zobj  !< local variable
  type(posvecl_t) :: zext  !< local variable
  type(posvecl_t) :: zextt  !< local variable
  type(posvecl_t) :: zorig  !< local variable
  type(posvecl_t) :: zorigt  !< local variable
  real(kr4) :: zsum  !< local variable
  real(kr4) :: zvol  !< local variable
  integer(ki2), dimension(3) :: iext  !< local variable
  integer(ki2), dimension(3) :: iexp  !< local variable
  integer(ki2) :: iz  !< local variable
  real(kr4), dimension(3) :: zcen  !< local variable
  real(kr4) :: zerg  !< local variable
  real(kr4) :: zbeta  !< local variable
  real(kr4) :: znorm  !< local variable
  integer(ki4) :: iupnode  !< local variable
  integer(ki4) :: inode  !< local variable
  integer(ki4) :: in1  !< local variable

  !      type(posvecl_t) :: btree_floatvec

  if (knode<=0) then
     !        call log_error(m_name,s_name,1,error_warning,'Empty node found in binary tree')
     pres=0.
  else
     ! calculate volume

     ! extent
     iext=btree%exten(:,btree%desc(2,knode))
     ! convert
     iexp=2**iext
     zext=btree_floatvec(iexp)
     zextt=position_invqtfm(zext,self%quantfm)
     iz=0
     iexp=(/iz,iz,iz/)
     zorig=btree_floatvec(iexp)
     zorigt=position_invqtfm(zorig,self%quantfm)
     zvol=(zextt%posvec(1)-zorigt%posvec(1))* &
 &   (zextt%posvec(2)-zorigt%posvec(2))* &
 &   (zextt%posvec(3)-zorigt%posvec(3))

     diagnostic_type: select case (kchar)
     case('density')
        if (zvol>0.) then
           iadr=btree%pter(2,knode)
           in=btree%objectls%list(iadr,2)
           if (in==0) then
              call log_error(m_name,s_name,2,error_warning,'Empty node found in binary tree')
              !! empty list
              pres=0.
           else
              zsum=0.
              do jj=1,in
                 iaobj=btree%objectls%list(iadr+jj,2)
                 zobj=self%obj(iaobj)%weight
                 zsum=zsum+zobj
              end do
              pres=zsum/zvol
           end if
        else

           call log_error(m_name,s_name,3,error_warning,'Zero volume node found in binary tree')
           pres=0.
        end if

     case('adaptive density')
        if (zvol>0.) then
           iadr=btree%pter(2,knode)
           in=btree%objectls%list(iadr,2)
           if (in==0) then
              call log_error(m_name,s_name,2,error_warning,'Empty node found in binary tree')
              !! empty list
              pres=0.
           else
              ! base node
              zsum=0.
              do jj=1,in
                 iaobj=btree%objectls%list(iadr+jj,2)
                 zobj=self%obj(iaobj)%weight
                 zsum=zsum+zobj
              end do
              pres=zsum
              ! adjacent node
              iupnode=btree%pter(1,knode)
              inode=btree%pter(2,iupnode)
              if (inode==knode) inode=btree%pter(3,iupnode)
              iadr=btree%pter(2,inode)
              in1=btree%objectls%list(iadr,2)
              if (in1/=0.AND.in1/=in) then
                 ! add contrib from adjacent node (note same vol guaranteed)
                 do jj=1,in1
                    iaobj=btree%objectls%list(iadr+jj,2)
                    zobj=self%obj(iaobj)%weight
                    zsum=zsum+zobj
                 end do
                 pres=0.5*zsum
              end if
              pres=pres/zvol
           end if
        else

           call log_error(m_name,s_name,3,error_warning,'Zero volume node found in binary tree')
           pres=0.
        end if

     case('flux density')
        if (zvol>0.) then
           zcen=0.5*(zextt%posvec+zorigt%posvec)
           zerg=sqrt( max(0.,dot_product(zcen,zcen)) )
           ! do not include c=slite in calc of vel
           zbeta = sqrt(zerg*(zerg+2.*const_rmass))/(zerg+const_rmass)
           iadr=btree%pter(2,knode)
           in=btree%objectls%list(iadr,2)
           if (in==0) then
              call log_error(m_name,s_name,4,error_warning,'Empty node found in binary tree')
              !! empty list
              pres=0.
           else
              zsum=0.
              do jj=1,in
                 iaobj=btree%objectls%list(iadr+jj,2)
                 zobj=self%obj(iaobj)%weight
                 zsum=zsum+zobj
              end do
              znorm=zbeta*zerg**2
              pres=znorm*zsum/zvol
           end if
        else

           call log_error(m_name,s_name,5,error_warning,'Zero volume node found in binary tree')
           pres=0.
        end if
     end select diagnostic_type
  end if

end subroutine geobjlist_querynode

function indict(ndim,dict,word)
  integer(ki4) :: indict !< dictionary index
  integer(ki4), intent(in) :: ndim !< dim of dict
  integer(ki4), dimension(2,ndim), intent(in) :: dict !< dictionary
  integer(ki4), intent(in) :: word !< dictionary
  integer(ki4) :: ji !< loop variable
  integer(ki4) :: ki !< loop variable
  indict=0
  !write(*,*) 'ndim=',ndim
  !write(*,*) 'word=',word
  !write(*,*) 'dict=',dict
  do ji=1,ndim
     do ki=1,2
        if (word==dict(ki,ji)) then
           indict=ji
           return
        end if
     end do
  end do
end function indict

subroutine misc_countnos(bigbuf,kfmt)
  character(len=*),intent(in) :: bigbuf !< buffer for input
  integer(ki4), intent(out) :: kfmt !< format of buffer - number of 3-vectors
  character(len=132) :: ibuf !< buffer for input/output
  integer(ki4) :: ilen !< length of string
  integer(ki4) :: iblan !< number of bank substrings
  integer(ki4) :: isw !< switch on if last character was not blank
  integer(ki4) :: ji !< loop variable
  iblan=0
  ibuf=adjustl(bigbuf)
  ilen=len_trim(ibuf)
  isw=1
  do ji=1,ilen
     if (ibuf(ji:ji)==' ') then
        if (isw/=0) then
           iblan=iblan+1
           isw=0
        end if
     else
        isw=1
     end if
  end do
  kfmt=(iblan+1)/3
end subroutine misc_countnos

end module geobjlist_m
