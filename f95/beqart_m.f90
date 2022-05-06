!> @addtogroup groupname4
!> @{
module beqart_m
!> @}
  use const_kind_m
  use const_numphys_h
  use position_m
  use posang_h
  use log_m
  use misc_m
  use fmesh_h
  use beqart_h
  use beq_h
  use fmesh_m
  use vfile_m
  use spl2d_m
  use posang_m

  implicit none
  private

! public subroutines
  public :: &
  beqart_initread,  & !< open file
  beqart_initwrite, & !< open new output file
  beqart_read, &  !< read in object
  beqart_scale, &  !< scale object
  beqart_interp, &  !< interpolate object
  beqart_tfm, &  !< transform field
  beqart_sptfm, &  !< special transform field
  beqart_write, &  !< write out object
  beqart_delete, & !< delete object
  beqart_closeread, & !< close read file
  beqart_closewrite !< close write file

! private variables
  character(*), parameter :: m_name='beqart_m' !< module name
  integer(ki4)  :: status   !< error status
  integer(ki4),save  :: ninca=5      !< control file unit number
  integer(ki4),save  :: noutca=6      !< output file unit number
  character(len=80), save :: inputfile !< control file name
  character(len=80), save :: outputfile !< output file name
  character(len=80) :: ibuff !< buffer
  integer(ki4)  :: ilog      !< for namelist dump after error
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer(ki4) :: ij !< loop counter
  real(kr8), dimension(3) :: b  !< vector for B field
  logical, parameter :: lbrief=.TRUE. !< input file has no labels
!! logical :: unitused !< flag to test unit is available

  contains
!---------------------------------------------------------------------
!> open file
subroutine beqart_initread(file,kin)

  !! arguments
  character(*), intent(in) :: file !< file name
  integer(ki4), intent(out),optional :: kin   !< input channel for object data structure
  !! local
  character(*), parameter :: s_name='beqart_initread' !< subroutine name


  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kin=i exit end if end do ninca=i

  !! open file
  inputfile=trim(file)
  call log_value("beqart data file",trim(inputfile))
  call misc_getfileunit(ninca)
  open(unit=ninca,file=inputfile,status='OLD',iostat=status)
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open beqart file, ",a)',inputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open beqart data file')
     stop
  end if
  if (present(kin)) kin=ninca

end  subroutine beqart_initread
!---------------------------------------------------------------------
!> open new file
subroutine beqart_initwrite(fileroot,kout)

  !! arguments
  character(*), intent(in) :: fileroot !< file root
  integer(ki4), intent(out),optional :: kout   !< output channel for object data structure
  !! local
  character(*), parameter :: s_name='beqart_initwrite' !< subroutine name

  !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then kout=i exit end if end do noutca=i

  !! open file
  outputfile=trim(fileroot)//"_beqart.out"
  call log_value("beqart output data file",trim(outputfile))
  call misc_getfileunit(noutca)
  open(unit=noutca,file=outputfile,status='NEW',iostat=status)
  if(status/=0)then
     open(unit=noutca,file=outputfile,status='REPLACE',iostat=status)
  end if
  if(status/=0)then
     !! error opening file
     print '("Fatal error: Unable to open output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot open output data file')
     stop
  end if
  if (present(kout)) kout=noutca

end subroutine beqart_initwrite
!---------------------------------------------------------------------
!> read beqart data
subroutine beqart_read(self,infile,kformat,kin)

  !! arguments
  type(beqart_t), intent(out) :: self   !< beqart data structure
  character(len=80),intent(in) :: infile !< name of input file
  character(*),intent(in) :: kformat !< format of input file
  integer(ki4), intent(in), optional :: kin   !< output channel for beqart data structure

  !! local
  character(*), parameter :: s_name='beqart_read' !< subroutine name
  integer(ki4) :: nin   !< output channel for beqart data structure
  character(len=80) :: icsuf  !< eqdsk file suffix
  integer(ki4) :: ierr !< error flag
  integer(ki4) :: ilent !< local variable

  integer(ki4) :: nxin   !< number of field values in ION x direction
  integer(ki4) :: nyin   !< number of field values in ION y direction
  integer(ki4) :: nzin   !< number of field values in ION z direction
  integer(ki4) :: ixfp   !< number+1 of field values in beam x direction
  integer(ki4) :: iyfp   !< number+1 of field values in beam y direction
  integer(ki4) :: izfp   !< number+1 of field values in beam z direction

  real(kr4), dimension(:), allocatable :: work !< real work for vector list
  integer(ki4) :: ip   !< size of vector list data
  integer(ki4), dimension(3) :: iadim   !< dimensions of vector list
  integer(ki4) :: ioff   !< offset in loop over vector list data
  integer(ki4) :: iopt   !< control vectorread termination on error

  if(present(kin).AND.kin/=0) then
     !! assume unit already open and reading infile
     ninca=kin
  else if(present(kin).AND.kin==0) then
     !! assume unit already open and know unit
  else
     !! get file unit do i=99,1,-1 inquire(i,opened=unitused) if(.not.unitused)then ninca=i exit end if end do

     !! open file
     call misc_getfileunit(ninca)
     open(unit=ninca,file=infile,status='OLD',form='FORMATTED',iostat=status)
     if(status/=0) then
        !! error opening file
        call log_error(m_name,s_name,1,error_fatal,'Error opening data structure file')
     else
        call log_error(m_name,s_name,2,log_info,'data structure file opened')
     end if
  end if

  ! set file type depending on suffix
  call misc_fsuffixget(infile,icsuf,ierr)
  if (ierr/=0) then
     call log_error(m_name,s_name,3,ierr,'Beq field data file has no suffix')
  end if
  !(icsuf=kformat

  file_suffix: select case (icsuf)
  case ('vtk')
     iopt=1
     call vfile_rvectorread(work,ip,iadim,infile,'Bcart',ninca,iopt)
     ! reformat work as bx,by and bz
     self%fmesh%nxf=iadim(1)
     self%fmesh%nyf=iadim(2)
     self%fmesh%nzf=iadim(3)

     if (self%fmesh%nxf*self%fmesh%nyf*self%fmesh%nzf>0) &
 &   allocate(self%bx(self%fmesh%nxf,self%fmesh%nyf,self%fmesh%nzf), &
     self%by(self%fmesh%nxf,self%fmesh%nyf,self%fmesh%nzf), &
     self%bz(self%fmesh%nxf,self%fmesh%nyf,self%fmesh%nzf), &
     stat=status)
     call log_alloc_check(m_name,s_name,9,status)

     ioff=0
     do k=1,self%fmesh%nzf
        do j=1,self%fmesh%nyf
           do i=1,self%fmesh%nxf
              self%bx(i,j,k)=work(ioff+1)
              self%by(i,j,k)=work(ioff+2)
              self%bz(i,j,k)=work(ioff+3)
              ioff=ioff+3
           end do
        end do
     end do

  case default

     ! BTOR format uses ION coordinates, order
     ! BX(1,1,1), BY(1,1,1), BZ(1,1,1)
     ! BX(1,2,1), BY(1,2,1), BZ(1,2,1)
     ! .....
     ! BX(2,NY,1), BY(2,NY,1), BZ(2,NY,1)
     ! BX(2,1,1), BY(2,1,1), BZ(2,1,1)
     ! BX(2,2,1), BY(2,2,1), BZ(2,2,1)
     ! .....
     ! BX(2,NY,1), BY(2,NY,1), BZ(2,NY,1)
     ! ...
     ! ...
     ! BX(NX,1,1), BY(NX,1,1), BZ(NX,1,1)
     ! BX(NX,2,1), BY(NX,2,1), BZ(NX,2,1)
     ! .....
     ! BX(NX,NY,1), BY(NX,NY,1), BZ(NX,NY,1)
     ! so fastest varying corresponds to Y ->-xb
     ! so next fastest varying corresponds to X ->-zb
     ! so slowest fastest varying corresponds to Z ->-yb
     ! Except * this cannot be right because then |B| in original
     ! fieldi.dat is increasing with
     ! increasing X, ie. decreasing with distance along beam to torus,
     ! so assume X -> zb
     ! *Also the vertical component of ITER field should be positive,
     ! so negated b(3) also
     !
     ! 'Original' 3D MAST format has order
     ! BX(1,1,1), BY(1,1,1), BZ(1,1,1)
     ! BX(1,1,2), BY(1,1,2), BZ(1,1,2)
     ! .....
     ! BX(1,1,NZ), BY(1,1,NZ), BZ(1,1,NZ)
     ! BX(2,1,1), BY(2,1,1), BZ(2,1,1)
     ! BX(2,1,2), BY(2,1,2), BZ(2,1,2)
     ! .....
     ! BX(2,1,NZ), BY(2,1,NZ), BZ(2,1,NZ)
     ! ...
     ! ...
     ! BX(NX,1,1), BY(NX,1,1), BZ(NX,1,1)
     ! BX(NX,1,2), BY(NX,1,2), BZ(NX,1,2)
     ! .....
     ! BX(NX,1,NZ), BY(NX,1,NZ), BZ(NX,1,NZ)

     if (.NOT.lbrief) read(ninca,*,iostat=status) ibuff
     if (.NOT.lbrief) call log_read_check(m_name,s_name,3,status)
     read(ninca,*,iostat=status) nxin
     call log_read_check(m_name,s_name,4,status)
     if (.NOT.lbrief) read(ninca,*,iostat=status) ibuff
     if (.NOT.lbrief) call log_read_check(m_name,s_name,5,status)
     read(ninca,*,iostat=status) nyin
     call log_read_check(m_name,s_name,6,status)
     if (.NOT.lbrief) read(ninca,*,iostat=status) ibuff
     if (.NOT.lbrief) call log_read_check(m_name,s_name,7,status)
     read(ninca,*,iostat=status) nzin
     call log_read_check(m_name,s_name,8,status)

     format1: select case (kformat)
     case('btor','special')
        self%fmesh%nzf=nxin
        self%fmesh%nxf=nyin
        self%fmesh%nyf=nzin
     case default
        self%fmesh%nxf=nxin
        self%fmesh%nyf=nyin
        self%fmesh%nzf=nzin
     end select format1

     ixfp=self%fmesh%nxf+1
     iyfp=self%fmesh%nyf+1
     izfp=self%fmesh%nzf+1
     if (self%fmesh%nxf*self%fmesh%nyf*self%fmesh%nzf>0) &
 &   allocate(self%bx(self%fmesh%nxf,self%fmesh%nyf,self%fmesh%nzf), &
     self%by(self%fmesh%nxf,self%fmesh%nyf,self%fmesh%nzf), &
     self%bz(self%fmesh%nxf,self%fmesh%nyf,self%fmesh%nzf), &
     stat=status)
     call log_alloc_check(m_name,s_name,9,status)

     if (.NOT.lbrief) read(kin,*,iostat=status) ibuff
     if (.NOT.lbrief) call log_read_check(m_name,s_name,10,status)
     format2: select case (kformat)
        ! remember ordering of loops controls order in which file read
     case('magnet')
        do k=1,self%fmesh%nxf
           do j=1,self%fmesh%nzf
              do i=1,self%fmesh%nyf
                 read(ninca,*,iostat=status) b
                 self%bx(ixfp-k,i,izfp-j)=b(3)
                 self%by(ixfp-k,i,izfp-j)=b(2)
                 self%bz(ixfp-k,i,izfp-j)=-b(1)
                 call log_read_check(m_name,s_name,11,status)
              end do
           end do
        end do
     case('btor')
        do k=1,self%fmesh%nyf
           do j=1,self%fmesh%nzf
              do i=1,self%fmesh%nxf
                 read(ninca,*,iostat=status) b
                 self%bx(ixfp-i,k,izfp-j)=-b(2)
                 self%by(ixfp-i,k,izfp-j)=b(3)
                 self%bz(ixfp-i,k,izfp-j)=-b(1)
                 call log_read_check(m_name,s_name,12,status)
              end do
           end do
        end do
     case('special')
        do k=1,self%fmesh%nyf
           do j=1,self%fmesh%nzf
              do i=1,self%fmesh%nxf
                 read(ninca,*,iostat=status) b
                 self%bx(ixfp-i,k,j)=-b(2)
                 self%by(ixfp-i,k,j)=-b(3)
                 self%bz(ixfp-i,k,j)=-b(1)
                 call log_read_check(m_name,s_name,15,status)
              end do
           end do
        end do
     case('apb')
        do j=1,self%fmesh%nyf
           do i=1,self%fmesh%nxf
              do k=1,self%fmesh%nzf
                 read(ninca,*,iostat=status) b
                 self%bx(i,j,izfp-k)=b(1)
                 self%by(i,j,izfp-k)=b(2)
                 self%bz(i,j,izfp-k)=b(3)
                 call log_read_check(m_name,s_name,20,status)
              end do
           end do
        end do
     case default
        do k=1,self%fmesh%nzf
           do j=1,self%fmesh%nyf
              do i=1,self%fmesh%nxf
                 read(ninca,*,iostat=status) b
                 self%bx(i,j,k)=b(1)
                 self%by(i,j,k)=b(2)
                 self%bz(i,j,k)=b(3)
                 call log_read_check(m_name,s_name,20,status)
              end do
           end do
        end do
     end select format2
     if (.NOT.lbrief) read(ninca,*,iostat=status) ibuff
     if (.NOT.lbrief) call log_read_check(m_name,s_name,25,status)

  end select file_suffix

end subroutine beqart_read
!---------------------------------------------------------------------
!> write beqart data
subroutine beqart_scale(self,pscale)

  !! arguments
  type(beqart_t), intent(inout) :: self   !< beqart data structure
  real(kr8), intent(in), dimension(3) :: pscale   !< scale factor for beqart data structure

  !! local
  character(*), parameter :: s_name='beqart_scale' !< subroutine name

  do k=1,self%fmesh%nzf
     do j=1,self%fmesh%nyf
        do i=1,self%fmesh%nxf
           self%bx(i,j,k)=pscale(1)*self%bx(i,j,k)
           self%by(i,j,k)=pscale(2)*self%by(i,j,k)
           self%bz(i,j,k)=pscale(3)*self%bz(i,j,k)
        end do
     end do
  end do

end subroutine beqart_scale
!---------------------------------------------------------------------
!> interpolate beqart
subroutine beqart_interp(self,pos,b)

  !! arguments
  type(beqart_t), intent(in) :: self   !< beqart data structure
  type(posvecl_t), intent(in) :: pos !< position vector in beam coordinates
  real(kr8), dimension(3), intent(out) :: b !< B-field

  !! local
  character(*), parameter :: s_name='beqart_interp' !< subroutine name
  real(kr8) :: zxof !< x coordinate offset
  real(kr8) :: zyof !< y coordinate offset
  real(kr8) :: zzof !< z coordinate offset
  real(kr8) :: zx !< x real quantised coordinate
  real(kr8) :: zy !< y real quantised coordinate
  real(kr8) :: zz !< z real quantised coordinate
  real(kr8) :: zfx !< x fractional coordinate
  real(kr8) :: zfy !< y fractional coordinate
  real(kr8) :: zfz !< z fractional coordinate
  real(kr8) :: zmx !<  complementary x fractional coordinate
  real(kr8) :: zmy !<  complementary y fractional coordinate
  real(kr8) :: zmz !<  complementary z fractional coordinate
  integer(ki4) :: ix !< x index
  integer(ki4) :: iy !< y index
  integer(ki4) :: iz !< z index

  !convert the current position into array subscripts for the magnetic field
  zxof=pos%posvec(1)-self%fmesh%x0f
  zyof=pos%posvec(2)-self%fmesh%y0f
  zzof=pos%posvec(3)-self%fmesh%z0f
  zx=zxof/self%fmesh%dxf
  zy=zyof/self%fmesh%dyf
  zz=zzof/self%fmesh%dzf
  ix=int(zx)
  iy=int(zy)
  iz=int(zz)

  ! ensure always get an answer on grid
  ix=min(max(0,ix),self%fmesh%nxf-2)
  iy=min(max(0,iy),self%fmesh%nyf-2)
  iz=min(max(0,iz),self%fmesh%nzf-2)

  ! use linear interpolation to set the vector B
  ! fractions of unity
  zfx=zx-real(ix)
  zfy=zy-real(iy)
  zfz=zz-real(iz)
  zmx=1-zfx
  zmy=1-zfy
  zmz=1-zfz
  ix=ix+1
  iy=iy+1
  iz=iz+1
  ! invoke formula
  b(1)=zmx*(zmy*(zmz*self%bx(ix,iy,iz)+zfz*self%bx(ix,iy,iz+1)) &
  +zfy*(zmz*self%bx(ix,iy+1,iz)+zfz*self%bx(ix,iy+1,iz+1))) &
  +zfx*(zmy*(zmz*self%bx(ix+1,iy,iz)+zfz*self%bx(ix+1,iy,iz+1)) &
  +zfy*(zmz*self%bx(ix+1,iy+1,iz)+zfz*self%bx(ix+1,iy+1,iz+1)))
  b(2)=zmx*(zmy*(zmz*self%by(ix,iy,iz)+zfz*self%by(ix,iy,iz+1)) &
  +zfy*(zmz*self%by(ix,iy+1,iz)+zfz*self%by(ix,iy+1,iz+1))) &
  +zfx*(zmy*(zmz*self%by(ix+1,iy,iz)+zfz*self%by(ix+1,iy,iz+1)) &
  +zfy*(zmz*self%by(ix+1,iy+1,iz)+zfz*self%by(ix+1,iy+1,iz+1)))
  b(3)=zmx*(zmy*(zmz*self%bz(ix,iy,iz)+zfz*self%bz(ix,iy,iz+1)) &
  +zfy*(zmz*self%bz(ix,iy+1,iz)+zfz*self%bz(ix,iy+1,iz+1))) &
  +zfx*(zmy*(zmz*self%bz(ix+1,iy,iz)+zfz*self%bz(ix+1,iy,iz+1)) &
  +zfy*(zmz*self%bz(ix+1,iy+1,iz)+zfz*self%bz(ix+1,iy+1,iz+1)))

end subroutine beqart_interp
!---------------------------------------------------------------------
!> transform field data
subroutine beqart_tfm(self,kopt,kunits)

  !! arguments
  type(beqart_t), intent(inout) :: self   !< beqart data structure
  integer(ki4), intent(in)  :: kopt   !< type of transform
  integer(ki4), intent(in) :: kunits !< units of result (-3 for millimetres)

  !! local
  character(*), parameter :: s_name='beqart_tfm' !< subroutine name
  type(posang_t)  :: zposang !< position and vector for ripple

  do k=1,self%fmesh%nzf
     zposang%pos(3)=self%fmesh%z0f+(k-1)*self%fmesh%dzf
     do j=1,self%fmesh%nyf
        zposang%pos(2)=self%fmesh%y0f+(j-1)*self%fmesh%dyf
        do i=1,self%fmesh%nxf
           zposang%pos(1)=self%fmesh%x0f+(i-1)*self%fmesh%dxf
           zposang%vec(1)=self%bx(i,j,k)
           zposang%vec(2)=self%by(i,j,k)
           zposang%vec(3)=self%bz(i,j,k)
           zposang%opt=kopt ! only 0,16,32 valid values
           ! converts cartesian to polar-toroidal coordinates
           call posang_invtfm(zposang,kunits)
           self%bx(i,j,k)=zposang%vec(1)
           self%by(i,j,k)=zposang%vec(2)
           self%bz(i,j,k)=zposang%vec(3)
        end do
     end do
  end do

end subroutine beqart_tfm
!---------------------------------------------------------------------
!> special transform field data
subroutine beqart_sptfm(self,kpower)

  !! arguments
  type(beqart_t), intent(inout) :: self   !< beqart data structure
  integer(ki4), intent(in)  :: kpower   !< power of R

  !! local
  character(*), parameter :: s_name='beqart_sptfm' !< subroutine name
  real(kr8) :: zr    !<  \f$ R(x,y,z) \f$
  real(kr8) :: zfac    !<  \f$ R(x,y,z) \f$ factor
  real(kr8), dimension(3)  :: zvecf    !<  Vacuum field

  ! arrange so that field array stores F=R^(rpower)x(BR,BZ,BT/R=I/R^2)
  ! (rpower=1)
  !                                = (R.BR,R.BZ,I/R)
  do k=1,self%fmesh%nzf
     do j=1,self%fmesh%nyf
        do i=1,self%fmesh%nxf
           zr=self%fmesh%x0f+(i-1)*self%fmesh%dxf
           zfac=zr**kpower
           zvecf(1)=self%bx(i,j,k)
           zvecf(2)=self%by(i,j,k)
           zvecf(3)=self%bz(i,j,k)
           zvecf=zfac*zvecf
           ! normalised for faster field update
           self%bx(i,j,k)=zvecf(1)
           self%by(i,j,k)=zvecf(2)
           self%bz(i,j,k)=zvecf(3)/zr
        end do
     end do
  end do

end subroutine beqart_sptfm
!---------------------------------------------------------------------
!> write beqart data
subroutine beqart_write(self,kout)

  !! arguments
  type(beqart_t), intent(in) :: self   !< beqart data structure
  integer(ki4), intent(in), optional :: kout   !< output channel for beqart data structure

  !! local
  character(*), parameter :: s_name='beqart_write' !< subroutine name
  integer(ki4) :: iout   !< output channel for beqart data structure

  !! sort out unit
  if(present(kout)) then
     iout=kout
  else
     iout=noutca
  end if

  write(iout,*,iostat=status) 'nxf'
  call log_write_check(m_name,s_name,3,status)
  write(iout,*,iostat=status) self%fmesh%nxf
  call log_write_check(m_name,s_name,4,status)
  write(iout,*,iostat=status) 'nyf'
  call log_write_check(m_name,s_name,5,status)
  write(iout,*,iostat=status) self%fmesh%nyf
  call log_write_check(m_name,s_name,6,status)
  write(iout,*,iostat=status) 'nzf'
  call log_write_check(m_name,s_name,7,status)
  write(iout,*,iostat=status) self%fmesh%nzf
  call log_write_check(m_name,s_name,8,status)

  write(iout,*,iostat=status) 'beq --- 3d cartesian'
  call log_write_check(m_name,s_name,10,status)
  do k=1,self%fmesh%nzf
     do j=1,self%fmesh%nxf
        do i=1,self%fmesh%nyf
           b(1)=self%bx(i,j,k)
           b(2)=self%by(i,j,k)
           b(3)=self%bz(i,j,k)
           write(iout,*,iostat=status) b
           call log_write_check(m_name,s_name,11,status)
        end do
     end do
  end do
  write(iout,*,iostat=status) 'end of 3d cartesian'
  call log_write_check(m_name,s_name,12,status)

end subroutine beqart_write
!---------------------------------------------------------------------
!> delete object
subroutine beqart_delete(self)

  !! arguments
  type(beqart_t), intent(inout) :: self !< type containing analytic equil parameters
  !! local
  character(*), parameter :: s_name='beqart_delete' !< subroutine name

  if (self%fmesh%nxf>0) deallocate(self%bx, self%by, self%bz)

end subroutine beqart_delete
!---------------------------------------------------------------------
!> close read file
subroutine beqart_closeread

  !! local
  character(*), parameter :: s_name='beqart_closeread' !< subroutine name

  !! close file
  close(unit=ninca,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close input file, ",a)',inputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close input data file')
     stop
  end if

end subroutine beqart_closeread

!---------------------------------------------------------------------
!> close write file
subroutine beqart_closewrite

  !! local
  character(*), parameter :: s_name='beqart_closewrite' !< subroutine name

  !! close file
  close(unit=noutca,iostat=status)
  if(status/=0)then
     !! error closing file
     print '("Fatal error: Unable to close output file, ",a)',outputfile
     call log_error(m_name,s_name,1,error_fatal,'Cannot close output data file')
     stop
  end if

end subroutine beqart_closewrite

end module beqart_m
