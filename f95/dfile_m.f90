!> @addtogroup groupname4
!> @{
module dfile_m
!> @}
  use const_kind_m
  use log_m
  use misc_m
  use datline_h
  use datline_m

  implicit none
  private

! public subroutines
  public :: &
 &dfile_init, & !< initialise dat file
 &dfile_initdatfile, & !< open new unit and insert gnu dat header line
 &dfile_close !< close dat file

! public types
!integer(ki4), parameter, public :: dfile_made_up_data=1

! private variables
  character(*), parameter :: m_name='dfile_m' !< module name
  character(len=80) :: ibuf1 !< buffer for input/output
  character(len=80) :: ibuf2 !< buffer for input/output
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  integer(ki4) :: l !< loop counter
  integer, save :: nunit  !< file unit for read/write
  integer :: status  !< status flag
  logical :: iltest !< logical flag

  contains
!---------------------------------------------------------------------
!> initialise dat file
subroutine dfile_init(fileroot,kunit,kwrite)

  !! arguments
  character(len=*), intent(in) :: fileroot !< file name root
  integer, intent(inout) :: kunit   !< unit number
  integer(ki4), intent(in), optional :: kwrite   !< if unity, open for writing

  !! local
  character(*), parameter :: s_name='dfile_init' !< subroutine name
  logical :: ilnew=.false. !< open new file

  !! open file
  if(present(kwrite)) then
     ilnew=(kwrite/=0)
  end if
  if (ilnew) then
     call misc_getfileunit(kunit)
     open(unit=kunit,file=trim(fileroot)//'.dat',status='NEW',iostat=status)
  else
     call misc_getfileunit(kunit)
     open(unit=kunit,file=trim(fileroot)//'.dat',status='OLD',form='FORMATTED',iostat=status)
  end if

  if(status/=0)then
     !! error opening file
     call log_error(m_name,s_name,2,error_fatal,'Error opening dat file')
  else
     call log_error(m_name,s_name,2,log_info,'dat file opened')
  end if

  nunit=kunit

end subroutine dfile_init
!---------------------------------------------------------------------
!> open new unit and insert gnu dat header line
subroutine dfile_initdatfile(fprint,descriptor,kunit)

  !! arguments
  character(len=*), intent(in) :: fprint !< file name root
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer, intent(inout) :: kunit   !< unit number

  !! local
  character(*), parameter :: s_name='dfile_initdatfile' !< subroutine name

  call misc_getfileunit(kunit)
  open(unit=kunit,file=trim(fprint)//'.dat')

  !! write header (or not)
  if (descriptor(1:1)=='#') then
     write(kunit,'(a)') descriptor
  else if (descriptor(1:4)=='null') then
     continue
  else
     write(kunit,'(''# '',a)') descriptor
  end if

  nunit=kunit

end subroutine dfile_initdatfile
!---------------------------------------------------------------------
!> close dat plot file on unit nunit
subroutine dfile_close

  !! local
  character(*), parameter :: s_name='dfile_close' !< subroutine name

  close(nunit,iostat=status)
  if(status/=0)then
     !! error closing file
     call log_error(m_name,s_name,1,error_fatal,'Error closing dat file')
  else
     call log_error(m_name,s_name,2,log_info,'dat file closed')
  end if

end subroutine dfile_close

end module dfile_m
