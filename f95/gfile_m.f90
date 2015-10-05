module gfile_m

  use const_kind_m
  use log_m

  implicit none
  private

! public subroutines
  public :: &
 &gfile_init, &
 &gfile_close
! public types


! private variables
  integer(ki4) :: i !< loop counter
  integer(ki4) :: j !< loop counter
  integer(ki4) :: k !< loop counter
  character(*), parameter :: m_name='gfile_m' !< module name
  integer(ki4) :: nprint  !< file unit for output
  integer(ki4) :: status  !< status flag

  contains
subroutine gfile_init(fprint,descriptor,kprint)

  !! arguments
  character(len=*), intent(in) :: fprint !< file name root
  character(len=*), intent(in) :: descriptor !< dataset descriptor
  integer(ki4), intent(inout) :: kprint   !< unit number

  !! local
  character(*), parameter :: s_name='gfile_init' !< subroutine name
  logical :: unitused !< flag to test unit is available

  !! open file

  do i=99,1,-1
     inquire(i,opened=unitused)
     if(.not.unitused)then
        kprint=i
        exit
     end if
  end do

  open(unit=kprint,file=trim(fprint)//'.gnu')

  !! write gnuplot header
  write(kprint,'(''#gnuplot data file'')')
  write(kprint,'(''# '',a)') descriptor
  write(kprint,'(''           '')')

  nprint=kprint

end subroutine gfile_init

!! close gnuplot print file on unit nprint
subroutine gfile_close

  close(nprint)

end subroutine gfile_close
end module gfile_m
