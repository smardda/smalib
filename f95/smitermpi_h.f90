!> @addtogroup groupname3
!> @{
module smitermpi_h
!> @}
  use const_kind_m
#ifdef WITH_MPI
  use mpi
#endif

  implicit none

!> type describing mpi environment
  type, public :: decomp_t
     integer(ki4) :: myrank = 0 !< Global Rank of this process
     integer(ki4) :: myrankgeoqshad = 0 !< Rank of this process in the geoqshad group
     integer(ki4) :: myrankgeoqres = 0 !< Rank of this process in the geoqres group
     integer(ki4) :: myrankhdsgen = 0 !< Rank of this process in the hdsgen group
     integer(ki4) :: myrankpowcal = 0 !< Rank of this process in the powcal group
     integer(ki4) :: myrankshadhds = 0 !< Rank of this process in the shadhds group
     integer(ki4) :: nproc = 1 !< Number of processes
     integer(ki4) :: nprocpowcal = 1 !< Number of processes in powcal group
     integer(ki4) :: nprocgeoqshad = 1 !< Number of processes in geoqshad powcal group
     integer(ki4) :: nprocgeoqres = 1 !< Number of processes in geoq res powcal group
     integer(ki4) :: nprochdsgen = 1 !< Number of processes in hdsgen powcal group
     integer(ki4) :: nprocshadhds = 1 !< Number of processes in shadow hdsgen group
     integer(ki4) :: firstpowe !< First track number
     integer(ki4) :: lastpowe !< Last and last track number
     integer(ki4) :: firstng !< First geometry index
     integer(ki4) :: lastng !< Last geometry index
     integer(ki4) :: mynpowe !< The local number of tracks
     integer(ki4) :: myng !< The local number of tracks x levels
   ! Arrays holding the size and displacements of tracks for each process
     integer(ki4),dimension(:),allocatable :: localsize !< local variable
     integer(ki4),dimension(:),allocatable :: localdisp !< local variable
     integer(ki4) :: commsmiter !< smiter communicator, currently=MPI_COMM_WORLD
     integer(ki4) :: commpowcal !< powcal communicator
     integer(ki4) :: commgeoqshad !< geoqshad and powcal communicator
     integer(ki4) :: commgeoqres !< geoqresults and powcal communicator
     integer(ki4) :: commhdsgen !< hdsgen and powcal communicator
     integer(ki4) :: commshadhds !< geoq shadow and hdsgen communicator
     integer(ki4) :: procgeoqshad = 0 !< The process assigned to run geoq shadow
     integer(ki4) :: procgeoqres = 0 !< The process assigned to run geoq results
     integer(ki4) :: prochdsgen = 0 !< The process assigned to run hdsgen
     logical :: powcal = .FALSE. !< If true this process runs the powcal step
     logical :: bcasthdsgen = .FALSE. !< If true then hdsgen expects data as broadcasts
     logical :: bcastpowcal = .FALSE. !< If true then powcal expects data as broadcasts
     logical :: readpowcalcontrol !< toggle to read the powcal pcontrol data
     logical :: bcastpowcalcontrol !< toggle to broadcast the powcal pcontrol data
     logical :: readpowcaldata !< toggle to read the powcal data
     logical :: bcastpowcaldata !< toggle to broadcast the powcal data
     logical :: bcastflag !< toggles powcal broadcasts
  end type decomp_t

! Global variables (yuck)
  public

  integer :: myrank_log = 0  ! Global variable to avoid changing all the log calls
  integer :: nproc_log = 1 ! global to keep log output in order for serial runs

end module smitermpi_h
