module bods_h

  use const_kind_m

  type, public :: bods_t
   !! bods data type
     integer(ki4), dimension(:), allocatable :: list !< list of bods
     integer(ki4)  :: nbod !< no of bods
     integer(ki4), dimension(:), allocatable :: indx !< index of bods
     integer(ki4)  :: nindx !< number in index
     integer(ki4)  :: maxindx !< allocated size of index array
     integer(ki4)  :: maxbodsf !< maximum number of bods in file
  end type bods_t

end module bods_h
