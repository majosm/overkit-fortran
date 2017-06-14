! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkGlobal

#ifdef f2003
  use, intrinsic iso_fortran_env, only : INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
#endif

  implicit none

  private

  ! API
  public :: ovk_rk
  public :: ovk_lk
  public :: ovk_bk
  public :: ovkCaseID
  public :: OVK_DEBUG
  public :: OVK_VERBOSE
  public :: OVK_NO_ERROR, OVK_IO_ERROR
  public :: OVK_NO_OVERLAP_PERIODIC, OVK_OVERLAP_PERIODIC
  public :: OVK_LITTLE_ENDIAN, OVK_BIG_ENDIAN
  public :: OVK_ALL_GRIDS
  public :: OVK_CONNECTION_NONE, OVK_CONNECTION_FRINGE, OVK_CONNECTION_FULL_GRID
  public :: OVK_INTERP_LINEAR, OVK_INTERP_CUBIC

  ! Internal
  public :: rk
  public :: lk
  public :: bk
  public :: INPUT_UNIT, OUTPUT_UNIT, ERROR_UNIT
  public :: MAX_ND
  public :: PATH_LENGTH
  public :: STRING_LENGTH
  public :: IntToString
  public :: LargeIntToString
  public :: TupleToString
  public :: CoordsToString

  integer, parameter :: ovk_rk = selected_real_kind(15, 307)
  integer, parameter :: ovk_lk = selected_int_kind(18)
  integer, parameter :: ovk_bk = selected_int_kind(1)

  character(len=256) :: ovkCaseID = ""

#ifdef OVERKIT_DEBUG
  logical, parameter :: OVK_DEBUG = .true.
#else
  logical, parameter :: OVK_DEBUG = .false.
#endif

#ifdef OVERKIT_VERBOSE
  logical, parameter :: OVK_VERBOSE = .true.
#else
  logical, parameter :: OVK_VERBOSE = .false.
#endif

  integer, parameter :: OVK_NO_ERROR = 0
  integer, parameter :: OVK_IO_ERROR = 1

  integer, parameter :: OVK_NO_OVERLAP_PERIODIC = 1
  integer, parameter :: OVK_OVERLAP_PERIODIC = 2

  integer, parameter :: OVK_LITTLE_ENDIAN = 1
  integer, parameter :: OVK_BIG_ENDIAN = 2

  integer, parameter :: OVK_ALL_GRIDS = -1

  integer, parameter :: OVK_CONNECTION_NONE = 0
  integer, parameter :: OVK_CONNECTION_FRINGE = 1
  integer, parameter :: OVK_CONNECTION_FULL_GRID = 2

  integer, parameter :: OVK_INTERP_LINEAR = 1
  integer, parameter :: OVK_INTERP_CUBIC = 2

#ifndef f2003
  integer, parameter :: INPUT_UNIT = 5
  integer, parameter :: OUTPUT_UNIT = 6
  integer, parameter :: ERROR_UNIT = 0
#endif

  integer, parameter :: rk = ovk_rk
  integer, parameter :: lk = ovk_lk
  integer, parameter :: bk = ovk_bk

  integer, parameter :: MAX_ND = 3

  integer, parameter :: PATH_LENGTH = 256

  ! Length of strings used in ToString function
  integer, parameter :: STRING_LENGTH = 256

contains

  pure elemental function IntToString(N) result(NString)

    integer, intent(in) :: N
    character(len=STRING_LENGTH) :: NString
  
    NString = LargeIntToString(int(N,kind=lk))

  end function IntToString

  pure elemental function LargeIntToString(N) result(NString)

    integer(lk), intent(in) :: N
    character(len=STRING_LENGTH) :: NString

    integer :: i, j
    integer :: NumDigits
    integer :: NumBeforeComma
    character(len=STRING_LENGTH) :: UnformattedNString

    write (UnformattedNString, '(i0)') N

    NumDigits = len_trim(UnformattedNString)
    NumBeforeComma = modulo(NumDigits-1, 3) + 1

    NString(:NumBeforeComma) = UnformattedNString(:NumBeforeComma)

    j = NumBeforeComma + 1
    do i = NumBeforeComma + 1, NumDigits
      if (modulo(i-NumBeforeComma-1, 3) == 0) then
        NString(j:j) = ','
        j = j + 1
      end if
      NString(j:j) = UnformattedNString(i:i)
      j = j + 1
    end do

    NString(j:) = ''

  end function LargeIntToString

  function TupleToString(Tuple) result(TupleString)

    integer, dimension(:), intent(in) :: Tuple
    character(len=STRING_LENGTH) :: TupleString

    integer :: i
    character(len=STRING_LENGTH), dimension(size(Tuple)) :: SeparatorStrings
    character(len=STRING_LENGTH) :: EntryString

    SeparatorStrings(:size(Tuple)-1) = ","
    SeparatorStrings(size(Tuple):) = ""

    TupleString = "("
    do i = 1, size(Tuple)
      write (EntryString, '(i0)') Tuple(i)
      TupleString = trim(TupleString) // trim(EntryString) // trim(SeparatorStrings(i))
    end do
    TupleString = trim(TupleString) // ")"

  end function TupleToString

  pure function CoordsToString(Coords) result(CoordsString)

    real(rk), dimension(:), intent(in) :: Coords
    character(len=STRING_LENGTH) :: CoordsString

    integer :: i
    character(len=STRING_LENGTH), dimension(size(Coords)) :: SeparatorStrings
    character(len=STRING_LENGTH) :: EntryString

    SeparatorStrings(:size(Coords)-1) = ","
    SeparatorStrings(size(Coords):) = ""

    CoordsString = "("
    do i = 1, size(Coords)
      write (EntryString, '(f10.4)') Coords(i)
      CoordsString = trim(CoordsString) // trim(EntryString) // trim(SeparatorStrings(i))
    end do
    CoordsString = trim(CoordsString) // ")"

  end function CoordsToString

end module ovkGlobal
