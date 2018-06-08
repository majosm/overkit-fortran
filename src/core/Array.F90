! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkArray

  use ovkGlobal
  implicit none

  private

  ! API
  public :: ovk_array_int
  public :: ovk_array_int_
  public :: ovk_array_large_int
  public :: ovk_array_large_int_
  public :: ovk_array_real
  public :: ovk_array_real_
  public :: ovk_array_logical
  public :: ovk_array_logical_
  public :: operator (==)
  public :: operator (/=)

  type ovk_array_int
    type(t_noconstruct) :: noconstruct
    integer(lk) :: n
    integer, dimension(:), allocatable :: values
  end type ovk_array_int

  type ovk_array_large_int
    type(t_noconstruct) :: noconstruct
    integer(lk) :: n
    integer(lk), dimension(:), allocatable :: values
  end type ovk_array_large_int

  type ovk_array_real
    type(t_noconstruct) :: noconstruct
    integer(lk) :: n
    real(rk), dimension(:), allocatable :: values
  end type ovk_array_real

  type ovk_array_logical
    type(t_noconstruct) :: noconstruct
    integer(lk) :: n
    logical(bk), dimension(:), allocatable :: values
  end type ovk_array_logical

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_array_int_
    module procedure ovk_array_int_Default
    module procedure ovk_array_int_Assigned_NoValues
    module procedure ovk_array_int_Assigned_Values_Scalar
    module procedure ovk_array_int_Assigned_Values_Array
  end interface ovk_array_int_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_array_large_int_
    module procedure ovk_array_large_int_Default
    module procedure ovk_array_large_int_Assigned_NoValues
    module procedure ovk_array_large_int_Assigned_Values_Scalar
    module procedure ovk_array_large_int_Assigned_Values_Array
  end interface ovk_array_large_int_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_array_real_
    module procedure ovk_array_real_Default
    module procedure ovk_array_real_Assigned_NoValues
    module procedure ovk_array_real_Assigned_Values_Scalar
    module procedure ovk_array_real_Assigned_Values_Array
  end interface ovk_array_real_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_array_logical_
    module procedure ovk_array_logical_Default
    module procedure ovk_array_logical_Assigned_NoValues
    module procedure ovk_array_logical_Assigned_Values_Scalar
    module procedure ovk_array_logical_Assigned_Values_Array
    module procedure ovk_array_logical_Assigned_Values_1Byte_Scalar
    module procedure ovk_array_logical_Assigned_Values_1Byte_Array
  end interface ovk_array_logical_

  interface operator (==)
    module procedure ovk_array_int_Equal
    module procedure ovk_array_large_int_Equal
    module procedure ovk_array_real_Equal
    module procedure ovk_array_logical_Equal
  end interface operator (==)

  interface operator (/=)
    module procedure ovk_array_int_NotEqual
    module procedure ovk_array_large_int_NotEqual
    module procedure ovk_array_real_NotEqual
    module procedure ovk_array_logical_NotEqual
  end interface operator (/=)

contains

  pure function ovk_array_int_Default() result(Array)

    type(ovk_array_int) :: Array

    Array%n = 0

  end function ovk_array_int_Default

  pure function ovk_array_large_int_Default() result(Array)

    type(ovk_array_large_int) :: Array

    Array%n = 0

  end function ovk_array_large_int_Default

  pure function ovk_array_real_Default() result(Array)

    type(ovk_array_real) :: Array

    Array%n = 0

  end function ovk_array_real_Default

  pure function ovk_array_logical_Default() result(Array)

    type(ovk_array_logical) :: Array

    Array%n = 0

  end function ovk_array_logical_Default

  pure function ovk_array_int_Assigned_NoValues(NumElements) result(Array)

    integer(lk), intent(in) :: NumElements
    type(ovk_array_int) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

  end function ovk_array_int_Assigned_NoValues

  pure function ovk_array_large_int_Assigned_NoValues(NumElements) result(Array)

    integer(lk), intent(in) :: NumElements
    type(ovk_array_large_int) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

  end function ovk_array_large_int_Assigned_NoValues

  pure function ovk_array_real_Assigned_NoValues(NumElements) result(Array)

    integer(lk), intent(in) :: NumElements
    type(ovk_array_real) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

  end function ovk_array_real_Assigned_NoValues

  pure function ovk_array_logical_Assigned_NoValues(NumElements) result(Array)

    integer(lk), intent(in) :: NumElements
    type(ovk_array_logical) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

  end function ovk_array_logical_Assigned_NoValues

  pure function ovk_array_int_Assigned_Values_Scalar(NumElements, Value) result(Array)

    integer(lk), intent(in) :: NumElements
    integer, intent(in) :: Value
    type(ovk_array_int) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Value

  end function ovk_array_int_Assigned_Values_Scalar

  pure function ovk_array_large_int_Assigned_Values_Scalar(NumElements, Value) result(Array)

    integer(lk), intent(in) :: NumElements
    integer(lk), intent(in) :: Value
    type(ovk_array_large_int) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Value

  end function ovk_array_large_int_Assigned_Values_Scalar

  pure function ovk_array_real_Assigned_Values_Scalar(NumElements, Value) result(Array)

    integer(lk), intent(in) :: NumElements
    real(rk), intent(in) :: Value
    type(ovk_array_real) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Value

  end function ovk_array_real_Assigned_Values_Scalar

  pure function ovk_array_logical_Assigned_Values_Scalar(NumElements, Value) result(Array)

    integer(lk), intent(in) :: NumElements
    logical, intent(in) :: Value
    type(ovk_array_logical) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Value

  end function ovk_array_logical_Assigned_Values_Scalar

  pure function ovk_array_logical_Assigned_Values_1Byte_Scalar(NumElements, Value) result(Array)

    integer(lk), intent(in) :: NumElements
    logical(bk), intent(in) :: Value
    type(ovk_array_logical) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Value

  end function ovk_array_logical_Assigned_Values_1Byte_Scalar

  pure function ovk_array_int_Assigned_Values_Array(NumElements, Values) result(Array)

    integer(lk), intent(in) :: NumElements
    integer, dimension(:), intent(in) :: Values
    type(ovk_array_int) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Values

  end function ovk_array_int_Assigned_Values_Array

  pure function ovk_array_large_int_Assigned_Values_Array(NumElements, Values) result(Array)

    integer(lk), intent(in) :: NumElements
    integer(lk), dimension(:), intent(in) :: Values
    type(ovk_array_large_int) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Values

  end function ovk_array_large_int_Assigned_Values_Array

  pure function ovk_array_real_Assigned_Values_Array(NumElements, Values) result(Array)

    integer(lk), intent(in) :: NumElements
    real(rk), dimension(:), intent(in) :: Values
    type(ovk_array_real) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Values

  end function ovk_array_real_Assigned_Values_Array

  pure function ovk_array_logical_Assigned_Values_Array(NumElements, Values) result(Array)

    integer(lk), intent(in) :: NumElements
    logical, dimension(:), intent(in) :: Values
    type(ovk_array_logical) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Values

  end function ovk_array_logical_Assigned_Values_Array

  pure function ovk_array_logical_Assigned_Values_1Byte_Array(NumElements, Values) result(Array)

    integer(lk), intent(in) :: NumElements
    logical(bk), dimension(:), intent(in) :: Values
    type(ovk_array_logical) :: Array

    Array%n = NumElements

    allocate(Array%values(Array%n))

    Array%values = Values

  end function ovk_array_logical_Assigned_Values_1Byte_Array

  pure function ovk_array_int_Equal(LeftArray, RightArray) result(Equal)

    type(ovk_array_int), intent(in) :: LeftArray, RightArray
    logical :: Equal

    Equal = LeftArray%n == RightArray%n

    if (Equal .and. LeftArray%n > 0) then
      Equal = Equal .and. all(LeftArray%values == RightArray%values)
    end if

  end function ovk_array_int_Equal

  pure function ovk_array_large_int_Equal(LeftArray, RightArray) result(Equal)

    type(ovk_array_large_int), intent(in) :: LeftArray, RightArray
    logical :: Equal

    Equal = LeftArray%n == RightArray%n

    if (Equal .and. LeftArray%n > 0) then
      Equal = Equal .and. all(LeftArray%values == RightArray%values)
    end if

  end function ovk_array_large_int_Equal

  pure function ovk_array_real_Equal(LeftArray, RightArray) result(Equal)

    type(ovk_array_real), intent(in) :: LeftArray, RightArray
    logical :: Equal

    Equal = LeftArray%n == RightArray%n

    if (Equal .and. LeftArray%n > 0) then
      Equal = Equal .and. all(LeftArray%values == RightArray%values)
    end if

  end function ovk_array_real_Equal

  pure function ovk_array_logical_Equal(LeftArray, RightArray) result(Equal)

    type(ovk_array_logical), intent(in) :: LeftArray, RightArray
    logical :: Equal

    Equal = LeftArray%n == RightArray%n

    if (Equal .and. LeftArray%n > 0) then
      Equal = Equal .and. all(LeftArray%values .eqv. RightArray%values)
    end if

  end function ovk_array_logical_Equal

  pure function ovk_array_int_NotEqual(LeftArray, RightArray) result(NotEqual)

    type(ovk_array_int), intent(in) :: LeftArray, RightArray
    logical :: NotEqual

    NotEqual = .not. ovk_array_int_Equal(LeftArray, RightArray)

  end function ovk_array_int_NotEqual

  pure function ovk_array_large_int_NotEqual(LeftArray, RightArray) result(NotEqual)

    type(ovk_array_large_int), intent(in) :: LeftArray, RightArray
    logical :: NotEqual

    NotEqual = .not. ovk_array_large_int_Equal(LeftArray, RightArray)

  end function ovk_array_large_int_NotEqual

  pure function ovk_array_real_NotEqual(LeftArray, RightArray) result(NotEqual)

    type(ovk_array_real), intent(in) :: LeftArray, RightArray
    logical :: NotEqual

    NotEqual = .not. ovk_array_real_Equal(LeftArray, RightArray)

  end function ovk_array_real_NotEqual

  pure function ovk_array_logical_NotEqual(LeftArray, RightArray) result(NotEqual)

    type(ovk_array_logical), intent(in) :: LeftArray, RightArray
    logical :: NotEqual

    NotEqual = .not. ovk_array_logical_Equal(LeftArray, RightArray)

  end function ovk_array_logical_NotEqual

end module ovkArray
