! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkField

  use ovkCart
  use ovkGlobal
  implicit none

  private

  ! API
  public :: ovk_field_int
  public :: ovk_field_int_
  public :: ovk_field_large_int
  public :: ovk_field_large_int_
  public :: ovk_field_real
  public :: ovk_field_real_
  public :: ovk_field_logical
  public :: ovk_field_logical_
  public :: operator (==)
  public :: operator (/=)
  public :: ovkExportField
  public :: ovkPrintField

  type ovk_field_int
    type(t_noconstruct) :: noconstruct
    type(ovk_cart) :: cart
    integer, dimension(:,:,:), allocatable :: values
  end type ovk_field_int

  type ovk_field_large_int
    type(t_noconstruct) :: noconstruct
    type(ovk_cart) :: cart
    integer(lk), dimension(:,:,:), allocatable :: values
  end type ovk_field_large_int

  type ovk_field_real
    type(t_noconstruct) :: noconstruct
    type(ovk_cart) :: cart
    real(rk), dimension(:,:,:), allocatable :: values
  end type ovk_field_real

  type ovk_field_logical
    type(t_noconstruct) :: noconstruct
    type(ovk_cart) :: cart
    logical(bk), dimension(:,:,:), allocatable :: values
  end type ovk_field_logical

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_field_int_
    module procedure ovk_field_int_Default
    module procedure ovk_field_int_Assigned_Empty
    module procedure ovk_field_int_Assigned_NoValues
    module procedure ovk_field_int_Assigned_Values_Scalar
    module procedure ovk_field_int_Assigned_Values_Rank2
    module procedure ovk_field_int_Assigned_Values_Rank3
  end interface ovk_field_int_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_field_large_int_
    module procedure ovk_field_large_int_Default
    module procedure ovk_field_large_int_Assigned_Empty
    module procedure ovk_field_large_int_Assigned_NoValues
    module procedure ovk_field_large_int_Assigned_Values_Scalar
    module procedure ovk_field_large_int_Assigned_Values_Rank2
    module procedure ovk_field_large_int_Assigned_Values_Rank3
  end interface ovk_field_large_int_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_field_real_
    module procedure ovk_field_real_Default
    module procedure ovk_field_real_Assigned_Empty
    module procedure ovk_field_real_Assigned_NoValues
    module procedure ovk_field_real_Assigned_Values_Scalar
    module procedure ovk_field_real_Assigned_Values_Rank2
    module procedure ovk_field_real_Assigned_Values_Rank3
  end interface ovk_field_real_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_field_logical_
    module procedure ovk_field_logical_Default
    module procedure ovk_field_logical_Assigned_Empty
    module procedure ovk_field_logical_Assigned_NoValues
    module procedure ovk_field_logical_Assigned_Values_Scalar
    module procedure ovk_field_logical_Assigned_Values_Rank2
    module procedure ovk_field_logical_Assigned_Values_Rank3
    module procedure ovk_field_logical_Assigned_Values_Scalar_1Byte
    module procedure ovk_field_logical_Assigned_Values_Rank2_1Byte
    module procedure ovk_field_logical_Assigned_Values_Rank3_1Byte
  end interface ovk_field_logical_

  interface operator (==)
    module procedure ovk_field_int_Equal
    module procedure ovk_field_large_int_Equal
    module procedure ovk_field_real_Equal
    module procedure ovk_field_logical_Equal
  end interface operator (==)

  interface operator (/=)
    module procedure ovk_field_int_NotEqual
    module procedure ovk_field_large_int_NotEqual
    module procedure ovk_field_real_NotEqual
    module procedure ovk_field_logical_NotEqual
  end interface operator (/=)

  interface ovkExportField
    module procedure ovkExportField_Integer
    module procedure ovkExportField_LargeInteger
    module procedure ovkExportField_Real
    module procedure ovkExportField_Logical
  end interface ovkExportField

  interface ovkPrintField
    module procedure ovkPrintField_Integer
    module procedure ovkPrintField_LargeInteger
    module procedure ovkPrintField_Real
    module procedure ovkPrintField_Logical
  end interface ovkPrintField

contains

  pure function ovk_field_int_Default() result(Field)

    type(ovk_field_int) :: Field

    Field%cart = ovk_cart_()

  end function ovk_field_int_Default

  pure function ovk_field_large_int_Default() result(Field)

    type(ovk_field_large_int) :: Field

    Field%cart = ovk_cart_()

  end function ovk_field_large_int_Default

  pure function ovk_field_real_Default() result(Field)

    type(ovk_field_real) :: Field

    Field%cart = ovk_cart_()

  end function ovk_field_real_Default

  pure function ovk_field_logical_Default() result(Field)

    type(ovk_field_logical) :: Field

    Field%cart = ovk_cart_()

  end function ovk_field_logical_Default

  pure function ovk_field_int_Assigned_Empty(NumDims) result(Field)

    integer, intent(in) :: NumDims
    type(ovk_field_int) :: Field

    Field%cart = ovk_cart_(NumDims)

  end function ovk_field_int_Assigned_Empty

  pure function ovk_field_large_int_Assigned_Empty(NumDims) result(Field)

    integer, intent(in) :: NumDims
    type(ovk_field_large_int) :: Field

    Field%cart = ovk_cart_(NumDims)

  end function ovk_field_large_int_Assigned_Empty

  pure function ovk_field_real_Assigned_Empty(NumDims) result(Field)

    integer, intent(in) :: NumDims
    type(ovk_field_real) :: Field

    Field%cart = ovk_cart_(NumDims)

  end function ovk_field_real_Assigned_Empty

  pure function ovk_field_logical_Assigned_Empty(NumDims) result(Field)

    integer, intent(in) :: NumDims
    type(ovk_field_logical) :: Field

    Field%cart = ovk_cart_(NumDims)

  end function ovk_field_logical_Assigned_Empty

  pure function ovk_field_int_Assigned_NoValues(Cart) result(Field)

    type(ovk_cart), intent(in) :: Cart
    type(ovk_field_int) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

  end function ovk_field_int_Assigned_NoValues

  pure function ovk_field_large_int_Assigned_NoValues(Cart) result(Field)

    type(ovk_cart), intent(in) :: Cart
    type(ovk_field_large_int) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

  end function ovk_field_large_int_Assigned_NoValues

  pure function ovk_field_real_Assigned_NoValues(Cart) result(Field)

    type(ovk_cart), intent(in) :: Cart
    type(ovk_field_real) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

  end function ovk_field_real_Assigned_NoValues

  pure function ovk_field_logical_Assigned_NoValues(Cart) result(Field)

    type(ovk_cart), intent(in) :: Cart
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

  end function ovk_field_logical_Assigned_NoValues

  pure function ovk_field_int_Assigned_Values_Scalar(Cart, Value) result(Field)

    type(ovk_cart), intent(in) :: Cart
    integer, intent(in) :: Value
    type(ovk_field_int) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Value

  end function ovk_field_int_Assigned_Values_Scalar

  pure function ovk_field_large_int_Assigned_Values_Scalar(Cart, Value) result(Field)

    type(ovk_cart), intent(in) :: Cart
    integer(lk), intent(in) :: Value
    type(ovk_field_large_int) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Value

  end function ovk_field_large_int_Assigned_Values_Scalar

  pure function ovk_field_real_Assigned_Values_Scalar(Cart, Value) result(Field)

    type(ovk_cart), intent(in) :: Cart
    real(rk), intent(in) :: Value
    type(ovk_field_real) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Value

  end function ovk_field_real_Assigned_Values_Scalar

  pure function ovk_field_logical_Assigned_Values_Scalar(Cart, Value) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical, intent(in) :: Value
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Value

  end function ovk_field_logical_Assigned_Values_Scalar

  pure function ovk_field_logical_Assigned_Values_Scalar_1Byte(Cart, Value) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical(bk), intent(in) :: Value
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Value

  end function ovk_field_logical_Assigned_Values_Scalar_1Byte

  pure function ovk_field_int_Assigned_Values_Rank2(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    integer, dimension(:,:), intent(in) :: Values
    type(ovk_field_int) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values(:,:,1) = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1)

  end function ovk_field_int_Assigned_Values_Rank2

  pure function ovk_field_large_int_Assigned_Values_Rank2(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    integer(lk), dimension(:,:), intent(in) :: Values
    type(ovk_field_large_int) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values(:,:,1) = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1)

  end function ovk_field_large_int_Assigned_Values_Rank2

  pure function ovk_field_real_Assigned_Values_Rank2(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    real(rk), dimension(:,:), intent(in) :: Values
    type(ovk_field_real) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values(:,:,1) = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1)

  end function ovk_field_real_Assigned_Values_Rank2

  pure function ovk_field_logical_Assigned_Values_Rank2(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical, dimension(:,:), intent(in) :: Values
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values(:,:,1) = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1)

  end function ovk_field_logical_Assigned_Values_Rank2

  pure function ovk_field_logical_Assigned_Values_Rank2_1Byte(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical(bk), dimension(:,:), intent(in) :: Values
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values(:,:,1) = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1)

  end function ovk_field_logical_Assigned_Values_Rank2_1Byte

  pure function ovk_field_int_Assigned_Values_Rank3(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    integer, dimension(:,:,:), intent(in) :: Values
    type(ovk_field_int) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1,:Field%cart%ie(3)-Field%cart%is(3)+1)

  end function ovk_field_int_Assigned_Values_Rank3

  pure function ovk_field_large_int_Assigned_Values_Rank3(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    integer(lk), dimension(:,:,:), intent(in) :: Values
    type(ovk_field_large_int) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1,:Field%cart%ie(3)-Field%cart%is(3)+1)

  end function ovk_field_large_int_Assigned_Values_Rank3

  pure function ovk_field_real_Assigned_Values_Rank3(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    real(rk), dimension(:,:,:), intent(in) :: Values
    type(ovk_field_real) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1,:Field%cart%ie(3)-Field%cart%is(3)+1)

  end function ovk_field_real_Assigned_Values_Rank3

  pure function ovk_field_logical_Assigned_Values_Rank3(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical, dimension(:,:,:), intent(in) :: Values
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1,:Field%cart%ie(3)-Field%cart%is(3)+1)

  end function ovk_field_logical_Assigned_Values_Rank3

  pure function ovk_field_logical_Assigned_Values_Rank3_1Byte(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical(bk), dimension(:,:,:), intent(in) :: Values
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1,:Field%cart%ie(3)-Field%cart%is(3)+1)

  end function ovk_field_logical_Assigned_Values_Rank3_1Byte

  pure function ovk_field_int_Equal(LeftField, RightField) result(Equal)

    type(ovk_field_int), intent(in) :: LeftField, RightField
    logical :: Equal

    Equal = LeftField%cart == RightField%cart

    if (Equal .and. all(LeftField%cart%ie >= LeftField%cart%is)) then
      Equal = Equal .and. all(LeftField%values == RightField%values)
    end if

  end function ovk_field_int_Equal

  pure function ovk_field_large_int_Equal(LeftField, RightField) result(Equal)

    type(ovk_field_large_int), intent(in) :: LeftField, RightField
    logical :: Equal

    Equal = LeftField%cart == RightField%cart

    if (Equal .and. all(LeftField%cart%ie >= LeftField%cart%is)) then
      Equal = Equal .and. all(LeftField%values == RightField%values)
    end if

  end function ovk_field_large_int_Equal

  pure function ovk_field_real_Equal(LeftField, RightField) result(Equal)

    type(ovk_field_real), intent(in) :: LeftField, RightField
    logical :: Equal

    Equal = LeftField%cart == RightField%cart

    if (Equal .and. all(LeftField%cart%ie >= LeftField%cart%is)) then
      Equal = Equal .and. all(LeftField%values == RightField%values)
    end if

  end function ovk_field_real_Equal

  pure function ovk_field_logical_Equal(LeftField, RightField) result(Equal)

    type(ovk_field_logical), intent(in) :: LeftField, RightField
    logical :: Equal

    Equal = LeftField%cart == RightField%cart

    if (Equal .and. all(LeftField%cart%ie >= LeftField%cart%is)) then
      Equal = Equal .and. all(LeftField%values .eqv. RightField%values)
    end if

  end function ovk_field_logical_Equal

  pure function ovk_field_int_NotEqual(LeftField, RightField) result(NotEqual)

    type(ovk_field_int), intent(in) :: LeftField, RightField
    logical :: NotEqual

    NotEqual = .not. ovk_field_int_Equal(LeftField, RightField)

  end function ovk_field_int_NotEqual

  pure function ovk_field_large_int_NotEqual(LeftField, RightField) result(NotEqual)

    type(ovk_field_large_int), intent(in) :: LeftField, RightField
    logical :: NotEqual

    NotEqual = .not. ovk_field_large_int_Equal(LeftField, RightField)

  end function ovk_field_large_int_NotEqual

  pure function ovk_field_real_NotEqual(LeftField, RightField) result(NotEqual)

    type(ovk_field_real), intent(in) :: LeftField, RightField
    logical :: NotEqual

    NotEqual = .not. ovk_field_real_Equal(LeftField, RightField)

  end function ovk_field_real_NotEqual

  pure function ovk_field_logical_NotEqual(LeftField, RightField) result(NotEqual)

    type(ovk_field_logical), intent(in) :: LeftField, RightField
    logical :: NotEqual

    NotEqual = .not. ovk_field_logical_Equal(LeftField, RightField)

  end function ovk_field_logical_NotEqual

  subroutine ovkExportField_Integer(Field, ExportCart, ExportedField)

    type(ovk_field_int), intent(in) :: Field
    type(ovk_cart), intent(in) :: ExportCart
    type(ovk_field_int), intent(out) :: ExportedField

    integer :: i, j, k
    integer, dimension(MAX_ND) :: Point

    if (OVK_DEBUG) then
      if (.not. ovkCartIsCompatible(ExportCart, Field%cart)) then
        write (ERROR_UNIT, '(a)') "ERROR: Export cart is incompatible with field."
        stop 1
      end if
    end if

    ExportedField = ovk_field_int_(ExportCart)

    do k = ExportCart%is(3), ExportCart%ie(3)
      do j = ExportCart%is(2), ExportCart%ie(2)
        do i = ExportCart%is(1), ExportCart%ie(1)
          Point = [i,j,k]
          Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
          ExportedField%values(i,j,k) = Field%values(Point(1),Point(2),Point(3))
        end do
      end do
    end do

  end subroutine ovkExportField_Integer

  subroutine ovkExportField_LargeInteger(Field, ExportCart, ExportedField)

    type(ovk_field_large_int), intent(in) :: Field
    type(ovk_cart), intent(in) :: ExportCart
    type(ovk_field_large_int), intent(out) :: ExportedField

    integer :: i, j, k
    integer, dimension(MAX_ND) :: Point

    if (OVK_DEBUG) then
      if (.not. ovkCartIsCompatible(ExportCart, Field%cart)) then
        write (ERROR_UNIT, '(a)') "ERROR: Export cart is incompatible with field."
        stop 1
      end if
    end if

    ExportedField = ovk_field_large_int_(ExportCart)

    do k = ExportCart%is(3), ExportCart%ie(3)
      do j = ExportCart%is(2), ExportCart%ie(2)
        do i = ExportCart%is(1), ExportCart%ie(1)
          Point = [i,j,k]
          Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
          ExportedField%values(i,j,k) = Field%values(Point(1),Point(2),Point(3))
        end do
      end do
    end do

  end subroutine ovkExportField_LargeInteger

  subroutine ovkExportField_Real(Field, ExportCart, ExportedField)

    type(ovk_field_real), intent(in) :: Field
    type(ovk_cart), intent(in) :: ExportCart
    type(ovk_field_real), intent(out) :: ExportedField

    integer :: i, j, k
    integer, dimension(MAX_ND) :: Point

    if (OVK_DEBUG) then
      if (.not. ovkCartIsCompatible(ExportCart, Field%cart)) then
        write (ERROR_UNIT, '(a)') "ERROR: Export cart is incompatible with field."
        stop 1
      end if
    end if

    ExportedField = ovk_field_real_(ExportCart)

    do k = ExportCart%is(3), ExportCart%ie(3)
      do j = ExportCart%is(2), ExportCart%ie(2)
        do i = ExportCart%is(1), ExportCart%ie(1)
          Point = [i,j,k]
          Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
          ExportedField%values(i,j,k) = Field%values(Point(1),Point(2),Point(3))
        end do
      end do
    end do

  end subroutine ovkExportField_Real

  subroutine ovkExportField_Logical(Field, ExportCart, ExportedField)

    type(ovk_field_logical), intent(in) :: Field
    type(ovk_cart), intent(in) :: ExportCart
    type(ovk_field_logical), intent(out) :: ExportedField

    integer :: i, j, k
    integer, dimension(MAX_ND) :: Point

    if (OVK_DEBUG) then
      if (.not. ovkCartIsCompatible(ExportCart, Field%cart)) then
        write (ERROR_UNIT, '(a)') "ERROR: Export cart is incompatible with field."
        stop 1
      end if
    end if

    ExportedField = ovk_field_logical_(ExportCart)

    do k = ExportCart%is(3), ExportCart%ie(3)
      do j = ExportCart%is(2), ExportCart%ie(2)
        do i = ExportCart%is(1), ExportCart%ie(1)
          Point = [i,j,k]
          Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
          ExportedField%values(i,j,k) = Field%values(Point(1),Point(2),Point(3))
        end do
      end do
    end do

  end subroutine ovkExportField_Logical

  subroutine ovkPrintField_Integer(Field, OutputUnit, StartIndex, EndIndex, NumDigits)

    type(ovk_field_int), intent(in) :: Field
    integer, intent(in) :: OutputUnit
    integer, dimension(Field%cart%nd), intent(in), optional :: StartIndex, EndIndex
    integer, intent(in), optional :: NumDigits

    integer, dimension(MAX_ND) :: StartIndex_, EndIndex_
    integer :: NumDigits_
    integer :: i, j
    integer :: IDir, JDir, KDir
    integer :: IStart, IEnd, JStart, JEnd
    integer :: KSlice
    character(len=16) :: FormatString
    integer, dimension(MAX_ND) :: Point

    if (OVK_DEBUG) then
      if (Field%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    if (present(StartIndex)) then
      StartIndex_(:Field%cart%nd) = StartIndex
      StartIndex_(Field%cart%nd+1:) = 1
    else
      StartIndex_ = Field%cart%is
    end if

    if (present(EndIndex)) then
      EndIndex_(:Field%cart%nd) = EndIndex
      EndIndex_(Field%cart%nd+1:) = 1
    else
      EndIndex_ = Field%cart%ie
      EndIndex_(MAX_ND) = Field%cart%is(MAX_ND)
    end if

    if (present(NumDigits)) then
      NumDigits_ = NumDigits
    else
      NumDigits_ = 4
    end if

    call FindSlice(StartIndex_, EndIndex_, IDir, JDir, KDir, IStart, IEnd, JStart, JEnd, KSlice)

    write (FormatString, '(a,i0,a)') "(a,i", NumDigits_, ",a)"

    call PrintRowHeaderFooter(OutputUnit, IStart, IEnd, NumDigits_)

    do j = JEnd, JStart, -1

      write (OutputUnit, '(a)', advance='no') " | "

      do i = IStart, IEnd
        Point(IDir) = i
        Point(JDir) = j
        Point(KDir) = KSlice
        Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
        write (OutputUnit, FormatString, advance='no') " ", Field%values(Point(1),Point(2), &
          Point(3)), " "
      end do

      write (OutputUnit, '(a)') " | "

    end do

    call PrintRowHeaderFooter(OutputUnit, IStart, IEnd, NumDigits_)

    write (OutputUnit, '(a)') ""

  end subroutine ovkPrintField_Integer

  subroutine ovkPrintField_LargeInteger(Field, OutputUnit, StartIndex, EndIndex, NumDigits)

    type(ovk_field_large_int), intent(in) :: Field
    integer, intent(in) :: OutputUnit
    integer, dimension(Field%cart%nd), intent(in), optional :: StartIndex, EndIndex
    integer, intent(in), optional :: NumDigits

    integer, dimension(MAX_ND) :: StartIndex_, EndIndex_
    integer :: NumDigits_
    integer :: i, j
    integer :: IDir, JDir, KDir
    integer :: IStart, IEnd, JStart, JEnd
    integer :: KSlice
    character(len=16) :: FormatString
    integer, dimension(MAX_ND) :: Point

    if (OVK_DEBUG) then
      if (Field%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    if (present(StartIndex)) then
      StartIndex_(:Field%cart%nd) = StartIndex
      StartIndex_(Field%cart%nd+1:) = 1
    else
      StartIndex_ = Field%cart%is
    end if

    if (present(EndIndex)) then
      EndIndex_(:Field%cart%nd) = EndIndex
      EndIndex_(Field%cart%nd+1:) = 1
    else
      EndIndex_ = Field%cart%ie
      EndIndex_(MAX_ND) = Field%cart%is(MAX_ND)
    end if

    if (present(NumDigits)) then
      NumDigits_ = NumDigits
    else
      NumDigits_ = 4
    end if

    call FindSlice(StartIndex_, EndIndex_, IDir, JDir, KDir, IStart, IEnd, JStart, JEnd, KSlice)

    write (FormatString, '(a,i0,a)') "(a,i", NumDigits_, ",a)"

    call PrintRowHeaderFooter(OutputUnit, IStart, IEnd, NumDigits_)

    do j = JEnd, JStart, -1

      write (OutputUnit, '(a)', advance='no') " | "

      do i = IStart, IEnd
        Point(IDir) = i
        Point(JDir) = j
        Point(KDir) = KSlice
        Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
        write (OutputUnit, FormatString, advance='no') " ", Field%values(Point(1),Point(2), &
          Point(3)), " "
      end do

      write (OutputUnit, '(a)') " | "

    end do

    call PrintRowHeaderFooter(OutputUnit, IStart, IEnd, NumDigits_)

    write (OutputUnit, '(a)') ""

  end subroutine ovkPrintField_LargeInteger

  subroutine ovkPrintField_Real(Field, OutputUnit, StartIndex, EndIndex, NumDigits, &
    NumFractionalDigits, NumExponentDigits)

    type(ovk_field_real), intent(in) :: Field
    integer, intent(in) :: OutputUnit
    integer, dimension(Field%cart%nd), intent(in), optional :: StartIndex, EndIndex
    integer, intent(in), optional :: NumDigits
    integer, intent(in), optional :: NumFractionalDigits
    integer, intent(in), optional :: NumExponentDigits

    integer, dimension(MAX_ND) :: StartIndex_, EndIndex_
    integer :: NumDigits_
    integer :: NumFractionalDigits_
    integer :: NumExponentDigits_
    integer :: i, j
    integer :: IDir, JDir, KDir
    integer :: IStart, IEnd, JStart, JEnd
    integer :: KSlice
    character(len=16) :: FormatString
    character(len=256) :: TestString
    integer :: ElementLength
    integer, dimension(MAX_ND) :: Point

    if (OVK_DEBUG) then
      if (Field%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    if (present(StartIndex)) then
      StartIndex_(:Field%cart%nd) = StartIndex
      StartIndex_(Field%cart%nd+1:) = 1
    else
      StartIndex_ = Field%cart%is
    end if

    if (present(EndIndex)) then
      EndIndex_(:Field%cart%nd) = EndIndex
      EndIndex_(Field%cart%nd+1:) = 1
    else
      EndIndex_ = Field%cart%ie
      EndIndex_(MAX_ND) = Field%cart%is(MAX_ND)
    end if

    if (present(NumDigits)) then
      NumDigits_ = NumDigits
    else
      NumDigits_ = 6
    end if

    if (present(NumFractionalDigits)) then
      NumFractionalDigits_ = NumFractionalDigits
    else
      NumFractionalDigits_ = 2
    end if

    if (present(NumExponentDigits)) then
      NumExponentDigits_ = NumExponentDigits
    else
      NumExponentDigits_ = 0
    end if

    call FindSlice(StartIndex_, EndIndex_, IDir, JDir, KDir, IStart, IEnd, JStart, JEnd, KSlice)

    if (NumExponentDigits_ > 0) then
      write (FormatString, '(a,3(i0,a))') "(a,e", NumDigits_, ".",  NumFractionalDigits_, "e", &
        NumExponentDigits_, ",a)"
    else
      write (FormatString, '(a,2(i0,a))') "(a,f", NumDigits_, ".", NumFractionalDigits_, ",a)"
    end if

    write (TestString, FormatString) " ", 0._rk, " "
    ElementLength = len_trim(TestString)-1

    call PrintRowHeaderFooter(OutputUnit, IStart, IEnd, ElementLength)

    do j = JEnd, JStart, -1

      write (OutputUnit, '(a)', advance='no') " | "

      do i = IStart, IEnd
        Point(IDir) = i
        Point(JDir) = j
        Point(KDir) = KSlice
        Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
        write (OutputUnit, FormatString, advance='no') " ", Field%values(Point(1),Point(2), &
          Point(3)), " "
      end do

      write (OutputUnit, '(a)') " | "

    end do

    call PrintRowHeaderFooter(OutputUnit, IStart, IEnd, ElementLength)

    write (OutputUnit, '(a)') ""

  end subroutine ovkPrintField_Real

  subroutine ovkPrintField_Logical(Field, OutputUnit, StartIndex, EndIndex)

    type(ovk_field_logical), intent(in) :: Field
    integer, intent(in) :: OutputUnit
    integer, dimension(Field%cart%nd), intent(in), optional :: StartIndex, EndIndex

    integer, dimension(MAX_ND) :: StartIndex_, EndIndex_
    integer :: i, j
    integer :: IDir, JDir, KDir
    integer :: IStart, IEnd, JStart, JEnd
    integer :: KSlice
    integer, dimension(MAX_ND) :: Point

    if (OVK_DEBUG) then
      if (Field%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    if (present(StartIndex)) then
      StartIndex_(:Field%cart%nd) = StartIndex
      StartIndex_(Field%cart%nd+1:) = 1
    else
      StartIndex_ = Field%cart%is
    end if

    if (present(EndIndex)) then
      EndIndex_(:Field%cart%nd) = EndIndex
      EndIndex_(Field%cart%nd+1:) = 1
    else
      EndIndex_ = Field%cart%ie
      EndIndex_(MAX_ND) = Field%cart%is(MAX_ND)
    end if

    call FindSlice(StartIndex_, EndIndex_, IDir, JDir, KDir, IStart, IEnd, JStart, JEnd, KSlice)

    call PrintRowHeaderFooter(OutputUnit, IStart, IEnd, 1)

    do j = JEnd, JStart, -1

      write (OutputUnit, '(a)', advance='no') " | "

      do i = IStart, IEnd
        Point(IDir) = i
        Point(JDir) = j
        Point(KDir) = KSlice
        Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
        write (OutputUnit, '(a)', advance='no') merge(" X ", "   ", Field%values(Point(1), &
          Point(2),Point(3)))
      end do

      write (OutputUnit, '(a)') " | "

    end do

    call PrintRowHeaderFooter(OutputUnit, IStart, IEnd, 1)

    write (OutputUnit, '(a)') ""

  end subroutine ovkPrintField_Logical

  subroutine FindSlice(StartIndex, EndIndex, IDir, JDir, KDir, IStart, IEnd, JStart, JEnd, KSlice)

    integer, dimension(MAX_ND), intent(in) :: StartIndex, EndIndex
    integer, intent(out) :: IDir, JDir, KDir
    integer, intent(out) :: IStart, IEnd, JStart, JEnd
    integer, intent(out) :: KSlice

    integer :: d

    do d = 1, MAX_ND
      if (EndIndex(d) == StartIndex(d)) then
        exit
      end if
    end do
    KDir = d

    if (OVK_DEBUG) then
      if (KDir > MAX_ND) then
        write (ERROR_UNIT, '(a)') "ERROR: Range to print must be two-dimensional."
        stop 1
      end if
    end if

    IDir = modulo(KDir,MAX_ND)+1
    JDir = modulo(KDir+1,MAX_ND)+1

    IStart = StartIndex(IDir)
    IEnd = EndIndex(IDir)

    JStart = StartIndex(JDir)
    JEnd = EndIndex(JDir)

    KSlice = StartIndex(KDir)

  end subroutine FindSlice

  subroutine PrintRowHeaderFooter(OutputUnit, IStart, IEnd, ElementLength)

    integer, intent(in) :: OutputUnit
    integer, intent(in) :: IStart, IEnd
    integer, intent(in) :: ElementLength

    integer :: i, l

    write (OutputUnit, '(a)', advance='no') "   "
    do i = IStart, IEnd
      write (OutputUnit, '(a)', advance='no') " "
      do l = 1, ElementLength
        write (OutputUnit, '(a)', advance='no') "-"
      end do
      write (OutputUnit, '(a)', advance='no') " "
    end do
    write (OutputUnit, '(a)') "   "

  end subroutine PrintRowHeaderFooter

end module ovkField
