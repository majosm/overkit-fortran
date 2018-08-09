! Copyright (c) 2018 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkField

  use ovkCart
  use ovkGlobal
  implicit none

  private

  ! Public API
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
  public :: ovkFieldPeriodicFill
  public :: ovkGetFieldPatch
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
    module procedure ovk_field_logical_Assigned_Values_1Byte_Scalar
    module procedure ovk_field_logical_Assigned_Values_1Byte_Rank2
    module procedure ovk_field_logical_Assigned_Values_1Byte_Rank3
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

  interface ovkFieldPeriodicFill
    module procedure ovkFieldPeriodicFill_Integer
    module procedure ovkFieldPeriodicFill_LargeInteger
    module procedure ovkFieldPeriodicFill_Real
    module procedure ovkFieldPeriodicFill_Logical
  end interface ovkFieldPeriodicFill

  interface ovkGetFieldPatch
    module procedure ovkGetFieldPatch_Integer_Rank2
    module procedure ovkGetFieldPatch_Integer_Rank3
    module procedure ovkGetFieldPatch_LargeInteger_Rank2
    module procedure ovkGetFieldPatch_LargeInteger_Rank3
    module procedure ovkGetFieldPatch_Real_Rank2
    module procedure ovkGetFieldPatch_Real_Rank3
    module procedure ovkGetFieldPatch_Logical_Rank2
    module procedure ovkGetFieldPatch_Logical_Rank3
    module procedure ovkGetFieldPatch_Logical1Byte_Rank2
    module procedure ovkGetFieldPatch_Logical1Byte_Rank3
  end interface ovkGetFieldPatch

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

  pure function ovk_field_logical_Assigned_Values_1Byte_Scalar(Cart, Value) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical(bk), intent(in) :: Value
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Value

  end function ovk_field_logical_Assigned_Values_1Byte_Scalar

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

  pure function ovk_field_logical_Assigned_Values_1Byte_Rank2(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical(bk), dimension(:,:), intent(in) :: Values
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values(:,:,1) = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1)

  end function ovk_field_logical_Assigned_Values_1Byte_Rank2

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

  pure function ovk_field_logical_Assigned_Values_1Byte_Rank3(Cart, Values) result(Field)

    type(ovk_cart), intent(in) :: Cart
    logical(bk), dimension(:,:,:), intent(in) :: Values
    type(ovk_field_logical) :: Field

    Field%cart = Cart

    allocate(Field%values(Field%cart%is(1):Field%cart%ie(1),Field%cart%is(2):Field%cart%ie(2), &
      Field%cart%is(3):Field%cart%ie(3)))

    Field%values = Values(:Field%cart%ie(1)-Field%cart%is(1)+1, &
      :Field%cart%ie(2)-Field%cart%is(2)+1,:Field%cart%ie(3)-Field%cart%is(3)+1)

  end function ovk_field_logical_Assigned_Values_1Byte_Rank3

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

  subroutine ovkFieldPeriodicFill_Integer(Field, PrincipalCart)

    type(ovk_field_int), intent(inout) :: Field
    type(ovk_cart), intent(in) :: PrincipalCart

    integer :: i, j, k
    type(ovk_cart) :: ExtendedCart
    integer, dimension(MAX_DIMS) :: Point, AdjustedPoint

    ExtendedCart = Field%cart

    if (ovkCartIsCompatible(ExtendedCart, PrincipalCart)) then
      if (ExtendedCart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        if (ExtendedCart%periodic(1)) then
          Field%values(ExtendedCart%ie(1),:,:) = Field%values(ExtendedCart%is(1),:,:)
        end if
        if (ExtendedCart%periodic(2)) then
          Field%values(:,ExtendedCart%ie(2),:) = Field%values(:,ExtendedCart%is(2),:)
        end if
        if (ExtendedCart%periodic(3)) then
          Field%values(:,:,ExtendedCart%ie(3)) = Field%values(:,:,ExtendedCart%is(3))
        end if
      end if
    else
      AdjustedPoint = 1
      do k = ExtendedCart%is(3), ExtendedCart%ie(3)
        do j = ExtendedCart%is(2), ExtendedCart%ie(2)
          do i = ExtendedCart%is(1), ExtendedCart%ie(1)
            Point = [i,j,k]
            if (.not. ovkCartContains(PrincipalCart, Point)) then
              AdjustedPoint(:PrincipalCart%nd) = ovkCartPeriodicAdjust(PrincipalCart, Point)
              if (ovkCartContains(PrincipalCart, AdjustedPoint)) then
                Field%values(Point(1),Point(2),Point(3)) = Field%values(AdjustedPoint(1), &
                  AdjustedPoint(2),AdjustedPoint(3))
              end if
            end if
          end do
        end do
      end do
    end if

  end subroutine ovkFieldPeriodicFill_Integer

  subroutine ovkFieldPeriodicFill_LargeInteger(Field, PrincipalCart)

    type(ovk_field_large_int), intent(inout) :: Field
    type(ovk_cart), intent(in) :: PrincipalCart

    integer :: i, j, k
    type(ovk_cart) :: ExtendedCart
    integer, dimension(MAX_DIMS) :: Point, AdjustedPoint

    ExtendedCart = Field%cart

    if (ovkCartIsCompatible(ExtendedCart, PrincipalCart)) then
      if (ExtendedCart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        if (ExtendedCart%periodic(1)) then
          Field%values(ExtendedCart%ie(1),:,:) = Field%values(ExtendedCart%is(1),:,:)
        end if
        if (ExtendedCart%periodic(2)) then
          Field%values(:,ExtendedCart%ie(2),:) = Field%values(:,ExtendedCart%is(2),:)
        end if
        if (ExtendedCart%periodic(3)) then
          Field%values(:,:,ExtendedCart%ie(3)) = Field%values(:,:,ExtendedCart%is(3))
        end if
      end if
    else
      AdjustedPoint = 1
      do k = ExtendedCart%is(3), ExtendedCart%ie(3)
        do j = ExtendedCart%is(2), ExtendedCart%ie(2)
          do i = ExtendedCart%is(1), ExtendedCart%ie(1)
            Point = [i,j,k]
            if (.not. ovkCartContains(PrincipalCart, Point)) then
              AdjustedPoint(:PrincipalCart%nd) = ovkCartPeriodicAdjust(PrincipalCart, Point)
              if (ovkCartContains(PrincipalCart, AdjustedPoint)) then
                Field%values(Point(1),Point(2),Point(3)) = Field%values(AdjustedPoint(1), &
                  AdjustedPoint(2),AdjustedPoint(3))
              end if
            end if
          end do
        end do
      end do
    end if

  end subroutine ovkFieldPeriodicFill_LargeInteger

  subroutine ovkFieldPeriodicFill_Real(Field, PrincipalCart)

    type(ovk_field_real), intent(inout) :: Field
    type(ovk_cart), intent(in) :: PrincipalCart

    integer :: i, j, k
    type(ovk_cart) :: ExtendedCart
    integer, dimension(MAX_DIMS) :: Point, AdjustedPoint

    ExtendedCart = Field%cart

    if (ovkCartIsCompatible(ExtendedCart, PrincipalCart)) then
      if (ExtendedCart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        if (ExtendedCart%periodic(1)) then
          Field%values(ExtendedCart%ie(1),:,:) = Field%values(ExtendedCart%is(1),:,:)
        end if
        if (ExtendedCart%periodic(2)) then
          Field%values(:,ExtendedCart%ie(2),:) = Field%values(:,ExtendedCart%is(2),:)
        end if
        if (ExtendedCart%periodic(3)) then
          Field%values(:,:,ExtendedCart%ie(3)) = Field%values(:,:,ExtendedCart%is(3))
        end if
      end if
    else
      AdjustedPoint = 1
      do k = ExtendedCart%is(3), ExtendedCart%ie(3)
        do j = ExtendedCart%is(2), ExtendedCart%ie(2)
          do i = ExtendedCart%is(1), ExtendedCart%ie(1)
            Point = [i,j,k]
            if (.not. ovkCartContains(PrincipalCart, Point)) then
              AdjustedPoint(:PrincipalCart%nd) = ovkCartPeriodicAdjust(PrincipalCart, Point)
              if (ovkCartContains(PrincipalCart, AdjustedPoint)) then
                Field%values(Point(1),Point(2),Point(3)) = Field%values(AdjustedPoint(1), &
                  AdjustedPoint(2),AdjustedPoint(3))
              end if
            end if
          end do
        end do
      end do
    end if

  end subroutine ovkFieldPeriodicFill_Real

  subroutine ovkFieldPeriodicFill_Logical(Field, PrincipalCart)

    type(ovk_field_logical), intent(inout) :: Field
    type(ovk_cart), intent(in) :: PrincipalCart

    integer :: i, j, k
    type(ovk_cart) :: ExtendedCart
    integer, dimension(MAX_DIMS) :: Point, AdjustedPoint

    ExtendedCart = Field%cart

    if (ovkCartIsCompatible(ExtendedCart, PrincipalCart)) then
      if (ExtendedCart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        if (ExtendedCart%periodic(1)) then
          Field%values(ExtendedCart%ie(1),:,:) = Field%values(ExtendedCart%is(1),:,:)
        end if
        if (ExtendedCart%periodic(2)) then
          Field%values(:,ExtendedCart%ie(2),:) = Field%values(:,ExtendedCart%is(2),:)
        end if
        if (ExtendedCart%periodic(3)) then
          Field%values(:,:,ExtendedCart%ie(3)) = Field%values(:,:,ExtendedCart%is(3))
        end if
      end if
    else
      AdjustedPoint = 1
      do k = ExtendedCart%is(3), ExtendedCart%ie(3)
        do j = ExtendedCart%is(2), ExtendedCart%ie(2)
          do i = ExtendedCart%is(1), ExtendedCart%ie(1)
            Point = [i,j,k]
            if (.not. ovkCartContains(PrincipalCart, Point)) then
              AdjustedPoint(:PrincipalCart%nd) = ovkCartPeriodicAdjust(PrincipalCart, Point)
              if (ovkCartContains(PrincipalCart, AdjustedPoint)) then
                Field%values(Point(1),Point(2),Point(3)) = Field%values(AdjustedPoint(1), &
                  AdjustedPoint(2),AdjustedPoint(3))
              end if
            end if
          end do
        end do
      end do
    end if

  end subroutine ovkFieldPeriodicFill_Logical

  subroutine ovkGetFieldPatch_Integer_Rank2(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_int), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    integer, dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2)), &
      intent(out) :: PatchData

    integer :: i, j
    logical :: AwayFromEdge
    integer, dimension(2) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          PatchData(i,j) = Field%values(i,j,1)
        end do
      end do
    else
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          Point = [i,j]
          Point = ovkCartPeriodicAdjust(Field%cart, Point)
          PatchData(i,j) = Field%values(Point(1),Point(2),1)
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_Integer_Rank2

  subroutine ovkGetFieldPatch_LargeInteger_Rank2(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_large_int), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    integer(lk), dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2)), &
      intent(out) :: PatchData

    integer :: i, j
    logical :: AwayFromEdge
    integer, dimension(2) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          PatchData(i,j) = Field%values(i,j,1)
        end do
      end do
    else
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          Point = [i,j]
          Point = ovkCartPeriodicAdjust(Field%cart, Point)
          PatchData(i,j) = Field%values(Point(1),Point(2),1)
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_LargeInteger_Rank2

  subroutine ovkGetFieldPatch_Real_Rank2(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_real), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    real(rk), dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2)), &
      intent(out) :: PatchData

    integer :: i, j
    logical :: AwayFromEdge
    integer, dimension(2) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          PatchData(i,j) = Field%values(i,j,1)
        end do
      end do
    else
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          Point = [i,j]
          Point = ovkCartPeriodicAdjust(Field%cart, Point)
          PatchData(i,j) = Field%values(Point(1),Point(2),1)
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_Real_Rank2

  subroutine ovkGetFieldPatch_Logical_Rank2(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_logical), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    logical, dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2)), &
      intent(out) :: PatchData

    integer :: i, j
    logical :: AwayFromEdge
    integer, dimension(2) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          PatchData(i,j) = Field%values(i,j,1)
        end do
      end do
    else
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          Point = [i,j]
          Point = ovkCartPeriodicAdjust(Field%cart, Point)
          PatchData(i,j) = Field%values(Point(1),Point(2),1)
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_Logical_Rank2

  subroutine ovkGetFieldPatch_Logical1Byte_Rank2(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_logical), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    logical(bk), dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2)), &
      intent(out) :: PatchData

    integer :: i, j
    logical :: AwayFromEdge
    integer, dimension(2) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          PatchData(i,j) = Field%values(i,j,1)
        end do
      end do
    else
      do j = PatchBegin(2), PatchEnd(2)
        do i = PatchBegin(1), PatchEnd(1)
          Point = [i,j]
          Point = ovkCartPeriodicAdjust(Field%cart, Point)
          PatchData(i,j) = Field%values(Point(1),Point(2),1)
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_Logical1Byte_Rank2

  subroutine ovkGetFieldPatch_Integer_Rank3(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_int), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    integer, dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2), &
      PatchBegin(3):PatchEnd(3)), intent(out) :: PatchData

    integer :: i, j, k
    logical :: AwayFromEdge
    integer, dimension(3) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            PatchData(i,j,k) = Field%values(i,j,k)
          end do
        end do
      end do
    else
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            Point = [i,j,k]
            Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
            PatchData(i,j,k) = Field%values(Point(1),Point(2),Point(3))
          end do
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_Integer_Rank3

  subroutine ovkGetFieldPatch_LargeInteger_Rank3(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_large_int), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    integer(lk), dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2), &
      PatchBegin(3):PatchEnd(3)), intent(out) :: PatchData

    integer :: i, j, k
    logical :: AwayFromEdge
    integer, dimension(3) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            PatchData(i,j,k) = Field%values(i,j,k)
          end do
        end do
      end do
    else
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            Point = [i,j,k]
            Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
            PatchData(i,j,k) = Field%values(Point(1),Point(2),Point(3))
          end do
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_LargeInteger_Rank3

  subroutine ovkGetFieldPatch_Real_Rank3(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_real), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    real(rk), dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2), &
      PatchBegin(3):PatchEnd(3)), intent(out) :: PatchData

    integer :: i, j, k
    logical :: AwayFromEdge
    integer, dimension(3) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            PatchData(i,j,k) = Field%values(i,j,k)
          end do
        end do
      end do
    else
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            Point = [i,j,k]
            Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
            PatchData(i,j,k) = Field%values(Point(1),Point(2),Point(3))
          end do
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_Real_Rank3

  subroutine ovkGetFieldPatch_Logical_Rank3(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_logical), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    logical, dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2), &
      PatchBegin(3):PatchEnd(3)), intent(out) :: PatchData

    integer :: i, j, k
    logical :: AwayFromEdge
    integer, dimension(3) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            PatchData(i,j,k) = Field%values(i,j,k)
          end do
        end do
      end do
    else
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            Point = [i,j,k]
            Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
            PatchData(i,j,k) = Field%values(Point(1),Point(2),Point(3))
          end do
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_Logical_Rank3

  subroutine ovkGetFieldPatch_Logical1Byte_Rank3(Field, PatchBegin, PatchEnd, PatchData)

    type(ovk_field_logical), intent(in) :: Field
    integer, dimension(:), intent(in) :: PatchBegin, PatchEnd
    logical(bk), dimension(PatchBegin(1):PatchEnd(1),PatchBegin(2):PatchEnd(2), &
      PatchBegin(3):PatchEnd(3)), intent(out) :: PatchData

    integer :: i, j, k
    logical :: AwayFromEdge
    integer, dimension(3) :: Point

    AwayFromEdge = ovkCartContains(Field%cart, PatchBegin) .and. &
      ovkCartContains(Field%cart, PatchEnd)

    if (AwayFromEdge) then
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            PatchData(i,j,k) = Field%values(i,j,k)
          end do
        end do
      end do
    else
      do k = PatchBegin(3), PatchEnd(3)
        do j = PatchBegin(2), PatchEnd(2)
          do i = PatchBegin(1), PatchEnd(1)
            Point = [i,j,k]
            Point(:Field%cart%nd) = ovkCartPeriodicAdjust(Field%cart, Point)
            PatchData(i,j,k) = Field%values(Point(1),Point(2),Point(3))
          end do
        end do
      end do
    end if

  end subroutine ovkGetFieldPatch_Logical1Byte_Rank3

![COVERAGEIGNORE]
  subroutine ovkPrintField_Integer(Field, OutputUnit, StartIndex, EndIndex, NumDigits)

    type(ovk_field_int), intent(in) :: Field
    integer, intent(in) :: OutputUnit
    integer, dimension(Field%cart%nd), intent(in), optional :: StartIndex, EndIndex
    integer, intent(in), optional :: NumDigits

    integer, dimension(MAX_DIMS) :: StartIndex_, EndIndex_
    integer :: NumDigits_
    integer :: i, j
    integer :: IDir, JDir, KDir
    integer :: IStart, IEnd, JStart, JEnd
    integer :: KSlice
    character(len=16) :: FormatString
    integer, dimension(MAX_DIMS) :: Point

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
      EndIndex_(MAX_DIMS) = Field%cart%is(MAX_DIMS)
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

    integer, dimension(MAX_DIMS) :: StartIndex_, EndIndex_
    integer :: NumDigits_
    integer :: i, j
    integer :: IDir, JDir, KDir
    integer :: IStart, IEnd, JStart, JEnd
    integer :: KSlice
    character(len=16) :: FormatString
    integer, dimension(MAX_DIMS) :: Point

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
      EndIndex_(MAX_DIMS) = Field%cart%is(MAX_DIMS)
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

    integer, dimension(MAX_DIMS) :: StartIndex_, EndIndex_
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
    integer, dimension(MAX_DIMS) :: Point

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
      EndIndex_(MAX_DIMS) = Field%cart%is(MAX_DIMS)
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

    integer, dimension(MAX_DIMS) :: StartIndex_, EndIndex_
    integer :: i, j
    integer :: IDir, JDir, KDir
    integer :: IStart, IEnd, JStart, JEnd
    integer :: KSlice
    integer, dimension(MAX_DIMS) :: Point

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
      EndIndex_(MAX_DIMS) = Field%cart%is(MAX_DIMS)
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

    integer, dimension(MAX_DIMS), intent(in) :: StartIndex, EndIndex
    integer, intent(out) :: IDir, JDir, KDir
    integer, intent(out) :: IStart, IEnd, JStart, JEnd
    integer, intent(out) :: KSlice

    integer :: d

    do d = 1, MAX_DIMS
      if (EndIndex(d) == StartIndex(d)) then
        exit
      end if
    end do
    KDir = d

    if (OVK_DEBUG) then
      if (KDir > MAX_DIMS) then
        write (ERROR_UNIT, '(a)') "ERROR: Range to print must be two-dimensional."
        stop 1
      end if
    end if

    IDir = modulo(KDir,MAX_DIMS)+1
    JDir = modulo(KDir+1,MAX_DIMS)+1

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
![/COVERAGEIGNORE]

end module ovkField
