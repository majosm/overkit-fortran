! Copyright (c) 2018 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkCart

  use ovkGlobal
  implicit none

  private

  ! Public API
  public :: ovk_cart
  public :: ovk_cart_
  public :: operator (==)
  public :: operator (/=)
  public :: ovkCartIsEmpty
  public :: ovkCartSize
  public :: ovkCartCount
  public :: ovkCartTupleToIndex
  public :: ovkCartIndexToTuple
  public :: ovkCartPeriodicAdjust
  public :: ovkCartContains
  public :: ovkCartClamp
  public :: ovkCartIsCompatible
  public :: ovkCartConvertPeriodicStorage
  public :: ovkCartPointToCell

  type ovk_cart
    type(t_noconstruct) :: noconstruct
    integer :: nd
    integer, dimension(MAX_DIMS) :: is, ie
    logical, dimension(MAX_DIMS) :: periodic
    integer :: periodic_storage
  end type ovk_cart

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_cart_
    module procedure ovk_cart_Default
    module procedure ovk_cart_Assigned_Empty
    module procedure ovk_cart_Assigned_EmptyPeriodic
    module procedure ovk_cart_Assigned_EmptyPeriodicWithStorage
    module procedure ovk_cart_Assigned_NumPoints
    module procedure ovk_cart_Assigned_NumPointsPeriodic
    module procedure ovk_cart_Assigned_NumPointsPeriodicWithStorage
    module procedure ovk_cart_Assigned_StartEnd
    module procedure ovk_cart_Assigned_StartEndPeriodic
    module procedure ovk_cart_Assigned_StartEndPeriodicWithStorage
  end interface ovk_cart_

  interface operator (==)
    module procedure ovk_cart_Equal
  end interface operator (==)

  interface operator (/=)
    module procedure ovk_cart_NotEqual
  end interface operator (/=)

contains

  pure function ovk_cart_Default() result(Cart)

    type(ovk_cart) :: Cart

    Cart = ovk_cart_Assigned_Empty(2)

  end function ovk_cart_Default

  pure function ovk_cart_Assigned_Empty(NumDims) result(Cart)

    integer, intent(in) :: NumDims
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is = 1
    Cart%ie(:NumDims) = 0
    Cart%ie(NumDims+1:) = 1
    Cart%periodic = .false.
    Cart%periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE

  end function ovk_cart_Assigned_Empty

  pure function ovk_cart_Assigned_EmptyPeriodic(NumDims, Periodic) result(Cart)

    integer, intent(in) :: NumDims
    logical, dimension(NumDims), intent(in) :: Periodic
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is = 1
    Cart%ie(:NumDims) = 0
    Cart%ie(NumDims+1:) = 1
    Cart%periodic(:NumDims) = Periodic
    Cart%periodic(NumDims+1:) = .false.
    Cart%periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE

  end function ovk_cart_Assigned_EmptyPeriodic

  pure function ovk_cart_Assigned_EmptyPeriodicWithStorage(NumDims, Periodic, PeriodicStorage) &
    result(Cart)

    integer, intent(in) :: NumDims
    logical, dimension(NumDims), intent(in) :: Periodic
    integer, intent(in) :: PeriodicStorage
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is = 1
    Cart%ie(:NumDims) = 0
    Cart%ie(NumDims+1:) = 1
    Cart%periodic(:NumDims) = Periodic
    Cart%periodic(NumDims+1:) = .false.
    Cart%periodic_storage = PeriodicStorage

  end function ovk_cart_Assigned_EmptyPeriodicWithStorage

  pure function ovk_cart_Assigned_NumPoints(NumDims, NumPoints) result(Cart)

    integer, intent(in) :: NumDims
    integer, dimension(NumDims), intent(in) :: NumPoints
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is = 1
    Cart%ie(:NumDims) = NumPoints
    Cart%ie(NumDims+1:) = 1
    Cart%periodic = .false.
    Cart%periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE

  end function ovk_cart_Assigned_NumPoints

  pure function ovk_cart_Assigned_NumPointsPeriodic(NumDims, NumPoints, Periodic) result(Cart)

    integer, intent(in) :: NumDims
    integer, dimension(NumDims), intent(in) :: NumPoints
    logical, dimension(NumDims), intent(in) :: Periodic
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is = 1
    Cart%ie(:NumDims) = NumPoints
    Cart%ie(NumDims+1:) = 1
    Cart%periodic(:NumDims) = Periodic
    Cart%periodic(NumDims+1:) = .false.
    Cart%periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE

  end function ovk_cart_Assigned_NumPointsPeriodic

  pure function ovk_cart_Assigned_NumPointsPeriodicWithStorage(NumDims, NumPoints, Periodic, &
    PeriodicStorage) result(Cart)

    integer, intent(in) :: NumDims
    integer, dimension(NumDims), intent(in) :: NumPoints
    logical, dimension(NumDims), intent(in) :: Periodic
    integer, intent(in) :: PeriodicStorage
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is = 1
    Cart%ie(:NumDims) = NumPoints
    Cart%ie(NumDims+1:) = 1
    Cart%periodic(:NumDims) = Periodic
    Cart%periodic(NumDims+1:) = .false.
    Cart%periodic_storage = PeriodicStorage

  end function ovk_cart_Assigned_NumPointsPeriodicWithStorage

  pure function ovk_cart_Assigned_StartEnd(NumDims, StartIndex, EndIndex) result(Cart)

    integer, intent(in) :: NumDims
    integer, dimension(NumDims), intent(in) :: StartIndex, EndIndex
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is(:NumDims) = StartIndex
    Cart%is(NumDims+1:) = 1
    Cart%ie(:NumDims) = EndIndex
    Cart%ie(NumDims+1:) = 1
    Cart%periodic = .false.
    Cart%periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE

  end function ovk_cart_Assigned_StartEnd

  pure function ovk_cart_Assigned_StartEndPeriodic(NumDims, StartIndex, EndIndex, Periodic) &
    result(Cart)

    integer, intent(in) :: NumDims
    integer, dimension(NumDims), intent(in) :: StartIndex, EndIndex
    logical, dimension(NumDims), intent(in) :: Periodic
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is(:NumDims) = StartIndex
    Cart%is(NumDims+1:) = 1
    Cart%ie(:NumDims) = EndIndex
    Cart%ie(NumDims+1:) = 1
    Cart%periodic(:NumDims) = Periodic
    Cart%periodic(NumDims+1:) = .false.
    Cart%periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE

  end function ovk_cart_Assigned_StartEndPeriodic

  pure function ovk_cart_Assigned_StartEndPeriodicWithStorage(NumDims, StartIndex, EndIndex, &
    Periodic, PeriodicStorage) result(Cart)

    integer, intent(in) :: NumDims
    integer, dimension(NumDims), intent(in) :: StartIndex, EndIndex
    logical, dimension(NumDims), intent(in) :: Periodic
    integer, intent(in) :: PeriodicStorage
    type(ovk_cart) :: Cart

    Cart%nd = NumDims
    Cart%is(:NumDims) = StartIndex
    Cart%is(NumDims+1:) = 1
    Cart%ie(:NumDims) = EndIndex
    Cart%ie(NumDims+1:) = 1
    Cart%periodic(:NumDims) = Periodic
    Cart%periodic(NumDims+1:) = .false.
    Cart%periodic_storage = PeriodicStorage

  end function ovk_cart_Assigned_StartEndPeriodicWithStorage

  pure function ovk_cart_Equal(LeftCart, RightCart) result(Equal)

    type(ovk_cart), intent(in) :: LeftCart, RightCart
    logical :: Equal

    Equal = &
      LeftCart%nd == RightCart%nd .and. &
      all(LeftCart%is == RightCart%is) .and. &
      all(LeftCart%ie == RightCart%ie) .and. &
      all(LeftCart%periodic .eqv. RightCart%periodic) .and. &
      LeftCart%periodic_storage == RightCart%periodic_storage

  end function ovk_cart_Equal

  pure function ovk_cart_NotEqual(LeftCart, RightCart) result(NotEqual)

    type(ovk_cart), intent(in) :: LeftCart, RightCart
    logical :: NotEqual

    NotEqual = .not. ovk_cart_Equal(LeftCart, RightCart)

  end function ovk_cart_NotEqual

  pure function ovkCartIsEmpty(Cart) result(CartIsEmpty)

    type(ovk_cart), intent(in) :: Cart
    logical :: CartIsEmpty

    CartIsEmpty = any(Cart%ie(:Cart%nd) < Cart%is(:Cart%nd))

  end function ovkCartIsEmpty

  pure function ovkCartSize(Cart) result(NumPoints)

    type(ovk_cart), intent(in) :: Cart
    integer, dimension(Cart%nd) :: NumPoints

    NumPoints = Cart%ie(:Cart%nd) - Cart%is(:Cart%nd) + 1

  end function ovkCartSize

  pure function ovkCartCount(Cart) result(NumPointsTotal)

    type(ovk_cart), intent(in) :: Cart
    integer(lk) :: NumPointsTotal

    integer(lk), dimension(Cart%nd) :: NumPoints

    NumPoints = int(Cart%ie(:Cart%nd)-Cart%is(:Cart%nd)+1, kind=lk)

    NumPointsTotal = product(NumPoints)

  end function ovkCartCount

  pure function ovkCartTupleToIndex(Cart, Tuple) result(Ind)

    type(ovk_cart), intent(in) :: Cart
    integer, dimension(Cart%nd), intent(in) :: Tuple
    integer(lk) :: Ind

    integer :: d
    integer(lk), dimension(Cart%nd) :: N
    integer(lk) :: Stride

    N = Cart%ie(:Cart%nd) - Cart%is(:Cart%nd) + 1

    Stride = 1_lk
    Ind = 1_lk + Stride*(Tuple(1) - Cart%is(1))
    do d = 2, Cart%nd
      Stride = Stride*N(d-1)
      Ind = Ind + Stride*(Tuple(d) - Cart%is(d))
    end do

  end function ovkCartTupleToIndex

  pure function ovkCartIndexToTuple(Cart, Ind) result(Tuple)

    type(ovk_cart), intent(in) :: Cart
    integer(lk), intent(in) :: Ind
    integer, dimension(Cart%nd) :: Tuple

    integer :: d
    integer(lk), dimension(Cart%nd) :: N
    integer(lk), dimension(Cart%nd) :: Strides
    integer(lk) :: ReducedInd
    integer :: Offset

    N = Cart%ie(:Cart%nd) - Cart%is(:Cart%nd) + 1

    Strides(1) = 1
    do d = 2, Cart%nd
      Strides(d) = Strides(d-1)*N(d-1)
    end do

    ReducedInd = Ind-1
    do d = Cart%nd, 1, -1
      Offset = int(ReducedInd/Strides(d))
      Tuple(d) = Cart%is(d) + Offset
      ReducedInd = ReducedInd - Offset*Strides(d)
    end do

  end function ovkCartIndexToTuple

  pure function ovkCartPeriodicAdjust(Cart, Point) result(AdjustedPoint)

    type(ovk_cart), intent(in) :: Cart
    integer, dimension(Cart%nd), intent(in) :: Point
    integer, dimension(Cart%nd) :: AdjustedPoint

    integer, dimension(Cart%nd) :: NumPeriod

    NumPeriod = Cart%ie(:Cart%nd) - Cart%is(:Cart%nd) + 1
    NumPeriod = merge(NumPeriod-1, NumPeriod, Cart%periodic(:Cart%nd) .and. &
        Cart%periodic_storage == OVK_PERIODIC_STORAGE_DUPLICATED)

    AdjustedPoint = merge(modulo(Point-1, NumPeriod)+1, Point, Cart%periodic(:Cart%nd))

  end function ovkCartPeriodicAdjust

  pure function ovkCartContains(Cart, Point) result(CartContains)

    type(ovk_cart), intent(in) :: Cart
    integer, dimension(Cart%nd), intent(in) :: Point
    logical :: CartContains

    CartContains = all(Point >= Cart%is(:Cart%nd)) .and. all(Point <= Cart%ie(:Cart%nd))

  end function ovkCartContains

  pure function ovkCartClamp(Cart, Point) result(ClampedPoint)

    type(ovk_cart), intent(in) :: Cart
    integer, dimension(Cart%nd), intent(in) :: Point
    integer, dimension(Cart%nd) :: ClampedPoint

    ClampedPoint = min(max(Point, Cart%is(:Cart%nd)), Cart%ie(:Cart%nd))

  end function ovkCartClamp

  pure function ovkCartIsCompatible(Cart, OtherCart) result(CartIsCompatible)

    type(ovk_cart), intent(in) :: Cart, OtherCart
    logical :: CartIsCompatible

    type(ovk_cart) :: ConvertedCart

    ConvertedCart = ovkCartConvertPeriodicStorage(Cart, OtherCart%periodic_storage)

    CartIsCompatible = ConvertedCart == OtherCart

  end function ovkCartIsCompatible

  pure function ovkCartConvertPeriodicStorage(Cart, PeriodicStorage) result(ConvertedCart)

    type(ovk_cart), intent(in) :: Cart
    integer, intent(in) :: PeriodicStorage
    type(ovk_cart) :: ConvertedCart

    ConvertedCart = Cart
    ConvertedCart%periodic_storage = PeriodicStorage

    if (ConvertedCart%periodic_storage /= Cart%periodic_storage) then
      if (Cart%periodic_storage == OVK_PERIODIC_STORAGE_DUPLICATED) then
        ConvertedCart%ie = merge(ConvertedCart%ie-1, ConvertedCart%ie, ConvertedCart%periodic)
      else
        ConvertedCart%ie = merge(ConvertedCart%ie+1, ConvertedCart%ie, ConvertedCart%periodic)
      end if
    end if

  end function ovkCartConvertPeriodicStorage

  pure function ovkCartPointToCell(Cart) result(CellCart)

    type(ovk_cart), intent(in) :: Cart
    type(ovk_cart) :: CellCart

    CellCart = Cart

    if (Cart%periodic_storage == OVK_PERIODIC_STORAGE_UNIQUE) then
      CellCart%ie(:Cart%nd) = merge(Cart%ie(:Cart%nd), Cart%ie(:Cart%nd)-1, Cart%periodic(:Cart%nd))
    else
      CellCart%ie(:Cart%nd) = Cart%ie(:Cart%nd)-1
    end if

    CellCart%periodic_storage = OVK_PERIODIC_STORAGE_UNIQUE

  end function ovkCartPointToCell

end module ovkCart
