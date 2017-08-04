! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkBoundingBox

  use ovkGlobal
  implicit none

  private

  ! API
  public :: ovk_bbox
  public :: ovk_bbox_
  public :: operator (==)
  public :: operator (/=)
  public :: ovkBBOverlaps
  public :: ovkBBContains
  public :: ovkBBContainsPoint
  public :: ovkBBIsEmpty
  public :: ovkBBSize
  public :: ovkBBMove
  public :: ovkBBGrow
  public :: ovkBBScale
  public :: ovkBBExtend
  public :: ovkBBUnion
  public :: ovkBBIntersect
  public :: ovkBBFromPoints

  type ovk_bbox
    type(t_noconstruct) :: noconstruct
    integer :: nd
    real(rk), dimension(MAX_ND) :: b, e
  end type ovk_bbox

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_bbox_
    module procedure ovk_bbox_Default
    module procedure ovk_bbox_Assigned_Empty
    module procedure ovk_bbox_Assigned_BeginEnd
  end interface ovk_bbox_

  interface operator (==)
    module procedure ovk_bbox_Equal
  end interface operator (==)

  interface operator (/=)
    module procedure ovk_bbox_NotEqual
  end interface operator (/=)

  interface ovkBBFromPoints
    module procedure ovkBBFromPoints_Rank2
    module procedure ovkBBFromPoints_Rank3
    module procedure ovkBBFromPoints_Rank4
  end interface ovkBBFromPoints

contains

  pure function ovk_bbox_Default() result(BBox)

    type(ovk_bbox) :: BBox

    BBox = ovk_bbox_Assigned_Empty(2)

  end function ovk_bbox_Default

  pure function ovk_bbox_Assigned_Empty(NumDims) result(BBox)

    integer, intent(in) :: NumDims
    type(ovk_bbox) :: BBox

    BBox%nd = NumDims
    BBox%b = 0._rk
    BBox%e(:NumDims) = -1._rk
    BBox%e(NumDims+1:) = 0._rk

  end function ovk_bbox_Assigned_Empty

  pure function ovk_bbox_Assigned_BeginEnd(NumDims, B, E) result(BBox)

    integer, intent(in) :: NumDims
    real(rk), dimension(NumDims), intent(in) :: B, E
    type(ovk_bbox) :: BBox

    BBox%nd = NumDims
    BBox%b(:BBox%nd) = B
    BBox%b(BBox%nd+1:) = 0._rk
    BBox%e(:BBox%nd) = E
    BBox%e(BBox%nd+1:) = 0._rk

  end function ovk_bbox_Assigned_BeginEnd

  pure function ovk_bbox_Equal(LeftBBox, RightBBox) result(Equal)

    type(ovk_bbox), intent(in) :: LeftBBox, RightBBox
    logical :: Equal

    Equal = &
      LeftBBox%nd == RightBBox%nd .and. &
      all(LeftBBox%b == RightBBox%b) .and. &
      all(LeftBBox%e == RightBBox%e)

  end function ovk_bbox_Equal

  pure function ovk_bbox_NotEqual(LeftBBox, RightBBox) result(NotEqual)

    type(ovk_bbox), intent(in) :: LeftBBox, RightBBox
    logical :: NotEqual

    NotEqual = .not. ovk_bbox_Equal(LeftBBox, RightBBox)

  end function ovk_bbox_NotEqual

  pure function ovkBBOverlaps(LeftBBox, RightBBox) result(BBOverlaps)

    type(ovk_bbox), intent(in) :: LeftBBox, RightBBox
    logical :: BBOverlaps

    BBOverlaps = &
      all(LeftBBox%e(:LeftBBox%nd) >= LeftBBox%b(:LeftBBox%nd)) .and. &
      all(RightBBox%e(:RightBBox%nd) >= RightBBox%b(:RightBBox%nd)) .and. &
      all(RightBBox%e(:RightBBox%nd) >= LeftBBox%b(:LeftBBox%nd)) .and. &
      all(LeftBBox%e(:LeftBBox%nd) >= RightBBox%b(:RightBBox%nd))

  end function ovkBBOverlaps

  pure function ovkBBContains(LeftBBox, RightBBox) result(BBContains)

    type(ovk_bbox), intent(in) :: LeftBBox, RightBBox
    logical :: BBContains

    BBContains = &
      all(LeftBBox%e(:LeftBBox%nd) >= LeftBBox%b(:LeftBBox%nd)) .and. &
      all(RightBBox%e(:RightBBox%nd) >= RightBBox%b(:RightBBox%nd)) .and. &
      all(RightBBox%b(:RightBBox%nd) >= LeftBBox%b(:LeftBBox%nd)) .and. &
      all(LeftBBox%e(:LeftBBox%nd) >= RightBBox%e(:RightBBox%nd))

  end function ovkBBContains

  pure function ovkBBContainsPoint(BBox, Point) result(BBContainsPoint)

    type(ovk_bbox), intent(in) :: BBox
    real(rk), dimension(BBox%nd), intent(in) :: Point
    logical :: BBContainsPoint

    BBContainsPoint = all(Point >= BBox%b(:BBox%nd) .and. Point <= BBox%e(:BBox%nd))

  end function ovkBBContainsPoint

  pure function ovkBBIsEmpty(BBox) result(BBIsEmpty)

    type(ovk_bbox), intent(in) :: BBox
    logical :: BBIsEmpty

    BBIsEmpty = any(BBox%e(:BBox%nd) < BBox%b(:BBox%nd))

  end function ovkBBIsEmpty

  pure function ovkBBSize(BBox) result(BBSize)

    type(ovk_bbox), intent(in) :: BBox
    real(rk), dimension(BBox%nd) :: BBSize

    BBSize = max(BBox%e(:BBox%nd) - BBox%b(:BBox%nd), 0._rk)

  end function ovkBBSize

  pure function ovkBBMove(BBox, Amount) result(BBMove)

    type(ovk_bbox), intent(in) :: BBox
    real(rk), dimension(BBox%nd), intent(in) :: Amount
    type(ovk_bbox) :: BBMove

    BBMove%nd = BBox%nd
    BBMove%b(:BBox%nd) = BBox%b(:BBox%nd) + Amount
    BBMove%b(BBox%nd+1:) = BBox%b(BBox%nd+1:)
    BBMove%e(:BBox%nd) = BBox%e(:BBox%nd) + Amount
    BBMove%e(BBox%nd+1:) = BBox%e(BBox%nd+1:)

  end function ovkBBMove

  pure function ovkBBGrow(BBox, Amount) result(BBGrow)

    type(ovk_bbox), intent(in) :: BBox
    real(rk), intent(in) :: Amount
    type(ovk_bbox) :: BBGrow

    BBGrow%nd = BBox%nd
    BBGrow%b(:BBox%nd) = BBox%b(:BBox%nd) - Amount
    BBGrow%b(BBox%nd+1:) = BBox%b(BBox%nd+1:)
    BBGrow%e(:BBox%nd) = BBox%e(:BBox%nd) + Amount
    BBGrow%e(BBox%nd+1:) = BBox%e(BBox%nd+1:)

  end function ovkBBGrow

  pure function ovkBBScale(BBox, Factor) result(BBScale)

    type(ovk_bbox), intent(in) :: BBox
    real(rk), intent(in) :: Factor
    type(ovk_bbox) :: BBScale

    real(rk), dimension(BBox%nd) :: BBCenter, BBHalfSize

    BBCenter = 0.5_rk * (BBox%b(:BBox%nd) + BBox%e(:BBox%nd))
    BBHalfSize = 0.5_rk * Factor * (BBox%e(:BBox%nd) - BBox%b(:BBox%nd))

    BBScale%nd = BBox%nd
    BBScale%b(:BBox%nd) = BBCenter - BBHalfSize
    BBScale%b(BBox%nd+1:) = BBox%b(BBox%nd+1:)
    BBScale%e(:BBox%nd) = BBCenter + BBHalfSize
    BBScale%e(BBox%nd+1:) = BBox%e(BBox%nd+1:)

  end function ovkBBScale

  pure function ovkBBExtend(BBox, Point) result(BBExtend)

    type(ovk_bbox), intent(in) :: BBox
    real(rk), dimension(BBox%nd), intent(in) :: Point
    type(ovk_bbox) :: BBExtend

    BBExtend%nd = BBox%nd
    BBExtend%b(:BBox%nd) = min(BBox%b(:BBox%nd), Point)
    BBExtend%b(BBox%nd+1:) = BBox%b(BBox%nd+1:)
    BBExtend%e(:BBox%nd) = max(BBox%e(:BBox%nd), Point)
    BBExtend%e(BBox%nd+1:) = BBox%e(BBox%nd+1:)

  end function ovkBBExtend

  pure function ovkBBUnion(LeftBBox, RightBBox) result(BBUnion)

    type(ovk_bbox), intent(in) :: LeftBBox
    type(ovk_bbox), intent(in) :: RightBBox
    type(ovk_bbox) :: BBUnion

    if (ovkBBIsEmpty(LeftBBox)) then
      BBUnion = RightBBox
    else if (ovkBBIsEmpty(RightBBox)) then
      BBUnion = LeftBBox
    else
      BBUnion%nd = LeftBBox%nd
      BBUnion%b(:LeftBBox%nd) = min(LeftBBox%b(:LeftBBox%nd), RightBBox%b(:RightBBox%nd))
      BBUnion%b(LeftBBox%nd+1:) = LeftBBox%b(LeftBBox%nd+1:)
      BBUnion%e(:LeftBBox%nd) = max(LeftBBox%e(:LeftBBox%nd), RightBBox%e(:RightBBox%nd))
      BBUnion%e(LeftBBox%nd+1:) = LeftBBox%e(LeftBBox%nd+1:)
    end if

  end function ovkBBUnion

  pure function ovkBBIntersect(LeftBBox, RightBBox) result(BBIntersect)

    type(ovk_bbox), intent(in) :: LeftBBox
    type(ovk_bbox), intent(in) :: RightBBox
    type(ovk_bbox) :: BBIntersect

    BBIntersect%nd = LeftBBox%nd
    BBIntersect%b(:LeftBBox%nd) = max(LeftBBox%b(:LeftBBox%nd), RightBBox%b(:RightBBox%nd))
    BBIntersect%b(LeftBBox%nd+1:) = LeftBBox%b(LeftBBox%nd+1:)
    BBIntersect%e(:LeftBBox%nd) = min(LeftBBox%e(:LeftBBox%nd), RightBBox%e(:RightBBox%nd))
    BBIntersect%e(LeftBBox%nd+1:) = LeftBBox%e(LeftBBox%nd+1:)

  end function ovkBBIntersect

  pure function ovkBBFromPoints_Rank2(Points) result(BBFromPoints)

    real(rk), dimension(:,:), intent(in) :: Points
    type(ovk_bbox) :: BBFromPoints

    integer :: i
    integer :: NumDims

    NumDims = size(Points,1)
    BBFromPoints%nd = NumDims
    if (size(Points,2) > 0) then
      BBFromPoints%b(:NumDims) = Points(:,1)
      BBFromPoints%b(NumDims+1:) = 0._rk
      BBFromPoints%e(:NumDims) = Points(:,1)
      BBFromPoints%e(NumDims+1:) = 0._rk
      do i = 1, size(Points,2)
        BBFromPoints%b(:NumDims) = min(BBFromPoints%b(:NumDims), Points(:,i))
        BBFromPoints%e(:NumDims) = max(BBFromPoints%e(:NumDims), Points(:,i))
      end do
    else
      BBFromPoints%b(:NumDims) = 0._rk
      BBFromPoints%b(NumDims+1:) = 0._rk
      BBFromPoints%e(:NumDims) = -1._rk
      BBFromPoints%e(NumDims+1:) = 0._rk
    end if

  end function ovkBBFromPoints_Rank2

  pure function ovkBBFromPoints_Rank3(Points) result(BBFromPoints)

    real(rk), dimension(:,:,:), intent(in) :: Points
    type(ovk_bbox) :: BBFromPoints

    integer :: i, j
    integer :: NumDims

    NumDims = size(Points,1)
    BBFromPoints%nd = NumDims
    if (size(Points,2) > 0 .and. size(Points,3) > 0) then
      BBFromPoints%b(:NumDims) = Points(:,1,1)
      BBFromPoints%b(NumDims+1:) = 0._rk
      BBFromPoints%e(:NumDims) = Points(:,1,1)
      BBFromPoints%e(NumDims+1:) = 0._rk
      do j = 1, size(Points,3)
        do i = 1, size(Points,2)
          BBFromPoints%b(:NumDims) = min(BBFromPoints%b(:NumDims), Points(:,i,j))
          BBFromPoints%e(:NumDims) = max(BBFromPoints%e(:NumDims), Points(:,i,j))
        end do
      end do
    else
      BBFromPoints%b(:NumDims) = 0._rk
      BBFromPoints%b(NumDims+1:) = 0._rk
      BBFromPoints%e(:NumDims) = -1._rk
      BBFromPoints%e(NumDims+1:) = 0._rk
    end if

  end function ovkBBFromPoints_Rank3

  pure function ovkBBFromPoints_Rank4(Points) result(BBFromPoints)

    real(rk), dimension(:,:,:,:), intent(in) :: Points
    type(ovk_bbox) :: BBFromPoints

    integer :: i, j, k
    integer :: NumDims

    NumDims = size(Points,1)
    BBFromPoints%nd = NumDims
    if (size(Points,2) > 0 .and. size(Points,3) > 0 .and. size(Points,4) > 0) then
      BBFromPoints%b(:NumDims) = Points(:,1,1,1)
      BBFromPoints%b(NumDims+1:) = 0._rk
      BBFromPoints%e(:NumDims) = Points(:,1,1,1)
      BBFromPoints%e(NumDims+1:) = 0._rk
      do k = 1, size(Points,4)
        do j = 1, size(Points,3)
          do i = 1, size(Points,2)
            BBFromPoints%b(:NumDims) = min(BBFromPoints%b(:NumDims), Points(:,i,j,k))
            BBFromPoints%e(:NumDims) = max(BBFromPoints%e(:NumDims), Points(:,i,j,k))
          end do
        end do
      end do
    else
      BBFromPoints%b(:NumDims) = 0._rk
      BBFromPoints%b(NumDims+1:) = 0._rk
      BBFromPoints%e(:NumDims) = -1._rk
      BBFromPoints%e(NumDims+1:) = 0._rk
    end if

  end function ovkBBFromPoints_Rank4

end module ovkBoundingBox
