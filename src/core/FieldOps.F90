! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkFieldOps

  use ovkCart
  use ovkField
  use ovkGlobal
  implicit none

  private

  ! API
  public :: ovkDetectEdge
  public :: ovkDilate
  public :: ovkErode
  public :: ovkConnectedComponents
  public :: ovkFloodFill
  public :: ovkThreshold
  public :: ovkDistanceField
  public :: ovkCountMask
  public :: ovkFilterState
  public :: OVK_INNER_EDGE, OVK_OUTER_EDGE

  integer, parameter :: OVK_INNER_EDGE = 1
  integer, parameter :: OVK_OUTER_EDGE = 2

  interface ovkThreshold
    module procedure ovkThreshold_Integer
    module procedure ovkThreshold_LargeInteger
    module procedure ovkThreshold_Real
  end interface ovkThreshold

contains

  subroutine ovkDetectEdge(Mask, EdgeType, BoundaryValue, IncludeBoundary, EdgeMask)

    type(ovk_field_logical), intent(in) :: Mask
    integer, intent(in) :: EdgeType
    integer, intent(in) :: BoundaryValue
    logical, intent(in) :: IncludeBoundary
    type(ovk_field_logical), intent(out) :: EdgeMask

    integer :: i, j, k, m, n, o
    type(ovk_cart) :: Cart
    type(ovk_cart) :: EdgeCart
    logical :: EdgeValue
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: MirrorPoint
    logical :: Value
    integer, dimension(MAX_ND) :: NeighborLower, NeighborUpper
    logical :: AwayFromCartEdge
    integer, dimension(MAX_ND) :: Neighbor
    logical :: NeighborValue

    Cart = Mask%cart
    EdgeCart = Mask%cart

    if (IncludeBoundary) then
      EdgeCart%is(:EdgeCart%nd) = EdgeCart%is(:EdgeCart%nd) - merge(0, 1, &
        EdgeCart%periodic(:EdgeCart%nd))
      EdgeCart%ie(:EdgeCart%nd) = EdgeCart%ie(:EdgeCart%nd) + merge(0, 1, &
        EdgeCart%periodic(:EdgeCart%nd))
    end if

    EdgeMask = ovk_field_logical_(EdgeCart, .false.)

    ! Points on inner edge will have Mask == .true., points on outer edge will have Mask == .false.
    EdgeValue = EdgeType == OVK_INNER_EDGE

    do k = EdgeCart%is(3), EdgeCart%ie(3)
      do j = EdgeCart%is(2), EdgeCart%ie(2)
        do i = EdgeCart%is(1), EdgeCart%ie(1)
          Point = [i,j,k]
          if (ovkCartContains(Cart, Point)) then
            Value = Mask%values(i,j,k)
          else
            select case (BoundaryValue)
            case (OVK_TRUE)
              Value = .true.
            case (OVK_FALSE)
              Value = .false.
            case (OVK_MIRROR)
              MirrorPoint = Point + max(Cart%is-Point,0) + min(Cart%ie-Point,0)
              Value = Mask%values(MirrorPoint(1),MirrorPoint(2),MirrorPoint(3))
            end select
          end if
          if (Value .eqv. EdgeValue) then
            NeighborLower(:Cart%nd) = Point(:Cart%nd)-1
            NeighborLower(Cart%nd+1:) = Point(Cart%nd+1:)
            NeighborUpper(:Cart%nd) = Point(:Cart%nd)+1
            NeighborUpper(Cart%nd+1:) = Point(Cart%nd+1:)
            AwayFromCartEdge = ovkCartContains(Cart, NeighborLower) .and. &
              ovkCartContains(Cart, NeighborUpper)
            if (AwayFromCartEdge) then
              L1: &
              do o = NeighborLower(3), NeighborUpper(3)
                do n = NeighborLower(2), NeighborUpper(2)
                  do m = NeighborLower(1), NeighborUpper(1)
                    Neighbor = [m,n,o]
                    NeighborValue = Mask%values(Neighbor(1),Neighbor(2),Neighbor(3))
                    if (NeighborValue .neqv. Value) then
                      EdgeMask%values(i,j,k) = .true.
                      exit L1
                    end if
                  end do
                end do
              end do L1
            else
              L2: &
              do o = NeighborLower(3), NeighborUpper(3)
                do n = NeighborLower(2), NeighborUpper(2)
                  do m = NeighborLower(1), NeighborUpper(1)
                    Neighbor = [m,n,o]
                    Neighbor(:Cart%nd) = ovkCartPeriodicAdjust(Cart, Neighbor)
                    if (ovkCartContains(Cart, Neighbor)) then
                      NeighborValue = Mask%values(Neighbor(1),Neighbor(2),Neighbor(3))
                    else
                      select case (BoundaryValue)
                      case (OVK_TRUE)
                        NeighborValue = .true.
                      case (OVK_FALSE)
                        NeighborValue = .false.
                      case (OVK_MIRROR)
                        MirrorPoint = Neighbor + max(Cart%is-Neighbor,0) + min(Cart%ie-Neighbor,0)
                        NeighborValue = Mask%values(MirrorPoint(1),MirrorPoint(2),MirrorPoint(3))
                      end select
                    end if
                    if (NeighborValue .neqv. Value) then
                      EdgeMask%values(i,j,k) = .true.
                      exit L2
                    end if
                  end do
                end do
              end do L2
            end if
          end if
        end do
      end do
    end do

  end subroutine ovkDetectEdge

  subroutine ovkDilate(Mask, Amount, BoundaryValue)

    type(ovk_field_logical), intent(inout) :: Mask
    integer, intent(in) :: Amount
    integer, intent(in) :: BoundaryValue

    call DilateErode(Mask, Amount, BoundaryValue)

  end subroutine ovkDilate

  subroutine ovkErode(Mask, Amount, BoundaryValue)

    type(ovk_field_logical), intent(inout) :: Mask
    integer, intent(in) :: Amount
    integer, intent(in) :: BoundaryValue

    call DilateErode(Mask, -Amount, BoundaryValue)

  end subroutine ovkErode

  subroutine DilateErode(Mask, Amount, BoundaryValue)

    type(ovk_field_logical), intent(inout) :: Mask
    integer, intent(in) :: Amount
    integer, intent(in) :: BoundaryValue

    integer :: i, j, k, m, n, o
    type(ovk_cart) :: Cart
    type(ovk_cart) :: EdgeCart
    integer :: FillDistance
    logical :: FillValue
    integer :: EdgeType
    integer :: EdgeBoundaryValue
    type(ovk_field_logical) :: EdgeMask
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: FillLower, FillUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: FillPoint

    if (Amount == 0) return

    if (Amount > 0) then
      FillDistance = Amount
      FillValue = .true.
      EdgeType = OVK_INNER_EDGE
      EdgeBoundaryValue = merge(OVK_TRUE, OVK_MIRROR, BoundaryValue == OVK_TRUE)
    else
      FillDistance = -Amount
      FillValue = .false.
      EdgeType = OVK_OUTER_EDGE
      EdgeBoundaryValue = merge(OVK_FALSE, OVK_MIRROR, BoundaryValue == OVK_FALSE)
    end if

    call ovkDetectEdge(Mask, EdgeType, EdgeBoundaryValue, .true., EdgeMask)

    Cart = Mask%cart
    EdgeCart = EdgeMask%cart

    do k = EdgeCart%is(3), EdgeCart%ie(3)
      do j = EdgeCart%is(2), EdgeCart%ie(2)
        do i = EdgeCart%is(1), EdgeCart%ie(1)
          if (EdgeMask%values(i,j,k)) then
            Point = [i,j,k]
            FillLower(:Cart%nd) = Point(:Cart%nd)-FillDistance
            FillLower(Cart%nd+1:) = Point(Cart%nd+1:)
            FillUpper(:Cart%nd) = Point(:Cart%nd)+FillDistance
            FillUpper(Cart%nd+1:) = Point(Cart%nd+1:)
            AwayFromEdge = ovkCartContains(Cart, FillLower) .and. ovkCartContains(Cart, FillUpper)
            if (AwayFromEdge) then
              do o = FillLower(3), FillUpper(3)
                do n = FillLower(2), FillUpper(2)
                  do m = FillLower(1), FillUpper(1)
                    FillPoint = [m,n,o]
                    Mask%values(FillPoint(1),FillPoint(2),FillPoint(3)) = FillValue
                  end do
                end do
              end do
            else
              do o = FillLower(3), FillUpper(3)
                do n = FillLower(2), FillUpper(2)
                  do m = FillLower(1), FillUpper(1)
                    FillPoint = [m,n,o]
                    FillPoint(:Cart%nd) = ovkCartPeriodicAdjust(Cart, FillPoint)
                    if (ovkCartContains(Cart, FillPoint)) then
                      Mask%values(FillPoint(1),FillPoint(2),FillPoint(3)) = FillValue
                    end if
                  end do
                end do
              end do
            end if
          end if
        end do
      end do
    end do

    call ovkFieldPeriodicFill(Mask)

  end subroutine DilateErode

  subroutine ovkConnectedComponents(Mask, NumComponents, ComponentLabels)

    type(ovk_field_logical), intent(in) :: Mask
    integer, intent(out) :: NumComponents
    type(ovk_field_int), intent(out) :: ComponentLabels

    integer :: i, j, k, m, n, o
    type(ovk_cart) :: Cart
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: NeighborLower, NeighborUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Neighbor
    integer :: Label, NeighborLabel
    integer, dimension(:), allocatable :: ReducedComponentLabel

    Cart = Mask%cart

    NumComponents = 0
    ComponentLabels = ovk_field_int_(Cart, 0)

    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          if (Mask%values(i,j,k)) then
            Point = [i,j,k]
            NeighborLower(:Cart%nd) = Point(:Cart%nd)-1
            NeighborLower(Cart%nd+1:) = Point(Cart%nd+1:)
            NeighborUpper(:Cart%nd) = Point(:Cart%nd)+1
            NeighborUpper(Cart%nd+1:) = Point(Cart%nd+1:)
            AwayFromEdge = ovkCartContains(Cart, NeighborLower) .and. &
              ovkCartContains(Cart, NeighborUpper)
            if (AwayFromEdge) then
              L1: &
              do o = NeighborLower(3), NeighborUpper(3)
                do n = NeighborLower(2), NeighborUpper(2)
                  do m = NeighborLower(1), NeighborUpper(1)
                    Neighbor = [m,n,o]
                    NeighborLabel = ComponentLabels%values(Neighbor(1),Neighbor(2),Neighbor(3))
                    if (NeighborLabel > 0) then
                      ComponentLabels%values(i,j,k) = NeighborLabel
                      exit L1
                    end if
                  end do
                end do
              end do L1
            else
              L2: &
              do o = NeighborLower(3), NeighborUpper(3)
                do n = NeighborLower(2), NeighborUpper(2)
                  do m = NeighborLower(1), NeighborUpper(1)
                    Neighbor = [m,n,o]
                    Neighbor(:Cart%nd) = ovkCartPeriodicAdjust(Cart, Neighbor)
                    if (ovkCartContains(Cart, Neighbor)) then
                      NeighborLabel = ComponentLabels%values(Neighbor(1),Neighbor(2),Neighbor(3))
                      ComponentLabels%values(i,j,k) = NeighborLabel
                      exit L2
                    end if
                  end do
                end do
              end do L2
            end if
            if (ComponentLabels%values(i,j,k) == 0) then
              NumComponents = NumComponents + 1
              ComponentLabels%values(i,j,k) = NumComponents
            end if
          end if
        end do
      end do
    end do

    allocate(ReducedComponentLabel(NumComponents))

    do i = 1, NumComponents
      ReducedComponentLabel(i) = i
    end do

    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          if (Mask%values(i,j,k)) then
            Point = [i,j,k]
            Label = ComponentLabels%values(i,j,k)
            NeighborLower(:Cart%nd) = Point(:Cart%nd)-1
            NeighborLower(Cart%nd+1:) = Point(Cart%nd+1:)
            NeighborUpper(:Cart%nd) = Point(:Cart%nd)+1
            NeighborUpper(Cart%nd+1:) = Point(Cart%nd+1:)
            AwayFromEdge = ovkCartContains(Cart, NeighborLower) .and. &
              ovkCartContains(Cart, NeighborUpper)
            if (AwayFromEdge) then
              do o = NeighborLower(3), NeighborUpper(3)
                do n = NeighborLower(2), NeighborUpper(2)
                  do m = NeighborLower(1), NeighborUpper(1)
                    Neighbor = [m,n,o]
                    NeighborLabel = ComponentLabels%values(Neighbor(1),Neighbor(2),Neighbor(3))
                    if (NeighborLabel > 0) then
                      ReducedComponentLabel(Label) = min(ReducedComponentLabel(Label), &
                        ReducedComponentLabel(NeighborLabel))
                    end if
                  end do
                end do
              end do
            else
              do o = NeighborLower(3), NeighborUpper(3)
                do n = NeighborLower(2), NeighborUpper(2)
                  do m = NeighborLower(1), NeighborUpper(1)
                    Neighbor = [m,n,o]
                    Neighbor(:Cart%nd) = ovkCartPeriodicAdjust(Cart, Neighbor)
                    if (ovkCartContains(Cart, Neighbor)) then
                      NeighborLabel = ComponentLabels%values(Neighbor(1),Neighbor(2),Neighbor(3))
                      if (NeighborLabel > 0) then
                        ReducedComponentLabel(Label) = min(ReducedComponentLabel(Label), &
                          ReducedComponentLabel(NeighborLabel))
                      end if
                    end if
                  end do
                end do
              end do
            end if
          end if
        end do
      end do
    end do

    do i = 1, NumComponents
      j = i
      do while (ReducedComponentLabel(j) /= j)
        j = ReducedComponentLabel(j)
      end do
      ReducedComponentLabel(i) = j
    end do

    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          Label = ComponentLabels%values(i,j,k)
          if (Label > 0) then
            ComponentLabels%values(i,j,k) = ReducedComponentLabel(Label)
          end if
        end do
      end do
    end do

    NumComponents = maxval(ReducedComponentLabel)

  end subroutine ovkConnectedComponents

  subroutine ovkFloodFill(Mask, BarrierMask)

    type(ovk_field_logical), intent(inout) :: Mask
    type(ovk_field_logical), intent(in) :: BarrierMask

    integer :: i, j, k
    type(ovk_field_logical) :: NonBarrierMask
    integer :: NumComponents
    type(ovk_field_int) :: ComponentLabels
    logical, dimension(:), allocatable :: IsFillComponent
    integer :: Label

    NonBarrierMask = ovk_field_logical_(Mask%cart)
    NonBarrierMask%values = .not. BarrierMask%values

    call ovkConnectedComponents(NonBarrierMask, NumComponents, ComponentLabels)

    allocate(IsFillComponent(NumComponents))
    IsFillComponent = .false.

    do k = Mask%cart%is(3), Mask%cart%ie(3)
      do j = Mask%cart%is(2), Mask%cart%ie(2)
        do i = Mask%cart%is(1), Mask%cart%ie(1)
          if (Mask%values(i,j,k) .and. .not. BarrierMask%values(i,j,k)) then
            Label = ComponentLabels%values(i,j,k)
            IsFillComponent(Label) = .true.
          end if
        end do
      end do
    end do

    do k = Mask%cart%is(3), Mask%cart%ie(3)
      do j = Mask%cart%is(2), Mask%cart%ie(2)
        do i = Mask%cart%is(1), Mask%cart%ie(1)
          if (.not. BarrierMask%values(i,j,k)) then
            Label = ComponentLabels%values(i,j,k)
            if (IsFillComponent(Label)) then
              Mask%values(i,j,k) = .true.
            end if
          end if
        end do
      end do
    end do

  end subroutine ovkFloodFill

  subroutine ovkThreshold_Integer(Field, ThresholdMask, Lower, Upper)

    type(ovk_field_int), intent(in) :: Field
    type(ovk_field_logical), intent(out) :: ThresholdMask
    integer, intent(in), optional :: Lower, Upper

    integer :: Lower_, Upper_
    integer :: i, j, k
    real(rk) :: Value

    if (present(Lower)) then
      Lower_ = Lower
    else
      Lower_ = -huge(0)
    end if

    if (present(Upper)) then
      Upper_ = Upper
    else
      Upper_ = huge(0)
    end if

    ThresholdMask = ovk_field_logical_(Field%cart, .false.)

    do k = Field%cart%is(3), Field%cart%ie(3)
      do j = Field%cart%is(2), Field%cart%ie(2)
        do i = Field%cart%is(1), Field%cart%ie(1)
          Value = Field%values(i,j,k)
          ThresholdMask%values(i,j,k) = Value >= Lower_ .and. Value <= Upper_
        end do
      end do
    end do

  end subroutine ovkThreshold_Integer

  subroutine ovkThreshold_LargeInteger(Field, ThresholdMask, Lower, Upper)

    type(ovk_field_large_int), intent(in) :: Field
    type(ovk_field_logical), intent(out) :: ThresholdMask
    integer(lk), intent(in), optional :: Lower, Upper

    integer(lk) :: Lower_, Upper_
    integer :: i, j, k
    real(rk) :: Value

    if (present(Lower)) then
      Lower_ = Lower
    else
      Lower_ = -huge(0_lk)
    end if

    if (present(Upper)) then
      Upper_ = Upper
    else
      Upper_ = huge(0_lk)
    end if

    ThresholdMask = ovk_field_logical_(Field%cart, .false.)

    do k = Field%cart%is(3), Field%cart%ie(3)
      do j = Field%cart%is(2), Field%cart%ie(2)
        do i = Field%cart%is(1), Field%cart%ie(1)
          Value = Field%values(i,j,k)
          ThresholdMask%values(i,j,k) = Value >= Lower_ .and. Value <= Upper_
        end do
      end do
    end do

  end subroutine ovkThreshold_LargeInteger

  subroutine ovkThreshold_Real(Field, ThresholdMask, Lower, Upper)

    type(ovk_field_real), intent(in) :: Field
    type(ovk_field_logical), intent(out) :: ThresholdMask
    real(rk), intent(in), optional :: Lower, Upper

    real(rk) :: Lower_, Upper_
    integer :: i, j, k
    real(rk) :: Value

    if (present(Lower)) then
      Lower_ = Lower
    else
      Lower_ = -huge(0._rk)
    end if

    if (present(Upper)) then
      Upper_ = Upper
    else
      Upper_ = huge(0._rk)
    end if

    ThresholdMask = ovk_field_logical_(Field%cart, .false.)

    do k = Field%cart%is(3), Field%cart%ie(3)
      do j = Field%cart%is(2), Field%cart%ie(2)
        do i = Field%cart%is(1), Field%cart%ie(1)
          Value = Field%values(i,j,k)
          ThresholdMask%values(i,j,k) = Value >= Lower_ .and. Value <= Upper_
        end do
      end do
    end do

  end subroutine ovkThreshold_Real

  subroutine ovkDistanceField(Mask, MaxDistance, BoundaryValue, Distances)

    type(ovk_field_logical), intent(in) :: Mask
    integer, intent(in) :: MaxDistance
    integer, intent(in) :: BoundaryValue
    type(ovk_field_int), intent(out) :: Distances

    integer :: i, j, k, d
    type(ovk_field_logical) :: NonMask
    integer :: NonMaskBoundaryValue
    type(ovk_field_logical) :: CoverMask
    type(ovk_field_logical) :: CoverEdgeMask

    Distances = ovk_field_int_(Mask%cart)

    NonMask = ovk_field_logical_(Mask%cart)
    NonMask%values = .not. Mask%values

    select case (BoundaryValue)
    case (OVK_TRUE)
      NonMaskBoundaryValue = OVK_FALSE
    case (OVK_FALSE)
      NonMaskBoundaryValue = OVK_TRUE
    case (OVK_MIRROR)
      NonMaskBoundaryValue = OVK_MIRROR
    end select

    call ovkDetectEdge(NonMask, OVK_OUTER_EDGE, NonMaskBoundaryValue, .true., CoverMask)

    do k = Mask%cart%is(3), Mask%cart%ie(3)
      do j = Mask%cart%is(2), Mask%cart%ie(2)
        do i = Mask%cart%is(1), Mask%cart%ie(1)
          if (Mask%values(i,j,k)) then
            if (CoverMask%values(i,j,k)) then
              Distances%values(i,j,k) = 0
            else
              Distances%values(i,j,k) = -MaxDistance
            end if
          else
            Distances%values(i,j,k) = MaxDistance
          end if
        end do
      end do
    end do

    do d = 1, MaxDistance-1
      call ovkDetectEdge(CoverMask, OVK_OUTER_EDGE, OVK_FALSE, .false., CoverEdgeMask)
      do k = Mask%cart%is(3), Mask%cart%ie(3)
        do j = Mask%cart%is(2), Mask%cart%ie(2)
          do i = Mask%cart%is(1), Mask%cart%ie(1)
            if (CoverEdgeMask%values(i,j,k)) then
              if (Mask%values(i,j,k)) then
                Distances%values(i,j,k) = -d
              else
                Distances%values(i,j,k) = d
              end if
            end if
          end do
        end do
      end do
      CoverMask%values = CoverMask%values .or. CoverEdgeMask%values
    end do

  end subroutine ovkDistanceField

  function ovkCountMask(Mask) result(NumTrue)

    type(ovk_field_logical), intent(in) :: Mask
    integer(lk) :: NumTrue

    integer :: i, j, k

    NumTrue = 0_lk
    do k = Mask%cart%is(3), Mask%cart%ie(3)
      do j = Mask%cart%is(2), Mask%cart%ie(2)
        do i = Mask%cart%is(1), Mask%cart%ie(1)
          if (Mask%values(i,j,k)) then
            NumTrue = NumTrue + 1_lk
          end if
        end do
      end do
    end do

  end function ovkCountMask

  subroutine ovkFilterState(State, FilterFlags, FilterCondition, FilterMask)

    type(ovk_field_int), intent(in) :: State
    integer, intent(in) :: FilterFlags
    integer, intent(in) :: FilterCondition
    type(ovk_field_logical), intent(out) :: FilterMask

    integer :: i, j, k

    FilterMask = ovk_field_logical_(State%cart)

    select case (FilterCondition)
    case (OVK_NONE)
      do k = State%cart%is(3), State%cart%ie(3)
        do j = State%cart%is(2), State%cart%ie(2)
          do i = State%cart%is(1), State%cart%ie(1)
            FilterMask%values(i,j,k) = iand(State%values(i,j,k), FilterFlags) == 0
          end do
        end do
      end do
    case (OVK_ANY)
      do k = State%cart%is(3), State%cart%ie(3)
        do j = State%cart%is(2), State%cart%ie(2)
          do i = State%cart%is(1), State%cart%ie(1)
            FilterMask%values(i,j,k) = iand(State%values(i,j,k), FilterFlags) /= 0
          end do
        end do
      end do
    case (OVK_ALL)
      do k = State%cart%is(3), State%cart%ie(3)
        do j = State%cart%is(2), State%cart%ie(2)
          do i = State%cart%is(1), State%cart%ie(1)
            FilterMask%values(i,j,k) = iand(State%values(i,j,k), FilterFlags) == FilterFlags
          end do
        end do
      end do
    end select

  end subroutine ovkFilterState

end module ovkFieldOps
