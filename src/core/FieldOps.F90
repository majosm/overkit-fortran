! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkFieldOps

  use ovkCart
  use ovkField
  use ovkGlobal
  implicit none

  private

  ! API
  public :: ovkFindMaskEdge
  public :: ovkGrowMask
  public :: ovkConnectedComponents
  public :: ovkFillMask
  public :: ovkDistanceField
  public :: ovkGenerateNearEdgeMask
  public :: ovkGenerateThresholdMask
  public :: ovkCountMask
  public :: OVK_EDGE_TYPE_INNER, OVK_EDGE_TYPE_OUTER

  integer, parameter :: OVK_EDGE_TYPE_INNER = 1
  integer, parameter :: OVK_EDGE_TYPE_OUTER = 2

contains

  subroutine ovkFindMaskEdge(Mask, EdgeType, OuterValue, EdgeMask)

    type(ovk_field_logical), intent(in) :: Mask
    integer, intent(in) :: EdgeType
    integer, intent(in) :: OuterValue
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

    if (OVK_DEBUG) then
      if (Mask%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    Cart = Mask%cart

    EdgeCart = Mask%cart
    if (EdgeType == OVK_EDGE_TYPE_OUTER) then
      EdgeCart%is(:EdgeCart%nd) = EdgeCart%is(:EdgeCart%nd) - merge(0, 1, &
        EdgeCart%periodic(:EdgeCart%nd))
      EdgeCart%ie(:EdgeCart%nd) = EdgeCart%ie(:EdgeCart%nd) + merge(0, 1, &
        EdgeCart%periodic(:EdgeCart%nd))
    end if

    EdgeMask = ovk_field_logical_(EdgeCart, .false.)

    ! Points on inner edge will have Mask == .true., points on outer edge will have Mask == .false.
    EdgeValue = EdgeType == OVK_EDGE_TYPE_INNER

    do k = EdgeCart%is(3), EdgeCart%ie(3)
      do j = EdgeCart%is(2), EdgeCart%ie(2)
        do i = EdgeCart%is(1), EdgeCart%ie(1)
          Point = [i,j,k]
          if (ovkCartContains(Cart, Point)) then
            Value = Mask%values(i,j,k)
          else
            select case (OuterValue)
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
                      select case (OuterValue)
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

  end subroutine ovkFindMaskEdge

  subroutine ovkGrowMask(Mask, Amount)

    type(ovk_field_logical), intent(inout) :: Mask
    integer, intent(in) :: Amount

    integer :: i, j, k, m, n, o
    type(ovk_cart) :: Cart
    integer :: FillDistance
    logical :: FillValue
    integer :: EdgeType
    integer :: OuterValue
    type(ovk_field_logical) :: EdgeMask
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: FillLower, FillUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: FillPoint

    if (OVK_DEBUG) then
      if (Mask%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    if (Amount == 0) return

    FillDistance = abs(Amount)
    FillValue = Amount > 0

    EdgeType = merge(OVK_EDGE_TYPE_INNER, OVK_EDGE_TYPE_OUTER, Amount > 0)
    OuterValue = merge(OVK_TRUE, OVK_FALSE, FillValue)

    Cart = Mask%cart

    call ovkFindMaskEdge(Mask, EdgeType, OuterValue, EdgeMask)

    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          if (EdgeMask%values(i,j,k)) then
            Point = [i,j,k]
            FillLower(:Cart%nd) = Point(:Cart%nd)-FillDistance
            FillLower(Cart%nd+1:) = Point(Cart%nd+1:)
            FillUpper(:Cart%nd) = Point(:Cart%nd)+FillDistance
            FillUpper(Cart%nd+1:) = Point(Cart%nd+1:)
            AwayFromEdge = ovkCartContains(Cart, FillLower) .and. &
              ovkCartContains(Cart, FillUpper)
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

  end subroutine ovkGrowMask

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

    if (OVK_DEBUG) then
      if (Mask%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

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

  subroutine ovkFillMask(Mask, BarrierMask)

    type(ovk_field_logical), intent(inout) :: Mask
    type(ovk_field_logical), intent(in) :: BarrierMask

    integer :: i, j, k
    type(ovk_field_logical) :: NonBarrierMask
    integer :: NumComponents
    type(ovk_field_int) :: ComponentLabels
    logical, dimension(:), allocatable :: IsFillComponent
    integer :: Label

    if (OVK_DEBUG) then
      if (Mask%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

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

  end subroutine ovkFillMask

  subroutine ovkDistanceField(Mask, MaxDistance, OuterValue, Distances)

    type(ovk_field_logical), intent(in) :: Mask
    integer, intent(in) :: MaxDistance
    integer, intent(in) :: OuterValue
    type(ovk_field_int), intent(out) :: Distances

    integer :: i, j, k, d
    type(ovk_field_logical) :: NonMask
    integer :: NonMaskOuterValue
    type(ovk_field_logical) :: CoverMask
    type(ovk_field_logical) :: PrevCoverMask

    if (OVK_DEBUG) then
      if (Mask%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    ! Doing this the inefficient way for now -- will eventually need to use a better algorithm

    Distances = ovk_field_int_(Mask%cart)

    NonMask = ovk_field_logical_(Mask%cart)
    NonMask%values = .not. Mask%values

    select case (OuterValue)
    case (OVK_TRUE)
      NonMaskOuterValue = OVK_FALSE
    case (OVK_FALSE)
      NonMaskOuterValue = OVK_TRUE
    case (OVK_MIRROR)
      NonMaskOuterValue = OVK_MIRROR
    end select

    call ovkFindMaskEdge(NonMask, OVK_EDGE_TYPE_OUTER, NonMaskOuterValue, CoverMask)

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

    PrevCoverMask = ovk_field_logical_(CoverMask%cart)
    do d = 1, MaxDistance-1
      PrevCoverMask%values = CoverMask%values
      call ovkGrowMask(CoverMask, 1)
      do k = Mask%cart%is(3), Mask%cart%ie(3)
        do j = Mask%cart%is(2), Mask%cart%ie(2)
          do i = Mask%cart%is(1), Mask%cart%ie(1)
            if (CoverMask%values(i,j,k) .and. .not. PrevCoverMask%values(i,j,k)) then
              if (Mask%values(i,j,k)) then
                Distances%values(i,j,k) = -d
              else
                Distances%values(i,j,k) = d
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine ovkDistanceField

  subroutine ovkGenerateNearEdgeMask(Mask, EdgeDistance, OuterValue, NearEdgeMask, IgnoredEdges)

    type(ovk_field_logical), intent(in) :: Mask
    integer, intent(in) :: EdgeDistance
    integer, intent(in) :: OuterValue
    type(ovk_field_logical), intent(out) :: NearEdgeMask
    type(ovk_field_logical), intent(in), optional :: IgnoredEdges

    type(ovk_field_logical) :: EdgeMask
    type(ovk_field_logical) :: NonSubsetMask
    type(ovk_field_logical) :: NonSubsetEdgeMask
    type(ovk_field_logical) :: SubsetMask
    type(ovk_field_logical) :: SubsetEdgeMask

    if (OVK_DEBUG) then
      if (Mask%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    call ovkFindMaskEdge(Mask, OVK_EDGE_TYPE_OUTER, OuterValue, EdgeMask)

    if (present(IgnoredEdges)) then

      NonSubsetMask = IgnoredEdges
      call ovkFindMaskEdge(NonSubsetMask, OVK_EDGE_TYPE_OUTER, OuterValue, NonSubsetEdgeMask)

      SubsetMask = ovk_field_logical_(Mask%cart)
      SubsetMask%values = Mask%values .and. .not. NonSubsetMask%values

      call ovkFindMaskEdge(SubsetMask, OVK_EDGE_TYPE_OUTER, OuterValue, SubsetEdgeMask)

      EdgeMask%values = EdgeMask%values .and. .not. (NonSubsetEdgeMask%values .and. .not. &
        SubsetEdgeMask%values)

    end if

    call ovkGrowMask(EdgeMask, EdgeDistance)

    NearEdgeMask = ovk_field_logical_(Mask%cart)
    NearEdgeMask%values = Mask%values .and. EdgeMask%values(Mask%cart%is(1):Mask%cart%ie(1), &
      Mask%cart%is(2):Mask%cart%ie(2),Mask%cart%is(3):Mask%cart%ie(3))

  end subroutine ovkGenerateNearEdgeMask

  subroutine ovkGenerateThresholdMask(Field, ThresholdMask, Lower, Upper, Subset)

    type(ovk_field_real), intent(in) :: Field
    type(ovk_field_logical), intent(out) :: ThresholdMask
    real(rk), intent(in), optional :: Lower, Upper
    type(ovk_field_logical), intent(in), optional :: Subset

    real(rk) :: Lower_, Upper_
    integer :: i, j, k
    logical :: IncludePoint
    real(rk) :: Value

    if (OVK_DEBUG) then
      if (Field%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

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
          if (present(Subset)) then
            IncludePoint = Subset%values(i,j,k)
          else
            IncludePoint = .true.
          end if
          if (IncludePoint) then
            Value = Field%values(i,j,k)
            ThresholdMask%values(i,j,k) = Value >= Lower_ .and. Value <= Upper_
          end if
        end do
      end do
    end do

  end subroutine ovkGenerateThresholdMask

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

end module ovkFieldOps