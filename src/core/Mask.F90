! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkMask

  use ovkCart
  use ovkField
  use ovkGlobal
  implicit none

  private

  ! API
  public :: ovkFindMaskEdge
  public :: ovkGrowMask
  public :: ovkGenerateExteriorMask
  public :: ovkGenerateNearEdgeMask
  public :: ovkCountMask
  public :: ovkMaskToIBlank
  public :: ovkPrintMask
  public :: OVK_EDGE_TYPE_INNER, OVK_EDGE_TYPE_OUTER

  integer, parameter :: OVK_EDGE_TYPE_INNER = 1
  integer, parameter :: OVK_EDGE_TYPE_OUTER = 2

contains

  subroutine ovkFindMaskEdge(Mask, EdgeType, EdgeMask, BoundaryValue)

    type(ovk_field_logical), intent(in) :: Mask
    integer, intent(in) :: EdgeType
    type(ovk_field_logical), intent(out) :: EdgeMask
    logical, intent(in), optional :: BoundaryValue

    logical :: BoundaryValue_
    integer :: i, j, k, m, n, o
    type(ovk_cart) :: Cart
    type(ovk_cart) :: EdgeCart
    logical :: InnerValue
    integer, dimension(MAX_ND) :: Point
    logical :: Value
    integer, dimension(MAX_ND) :: NeighborLower, NeighborUpper
    logical :: AwayFromBoundary
    integer, dimension(MAX_ND) :: Neighbor
    logical :: NeighborValue

    if (present(BoundaryValue)) then
      BoundaryValue_ = BoundaryValue
    else
      BoundaryValue_ = .false.
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

    InnerValue = EdgeType == OVK_EDGE_TYPE_INNER

    do k = EdgeCart%is(3), EdgeCart%ie(3)
      do j = EdgeCart%is(2), EdgeCart%ie(2)
        do i = EdgeCart%is(1), EdgeCart%ie(1)
          Point = [i,j,k]
          if (ovkCartContains(Cart, Point)) then
            Value = Mask%values(Point(1),Point(2),Point(3))
          else
            Value = BoundaryValue_
          end if
          if (Value .eqv. InnerValue) then
            NeighborLower(:Cart%nd) = Point(:Cart%nd)-1
            NeighborLower(Cart%nd+1:) = Point(Cart%nd+1:)
            NeighborUpper(:Cart%nd) = Point(:Cart%nd)+1
            NeighborUpper(Cart%nd+1:) = Point(Cart%nd+1:)
            AwayFromBoundary = ovkCartContains(Cart, NeighborLower) .and. &
              ovkCartContains(Cart, NeighborUpper)
            if (AwayFromBoundary) then
          L1: do o = NeighborLower(3), NeighborUpper(3)
                do n = NeighborLower(2), NeighborUpper(2)
                  do m = NeighborLower(1), NeighborUpper(1)
                    Neighbor = [m,n,o]
                    NeighborValue = Mask%values(Neighbor(1),Neighbor(2),Neighbor(3))
                    if (NeighborValue .neqv. InnerValue) then
                      EdgeMask%values(i,j,k) = .true.
                      exit L1
                    end if
                  end do
                end do
              end do L1
            else
          L2: do o = NeighborLower(3), NeighborUpper(3)
                do n = NeighborLower(2), NeighborUpper(2)
                  do m = NeighborLower(1), NeighborUpper(1)
                    Neighbor = [m,n,o]
                    Neighbor(:Cart%nd) = ovkCartPeriodicAdjust(Cart, Neighbor)
                    if (ovkCartContains(Cart, Neighbor)) then
                      NeighborValue = Mask%values(Neighbor(1),Neighbor(2),Neighbor(3))
                    else
                      NeighborValue = BoundaryValue_
                    end if
                    if (NeighborValue .neqv. InnerValue) then
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
    type(ovk_field_logical) :: EdgeMask
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: FillLower, FillUpper
    logical :: AwayFromBoundary
    integer, dimension(MAX_ND) :: FillPoint

    if (Amount == 0) return

    FillDistance = abs(Amount)
    FillValue = Amount > 0

    EdgeType = merge(OVK_EDGE_TYPE_INNER, OVK_EDGE_TYPE_OUTER, Amount > 0)

    Cart = Mask%cart

    call ovkFindMaskEdge(Mask, EdgeType, EdgeMask, BoundaryValue=FillValue)

    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          if (EdgeMask%values(i,j,k)) then
            Point = [i,j,k]
            FillLower(:Cart%nd) = Point(:Cart%nd)-FillDistance
            FillLower(Cart%nd+1:) = Point(Cart%nd+1:)
            FillUpper(:Cart%nd) = Point(:Cart%nd)+FillDistance
            FillUpper(Cart%nd+1:) = Point(Cart%nd+1:)
            AwayFromBoundary = ovkCartContains(Cart, FillLower) .and. &
              ovkCartContains(Cart, FillUpper)
            if (AwayFromBoundary) then
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

  subroutine ovkGenerateExteriorMask(Mask, ExteriorMask, EdgeMask)

    type(ovk_field_logical), intent(in) :: Mask
    type(ovk_field_logical), intent(out) :: ExteriorMask
    type(ovk_field_logical), intent(in), optional :: EdgeMask

    integer :: i, j, k, d, dOther1, dOther2
    type(ovk_field_int) :: State
    integer, dimension(MAX_ND) :: Point
    integer :: PrevState
    integer, parameter :: UNKNOWN = 1, INTERIOR = 2, EDGE = 3, EXTERIOR = 4

    if (present(EdgeMask)) then

      State = ovk_field_int_(Mask%cart)

      do k = Mask%cart%is(3), Mask%cart%ie(3)
        do j = Mask%cart%is(2), Mask%cart%ie(2)
          do i = Mask%cart%is(1), Mask%cart%ie(1)
            if (EdgeMask%values(i,j,k)) then
              State%values(i,j,k) = EDGE
            else if (Mask%values(i,j,k)) then
              State%values(i,j,k) = INTERIOR
            else
              State%values(i,j,k) = UNKNOWN
            end if
          end do
        end do
      end do

      do d = 1, MAX_ND
        dOther1 = modulo(d, MAX_ND) + 1
        dOther2 = modulo(d+1, MAX_ND) + 1
        do k = Mask%cart%is(dOther2), Mask%cart%ie(dOther2)
          do j = Mask%cart%is(dOther1), Mask%cart%ie(dOther1)
            Point(d) = Mask%cart%is(d)
            Point(dOther1) = j
            Point(dOther2) = k
            PrevState = State%values(Point(1),Point(2),Point(3))
            if (PrevState == UNKNOWN) then
              do i = Mask%cart%is(d)+1, Mask%cart%ie(d)
                Point(d) = i
                Point(dOther1) = j
                Point(dOther2) = k
                select case (State%values(Point(1),Point(2),Point(3)))
                case (INTERIOR)
                  PrevState = INTERIOR
                case (EDGE, EXTERIOR)
                  PrevState = EXTERIOR
                end select
                if (PrevState /= UNKNOWN) then
                  exit
                end if
              end do
            end if
            do i = Mask%cart%is(d), Mask%cart%ie(d)
              Point(d) = i
              Point(dOther1) = j
              Point(dOther2) = k
              if (State%values(Point(1),Point(2),Point(3)) == UNKNOWN) then
                select case (PrevState)
                case (INTERIOR)
                  State%values(Point(1),Point(2),Point(3)) = INTERIOR
                case (EDGE, EXTERIOR)
                  State%values(Point(1),Point(2),Point(3)) = EXTERIOR
                end select
              end if
              PrevState = State%values(Point(1),Point(2),Point(3))
            end do
          end do
        end do
      end do

      ExteriorMask = ovk_field_logical_(Mask%cart)
      ExteriorMask%values = State%values == EXTERIOR

    else

      ExteriorMask%values = .not. Mask%values

    end if

  end subroutine ovkGenerateExteriorMask

  subroutine ovkGenerateNearEdgeMask(Mask, EdgeType, EdgeDistance, NearEdgeMask)

    type(ovk_field_logical), intent(in) :: Mask
    integer, intent(in) :: EdgeType
    integer, intent(in) :: EdgeDistance
    type(ovk_field_logical), intent(out) :: NearEdgeMask

    type(ovk_field_logical) :: EdgeMask

    call ovkFindMaskEdge(Mask, EdgeType, EdgeMask)
    call ovkGrowMask(EdgeMask, EdgeDistance)

    NearEdgeMask = ovk_field_logical_(Mask%cart)
    NearEdgeMask%values = EdgeMask%values(Mask%cart%is(1):Mask%cart%ie(1), &
      Mask%cart%is(2):Mask%cart%ie(2),Mask%cart%is(3):Mask%cart%ie(3))

  end subroutine ovkGenerateNearEdgeMask

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

  subroutine ovkMaskToIBlank(Mask, IBlank, TrueValue, FalseValue)

    type(ovk_field_logical), intent(in) :: Mask
    type(ovk_field_int), intent(inout) :: IBlank
    integer, intent(in), optional :: TrueValue
    integer, intent(in), optional :: FalseValue

    if (present(TrueValue)) then
      IBlank%values = merge(TrueValue, IBlank%values, Mask%values)
    end if

    if (present(FalseValue)) then
      IBlank%values = merge(IBlank%values, FalseValue, Mask%values)
    end if

  end subroutine ovkMaskToIBlank

  subroutine ovkPrintMask(Mask, StartIndex, EndIndex)

    type(ovk_field_logical), intent(in) :: Mask
    integer, dimension(Mask%cart%nd), intent(in), optional :: StartIndex, EndIndex

    integer, dimension(MAX_ND) :: StartIndex_, EndIndex_
    integer :: i, j
    integer :: NormalDir
    integer :: IDir, JDir
    integer, dimension(2) :: StartSlice, EndSlice
    integer :: Slice
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: AdjustedPoint

    if (present(StartIndex)) then
      StartIndex_(:Mask%cart%nd) = StartIndex
      StartIndex_(Mask%cart%nd+1:) = 1
    else
      StartIndex_ = Mask%cart%is
    end if

    if (present(EndIndex)) then
      EndIndex_(:Mask%cart%nd) = EndIndex
      EndIndex_(Mask%cart%nd+1:) = 1
    else
      EndIndex_ = Mask%cart%ie
      EndIndex_(MAX_ND) = Mask%cart%is(MAX_ND)
    end if

    do i = 1, MAX_ND
      if (EndIndex_(i) == StartIndex_(i)) then
        exit
      end if
    end do
    NormalDir = i

    IDir = modulo(NormalDir,3) + 1
    JDir = modulo(NormalDir+1,3) + 1
    StartSlice = [StartIndex_(IDir),StartIndex_(JDir)]
    EndSlice = [EndIndex_(IDir),EndIndex_(JDir)]
    if (Mask%cart%nd == 3) then
      Slice = StartIndex_(NormalDir)
    else
      Slice = 1
    end if

    write (*, '(a)', advance='no') "   "
    do i = StartSlice(1), EndSlice(1)
      write (*, '(a)', advance='no') " - "
    end do
    write (*, '(a)') "   "

    do j = EndSlice(2), StartSlice(2), -1

      write (*, '(a)', advance='no') " | "

      do i = StartSlice(1), EndSlice(1)
        select case (NormalDir)
        case (1)
          Point = [Slice,i,j]
        case (2)
          Point = [j,Slice,i]
        case (3)
          Point = [i,j,Slice]
        end select
        AdjustedPoint(:Mask%cart%nd) = ovkCartPeriodicAdjust(Mask%cart, Point)
        AdjustedPoint(Mask%cart%nd+1:) = 1
        write(*, '(a)', advance='no') merge(" X ", "   ", Mask%values(AdjustedPoint(1), &
          AdjustedPoint(2), AdjustedPoint(3)))
      end do

      write (*, '(a)') " | "

    end do

    write (*, '(a)', advance='no') "   "
    do i = StartSlice(1), EndSlice(1)
      write (*, '(a)', advance='no') " - "
    end do
    write (*, '(a)') "   "
    write (*, '(a)') ""

  end subroutine ovkPrintMask

end module ovkMask
