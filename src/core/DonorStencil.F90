! Copyright (c) 2018 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkDonorStencil

  use ovkCart
  use ovkField
  use ovkFieldOps
  use ovkGeometryOps
  use ovkGlobal
  use ovkGrid
  use ovkLogger
  use ovkOverlap
  implicit none

  private

  ! Internal API
  public :: t_donor_stencil
  public :: t_donor_stencil_
  public :: CreateDonorStencil
  public :: DestroyDonorStencil
  public :: GetDonorStencilSize
  public :: FindDonors

  type t_donor_stencil
    type(t_noconstruct) :: noconstruct
    type(ovk_grid), pointer :: grid
    integer :: nd
    integer :: connection_type
    type(t_cubic), pointer :: cubic
  end type t_donor_stencil

  type t_cubic
    type(ovk_field_logical) :: shift_mask
    type(ovk_field_int), dimension(MAX_DIMS) :: shift_amounts
  end type t_cubic

contains

  pure function t_donor_stencil_() result(DonorStencil)

    type(t_donor_stencil) :: DonorStencil

    nullify(DonorStencil%grid)
    DonorStencil%nd = 2
    DonorStencil%connection_type = OVK_CONNECTION_NEAREST

    nullify(DonorStencil%cubic)

  end function t_donor_stencil_

  subroutine CreateDonorStencil(DonorStencil, Grid, ConnectionType)

    type(t_donor_stencil), intent(out) :: DonorStencil
    type(ovk_grid), pointer, intent(in) :: Grid
    integer, intent(in) :: ConnectionType

    DonorStencil%grid => Grid
    DonorStencil%nd = Grid%nd
    DonorStencil%connection_type = ConnectionType

    select case (ConnectionType)
    case (OVK_CONNECTION_NEAREST)
      ! No extra data
    case (OVK_CONNECTION_LINEAR)
      ! No extra data
    case (OVK_CONNECTION_CUBIC)
      call CreateDonorStencilCubic(DonorStencil)
    case default
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid connection type."
        stop 1
      end if
    end select

  end subroutine CreateDonorStencil

  subroutine CreateDonorStencilCubic(DonorStencil)

    type(t_donor_stencil), intent(inout) :: DonorStencil

    type(ovk_grid), pointer :: Grid
    integer :: NumDims
    type(t_cubic), pointer :: Cubic
    integer, dimension(MAX_DIMS) :: StencilCellLower, StencilCellUpper
    logical, dimension(:,:,:), allocatable :: StencilCellMask

    Grid => DonorStencil%grid
    NumDims = DonorStencil%nd

    allocate(DonorStencil%cubic)
    Cubic => DonorStencil%cubic

    StencilCellLower(:NumDims) = -1
    StencilCellLower(NumDims+1:) = 0
    StencilCellUpper(:NumDims) = 1
    StencilCellUpper(NumDims+1:) = 0

    allocate(StencilCellMask(StencilCellLower(1):StencilCellUpper(1), &
      StencilCellLower(2):StencilCellUpper(2),StencilCellLower(3):StencilCellUpper(3)))
    StencilCellMask = .true.

    call GenerateShifts(Grid, StencilCellLower, StencilCellUpper, StencilCellMask, &
      Cubic%shift_mask, Cubic%shift_amounts)

  end subroutine CreateDonorStencilCubic

  subroutine DestroyDonorStencil(DonorStencil)

    type(t_donor_stencil), intent(inout) :: DonorStencil

    select case (DonorStencil%connection_type)
    case (OVK_CONNECTION_NEAREST)
      ! No extra data
    case (OVK_CONNECTION_LINEAR)
      ! No extra data
    case (OVK_CONNECTION_CUBIC)
      call DestroyDonorStencilCubic(DonorStencil)
    end select

  end subroutine DestroyDonorStencil

  subroutine DestroyDonorStencilCubic(DonorStencil)

    type(t_donor_stencil), intent(inout) :: DonorStencil

    deallocate(DonorStencil%cubic)

  end subroutine DestroyDonorStencilCubic

  subroutine GetDonorStencilSize(DonorStencil, NumPoints)

    type(t_donor_stencil), intent(in) :: DonorStencil
    integer, dimension(DonorStencil%nd), intent(out) :: NumPoints

    select case (DonorStencil%connection_type)
    case (OVK_CONNECTION_NEAREST)
      NumPoints = 1
    case (OVK_CONNECTION_LINEAR)
      NumPoints = 2
    case (OVK_CONNECTION_CUBIC)
      NumPoints = 4
    end select

  end subroutine GetDonorStencilSize

  subroutine FindDonors(DonorGrid, ReceiverGrid, Overlap, DonorStencil, NumConnections, &
    ReceiverPoints, OverlapIndices, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_stencil), intent(in) :: DonorStencil
    integer(lk), intent(in) :: NumConnections
    integer, dimension(:,:), intent(in) :: ReceiverPoints
    integer(lk), dimension(:), intent(in) :: OverlapIndices
    integer, dimension(:,:,:), intent(out) :: DonorExtents
    real(rk), dimension(:,:), intent(out) :: DonorCoords
    real(rk), dimension(:,:,:), intent(out) :: DonorInterpCoefs

    select case (DonorStencil%connection_type)
    case (OVK_CONNECTION_NEAREST)
      call FindDonorsNearest(DonorGrid, ReceiverGrid, Overlap, DonorStencil, NumConnections, &
        ReceiverPoints, OverlapIndices, DonorExtents, DonorCoords, DonorInterpCoefs)
    case (OVK_CONNECTION_LINEAR)
      call FindDonorsLinear(DonorGrid, ReceiverGrid, Overlap, DonorStencil, NumConnections, &
        ReceiverPoints, OverlapIndices, DonorExtents, DonorCoords, DonorInterpCoefs)
    case (OVK_CONNECTION_CUBIC)
      call FindDonorsCubic(DonorGrid, ReceiverGrid, Overlap, DonorStencil, NumConnections, &
        ReceiverPoints, OverlapIndices, DonorExtents, DonorCoords, DonorInterpCoefs)
    end select

  end subroutine FindDonors

  subroutine FindDonorsNearest(DonorGrid, ReceiverGrid, Overlap, DonorStencil, NumConnections, &
    ReceiverPoints, OverlapIndices, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_stencil), intent(in) :: DonorStencil
    integer(lk), intent(in) :: NumConnections
    integer, dimension(:,:), intent(in) :: ReceiverPoints
    integer(lk), dimension(:), intent(in) :: OverlapIndices
    integer, dimension(:,:,:), intent(out) :: DonorExtents
    real(rk), dimension(:,:), intent(out) :: DonorCoords
    real(rk), dimension(:,:,:), intent(out) :: DonorInterpCoefs

    integer :: d
    integer(lk) :: l
    integer :: NumDims
    integer, dimension(MAX_DIMS) :: ReferenceCell
    real(rk), dimension(DonorGrid%nd) :: ReferenceCellCoords
    integer, dimension(MAX_DIMS) :: NearestPoint

    NumDims = DonorGrid%nd

    do l = 1_lk, NumConnections
      ReferenceCell = Overlap%cells(:,OverlapIndices(l))
      ReferenceCellCoords = Overlap%coords(:,OverlapIndices(l))
      NearestPoint(:NumDims) = ReferenceCell(:NumDims) + int(ReferenceCellCoords + 0.5_rk)
      NearestPoint(NumDims+1:) = 1
      NearestPoint(:NumDims) = ovkCartPeriodicAdjust(DonorGrid%cart, NearestPoint)
      DonorExtents(:NumDims,1,l) = NearestPoint(:NumDims)
      DonorExtents(:NumDims,2,l) = NearestPoint(:NumDims)
      DonorCoords(:,l) = 0._rk
      DonorInterpCoefs(:,:NumDims,l) = 0._rk
      do d = 1, NumDims
        DonorInterpCoefs(1,d,l) = 1._rk
      end do
    end do

  end subroutine FindDonorsNearest

  subroutine FindDonorsLinear(DonorGrid, ReceiverGrid, Overlap, DonorStencil, NumConnections, &
    ReceiverPoints, OverlapIndices, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_stencil), intent(in) :: DonorStencil
    integer(lk), intent(in) :: NumConnections
    integer, dimension(:,:), intent(in) :: ReceiverPoints
    integer(lk), dimension(:), intent(in) :: OverlapIndices
    integer, dimension(:,:,:), intent(out) :: DonorExtents
    real(rk), dimension(:,:), intent(out) :: DonorCoords
    real(rk), dimension(:,:,:), intent(out) :: DonorInterpCoefs

    integer :: d
    integer(lk) :: l
    integer :: NumDims
    integer, dimension(MAX_DIMS) :: ReferenceCell
    real(rk), dimension(DonorGrid%nd) :: ReferenceCellCoords

    NumDims = DonorGrid%nd

    do l = 1_lk, NumConnections
      ReferenceCell = Overlap%cells(:,OverlapIndices(l))
      ReferenceCellCoords = Overlap%coords(:,OverlapIndices(l))
      DonorExtents(:NumDims,1,l) = ReferenceCell(:NumDims)
      DonorExtents(:NumDims,2,l) = ReferenceCell(:NumDims) + 1
      DonorCoords(:,l) = ReferenceCellCoords
      DonorInterpCoefs(:,:NumDims,l) = 0._rk
      do d = 1, NumDims
        DonorInterpCoefs(:2,d,l) = ovkInterpBasisLinear(DonorCoords(d,l))
      end do
    end do

  end subroutine FindDonorsLinear

  subroutine FindDonorsCubic(DonorGrid, ReceiverGrid, Overlap, DonorStencil, NumConnections, &
    ReceiverPoints, OverlapIndices, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_stencil), intent(in) :: DonorStencil
    integer(lk), intent(in) :: NumConnections
    integer, dimension(:,:), intent(in) :: ReceiverPoints
    integer(lk), dimension(:), intent(in) :: OverlapIndices
    integer, dimension(:,:,:), intent(out) :: DonorExtents
    real(rk), dimension(:,:), intent(out) :: DonorCoords
    real(rk), dimension(:,:,:), intent(out) :: DonorInterpCoefs

    integer :: d
    integer(lk) :: l
    integer :: NumDims
    type(t_cubic), pointer :: Cubic
    integer, dimension(MAX_DIMS) :: ReferenceCell
    real(rk), dimension(DonorGrid%nd) :: ReferenceCellCoords
    integer, dimension(MAX_DIMS) :: StencilShift
    logical :: Success
    integer :: NumWarnings
    character(len=STRING_LENGTH) :: DonorCellString
    character(len=STRING_LENGTH) :: ReceiverPointString
    character(len=STRING_LENGTH) :: DonorGridIDString
    character(len=STRING_LENGTH) :: ReceiverGridIDString

    NumDims = DonorGrid%nd
    Cubic => DonorStencil%cubic

    NumWarnings = 0

    do l = 1, NumConnections
      ReferenceCell = Overlap%cells(:,OverlapIndices(l))
      ReferenceCellCoords = Overlap%coords(:,OverlapIndices(l))
      Success = .false.
      StencilShift = 0
      if (.not. Cubic%shift_mask%values(ReferenceCell(1),ReferenceCell(2),ReferenceCell(3))) then
        Success = .true.
      else
        do d = 1, NumDims
          StencilShift(d) = Cubic%shift_amounts(d)%values(ReferenceCell(1),ReferenceCell(2), &
            ReferenceCell(3))
        end do
        Success = any(StencilShift /= 0)
      end if
      if (Success) then
        DonorExtents(:NumDims,1,l) = ReferenceCell(:NumDims) + StencilShift(:NumDims) - 1
        DonorExtents(:NumDims,1,l) = ovkCartPeriodicAdjust(DonorGrid%cart, DonorExtents(:,1,l))
        DonorExtents(:NumDims,2,l) = DonorExtents(:NumDims,1,l) + 3
        DonorCoords(:,l) = ReferenceCellCoords - real(StencilShift(:NumDims),kind=rk)
        DonorInterpCoefs(:,:NumDims,l) = 0._rk
        do d = 1, NumDims
          DonorInterpCoefs(:,d,l) = ovkInterpBasisCubic(DonorCoords(d,l))
        end do
      else
        if (DonorGrid%logger%log_errors) then
          DonorCellString = TupleToString(DonorExtents(:NumDims,1,l))
          ReceiverPointString = TupleToString(ReceiverPoints(:NumDims,l))
          DonorGridIDString = IntToString(DonorGrid%id)
          ReceiverGridIDString = IntToString(ReceiverGrid%id)
          if (NumWarnings <= 100) then
            write (DonorGrid%logger%error_file, '(10a)') "WARNING: Could not use cubic ", &
              "interpolation for donor cell ", trim(DonorCellString), " of grid ", &
              trim(DonorGridIDString), " corresponding to receiver point ", &
              trim(ReceiverPointString), " of grid ", trim(ReceiverGridIDString), &
              "; using linear instead."
            if (NumWarnings == 100) then
              write (DonorGrid%logger%error_file, '(a)') "WARNING: Further warnings suppressed."
            end if
            NumWarnings = NumWarnings + 1
          end if
        end if
        DonorExtents(:NumDims,1,l) = ReferenceCell(:NumDims)
        DonorExtents(:NumDims,2,l) = ReferenceCell(:NumDims) + 1
        DonorCoords(:,l) = ReferenceCellCoords
        DonorInterpCoefs(:,:NumDims,l) = 0._rk
        do d = 1, NumDims
          DonorInterpCoefs(:2,d,l) = ovkInterpBasisLinear(DonorCoords(d,l))
        end do
      end if
    end do

  end subroutine FindDonorsCubic

  subroutine GenerateShifts(Grid, StencilCellLower, StencilCellUpper, StencilCellMask, &
    ShiftMask, ShiftAmounts)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(MAX_DIMS), intent(in) :: StencilCellLower, StencilCellUpper
    logical, dimension(StencilCellLower(1):StencilCellUpper(1),StencilCellLower(2): &
      StencilCellUpper(2),StencilCellLower(3):StencilCellUpper(3)), intent(in) :: StencilCellMask
    type(ovk_field_logical), intent(out) :: ShiftMask
    type(ovk_field_int), dimension(MAX_DIMS), intent(out) :: ShiftAmounts

    integer :: d, i, j, k, m, n, o
    integer :: NumDims
    type(ovk_cart) :: CellCart
    integer :: MaxStencilCellOffset
    type(ovk_field_int) :: CellEdgeDists
    type(ovk_field_int), dimension(:), allocatable :: CellEdgeDistGradients
    integer, dimension(MAX_DIMS) :: ReferenceCell
    integer, dimension(MAX_DIMS) :: CellStart, CellEnd
    logical :: AwayFromEdge
    integer, dimension(MAX_DIMS) :: Cell
    logical :: CellExists
    type(ovk_field_logical) :: FitsMask
    integer, dimension(MAX_DIMS) :: MinShiftedReferenceCell, MaxShiftedReferenceCell
    integer, dimension(MAX_DIMS) :: Shift
    integer, dimension(MAX_DIMS) :: ShiftedReferenceCell
    logical, dimension(-StencilCellUpper(1):-StencilCellLower(1),-StencilCellUpper(2): &
      -StencilCellLower(2),-StencilCellUpper(3):-StencilCellLower(3)) :: ShiftAllowed
    real(rk), dimension(Grid%nd) :: Gradient
    real(rk), dimension(Grid%nd) :: GradientDir
    real(rk), dimension(Grid%nd) :: Displacement
    real(rk), dimension(Grid%nd) :: DisplacementDir
    real(rk) :: DirectionQuality, DistanceQuality
    real(rk) :: ShiftQuality
    integer, dimension(MAX_DIMS) :: BestShift
    real(rk) :: BestShiftQuality

    NumDims = Grid%nd
    CellCart = Grid%cell_cart

    MaxStencilCellOffset = 0
    MaxStencilCellOffset = max(MaxStencilCellOffset, -minval(StencilCellLower))
    MaxStencilCellOffset = max(MaxStencilCellOffset, maxval(StencilCellUpper))

    call ComputeCellEdgeDistances(Grid, CellEdgeDists)

    allocate(CellEdgeDistGradients(NumDims))
    call ComputeCellEdgeDistanceGradients(Grid, CellEdgeDists, CellEdgeDistGradients)

    ShiftMask = ovk_field_logical_(Grid%cell_cart)

    do k = CellCart%is(3), CellCart%ie(3)
      do j = CellCart%is(2), CellCart%ie(2)
        do i = CellCart%is(1), CellCart%ie(1)
          if (CellEdgeDists%values(i,j,k) <= MaxStencilCellOffset) then
            ReferenceCell = [i,j,k]
            CellStart = ReferenceCell + StencilCellLower
            CellEnd = ReferenceCell + StencilCellUpper
            AwayFromEdge = ovkCartContains(CellCart, CellStart) .and. &
              ovkCartContains(CellCart, CellEnd)
            if (AwayFromEdge) then
          L1: do o = StencilCellLower(3), StencilCellUpper(3)
                do n = StencilCellLower(2), StencilCellUpper(2)
                  do m = StencilCellLower(1), StencilCellUpper(1)
                    if (StencilCellMask(m,n,o)) then
                      Cell = [ReferenceCell(1)+m,ReferenceCell(2)+n,ReferenceCell(3)+o]
                      CellExists = Grid%cell_mask%values(Cell(1),Cell(2),Cell(3))
                      if (.not. CellExists) then
                        ShiftMask%values(i,j,k) = .true.
                        exit L1
                      end if
                    end if
                  end do
                end do
              end do L1
            else
          L2: do o = StencilCellLower(3), StencilCellUpper(3)
                do n = StencilCellLower(2), StencilCellUpper(2)
                  do m = StencilCellLower(1), StencilCellUpper(1)
                    if (StencilCellMask(m,n,o)) then
                      Cell = [ReferenceCell(1)+m,ReferenceCell(2)+n,ReferenceCell(3)+o]
                      Cell(:NumDims) = ovkCartPeriodicAdjust(CellCart, Cell)
                      if (ovkCartContains(CellCart, Cell)) then
                        CellExists = Grid%cell_mask%values(Cell(1),Cell(2),Cell(3))
                      else
                        CellExists = .false.
                      end if
                      if (.not. CellExists) then
                        ShiftMask%values(i,j,k) = .true.
                        exit L2
                      end if
                    end if
                  end do
                end do
              end do L2
            end if
          else
            ShiftMask%values(i,j,k) = .false.
          end if
        end do
      end do
    end do

    FitsMask = ovk_field_logical_(CellCart)
    FitsMask%values = Grid%cell_mask%values .and. .not. ShiftMask%values

    do d = 1, MAX_DIMS
      ShiftAmounts(d) = ovk_field_int_(CellCart, 0)
    end do

    do k = CellCart%is(3), CellCart%ie(3)
      do j = CellCart%is(2), CellCart%ie(2)
        do i = CellCart%is(1), CellCart%ie(1)
          if (ShiftMask%values(i,j,k)) then
            ReferenceCell = [i,j,k]
            MinShiftedReferenceCell = ReferenceCell - StencilCellUpper
            MaxShiftedReferenceCell = ReferenceCell - StencilCellLower
            AwayFromEdge = ovkCartContains(CellCart, MinShiftedReferenceCell) .and. &
              ovkCartContains(CellCart, MaxShiftedReferenceCell)
            if (AwayFromEdge) then
              do o = StencilCellLower(3), StencilCellUpper(3)
                do n = StencilCellLower(2), StencilCellUpper(2)
                  do m = StencilCellLower(1), StencilCellUpper(1)
                    Shift = [-m,-n,-o]
                    ShiftedReferenceCell = ReferenceCell + Shift
                    ShiftAllowed(Shift(1),Shift(2),Shift(3)) = StencilCellMask(m,n,o) .and. &
                      FitsMask%values(ShiftedReferenceCell(1),ShiftedReferenceCell(2), &
                      ShiftedReferenceCell(3))
                  end do
                end do
              end do
            else
              do o = StencilCellLower(3), StencilCellUpper(3)
                do n = StencilCellLower(2), StencilCellUpper(2)
                  do m = StencilCellLower(1), StencilCellUpper(1)
                    Shift = [-m,-n,-o]
                    ShiftedReferenceCell = ReferenceCell + Shift
                    ShiftedReferenceCell(:NumDims) = ovkCartPeriodicAdjust(CellCart, &
                      ShiftedReferenceCell)
                    if (ovkCartContains(CellCart, ShiftedReferenceCell)) then
                      ShiftAllowed(Shift(1),Shift(2),Shift(3)) = StencilCellMask(m,n,o) .and. &
                        FitsMask%values(ShiftedReferenceCell(1),ShiftedReferenceCell(2), &
                        ShiftedReferenceCell(3))
                    else
                      ShiftAllowed(Shift(1),Shift(2),Shift(3)) = .false.
                    end if
                  end do
                end do
              end do
            end if
            do d = 1, NumDims
              Gradient(d) = real(CellEdgeDistGradients(d)%values(ReferenceCell(1), &
                ReferenceCell(2),ReferenceCell(3)), kind=rk)
            end do
            GradientDir = Gradient/max(sqrt(sum(Gradient**2)),1._rk)
            BestShift = 0
            BestShiftQuality = -huge(0._rk)
            do o = StencilCellLower(3), StencilCellUpper(3)
              do n = StencilCellLower(2), StencilCellUpper(2)
                do m = StencilCellLower(1), StencilCellUpper(1)
                  Shift = [-m,-n,-o]
                  if (ShiftAllowed(Shift(1),Shift(2),Shift(3))) then
                    Displacement = real(Shift(:NumDims), kind=rk)
                    DisplacementDir = Displacement/max(sqrt(sum(Displacement**2)),1._rk)
                    DirectionQuality = 0.5_rk*(1._rk + dot_product(GradientDir, DisplacementDir))
                    ! No particular justification for this coefficient other than to make
                    !   qual(dot-1,dist) = qual(dot,dist+1) + some small amount (e.g. 1.e-12)
                    ! so that 0 shift is preferred over distance 1 shift in gradient direction
                    DistanceQuality = (0.5_rk+1.e-12_rk)*(MaxStencilCellOffset-maxval(abs(Shift)))
                    ShiftQuality = DirectionQuality + DistanceQuality
                    if (ShiftQuality > BestShiftQuality) then
                      BestShift = Shift
                      BestShiftQuality = ShiftQuality
                    end if
                  end if
                end do
              end do
            end do
            do d = 1, MAX_DIMS
              ShiftAmounts(d)%values(i,j,k) = BestShift(d)
            end do
          end if
        end do
      end do
    end do

  end subroutine GenerateShifts

  subroutine ComputeCellEdgeDistances(Grid, CellEdgeDists)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_field_int), intent(out) :: CellEdgeDists

    integer :: NumDims
    type(ovk_cart) :: CellCart
    type(ovk_cart) :: PartiallyExtendedCellCart
    type(ovk_cart) :: ExtendedCellCart
    type(ovk_field_logical) :: NotMask
    type(ovk_field_int) :: PartialCellEdgeDists

    NumDims = Grid%nd
    CellCart = Grid%cell_cart

    ! Add extra ghost points in case we want to take the gradient
    ! Can't add them in the periodic directions yet
    PartiallyExtendedCellCart = CellCart
    PartiallyExtendedCellCart%is(:NumDims) = CellCart%is(:NumDims) - &
      merge(0, 1, CellCart%periodic(:NumDims))
    PartiallyExtendedCellCart%ie(:NumDims) = CellCart%ie(:NumDims) + &
      merge(0, 1, CellCart%periodic(:NumDims))

    NotMask = ovk_field_logical_(PartiallyExtendedCellCart, .true.)
    NotMask%values(CellCart%is(1):CellCart%ie(1),CellCart%is(2):CellCart%ie(2), &
      CellCart%is(3):CellCart%ie(3)) = .not. Grid%cell_mask%values

    call ovkDistanceField(NotMask, OVK_TRUE, PartialCellEdgeDists)

    ! Now add ghost points in the periodic directions
    ExtendedCellCart = CellCart
    ExtendedCellCart%is(:NumDims) = ExtendedCellCart%is(:NumDims) - 1
    ExtendedCellCart%ie(:NumDims) = ExtendedCellCart%ie(:NumDims) + 1
    ExtendedCellCart%periodic = .false.

    CellEdgeDists = ovk_field_int_(ExtendedCellCart)
    CellEdgeDists%values(PartiallyExtendedCellCart%is(1):PartiallyExtendedCellCart%ie(1), &
      PartiallyExtendedCellCart%is(2):PartiallyExtendedCellCart%ie(2), &
      PartiallyExtendedCellCart%is(3):PartiallyExtendedCellCart%ie(3)) = PartialCellEdgeDists%values
    call ovkFieldPeriodicFill(CellEdgeDists, PartiallyExtendedCellCart)

  end subroutine ComputeCellEdgeDistances

  subroutine ComputeCellEdgeDistanceGradients(Grid, CellEdgeDists, CellEdgeDistGradients)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_field_int), intent(in) :: CellEdgeDists
    type(ovk_field_int), dimension(:), intent(out) :: CellEdgeDistGradients

    integer :: d, i, j, k
    integer, dimension(MAX_DIMS) :: PrevOffset, NextOffset
    integer, dimension(MAX_DIMS) :: PrevCell, NextCell
    integer :: PrevDistance, NextDistance

    do d = 1, Grid%nd
      CellEdgeDistGradients(d) = ovk_field_int_(Grid%cell_cart, 0)
      PrevOffset = 0
      PrevOffset(d) = -1
      NextOffset = 0
      NextOffset(d) = 1
      do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
        do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
          do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
            PrevCell = [i+PrevOffset(1),j+PrevOffset(2),k+PrevOffset(3)]
            NextCell = [i+NextOffset(1),j+NextOffset(2),k+NextOffset(3)]
            PrevDistance = CellEdgeDists%values(PrevCell(1),PrevCell(2),PrevCell(3))
            NextDistance = CellEdgeDists%values(NextCell(1),NextCell(2),NextCell(3))
            CellEdgeDistGradients(d)%values(i,j,k) = NextDistance - PrevDistance
          end do
        end do
      end do
    end do

  end subroutine ComputeCellEdgeDistanceGradients

end module ovkDonorStencil
