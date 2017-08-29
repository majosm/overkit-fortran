! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkOverset

  use ovkAssembler
  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDomain
  use ovkDonorAccel
  use ovkDonors
  use ovkField
  use ovkFieldOps
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkInterp
  implicit none

  private

  ! API
  public :: ovkAssemble

contains

  subroutine ovkAssemble(Assembler)

    type(ovk_assembler), intent(inout) :: Assembler

    integer :: i, j, k, l, m, n, p, q, r
    integer :: ClockInitial, ClockFinal, ClockRate
    integer :: NumDims
    type(ovk_domain), pointer :: Domain
    integer :: NumGrids
    integer, dimension(:), allocatable :: IndexToID, IDToIndex
    type(ovk_grid), pointer :: Grid_m, Grid_n
    logical, dimension(:), allocatable :: InferBoundaries
    logical, dimension(:,:), allocatable :: Overlappable
    real(rk), dimension(:,:), allocatable :: OverlapTolerance
    logical, dimension(:,:), allocatable :: BoundaryHoleCutting
    logical, dimension(:,:), allocatable :: OverlapHoleCutting
    integer, dimension(:,:), allocatable :: ConnectionType
    logical, dimension(:,:), allocatable :: DisjointConnection
    integer, dimension(:), allocatable :: InterpScheme
    integer, dimension(:), allocatable :: FringeSize
    integer, dimension(:,:), allocatable :: FringePadding
    integer(lk) :: TotalPoints
    character(len=STRING_LENGTH) :: TotalPointsString
    character(len=STRING_LENGTH) :: iSString, iEString, jSString, jEString, kSString, kEString
    type(ovk_bbox) :: Bounds
    type(ovk_bbox), dimension(:), allocatable :: OverlapBounds
    real(rk), dimension(:), allocatable :: MaxOverlapTolerance
    integer(lk) :: PointCount
    type(ovk_donor_accel) :: DonorAccel
    type(ovk_donors), dimension(:,:), allocatable :: PairwiseDonors
    type(ovk_field_logical), dimension(:), allocatable :: OverlapMasks
    type(ovk_field_logical), dimension(:), allocatable :: CutMasks
    type(ovk_field_logical), dimension(:), allocatable :: FineMasks
    type(ovk_field_logical), dimension(:,:), allocatable :: PairwiseFineMasks
    type(ovk_field_logical) :: CoarseMask
    type(ovk_field_logical) :: PairwiseCoarseMask
    type(ovk_field_logical) :: CoarseEdgeMask
    type(ovk_field_logical) :: PaddingEdgeMask
    type(ovk_field_logical) :: NonGridMask
    type(ovk_field_logical) :: IgnoredEdgeMask
    type(ovk_field_logical), dimension(:), allocatable :: FringeMasks
    type(ovk_domain), pointer :: EditDomain
    type(ovk_grid), pointer :: EditGrid
    type(ovk_field_logical), pointer :: EditGridMask
    type(ovk_field_logical), pointer :: EditBoundaryMask
    type(ovk_field_logical), pointer :: EditInternalBoundaryMask
    type(ovk_field_logical), dimension(:), allocatable :: PaddingMasks
    integer :: PaddingAmount
    real(rk) :: CandidateCellDistance, BestCellDistance
    real(rk) :: CandidateCellDiff, BestCellDiff
    integer :: BestGrid
    type(ovk_field_logical) :: CoarseMask_m, CoarseMask_n
    type(ovk_field_logical) :: ReceiverMask_m, ReceiverMask_n
    type(ovk_field_logical) :: DonorMask
    type(ovk_field_logical) :: EdgeMask1
    type(ovk_field_logical) :: EdgeMask2
    type(ovk_field_logical) :: BoundaryMask
    type(ovk_field_logical) :: InteriorMask
    type(ovk_field_logical) :: ExteriorMask
    type(ovk_field_logical) :: ReceiverMask
    type(ovk_donors), dimension(:), allocatable :: Donors
    type(ovk_field_logical), dimension(:), allocatable :: OrphanMasks
    integer :: NumWarnings
    type(ovk_connectivity), pointer :: Connectivity
    type(ovk_grid_description) :: GridDescription
    integer :: StencilSize
    type(ovk_interp), pointer :: InterpData_n

    type(ovk_field_logical), dimension(:), allocatable :: ValidCellMasks
    type(ovk_field_logical) :: ValidCellInnerEdgeMask
    type(ovk_field_logical), dimension(:), allocatable :: ExpandableCellMasks
    type(ovk_field_int), dimension(:), allocatable :: CellQualities
    logical :: ExpandableCell
    integer, dimension(MAX_ND) :: ReceiverPoint
    integer, dimension(MAX_ND) :: DonorCell
    integer, dimension(MAX_ND) :: NeighborCellLower, NeighborCellUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: NeighborCell
    integer, dimension(MAX_ND) :: DonorCellShift
    integer :: BestCellQuality
    integer :: CellQuality
    integer, dimension(MAX_ND) :: ExpandedDonorCell
    real(rk), dimension(Assembler%properties%nd) :: ReceiverCoords
    real(rk), dimension(Assembler%properties%nd) :: DonorCellCoords
    real(rk), dimension(Assembler%properties%nd) :: ExpandedDonorCellCoords

    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    if (OVK_VERBOSE) then
      call system_clock(ClockInitial, ClockRate)
      write (*, '(a)') "Overset grid assembly started..."
    end if

    if (OVK_DEBUG) then
      if ( &
        Assembler%editor%properties_ref_count > 0 .or. &
        Assembler%editor%domain_ref_count > 0 .or. &
  !       Assembler%editor%overlap_ref_count > 0 .or. &
        Assembler%editor%connectivity_ref_count > 0 &
      ) then
        write (ERROR_UNIT, '(a)') "ERROR: Cannot perform assembly; assembler is still being edited."
        stop 1
      end if
      if (.not. Assembler%properties%manual_padding) then
        write (ERROR_UNIT, '(a)') "ERROR: Automatic padding is not currently supported."
        stop 1
      end if
    end if

    NumDims = Assembler%properties%nd
    Domain => Assembler%domain

    NumGrids = 0
    do n = 1, Assembler%properties%ngrids
      if (.not. ovkCartIsEmpty(Domain%grids(n)%cart)) then
        NumGrids = NumGrids + 1
      end if
    end do

    allocate(IndexToID(NumGrids))
    allocate(IDToIndex(Assembler%properties%ngrids))
    IndexToID = 0
    IDToIndex = 0
    m = 1
    do n = 1, Assembler%properties%ngrids
      if (.not. ovkCartIsEmpty(Domain%grids(n)%cart)) then
        IndexToID(m) = n
        IDToIndex(n) = m
        m = m + 1
      end if
    end do

    allocate(InferBoundaries(NumGrids))
    do n = 1, NumGrids
      q = IndexToID(n)
      InferBoundaries(n) = Assembler%properties%infer_boundaries(q)
    end do

    allocate(Overlappable(NumGrids,NumGrids))
    allocate(OverlapTolerance(NumGrids,NumGrids))
    do n = 1, NumGrids
      q = IndexToID(n)
      do m = 1, NumGrids
        p = IndexToID(m)
        Overlappable(m,n) = Assembler%properties%overlappable(p,q)
        OverlapTolerance(m,n) = Assembler%properties%overlap_tolerance(p,q)
      end do
    end do

    allocate(BoundaryHoleCutting(NumGrids,NumGrids))
    allocate(OverlapHoleCutting(NumGrids,NumGrids))
    do n = 1, NumGrids
      q = IndexToID(n)
      do m = 1, NumGrids
        p = IndexToID(m)
        BoundaryHoleCutting(m,n) = Assembler%properties%boundary_hole_cutting(p,q)
        OverlapHoleCutting(m,n) = Assembler%properties%overlap_hole_cutting(p,q)
      end do
    end do

    allocate(ConnectionType(NumGrids,NumGrids))
    allocate(DisjointConnection(NumGrids,NumGrids))
    allocate(InterpScheme(NumGrids))
    allocate(FringeSize(NumGrids))
    allocate(FringePadding(NumGrids,NumGrids))
    do n = 1, NumGrids
      q = IndexToID(n)
      do m = 1, NumGrids
        p = IndexToID(m)
        ConnectionType(m,n) = Assembler%properties%connection_type(p,q)
        DisjointConnection(m,n) = Assembler%properties%disjoint_connection(p,q)
        FringePadding(m,n) = Assembler%properties%fringe_padding(p,q)
      end do
      InterpScheme(n) = OVK_INTERP_LINEAR
      FringeSize(n) = 0
      do m = 1, NumGrids
        p = IndexToID(m)
        if (ConnectionType(m,n) == OVK_CONNECTION_FRINGE) then
          ! Interp scheme is currently assumed to be constant for each receiver grid
          InterpScheme(n) = Assembler%properties%interp_scheme(p,q)
          ! Fringe size is currently assumed to be constant for each receiver grid
          FringeSize(n) = Assembler%properties%fringe_size(p,q)
          exit
        end if
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Grid info:"
      write (*, '(3a)') "* Dimension: ", trim(IntToString(NumDims)), "D"
      write (*, '(2a)') "* Number of grids: ", trim(IntToString(NumGrids))
      TotalPoints = 0_lk
      do n = 1, NumGrids
        Grid_n => Domain%grids(IndexToID(n))
        TotalPoints = TotalPoints + ovkCartCount(Grid_n%cart)
      end do
      TotalPointsString = LargeIntToString(TotalPoints)
      write (*, '(2a)') "* Total number of grid points: ", trim(TotalPointsString)
      do n = 1, NumGrids
        Grid_n => Domain%grids(IndexToID(n))
        TotalPointsString = LargeIntToString(ovkCartCount(Grid_n%cart))
        iSString = IntToString(Grid_n%cart%is(1))
        iEString = IntToString(Grid_n%cart%ie(1))
        jSString = IntToString(Grid_n%cart%is(2))
        jEString = IntToString(Grid_n%cart%ie(2))
        kSString = IntToString(Grid_n%cart%is(3))
        kEString = IntToString(Grid_n%cart%ie(3))
        write (*, '(3a)', advance="no") "* Grid ", trim(IntToString(Grid_n%properties%id)), ": "
        write (*, '(2a)', advance="no") trim(TotalPointsString), " points "
        select case (NumDims)
        case (2)
          write (*, '(9a)', advance="no") "(i=", trim(iSString), ":", trim(iEString), &
            ", j=", trim(jSString), ":", trim(jEString), ")"
        case (3)
          write (*, '(13a)', advance="no") "(i=", trim(iSString), ":", trim(iEString), &
            ", j=", trim(jSString), ":", trim(jEString), ", k=", trim(kSString), ":", &
            trim(kEString), ")"
        end select
        write (*, '(a)') ""
      end do
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Searching for candidate donor/receiver pairs..."
    end if

    allocate(OverlapBounds(NumGrids))
    allocate(MaxOverlapTolerance(NumGrids))
    do m = 1, NumGrids
      OverlapBounds(m) = ovk_bbox_(NumDims)
      MaxOverlapTolerance(m) = 0._rk
      do n = 1, NumGrids
        if (Overlappable(m,n)) then
          Grid_m => Domain%grids(IndexToID(m))
          Grid_n => Domain%grids(IndexToID(n))
          Bounds = ovkBBIntersect(Grid_m%bounds, Grid_n%bounds)
          if (.not. ovkBBIsEmpty(Bounds)) then
            OverlapBounds(m) = ovkBBUnion(OverlapBounds(m), Bounds)
            MaxOverlapTolerance(m) = max(MaxOverlapTolerance(m), OverlapTolerance(m,n))
          else
            Overlappable(m,n) = .false.
            BoundaryHoleCutting(m,n) = .false.
            OverlapHoleCutting(m,n) = .false.
            ConnectionType(m,n) = OVK_CONNECTION_NONE
            Assembler%properties%overlappable(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%boundary_hole_cutting(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%overlap_hole_cutting(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%connection_type(IndexToID(m),IndexToID(n)) = OVK_CONNECTION_NONE
          end if
        end if
      end do
    end do

    allocate(PairwiseDonors(NumGrids,NumGrids))

    do m = 1, NumGrids
      Grid_m => Domain%grids(IndexToID(m))
      if (any(Overlappable(m,:))) then
        if (OVK_VERBOSE) then
          write (*, '(3a)') "* Generating donor search accelerator on grid ", &
            trim(IntToString(Grid_m%properties%id)), "..."
        end if
        call ovkGenerateDonorAccel(Grid_m, DonorAccel, Bounds=OverlapBounds(m), &
          OverlapTolerance=MaxOverlapTolerance(m))
      end if
      do n = 1, NumGrids
        Grid_n => Domain%grids(IndexToID(n))
        if (Overlappable(m,n)) then
          call ovkFindDonors(Grid_m, Grid_n, DonorAccel, PairwiseDonors(m,n), &
            OverlapTolerance=OverlapTolerance(m,n))
          if (OVK_VERBOSE) then
            PointCount = ovkCountMask(PairwiseDonors(m,n)%valid_mask)
            if (PointCount > 0) then
              write (*, '(7a)') "* ", trim(LargeIntToString(PointCount)), &
                " candidate donors from grid ", trim(IntToString(Grid_m%properties%id)), &
                " to grid ", trim(IntToString(Grid_n%properties%id)), " found."
            end if
          end if
        else
          call ovkMakeDonors(PairwiseDonors(m,n), Grid_n%cart)
          PairwiseDonors(m,n)%valid_mask%values = .false.
        end if
      end do
      if (any(Overlappable(m,:))) then
        call ovkDestroyDonorAccel(DonorAccel)
      end if
    end do

    deallocate(OverlapBounds)
    deallocate(MaxOverlapTolerance)

    allocate(OverlapMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      OverlapMasks(n) = ovk_field_logical_(Grid_n%cart, .false.)
      do m = 1, NumGrids
        if (Overlappable(m,n)) then
          OverlapMasks(n)%values = OverlapMasks(n)%values .or. PairwiseDonors(m,n)%valid_mask%values
        end if
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished searching for candidate donor/receiver pairs."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Inferring domain boundaries in non-overlapping regions..."
    end if

    call ovkEditAssemblerDomain(Assembler, EditDomain)

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      if (OVK_VERBOSE) then
        PointCount = 0_lk
      end if
      if (InferBoundaries(n)) then
        call ovkDetectEdge(Grid_n%grid_mask, OVK_INNER_EDGE, OVK_FALSE, .false., EdgeMask1)
        EdgeMask1%values = EdgeMask1%values .and. .not. OverlapMasks(n)%values
        call ovkEditDomainGrid(EditDomain, IndexToID(n), EditGrid)
        call ovkEditGridBoundaryMask(EditGrid, EditBoundaryMask)
        EditBoundaryMask%values = EditBoundaryMask%values .or. EdgeMask1%values
        call ovkReleaseGridBoundaryMask(EditGrid, EditBoundaryMask)
        call ovkReleaseDomainGrid(EditDomain, EditGrid)
        if (OVK_VERBOSE) then
          PointCount = ovkCountMask(EdgeMask1)
        end if
      end if
      if (OVK_VERBOSE) then
        if (PointCount > 0) then
          write (*, '(5a)') "* ", trim(LargeIntToString(PointCount)), &
            " points marked as boundaries on grid ", trim(IntToString(Grid_n%properties%id)), "."
        end if
      end if
    end do

    call ovkReleaseAssemblerDomain(Assembler, EditDomain)

    do m = 1, NumGrids
      do n = 1, NumGrids
        if (Overlappable(m,n)) then
          Grid_m => Domain%grids(IndexToID(m))
          Grid_n => Domain%grids(IndexToID(n))
          Bounds = ovkBBIntersect(Grid_m%bounds, Grid_n%bounds)
          if (ovkBBIsEmpty(Bounds)) then
            Overlappable(m,n) = .false.
            BoundaryHoleCutting(m,n) = .false.
            OverlapHoleCutting(m,n) = .false.
            ConnectionType(m,n) = OVK_CONNECTION_NONE
            Assembler%properties%overlappable(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%boundary_hole_cutting(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%overlap_hole_cutting(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%connection_type(IndexToID(m),IndexToID(n)) = OVK_CONNECTION_NONE
          end if
        end if
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished inferring domain boundaries in non-overlapping regions."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Cutting boundary holes..."
    end if

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      if (OVK_VERBOSE) then
        PointCount = 0_lk
      end if
      if (any(BoundaryHoleCutting(:,n))) then
        BoundaryMask = ovk_field_logical_(Grid_n%cart, .false.)
        InteriorMask = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          Grid_m => Domain%grids(IndexToID(m))
          if (BoundaryHoleCutting(m,n)) then
            call ovkDetectEdge(PairwiseDonors(m,n)%valid_mask, OVK_OUTER_EDGE, OVK_MIRROR, &
              .false., EdgeMask1)
            call ovkDetectEdge(PairwiseDonors(n,m)%valid_mask, OVK_INNER_EDGE, OVK_FALSE, &
              .false., EdgeMask2)
            call ovkGenerateDonorMask(Grid_n, Grid_m, PairwiseDonors(n,m), DonorMask, &
              ReceiverSubset=EdgeMask2)
            EdgeMask1%values = EdgeMask1%values .and. .not. DonorMask%values
            i = 0
            do while (any(EdgeMask1%values))
              call ovkDilate(DonorMask, 1, OVK_FALSE)
              EdgeMask1%values = EdgeMask1%values .and. .not. DonorMask%values
              i = i + 1
            end do
            call ovkGenerateDonorMask(Grid_n, Grid_m, PairwiseDonors(n,m), DonorMask, &
              ReceiverSubset=Grid_m%boundary_mask)
            do j = 1, i
              call ovkDilate(DonorMask, 1, OVK_FALSE)
            end do
            DonorMask%values = DonorMask%values .and. .not. PairwiseDonors(m,n)%valid_mask%values
            call ovkDetectEdge(DonorMask, OVK_OUTER_EDGE, OVK_FALSE, .false., EdgeMask1)
            EdgeMask1%values = EdgeMask1%values .and. PairwiseDonors(m,n)%valid_mask%values
            BoundaryMask%values = BoundaryMask%values .or. EdgeMask1%values
            call ovkDetectEdge(EdgeMask1, OVK_OUTER_EDGE, OVK_FALSE, .false., EdgeMask2)
            InteriorMask%values = InteriorMask%values .or. (EdgeMask2%values .and. &
              PairwiseDonors(m,n)%valid_mask%values)
          end if
        end do
        if (any(InteriorMask%values)) then
          BoundaryMask%values = BoundaryMask%values .or. Grid_n%boundary_mask%values
          call ovkFloodFill(InteriorMask, BoundaryMask)
          call ovkDetectEdge(InteriorMask, OVK_OUTER_EDGE, OVK_FALSE, .false., EdgeMask1)
          BoundaryMask%values = BoundaryMask%values .and. EdgeMask1%values
          InteriorMask%values = InteriorMask%values .or. BoundaryMask%values
          call ovkDetectEdge(InteriorMask, OVK_INNER_EDGE, OVK_FALSE, .false., BoundaryMask)
          call ovkDetectEdge(BoundaryMask, OVK_OUTER_EDGE, OVK_FALSE, .false., EdgeMask1)
          ExteriorMask = ovk_field_logical_(Grid_n%cart)
          ExteriorMask%values = EdgeMask1%values .and. .not. InteriorMask%values
          call ovkFloodFill(ExteriorMask, BoundaryMask)
          Grid_n%grid_mask%values = Grid_n%grid_mask%values .and. .not. ExteriorMask%values
          Grid_n%boundary_mask%values = Grid_n%boundary_mask%values .and. .not. ExteriorMask%values
          do m = 1, NumGrids
            if (Overlappable(m,n)) then
              PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
                .not. ExteriorMask%values
            end if
          end do
          do m = 1, NumGrids
            Grid_m => Domain%grids(IndexToID(m))
            if (Overlappable(n,m)) then
              call ovkGenerateReceiverMask(Grid_m, Grid_n, PairwiseDonors(n,m), ReceiverMask, &
                DonorSubset=ExteriorMask)
              PairwiseDonors(n,m)%valid_mask%values = PairwiseDonors(n,m)%valid_mask%values .and. &
                .not. ReceiverMask%values
            end if
          end do
          if (OVK_VERBOSE) then
            PointCount = ovkCountMask(ExteriorMask)
          end if
        end if
      end if
      if (OVK_VERBOSE) then
        if (PointCount > 0) then
          write (*, '(5a)') "* ", trim(LargeIntToString(PointCount)), &
            " points removed from grid ", trim(IntToString(Grid_n%properties%id)), "."
        end if
      end if
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished cutting boundary holes."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Cutting overlap holes..."
    end if

    allocate(FineMasks(NumGrids))
    allocate(PairwiseFineMasks(NumGrids,NumGrids))
    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      FineMasks(n) = Grid_n%grid_mask
      do m = 1, NumGrids
        PairwiseFineMasks(m,n) = Grid_n%grid_mask
      end do
    end do

    ! Cut out coarse points in overlapping regions
    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      do m = n+1, NumGrids
        Grid_m => Domain%grids(IndexToID(m))
        if (OverlapHoleCutting(m,n) .and. OverlapHoleCutting(n,m)) then
          ! Find the coarse points
          call ovkThreshold(PairwiseDonors(n,m)%res_diff, CoarseMask_m, Upper=0._rk)
          call ovkThreshold(PairwiseDonors(m,n)%res_diff, CoarseMask_n, Upper=0._rk)
          CoarseMask_m%values = CoarseMask_m%values .and. PairwiseDonors(n,m)%valid_mask%values
          CoarseMask_n%values = CoarseMask_n%values .and. PairwiseDonors(m,n)%valid_mask%values
          ! If any coarse points on grid m receive from coarse points on grid n, change them
          ! to fine
          if (any(CoarseMask_m%values)) then
            call ovkGenerateReceiverMask(Grid_m, Grid_n, PairwiseDonors(n,m), &
              ReceiverMask_m, ReceiverSubset=CoarseMask_m, DonorSubset=CoarseMask_n)
            FineMasks(m)%values = FineMasks(m)%values .or. ReceiverMask_m%values
          end if
          if (any(CoarseMask_n%values)) then
            call ovkGenerateReceiverMask(Grid_n, Grid_m, PairwiseDonors(m,n), &
              ReceiverMask_n, ReceiverSubset=CoarseMask_n, DonorSubset=CoarseMask_m)
            FineMasks(n)%values = FineMasks(n)%values .or. ReceiverMask_n%values
          end if
          FineMasks(m)%values = FineMasks(m)%values .and. .not. CoarseMask_m%values
          PairwiseFineMasks(n,m)%values = PairwiseFineMasks(n,m)%values .and. .not. &
            CoarseMask_m%values
          FineMasks(n)%values = FineMasks(n)%values .and. .not. CoarseMask_n%values
          PairwiseFineMasks(m,n)%values = PairwiseFineMasks(m,n)%values .and. .not. &
            CoarseMask_n%values
        else if (OverlapHoleCutting(m,n)) then
          ! Behave as if grid m is finer everywhere
          FineMasks(n)%values = FineMasks(n)%values .and. .not. &
            PairwiseDonors(m,n)%valid_mask%values
          PairwiseFineMasks(m,n)%values = PairwiseFineMasks(m,n)%values .and. .not. &
            PairwiseDonors(m,n)%valid_mask%values
        else if (OverlapHoleCutting(n,m)) then
          ! Behave as if grid n is finer everywhere
          FineMasks(m)%values = FineMasks(m)%values .and. .not. &
            PairwiseDonors(n,m)%valid_mask%values
          PairwiseFineMasks(n,m)%values = PairwiseFineMasks(n,m)%values .and. .not. &
            PairwiseDonors(n,m)%valid_mask%values
        end if
      end do
    end do

    allocate(PaddingMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      PaddingMasks(n) = FineMasks(n)
      CoarseMask = ovk_field_logical_(Grid_n%cart)
      CoarseMask%values = .not. FineMasks(n)%values
      call ovkDetectEdge(CoarseMask, OVK_OUTER_EDGE, OVK_MIRROR, .false., CoarseEdgeMask)
      do m = 1, NumGrids
        if (OverlapHoleCutting(m,n) .and. ConnectionType(m,n) == OVK_CONNECTION_FRINGE) then
          PairwiseCoarseMask = ovk_field_logical_(Grid_n%cart)
          PairwiseCoarseMask%values = .not. PairwiseFineMasks(m,n)%values
          call ovkDetectEdge(PairwiseCoarseMask, OVK_OUTER_EDGE, OVK_MIRROR, .false., &
            PaddingEdgeMask)
          PaddingEdgeMask%values = PaddingEdgeMask%values .and. CoarseEdgeMask%values
          call ovkDilate(PaddingEdgeMask, FringePadding(m,n), OVK_FALSE)
          PaddingMasks(n)%values = PaddingMasks(n)%values .or. &
            (PairwiseDonors(m,n)%valid_mask%values .and. PaddingEdgeMask%values)
        end if
      end do
    end do

    allocate(CutMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      CutMasks(n) = ovk_field_logical_(Grid_n%cart)
      CutMasks(n)%values = Grid_n%grid_mask%values .and. .not. PaddingMasks(n)%values
    end do

    call ovkEditAssemblerDomain(Assembler, EditDomain)

    do n = 1, NumGrids
      call ovkEditDomainGrid(EditDomain, IndexToID(n), EditGrid)
      call ovkEditGridMask(EditGrid, EditGridMask)
      call ovkEditGridBoundaryMask(EditGrid, EditBoundaryMask)
      call ovkEditGridInternalBoundaryMask(EditGrid, EditInternalBoundaryMask)
      EditGridMask%values = EditGridMask%values .and. .not. CutMasks(n)%values
      EditBoundaryMask%values = EditBoundaryMask%values .and. .not. CutMasks(n)%values
      EditInternalBoundaryMask%values = EditInternalBoundaryMask%values .and. .not. &
        CutMasks(n)%values
      call ovkReleaseGridMask(EditGrid, EditGridMask)
      call ovkReleaseGridBoundaryMask(EditGrid, EditBoundaryMask)
      call ovkReleaseGridInternalBoundaryMask(EditGrid, EditInternalBoundaryMask)
      call ovkReleaseDomainGrid(EditDomain, EditGrid)
      if (OVK_VERBOSE) then
        Grid_n => Domain%grids(IndexToID(n))
        PointCount = ovkCountMask(CutMasks(n))
        if (PointCount > 0) then
          write (*, '(5a)') "* ", trim(LargeIntToString(PointCount)), &
            " points removed from grid ", trim(IntToString(Grid_n%properties%id)), "."
        end if
      end if
    end do

    call ovkReleaseAssemblerDomain(Assembler, EditDomain)

    do m = 1, NumGrids
      do n = 1, NumGrids
        if (Overlappable(m,n)) then
          Grid_m => Domain%grids(IndexToID(m))
          Grid_n => Domain%grids(IndexToID(n))
          Bounds = ovkBBIntersect(Grid_m%bounds, Grid_n%bounds)
          if (ovkBBIsEmpty(Bounds)) then
            Overlappable(m,n) = .false.
            BoundaryHoleCutting(m,n) = .false.
            OverlapHoleCutting(m,n) = .false.
            ConnectionType(m,n) = OVK_CONNECTION_NONE
            Assembler%properties%overlappable(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%boundary_hole_cutting(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%overlap_hole_cutting(IndexToID(m),IndexToID(n)) = .false.
            Assembler%properties%connection_type(IndexToID(m),IndexToID(n)) = OVK_CONNECTION_NONE
          end if
        end if
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished cutting overlap holes."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Updating overlap information after hole cutting..."
    end if

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      do m = 1, NumGrids
        if (Overlappable(m,n)) then
          PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
            Grid_n%grid_mask%values
        end if
      end do
      do m = 1, NumGrids
        Grid_m => Domain%grids(IndexToID(m))
        if (Overlappable(m,n)) then
          NonGridMask = ovk_field_logical_(Grid_m%cart)
          NonGridMask%values = .not. Grid_m%grid_mask%values
          call ovkGenerateReceiverMask(Grid_n, Grid_m, PairwiseDonors(m,n), ReceiverMask, &
            DonorSubset=NonGridMask)
          PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
            .not. ReceiverMask%values
          do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
            do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
              do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
                if (PairwiseDonors(m,n)%valid_mask%values(i,j,k)) then
                  do l = 1, NumDims
                    DonorCell(l) = PairwiseDonors(m,n)%cells(l)%values(i,j,k)
                  end do
                  DonorCell(NumDims+1:) = 1
                  PairwiseDonors(m,n)%edge_dist%values(i,j,k) = Grid_m%cell_edge_dist%values( &
                    DonorCell(1),DonorCell(2),DonorCell(3))
                end if
              end do
            end do
          end do
        end if
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished updating overlap information."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Choosing best donor/receiver pairs..."
    end if

    allocate(FringeMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      FringeMasks(n) = ovk_field_logical_(Grid_n%cart, .false.)
      IgnoredEdgeMask = ovk_field_logical_(Grid_n%cart)
      IgnoredEdgeMask%values = Grid_n%boundary_mask%values .or. Grid_n%internal_boundary_mask%values
      call GenerateFringe(Grid_n%grid_mask, FringeSize(n), FringeMasks(n), &
        IgnoredEdges=IgnoredEdgeMask)
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        select case (ConnectionType(m,n))
        case (OVK_CONNECTION_FRINGE)
          PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
            FringeMasks(n)%values
        case (OVK_CONNECTION_NONE)
          PairwiseDonors(m,n)%valid_mask%values = .false.
        end select
      end do
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
        do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
          do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
            BestCellDistance = -huge(0._rk)
            BestCellDiff = huge(0._rk)
            BestGrid = 0
            do m = 1, NumGrids
              if (ConnectionType(m,n) /= OVK_CONNECTION_NONE .and. &
                PairwiseDonors(m,n)%valid_mask%values(i,j,k)) then
                if (DisjointConnection(n,m)) then
                  PaddingAmount = FringeSize(m) + FringePadding(n,m)
                else
                  PaddingAmount = FringePadding(n,m)
                end if
                CandidateCellDistance = real(PairwiseDonors(m,n)%edge_dist%values(i,j,k), &
                  kind=rk)/real(max(PaddingAmount,1),kind=rk)
                CandidateCellDiff = PairwiseDonors(m,n)%res_diff%values(i,j,k)
                if (CandidateCellDistance > BestCellDistance .or. (CandidateCellDistance >= &
                  1._rk-TOLERANCE .and. CandidateCellDiff < BestCellDiff)) then
                  if (BestGrid /= 0) then
                    PairwiseDonors(BestGrid,n)%valid_mask%values(i,j,k) = .false.
                  end if
                  BestCellDistance = CandidateCellDistance
                  BestCellDiff = CandidateCellDiff
                  BestGrid = m
                else
                  PairwiseDonors(m,n)%valid_mask%values(i,j,k) = .false.
                end if
              end if
            end do
          end do
        end do
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished choosing best donor/receiver pairs."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Merging remaining candidate donor/receiver pairs into final set..."
    end if

    allocate(Donors(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      call ovkMergeDonors(PairwiseDonors(:,n), Donors(n))
      if (OVK_VERBOSE) then
        PointCount = ovkCountMask(Donors(n)%valid_mask)
        if (PointCount > 0) then
          write (*, '(5a)') "* ", trim(LargeIntToString(PointCount)), &
            " receiver points on grid ", trim(IntToString(Grid_n%properties%id)), "."
        end if
      end if
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        call ovkDestroyDonors(PairwiseDonors(m,n))
      end do
    end do
    deallocate(PairwiseDonors)

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished merging candidate donor/receiver pairs."
    end if

    ! TODO: Need to update this to work with per-pair interp scheme parameter
    if (any(InterpScheme == OVK_INTERP_CUBIC)) then

      if (OVK_VERBOSE) then
        write (*, '(a)') "Expanding donor cells for cubic interpolation..."
      end if

      allocate(ValidCellMasks(NumGrids))
      allocate(ExpandableCellMasks(NumGrids))
      allocate(CellQualities(NumGrids))

      do n = 1, NumGrids
        Grid_n => Domain%grids(IndexToID(n))
        call GenerateValidCellMask(Grid_n, ValidCellMasks(n))
        call ovkDetectEdge(ValidCellMasks(n), OVK_INNER_EDGE, OVK_FALSE, .false., &
          ValidCellInnerEdgeMask)
        ExpandableCellMasks(n) = ovk_field_logical_(Grid_n%cell_cart)
        ExpandableCellMasks(n)%values = ValidCellMasks(n)%values .and. .not. &
          ValidCellInnerEdgeMask%values
        call CountNeighbors(ExpandableCellMasks(n), CellQualities(n))
      end do

      do n = 1, NumGrids
        Grid_n => Domain%grids(IndexToID(n))
        if (InterpScheme(n) /= OVK_INTERP_CUBIC) cycle
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              ReceiverPoint = [i,j,k]
              if (Donors(n)%valid_mask%values(i,j,k)) then
                m = IDToIndex(Donors(n)%grid_ids%values(i,j,k))
                Grid_m => Domain%grids(IndexToID(m))
                do l = 1, Grid_m%cart%nd
                  DonorCell(l) = Donors(n)%cells(l)%values(i,j,k)
                  DonorCellCoords(l) = Donors(n)%cell_coords(l)%values(i,j,k)
                end do
                DonorCell(Grid_m%cell_cart%nd+1:) = 1
                if (ExpandableCellMasks(m)%values(DonorCell(1),DonorCell(2),DonorCell(3))) then
                  ExpandableCell = .true.
                else
                  BestCellQuality = 0
                  DonorCellShift = 0
                  NeighborCellLower(:Grid_m%cart%nd) = DonorCell(:Grid_m%cart%nd)-1
                  NeighborCellLower(Grid_m%cart%nd+1:) = DonorCell(Grid_m%cart%nd+1:)
                  NeighborCellUpper(:Grid_m%cart%nd) = DonorCell(:Grid_m%cart%nd)+1
                  NeighborCellUpper(Grid_m%cart%nd+1:) = DonorCell(Grid_m%cart%nd+1:)
                  AwayFromEdge = ovkCartContains(Grid_m%cell_cart, NeighborCellLower) .and. &
                    ovkCartContains(Grid_m%cell_cart, NeighborCellUpper)
                  if (AwayFromEdge) then
                    do r = NeighborCellLower(3), NeighborCellUpper(3)
                      do q = NeighborCellLower(2), NeighborCellUpper(2)
                        do p = NeighborCellLower(1), NeighborCellUpper(1)
                          NeighborCell = [p,q,r]
                          if (ExpandableCellMasks(m)%values(NeighborCell(1),NeighborCell(2), &
                            NeighborCell(3))) then
                            CellQuality = CellQualities(m)%values(NeighborCell(1),NeighborCell(2),NeighborCell(3))
                            if (CellQuality > BestCellQuality) then
                              DonorCellShift = [p-DonorCell(1),q-DonorCell(2),r-DonorCell(3)]
                              BestCellQuality = CellQuality
                            end if
                          end if
                        end do
                      end do
                    end do
                  else
                    do r = NeighborCellLower(3), NeighborCellUpper(3)
                      do q = NeighborCellLower(2), NeighborCellUpper(2)
                        do p = NeighborCellLower(1), NeighborCellUpper(1)
                          NeighborCell = [p,q,r]
                          NeighborCell(:Grid_m%cell_cart%nd) = ovkCartPeriodicAdjust( &
                            Grid_m%cell_cart, NeighborCell)
                          if (ovkCartContains(Grid_m%cell_cart, NeighborCell)) then
                            if (ExpandableCellMasks(m)%values(NeighborCell(1),NeighborCell(2), &
                              NeighborCell(3))) then
                              CellQuality = CellQualities(m)%values(NeighborCell(1),NeighborCell(2),NeighborCell(3))
                              if (CellQuality > BestCellQuality) then
                                DonorCellShift = [p-DonorCell(1),q-DonorCell(2),r-DonorCell(3)]
                                BestCellQuality = CellQuality
                              end if
                            end if
                          end if
                        end do
                      end do
                    end do
                  end if
                  if (BestCellQuality > 0) then
                    ExpandableCell = .true.
                    DonorCell = DonorCell + DonorCellShift
                    DonorCellCoords = DonorCellCoords - real(DonorCellShift(:Grid_m%cell_cart%nd),kind=rk)
                  else
                    ExpandableCell = .false.
                    if (OVK_VERBOSE) then
                      write (ERROR_UNIT, '(8a)') "WARNING: Could not use cubic interpolation for donor cell ", &
                        trim(TupleToString(DonorCell(:Grid_m%cart%nd))), " of grid ", trim(IntToString(Grid_m%properties%id)), &
                        " corresponding to receiver point ", trim(TupleToString(ReceiverPoint(:Grid_n%cart%nd))), &
                        " of grid ", trim(IntToString(Grid_n%properties%id))
                    end if
                  end if
                end if
                if (ExpandableCell) then
                  do l = 1, Grid_n%cart%nd
                    ReceiverCoords(l) = Grid_n%xyz(l)%values(ReceiverPoint(1),ReceiverPoint(2), &
                      ReceiverPoint(3))
                  end do
                  call ExpandDonorCell(Grid_m, DonorCell, DonorCellCoords, &
                    ReceiverCoords, ExpandedDonorCell, ExpandedDonorCellCoords)
                  ExpandedDonorCell(Grid_m%cart%nd+1:) = 1
                  do l = 1, Grid_m%cart%nd
                    Donors(n)%cells(l)%values(i,j,k) = ExpandedDonorCell(l)
                    Donors(n)%cell_coords(l)%values(i,j,k) = ExpandedDonorCellCoords(l)
                  end do
                  Donors(n)%cell_extents%values(i,j,k) = 4
                end if
              end if
            end do
          end do
        end do
      end do

      if (OVK_VERBOSE) then
        write (*, '(a)') "Finished expanding donor cells."
      end if

    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Searching for orphan points..."
    end if

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      do m = 1, NumGrids
        Grid_m => Domain%grids(IndexToID(m))
        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE .and. DisjointConnection(m,n)) then
          call ovkGenerateReceiverMask(Grid_n, Grid_m, Donors(n), ReceiverMask, &
            DonorSubset=Donors(m)%valid_mask)
          Donors(n)%valid_mask%values = Donors(n)%valid_mask%values .and. .not. &
            ReceiverMask%values
        end if
      end do
    end do

    allocate(OrphanMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      OrphanMasks(n) = ovk_field_logical_(Grid_n%cart)
      OrphanMasks(n)%values = FringeMasks(n)%values .and. .not. Donors(n)%valid_mask%values
      if (OVK_VERBOSE) then
        NumWarnings = 0
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              ReceiverPoint = [i,j,k]
              if (OrphanMasks(n)%values(i,j,k)) then
                if (NumWarnings <= 100) then
                  write (ERROR_UNIT, '(4a)') "WARNING: Orphan detected at point ", &
                    trim(TupleToString(ReceiverPoint(:Grid_n%cart%nd))), " of grid ", &
                    trim(IntToString(Grid_n%properties%id))
                  if (NumWarnings == 100) then
                    write (ERROR_UNIT, '(a)') "Further warnings suppressed."
                  end if
                  NumWarnings = NumWarnings + 1
                end if
              end if
            end do
          end do
        end do
        PointCount = ovkCountMask(OrphanMasks(n))
        if (PointCount > 0) then
          write (*, '(5a)') "* ", trim(LargeIntToString(PointCount)), &
            " orphan points found on grid ", trim(IntToString(Grid_n%properties%id)), "."
        end if
      end if
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished searching for orphan points."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Generating interpolation data..."
    end if

    call ovkEditAssemblerConnectivity(Assembler, Connectivity)

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      StencilSize = 2
      do m = 1, NumGrids
        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE .and. InterpScheme(n) == OVK_INTERP_CUBIC) &
          then
          StencilSize = 4
          exit
        end if
      end do
      call ovkGetGridDescription(Grid_n, GridDescription)
      call ovkCreateConnectivityInterpData(Connectivity, IndexToID(n), GridDescription, StencilSize)
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grids(IndexToID(n))
      InterpData_n => Connectivity%interp_data(IndexToID(n))
      if (any(ConnectionType(:,n) /= OVK_CONNECTION_NONE)) then
        call ovkFillInterpData(InterpData_n, Donors(n), OrphanMasks(n), InterpScheme=InterpScheme(n))
        if (OVK_VERBOSE) then
          write (*, '(3a)') "* Generated interpolation data for grid ", &
            trim(IntToString(Grid_n%properties%id)), "."
        end if
      end if
    end do

    call ovkReleaseAssemblerConnectivity(Assembler, Connectivity)

    do m = 1, NumGrids
      call ovkDestroyDonors(Donors(m))
    end do
    deallocate(Donors)

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished generating interpolation data."
    end if

    if (OVK_VERBOSE) then
      call system_clock(ClockFinal, ClockRate)
      write (*, '(a,f0.3,a)') "Overset grid assembly finished (time: ", &
        real(ClockFinal-ClockInitial,kind=rk)/real(ClockRate,kind=rk), " seconds)."
    end if

  end subroutine ovkAssemble

  subroutine GenerateValidCellMask(Grid, ValidCellMask)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_field_logical), intent(out) :: ValidCellMask

    integer :: i, j, k
    integer, dimension(MAX_ND) :: Cell
    logical, dimension(2**Grid%cart%nd) :: VertexGridMaskValues

    ValidCellMask = ovk_field_logical_(Grid%cell_cart)

    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          Cell = [i,j,k]
          call ovkGetCellVertexData(Grid, Cell, VertexGridMaskValues=VertexGridMaskValues)
          ValidCellMask%values(i,j,k) = all(VertexGridMaskValues)
        end do
      end do
    end do

  end subroutine GenerateValidCellMask

  subroutine CountNeighbors(Mask, NeighborCounts)

    type(ovk_field_logical), intent(in) :: Mask
    type(ovk_field_int), intent(out) :: NeighborCounts

    integer :: i, j, k, m, n, o
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: NeighborLower, NeighborUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Neighbor

    NeighborCounts = ovk_field_int_(Mask%cart, 0)

    do k = Mask%cart%is(3), Mask%cart%ie(3)
      do j = Mask%cart%is(2), Mask%cart%ie(2)
        do i = Mask%cart%is(1), Mask%cart%ie(1)
          Point = [i,j,k]
          NeighborLower(:Mask%cart%nd) = Point(:Mask%cart%nd)-1
          NeighborLower(Mask%cart%nd+1:) = Point(Mask%cart%nd+1:)
          NeighborUpper(:Mask%cart%nd) = Point(:Mask%cart%nd)+1
          NeighborUpper(Mask%cart%nd+1:) = Point(Mask%cart%nd+1:)
          AwayFromEdge = ovkCartContains(Mask%cart, NeighborLower) .and. &
            ovkCartContains(Mask%cart, NeighborUpper)
          if (AwayFromEdge) then
            do o = NeighborLower(3), NeighborUpper(3)
              do n = NeighborLower(2), NeighborUpper(2)
                do m = NeighborLower(1), NeighborUpper(1)
                  Neighbor = [m,n,o]
                  if (Mask%values(Neighbor(1),Neighbor(2),Neighbor(3))) then
                    NeighborCounts%values(i,j,k) = NeighborCounts%values(i,j,k) + 1
                  end if
                end do
              end do
            end do
          else
            do o = NeighborLower(3), NeighborUpper(3)
              do n = NeighborLower(2), NeighborUpper(2)
                do m = NeighborLower(1), NeighborUpper(1)
                  Neighbor = [m,n,o]
                  Neighbor(:Mask%cart%nd) = ovkCartPeriodicAdjust(Mask%cart, Neighbor)
                  if (ovkCartContains(Mask%cart, Neighbor)) then
                    if (Mask%values(Neighbor(1),Neighbor(2),Neighbor(3))) then
                      NeighborCounts%values(i,j,k) = NeighborCounts%values(i,j,k) + 1
                    end if
                  end if
                end do
              end do
            end do
          end if
        end do
      end do
    end do

  end subroutine CountNeighbors

  subroutine ExpandDonorCell(DonorGrid, DonorCell, DonorCellCoords, ReceiverCoords, &
    ExpandedDonorCell, ExpandedDonorCellCoords)

    type(ovk_grid), intent(in) :: DonorGrid
    integer, dimension(DonorGrid%cart%nd), intent(in) :: DonorCell
    real(rk), dimension(DonorGrid%cart%nd), intent(in) :: DonorCellCoords
    real(rk), dimension(DonorGrid%cart%nd), intent(in) :: ReceiverCoords
    integer, dimension(DonorGrid%cart%nd), intent(out) :: ExpandedDonorCell
    real(rk), dimension(DonorGrid%cart%nd), intent(out) :: ExpandedDonorCellCoords

    integer, dimension(MAX_ND) :: ExpandedDonorLower, ExpandedDonorUpper
    logical :: AwayFromEdge
    integer :: i, j
    integer, dimension(MAX_ND) :: Vertex
    integer, dimension(MAX_ND) :: AdjustedVertex
    real(rk), dimension(DonorGrid%cart%nd) :: PrincipalCoords
    real(rk), dimension(DonorGrid%cart%nd,4**DonorGrid%cart%nd) :: VertexCoords
    logical :: Success

    ExpandedDonorLower(:DonorGrid%cart%nd) = DonorCell-1
    ExpandedDonorLower(DonorGrid%cart%nd+1:) = 1

    ExpandedDonorUpper(:DonorGrid%cart%nd) = DonorCell+2
    ExpandedDonorUpper(DonorGrid%cart%nd+1:) = 1

    AwayFromEdge = ovkCartContains(DonorGrid%cart, ExpandedDonorLower) .and. &
      ovkCartContains(DonorGrid%cart, ExpandedDonorUpper)

    ExpandedDonorCell = ovkCartPeriodicAdjust(DonorGrid%cell_cart, ExpandedDonorLower)

    if (AwayFromEdge) then
      do i = 1, 4**DonorGrid%cart%nd
        Vertex = ExpandedDonorLower + [(modulo((i-1)/4**j,4),j=0,MAX_ND-1)]
        do j = 1, DonorGrid%cart%nd
          VertexCoords(j,i) = DonorGrid%xyz(j)%values(Vertex(1),Vertex(2),Vertex(3))
        end do
      end do
    else
      do i = 1, 4**DonorGrid%cart%nd
        Vertex = ExpandedDonorLower + [(modulo((i-1)/4**j,4),j=0,MAX_ND-1)]
        AdjustedVertex(:DonorGrid%cart%nd) = ovkCartPeriodicAdjust(DonorGrid%cart, Vertex)
        AdjustedVertex(DonorGrid%cart%nd+1:) = 1
        do j = 1, DonorGrid%cart%nd
          PrincipalCoords(j) = DonorGrid%xyz(j)%values(AdjustedVertex(1),AdjustedVertex(2), &
            AdjustedVertex(3))
        end do
        VertexCoords(:,i) = ovkPeriodicExtend(DonorGrid%cart, DonorGrid%properties%periodic_length, &
          Vertex, PrincipalCoords)
      end do
    end if

    Success = .true.

    select case (DonorGrid%properties%geometry_type)
    case (OVK_GRID_GEOMETRY_CARTESIAN)
      select case (DonorGrid%cart%nd)
      case (2)
        ExpandedDonorCellCoords = ovkRectangleIsoInverseCubic(VertexCoords, ReceiverCoords)
      case (3)
        ExpandedDonorCellCoords = ovkCuboidIsoInverseCubic(VertexCoords, ReceiverCoords)
      end select
    case (OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN)
      select case (DonorGrid%cart%nd)
      case (2)
        ExpandedDonorCellCoords = ovkOrientedRectangleIsoInverseCubic(VertexCoords, ReceiverCoords)
      case (3)
        ExpandedDonorCellCoords = ovkOrientedCuboidIsoInverseCubic(VertexCoords, ReceiverCoords)
      end select
    case default
      select case (DonorGrid%cart%nd)
      case (2)
        ExpandedDonorCellCoords = ovkQuadIsoInverseCubic(VertexCoords, ReceiverCoords, &
          Guess=DonorCellCoords, Success=Success)
      case (3)
        ExpandedDonorCellCoords = ovkHexahedronIsoInverseCubic(VertexCoords, ReceiverCoords, &
          Guess=DonorCellCoords, Success=Success)
      end select
    end select

    if (OVK_VERBOSE) then
      if (.not. Success) then
        write (ERROR_UNIT, '(6a)') "WARNING: Failed to compute isoparametric coordinates of ", &
          trim(CoordsToString(ReceiverCoords)), " in expanded donor cell surrounding cell ", &
          trim(TupleToString(DonorCell)), " of grid ", trim(IntToString(DonorGrid%properties%id))
      end if
    end if

  end subroutine ExpandDonorCell

  subroutine GenerateFringe(GridMask, FringeSize, FringeMask, IgnoredEdges)

    type(ovk_field_logical), intent(in) :: GridMask
    integer, intent(in) :: FringeSize
    type(ovk_field_logical), intent(out) :: FringeMask
    type(ovk_field_logical), intent(in), optional :: IgnoredEdges

    type(ovk_field_logical) :: EdgeMask
    type(ovk_field_logical) :: NonSubsetMask
    type(ovk_field_logical) :: NonSubsetEdgeMask
    type(ovk_field_logical) :: SubsetMask
    type(ovk_field_logical) :: SubsetEdgeMask

    if (OVK_DEBUG) then
      if (GridMask%cart%periodic_storage == OVK_OVERLAP_PERIODIC) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_OVERLAP_PERIODIC is not currently supported."
        stop 1
      end if
    end if

    call ovkDetectEdge(GridMask, OVK_OUTER_EDGE, OVK_FALSE, .true., EdgeMask)

    if (present(IgnoredEdges)) then

      NonSubsetMask = IgnoredEdges
      call ovkDetectEdge(NonSubsetMask, OVK_OUTER_EDGE, OVK_FALSE, .true., NonSubsetEdgeMask)

      SubsetMask = ovk_field_logical_(GridMask%cart)
      SubsetMask%values = GridMask%values .and. .not. NonSubsetMask%values

      call ovkDetectEdge(SubsetMask, OVK_OUTER_EDGE, OVK_FALSE, .true., SubsetEdgeMask)

      EdgeMask%values = EdgeMask%values .and. .not. (NonSubsetEdgeMask%values .and. .not. &
        SubsetEdgeMask%values)

    end if

    call ovkDilate(EdgeMask, FringeSize, OVK_FALSE)

    FringeMask = ovk_field_logical_(GridMask%cart)
    FringeMask%values = GridMask%values .and. EdgeMask%values(&
      GridMask%cart%is(1):GridMask%cart%ie(1),GridMask%cart%is(2):GridMask%cart%ie(2), &
      GridMask%cart%is(3):GridMask%cart%ie(3))

  end subroutine GenerateFringe

end module ovkOverset
