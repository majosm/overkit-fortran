! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkAssembly

  use ovkArray
  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDomain
  use ovkOverlap
  use ovkOverlapAccel
  use ovkField
  use ovkFieldOps
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkPLOT3D
  implicit none

  private

  ! API
  public :: ovkAssemble

  type t_reduced_domain_info
    integer :: ngrids
    integer, dimension(:), pointer :: index_to_id
    logical, dimension(:,:), pointer :: overlappable
    real(rk), dimension(:,:), pointer :: overlap_tolerance
    real(rk), dimension(:), pointer :: overlap_accel_quality_adjust
    logical, dimension(:), pointer :: infer_boundaries
    logical, dimension(:,:), pointer :: boundary_hole_cutting
    integer, dimension(:,:), pointer :: occludes
    integer, dimension(:,:), pointer :: connection_type
    integer, dimension(:,:), pointer :: interp_scheme
    integer, dimension(:), pointer :: fringe_size
    integer, dimension(:,:), pointer :: edge_padding
    logical, dimension(:,:), pointer :: overlap_minimization
  end type t_reduced_domain_info

contains

  subroutine ovkAssemble(Domain)

    type(ovk_domain), intent(inout) :: Domain

    integer :: ClockInitial, ClockFinal, ClockRate
    type(t_reduced_domain_info) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), allocatable :: OverlapResolutions
    type(ovk_field_int), dimension(:), allocatable :: DonorGridIDs

    if (Domain%logger%verbose) then
      call system_clock(ClockInitial, ClockRate)
      write (*, '(a)') "Overset grid assembly started..."
    end if

    if (OVK_DEBUG) then
      if ( &
        Domain%properties_edit_ref_count > 0 .or. &
        any(Domain%grid_edit_ref_counts > 0) .or. &
        any(Domain%connectivity_edit_ref_counts > 0) &
      ) then
        write (ERROR_UNIT, '(a)') "ERROR: Cannot perform assembly; domain is still being edited."
        stop 1
      end if
    end if

    call InitAssembly(Domain, ReducedDomainInfo)

    call CollideGrids(Domain, ReducedDomainInfo)

    call InferNonOverlappingBoundaries(Domain, ReducedDomainInfo)

    call CutBoundaryHoles(Domain, ReducedDomainInfo)

    allocate(OverlapResolutions(ReducedDomainInfo%ngrids,ReducedDomainInfo%ngrids))

    call ComputeOverlapResolutions(Domain, ReducedDomainInfo, OverlapResolutions)

    call DetectOccludedPoints(Domain, ReducedDomainInfo, OverlapResolutions)

    call LocateCandidateReceivers(Domain, ReducedDomainInfo)

    allocate(DonorGridIDs(ReducedDomainInfo%ngrids))

    call ChooseDonors(Domain, ReducedDomainInfo, OverlapResolutions, DonorGridIDs)

    deallocate(OverlapResolutions)

    call MinimizeOverlap(Domain, ReducedDomainInfo, DonorGridIDs)

    call FillConnectivity(Domain, ReducedDomainInfo, DonorGridIDs)

    call FinalizeAssembly(Domain, ReducedDomainInfo)

    if (Domain%logger%verbose) then
      call system_clock(ClockFinal, ClockRate)
      write (*, '(a,f0.3,a)') "Overset grid assembly finished (time: ", &
        real(ClockFinal-ClockInitial,kind=rk)/real(ClockRate,kind=rk), " seconds)."
    end if

  end subroutine ovkAssemble

  subroutine InitAssembly(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(out) :: ReducedDomainInfo

    integer :: m, n, p, q
    integer :: NumGrids
    logical :: GridExists
    logical :: GridIsOverlappable
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: Overlappable
    real(rk), dimension(:,:), pointer :: OverlapTolerance
    real(rk), dimension(:), pointer :: OverlapAccelQualityAdjust
    logical, dimension(:), pointer :: InferBoundaries
    logical, dimension(:,:), pointer :: BoundaryHoleCutting
    integer, dimension(:,:), pointer :: Occludes
    integer, dimension(:,:), pointer :: ConnectionType
    integer, dimension(:,:), pointer :: InterpScheme
    integer, dimension(:), pointer :: FringeSize
    integer, dimension(:,:), pointer :: EdgePadding
    logical, dimension(:,:), pointer :: OverlapMinimization

    call PrintDomainSummary(Domain)

    ! Ignore empty and non-overlapping grids
    NumGrids = 0
    do n = 1, Domain%properties%ngrids
      GridExists = .not. ovkCartIsEmpty(Domain%grid(n)%cart)
      GridIsOverlappable = any(Domain%properties%overlappable(:,n)) .or. &
        any(Domain%properties%overlappable(n,:))
      if (GridExists .and. GridIsOverlappable) then
        NumGrids = NumGrids + 1
      end if
    end do

    ReducedDomainInfo%ngrids = NumGrids

    allocate(ReducedDomainInfo%index_to_id(NumGrids))
    IndexToID => ReducedDomainInfo%index_to_id

    m = 1
    do n = 1, Domain%properties%ngrids
      GridExists = .not. ovkCartIsEmpty(Domain%grid(n)%cart)
      GridIsOverlappable = any(Domain%properties%overlappable(:,n)) .or. &
        any(Domain%properties%overlappable(n,:))
      if (GridExists .and. GridIsOverlappable) then
        IndexToID(m) = n
        m = m + 1
      end if
    end do

    allocate(ReducedDomainInfo%overlappable(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%overlap_tolerance(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%overlap_accel_quality_adjust(NumGrids))
    Overlappable => ReducedDomainInfo%overlappable
    OverlapTolerance => ReducedDomainInfo%overlap_tolerance
    OverlapAccelQualityAdjust => ReducedDomainInfo%overlap_accel_quality_adjust

    do n = 1, NumGrids
      q = IndexToID(n)
      do m = 1, NumGrids
        p = IndexToID(m)
        Overlappable(m,n) = Domain%properties%overlappable(p,q)
        OverlapTolerance(m,n) = Domain%properties%overlap_tolerance(p,q)
      end do
      OverlapAccelQualityAdjust(n) = Domain%properties%overlap_accel_quality_adjust(q)
    end do

    allocate(ReducedDomainInfo%infer_boundaries(NumGrids))
    allocate(ReducedDomainInfo%boundary_hole_cutting(NumGrids,NumGrids))
    InferBoundaries => ReducedDomainInfo%infer_boundaries
    BoundaryHoleCutting => ReducedDomainInfo%boundary_hole_cutting

    do n = 1, NumGrids
      q = IndexToID(n)
      InferBoundaries(n) = Domain%properties%infer_boundaries(q)
      do m = 1, NumGrids
        p = IndexToID(m)
        BoundaryHoleCutting(m,n) = Domain%properties%boundary_hole_cutting(p,q)
      end do
    end do

    allocate(ReducedDomainInfo%occludes(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%connection_type(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%interp_scheme(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%fringe_size(NumGrids))
    allocate(ReducedDomainInfo%edge_padding(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%overlap_minimization(NumGrids,NumGrids))
    Occludes => ReducedDomainInfo%occludes
    ConnectionType => ReducedDomainInfo%connection_type
    InterpScheme => ReducedDomainInfo%interp_scheme
    FringeSize => ReducedDomainInfo%fringe_size
    EdgePadding => ReducedDomainInfo%edge_padding
    OverlapMinimization => ReducedDomainInfo%overlap_minimization

    do n = 1, NumGrids
      q = IndexToID(n)
      FringeSize(n) = Domain%properties%fringe_size(q)
      do m = 1, NumGrids
        p = IndexToID(m)
        Occludes(m,n) = Domain%properties%occludes(p,q)
        ConnectionType(m,n) = Domain%properties%connection_type(p,q)
        InterpScheme(m,n) = Domain%properties%interp_scheme(p,q)
        EdgePadding(m,n) = Domain%properties%edge_padding(p,q)
        OverlapMinimization(m,n) = Domain%properties%overlap_minimization(p,q)
      end do
    end do

  end subroutine InitAssembly

  subroutine CollideGrids(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(in) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: Overlappable
    real(rk), dimension(:,:), pointer :: OverlapTolerance
    real(rk), dimension(:), pointer :: OverlapAccelQualityAdjust
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_bbox), dimension(:,:), allocatable :: Bounds
    type(ovk_bbox) :: AccelBounds
    real(rk) :: AccelMaxOverlapTolerance
    type(t_overlap_accel) :: OverlapAccel
    type(ovk_overlap), pointer :: Overlap

    if (Domain%logger%verbose) then
      write (*, '(a)') "Detecting overlap between grids..."
    end if

    NumDims = Domain%properties%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    Overlappable => ReducedDomainInfo%overlappable
    OverlapTolerance => ReducedDomainInfo%overlap_tolerance
    OverlapAccelQualityAdjust => ReducedDomainInfo%overlap_accel_quality_adjust

    allocate(Bounds(NumGrids,NumGrids))

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (Overlappable(m,n)) then
          Bounds(m,n) = ovkBBIntersect(Domain%grid(m)%bounds, Domain%grid(n)%bounds)
        else
          Bounds(m,n) = ovk_bbox_(NumDims)
        end if
      end do
    end do

    do m = 1, NumGrids
      Grid_m => Domain%grid(IndexToID(m))
      if (any(Overlappable(m,:))) then
        if (Domain%logger%verbose) then
          write (*, '(3a)') "* Generating overlap search accelerator on grid ", &
            trim(IntToString(Grid_m%properties%id)), "..."
        end if
        AccelBounds = ovk_bbox_(NumDims)
        AccelMaxOverlapTolerance = 0._rk
        do n = 1, NumGrids
          if (Overlappable(m,n)) then
            AccelBounds = ovkBBUnion(AccelBounds, Bounds(m,n))
            AccelMaxOverlapTolerance = max(AccelMaxOverlapTolerance, OverlapTolerance(m,n))
          end if
        end do
        call GenerateOverlapAccel(Grid_m, OverlapAccel, AccelBounds, AccelMaxOverlapTolerance, &
          OverlapAccelQualityAdjust(m))
        if (Domain%logger%verbose) then
          write (*, '(3a)') "* Finished generating overlap search accelerator on grid ", &
            trim(IntToString(Grid_m%properties%id)), "."
        end if
        do n = 1, NumGrids
          Grid_n => Domain%grid(IndexToID(n))
          if (Overlappable(m,n)) then
            Overlap => Domain%overlap(Grid_m%properties%id, Grid_n%properties%id)
            call DetectOverlap(Grid_m, Grid_n, OverlapAccel, Bounds(m,n), OverlapTolerance(m,n), &
              Overlap)
            if (Domain%logger%verbose) then
              if (Overlap%properties%noverlap > 0_lk) then
                write (*, '(7a)') "* Detected ", trim(LargeIntToString( &
                  Overlap%properties%noverlap)), " points on grid ", &
                  trim(IntToString(Grid_n%properties%id)), " overlapped by grid ", &
                  trim(IntToString(Grid_m%properties%id)), "."
              end if
            end if
          end if
        end do
        call DestroyOverlapAccel(OverlapAccel)
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished detecting overlap between grids."
    end if

  end subroutine CollideGrids

  subroutine InferNonOverlappingBoundaries(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, m, n
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:), pointer :: InferBoundaries
    type(ovk_field_logical) :: InferredBoundaryMask
    type(ovk_field_logical) :: InitialBoundaryMask
    type(ovk_grid), pointer :: Grid
    type(ovk_overlap), pointer :: Overlap
    type(ovk_field_int), pointer :: State
    integer(lk) :: NumBoundaryPoints

    if (Domain%logger%verbose) then
      write (*, '(a)') "Inferring domain boundaries in non-overlapping regions..."
    end if

    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    InferBoundaries => ReducedDomainInfo%infer_boundaries

    do n = 1, NumGrids
      Grid => Domain%grid(IndexToID(n))
      if (InferBoundaries(n)) then
        call ovkDetectEdge(Grid%mask, OVK_INNER_EDGE, OVK_FALSE, .false., InferredBoundaryMask)
        call ovkFilterGridState(Grid, OVK_STATE_DOMAIN_BOUNDARY, OVK_ALL, InitialBoundaryMask)
        InferredBoundaryMask%values = InferredBoundaryMask%values .and. .not. &
          InitialBoundaryMask%values
        do m = 1, NumGrids
          Overlap => Domain%overlap(IndexToID(m),IndexToID(n))
          if (OverlapExists(Overlap)) then
            InferredBoundaryMask%values = InferredBoundaryMask%values .and. .not. &
              Overlap%mask%values
          end if
        end do
        call ovkEditGridState(Grid, State)
        do k = Grid%cart%is(3), Grid%cart%ie(3)
          do j = Grid%cart%is(2), Grid%cart%ie(2)
            do i = Grid%cart%is(1), Grid%cart%ie(1)
              if (InferredBoundaryMask%values(i,j,k)) then
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_DOMAIN_BOUNDARY, &
                  OVK_STATE_INFERRED_DOMAIN_BOUNDARY))
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid, State)
        if (Domain%logger%verbose) then
          NumBoundaryPoints = ovkCountMask(InferredBoundaryMask)
          if (NumBoundaryPoints > 0_lk) then
            write (*, '(5a)') "* ", trim(LargeIntToString(NumBoundaryPoints)), &
              " points marked as boundaries on grid ", trim(IntToString(Grid%properties%id)), "."
          end if
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished inferring domain boundaries in non-overlapping regions."
    end if

  end subroutine InferNonOverlappingBoundaries

  subroutine CutBoundaryHoles(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: BoundaryHoleCutting
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap_mn, Overlap_nm
    type(ovk_field_logical), dimension(:), allocatable :: BoundaryHoleMasks
    type(ovk_field_logical) :: EdgeMask1, EdgeMask2
    type(ovk_field_logical) :: BoundaryMask
    type(ovk_field_logical) :: InteriorMask
    type(ovk_field_logical) :: SpuriousBoundaryMask
    type(ovk_field_logical) :: OverlappingMask
    integer(lk) :: NumRemoved
    logical, dimension(:), allocatable :: UpdateGrid
    type(ovk_field_int), pointer :: State

    if (Domain%logger%verbose) then
      write (*, '(a)') "Cutting boundary holes..."
    end if

    NumDims = Domain%properties%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    BoundaryHoleCutting => ReducedDomainInfo%boundary_hole_cutting

    allocate(BoundaryHoleMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(BoundaryHoleCutting(:,n))) then
        BoundaryMask = ovk_field_logical_(Grid_n%cart, .false.)
        InteriorMask = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          Grid_m => Domain%grid(IndexToID(m))
          Overlap_mn => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
          Overlap_nm => Domain%overlap(Grid_n%properties%id,Grid_m%properties%id)
          if (BoundaryHoleCutting(m,n)) then
            call ovkDetectEdge(Overlap_mn%mask, OVK_OUTER_EDGE, OVK_MIRROR, .false., EdgeMask1)
            call ovkDetectEdge(Overlap_nm%mask, OVK_INNER_EDGE, OVK_FALSE, .false., EdgeMask2)
            call ovkFindOverlappingPoints(Grid_n, Grid_m, Overlap_nm, EdgeMask2, OverlappingMask)
            EdgeMask1%values = EdgeMask1%values .and. .not. OverlappingMask%values
            i = 0
            ! Explicit conversion to logical in order to work around GCC 4.7 bug
            do while (logical(any(EdgeMask1%values)))
              call ovkDilate(OverlappingMask, 1, OVK_FALSE)
              EdgeMask1%values = EdgeMask1%values .and. .not. OverlappingMask%values
              i = i + 1
            end do
            call ovkFindOverlappingPoints(Grid_n, Grid_m, Overlap_nm, Grid_m%boundary_mask, &
              OverlappingMask)
            do j = 1, i
              call ovkDilate(OverlappingMask, 1, OVK_FALSE)
            end do
            OverlappingMask%values = OverlappingMask%values .and. .not. Overlap_mn%mask%values
            call ovkDetectEdge(OverlappingMask, OVK_OUTER_EDGE, OVK_FALSE, .false., EdgeMask1)
            EdgeMask1%values = EdgeMask1%values .and. Overlap_mn%mask%values
            BoundaryMask%values = BoundaryMask%values .or. EdgeMask1%values
            call ovkDetectEdge(EdgeMask1, OVK_OUTER_EDGE, OVK_FALSE, .false., EdgeMask2)
            InteriorMask%values = InteriorMask%values .or. (EdgeMask2%values .and. &
              Overlap_mn%mask%values)
          end if
        end do
        BoundaryHoleMasks(n) = ovk_field_logical_(Grid_n%cart, .false.)
        if (any(InteriorMask%values)) then
          BoundaryMask%values = BoundaryMask%values .or. Grid_n%boundary_mask%values
          call ovkFlood(InteriorMask, BoundaryMask)
          call ovkDetectEdge(Grid_n%mask, OVK_INNER_EDGE, OVK_FALSE, .false., EdgeMask1)
          call ovkDetectEdge(InteriorMask, OVK_OUTER_EDGE, OVK_FALSE, .false., EdgeMask2)
          SpuriousBoundaryMask = ovk_field_logical_(Grid_n%cart)
          SpuriousBoundaryMask%values = Grid_n%boundary_mask%values .and. EdgeMask1%values &
            .and. .not. EdgeMask2%values
          BoundaryMask%values = BoundaryMask%values .and. .not. SpuriousBoundaryMask%values
          BoundaryHoleMasks(n)%values = Grid_n%mask%values .and. .not. (InteriorMask%values .or. &
            BoundaryMask%values)
          if (Domain%logger%verbose) then
            NumRemoved = ovkCountMask(BoundaryHoleMasks(n))
            if (NumRemoved > 0_lk) then
              write (*, '(5a)') "* ", trim(LargeIntToString(NumRemoved)), &
                " points removed from grid ", trim(IntToString(Grid_n%properties%id)), "."
            end if
          end if
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished cutting boundary holes."
    end if

    if (Domain%logger%verbose) then
      write (*, '(a)') "Updating grids and overlap information after boundary hole cutting..."
    end if

    allocate(UpdateGrid(NumGrids))

    do n = 1, NumGrids
      UpdateGrid(n) = .false.
      if (any(BoundaryHoleCutting(:,n))) then
        if (any(BoundaryHoleMasks(n)%values)) then
          UpdateGrid(n) = .true.
        end if
      end if
    end do

    do n = 1, NumGrids
      if (UpdateGrid(n)) then
        Grid_n => Domain%grid(IndexToID(n))
        call ovkEditGridState(Grid_n, State)
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              if (BoundaryHoleMasks(n)%values(i,j,k)) then
                State%values(i,j,k) = iand(State%values(i,j,k), not(OVK_STATE_GRID))
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_HOLE, &
                  OVK_STATE_BOUNDARY_HOLE))
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid_n, State)
      end if
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done updating grid ", trim(IntToString(Grid_n%properties%id)), "."
      end if
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (UpdateGrid(m) .or. UpdateGrid(n)) then
          Grid_m => Domain%grid(IndexToID(m))
          Grid_n => Domain%grid(IndexToID(n))
          Overlap_mn => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
          if (OverlapExists(Overlap_mn)) then
            call UpdateOverlapAfterCut(Grid_m, Grid_n, Overlap_mn)
          end if
        end if
      end do
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done updating overlap information on grid ", &
          trim(IntToString(Grid_n%properties%id)), "."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished updating grids and overlap information."
    end if

  end subroutine CutBoundaryHoles

  subroutine ComputeOverlapResolutions(Domain, ReducedDomainInfo, OverlapResolutions)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), intent(out) :: OverlapResolutions

    integer :: m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID

    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap

    NumDims = Domain%properties%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = 1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        Overlap => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
        if (OverlapExists(Overlap)) then
          call ovkOverlapCollect(Grid_m, Overlap, OVK_COLLECT_INTERPOLATE, Grid_m%resolution, &
            OverlapResolutions(m,n))
        end if
      end do
    end do

  end subroutine ComputeOverlapResolutions

  subroutine DetectOccludedPoints(Domain, ReducedDomainInfo, OverlapResolutions)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), intent(in) :: OverlapResolutions

    integer :: i, j, k, m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:,:), pointer :: Occludes
    integer, dimension(:,:), pointer :: EdgePadding
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap_mn, Overlap_nm
    type(ovk_field_logical), dimension(:,:), allocatable :: PairwiseOcclusionMasks
    type(ovk_field_logical) :: OverlappedMask_m, OverlappedMask_n
    integer(lk) :: PointCount
    integer(lk), dimension(:), allocatable :: PointCounts
    type(ovk_field_logical), dimension(:), allocatable :: OcclusionMasks
    integer :: MinPadding
    type(ovk_field_logical) :: OcclusionOrHoleMask
    type(ovk_field_logical) :: OcclusionEdgeMask
    type(ovk_field_logical) :: PaddingMask
    type(ovk_field_int), pointer :: State

    NumDims = Domain%properties%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    Occludes => ReducedDomainInfo%occludes
    EdgePadding => ReducedDomainInfo%edge_padding

    if (Domain%logger%verbose) then
      write (*, '(a)') "Detecting pairwise occlusion..."
    end if

    allocate(PairwiseOcclusionMasks(NumGrids,NumGrids))
    allocate(OcclusionMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = n+1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        Overlap_mn => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
        Overlap_nm => Domain%overlap(Grid_n%properties%id,Grid_m%properties%id)
        select case (Occludes(m,n))
        case (OVK_AUTO)
          call FindCoarsePoints(Grid_n, Overlap_mn, OverlapResolutions(m,n), &
            PairwiseOcclusionMasks(m,n))
        case (OVK_TRUE)
          PairwiseOcclusionMasks(m,n) = Overlap_mn%mask
        case (OVK_FALSE)
          PairwiseOcclusionMasks(m,n) = ovk_field_logical_(NumDims)
        end select
        select case (Occludes(n,m))
        case (OVK_AUTO)
          call FindCoarsePoints(Grid_m, Overlap_nm, OverlapResolutions(n,m), &
            PairwiseOcclusionMasks(n,m))
        case (OVK_TRUE)
          PairwiseOcclusionMasks(n,m) = Overlap_nm%mask
        case (OVK_FALSE)
          PairwiseOcclusionMasks(n,m) = ovk_field_logical_(NumDims)
        end select
        if (Occludes(m,n) == OVK_AUTO .and. Occludes(n,m) == OVK_AUTO) then
          ! Exclude occluded points that are overlapped by occluded points
          call ovkFindOverlappedPoints(Grid_m, Grid_n, Overlap_mn, PairwiseOcclusionMasks(n,m), &
            OverlappedMask_n)
          PairwiseOcclusionMasks(m,n)%values = PairwiseOcclusionMasks(m,n)%values .and. .not. &
            OverlappedMask_n%values
          call ovkFindOverlappedPoints(Grid_n, Grid_m, Overlap_nm, PairwiseOcclusionMasks(m,n), &
            OverlappedMask_m)
          PairwiseOcclusionMasks(n,m)%values = PairwiseOcclusionMasks(n,m)%values .and. .not. &
            OverlappedMask_m%values
        end if
        if (Domain%logger%verbose) then
          PointCount = ovkCountMask(PairwiseOcclusionMasks(m,n))
          if (PointCount > 0_lk) then
            write (*, '(7a)') "* ", trim(LargeIntToString(PointCount)), " points on grid ", &
              trim(IntToString(Grid_n%properties%id)), " occluded by grid ", &
              trim(IntToString(Grid_m%properties%id)), "."
          end if
          PointCount = ovkCountMask(PairwiseOcclusionMasks(n,m))
          if (PointCount > 0_lk) then
            write (*, '(7a)') "* ", trim(LargeIntToString(PointCount)), " points on grid ", &
              trim(IntToString(Grid_m%properties%id)), " occluded by grid ", &
              trim(IntToString(Grid_n%properties%id)), "."
          end if
        end if
      end do
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished detecting pairwise occlusion."
    end if

    if (Domain%logger%verbose) then
      write (*, '(a)') "Applying edge padding..."
    end if

    if (Domain%logger%verbose) then
      allocate(PointCounts(NumGrids))
    end if

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(Occludes(:,n) /= OVK_FALSE)) then
        if (Domain%logger%verbose) then
          do m = 1, NumGrids
            if (Occludes(m,n) /= OVK_FALSE) then
              PointCounts(m) = ovkCountMask(PairwiseOcclusionMasks(m,n))
            end if
          end do
        end if
        MinPadding = huge(0)
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_FALSE) then
            MinPadding = min(MinPadding, EdgePadding(m,n))
          end if
        end do
        OcclusionOrHoleMask = ovk_field_logical_(Grid_n%cart)
        OcclusionOrHoleMask%values = .not. Grid_n%mask%values
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_FALSE) then
            OcclusionOrHoleMask%values = OcclusionOrHoleMask%values .or. &
              PairwiseOcclusionMasks(m,n)%values
          end if
        end do
        call ovkDetectEdge(OcclusionOrHoleMask, OVK_OUTER_EDGE, OVK_MIRROR, .false., &
          OcclusionEdgeMask)
        OcclusionMasks(n) = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_FALSE) then
            if (EdgePadding(m,n) > MinPadding) then
              OcclusionOrHoleMask%values = PairwiseOcclusionMasks(m,n)%values .or. .not. &
                Grid_n%mask%values
              call ovkDetectEdge(OcclusionOrHoleMask, OVK_OUTER_EDGE, OVK_MIRROR, .false., &
                PaddingMask)
              PaddingMask%values = PaddingMask%values .and. OcclusionEdgeMask%values
              call ovkDilate(PaddingMask, EdgePadding(m,n)-MinPadding, OVK_FALSE)
              PairwiseOcclusionMasks(m,n)%values = PairwiseOcclusionMasks(m,n)%values .and. .not. &
                PaddingMask%values
            end if
            OcclusionMasks(n)%values = OcclusionMasks(n)%values .or. &
              PairwiseOcclusionMasks(m,n)%values
          end if
        end do
        OcclusionOrHoleMask%values = OcclusionMasks(n)%values .or. .not. Grid_n%mask%values
        call ovkErode(OcclusionOrHoleMask, MinPadding, OVK_TRUE)
        OcclusionMasks(n)%values = OcclusionOrHoleMask%values .and. Grid_n%mask%values
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_FALSE) then
            PairwiseOcclusionMasks(m,n)%values = PairwiseOcclusionMasks(m,n)%values .and. &
              OcclusionMasks(n)%values
            if (Domain%logger%verbose) then
              Grid_m => Domain%grid(IndexToID(m))
              PointCount = PointCounts(m) - ovkCountMask(PairwiseOcclusionMasks(m,n))
              if (PointCount > 0_lk) then
                write (*, '(7a)') "* ", trim(LargeIntToString(PointCount)), " points on grid ", &
                  trim(IntToString(Grid_n%properties%id)), " marked as not occluded by grid ", &
                  trim(IntToString(Grid_m%properties%id)), " due to edge padding."
              end if
            end if
          end if
        end do
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished applying edge padding."
    end if

    if (Domain%logger%verbose) then
      write (*, '(a)') "Accumulating occlusion..."
    end if

    do n = 1, NumGrids
      if (any(Occludes(:,n) /= OVK_FALSE)) then
        Grid_n => Domain%grid(IndexToID(n))
        call ovkEditGridState(Grid_n, State)
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              if (OcclusionMasks(n)%values(i,j,k)) then
                State%values(i,j,k) = ior(State%values(i,j,k),OVK_STATE_OCCLUDED)
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid_n, State)
        if (Domain%logger%verbose) then
          PointCount = ovkCountMask(OcclusionMasks(n))
          if (PointCount > 0_lk) then
            write (*, '(5a)') "* ", trim(LargeIntToString(PointCount)), &
              " occluded points on grid ", trim(IntToString(Grid_n%properties%id)), "."
          end if
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished accumulating occlusion."
    end if

  end subroutine DetectOccludedPoints

  subroutine FindCoarsePoints(OverlappedGrid, Overlap, OverlapResolutions, CoarseMask)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_real), intent(in) :: OverlapResolutions
    type(ovk_field_logical), intent(out) :: CoarseMask

    integer :: i, j, k
    integer(lk) :: l

    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    CoarseMask = ovk_field_logical_(OverlappedGrid%cart, .false.)

    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (Overlap%mask%values(i,j,k)) then
            CoarseMask%values(i,j,k) = OverlappedGrid%resolution%values(i,j,k) < &
              (1._rk-TOLERANCE) * OverlapResolutions%values(l)
            l = l + 1_lk
          end if
        end do
      end do
    end do

  end subroutine FindCoarsePoints

  subroutine LocateCandidateReceivers(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:), pointer :: FringeSize
    type(ovk_grid), pointer :: Grid
    type(ovk_field_logical) :: EdgeMask
    type(ovk_field_logical) :: NonSubsetMask
    type(ovk_field_logical) :: NonSubsetEdgeMask
    type(ovk_field_logical) :: SubsetMask
    type(ovk_field_logical) :: SubsetEdgeMask
    type(ovk_field_logical) :: OcclusionMask
    type(ovk_field_logical) :: ReceiverMask
    type(ovk_field_int), pointer :: State
    integer(lk) :: NumReceivers

    NumDims = Domain%properties%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    FringeSize => ReducedDomainInfo%fringe_size

    if (Domain%logger%verbose) then
      write (*, '(a)') "Locating candidate receiver points..."
    end if

    do n = 1, NumGrids
      Grid => Domain%grid(IndexToID(n))
      ReceiverMask = ovk_field_logical_(Grid%cart, .false.)
      if (FringeSize(n) > 0) then
        call ovkDetectEdge(Grid%mask, OVK_OUTER_EDGE, OVK_FALSE, .true., EdgeMask)
        NonSubsetMask = ovk_field_logical_(Grid%cart)
        NonSubsetMask%values = Grid%boundary_mask%values .or. Grid%internal_boundary_mask%values
        call ovkDetectEdge(NonSubsetMask, OVK_OUTER_EDGE, OVK_FALSE, .true., NonSubsetEdgeMask)
        SubsetMask = ovk_field_logical_(Grid%cart)
        SubsetMask%values = Grid%mask%values .and. .not. NonSubsetMask%values
        call ovkDetectEdge(SubsetMask, OVK_OUTER_EDGE, OVK_FALSE, .true., SubsetEdgeMask)
        EdgeMask%values = EdgeMask%values .and. .not. (NonSubsetEdgeMask%values .and. .not. &
          SubsetEdgeMask%values)
        call ovkDilate(EdgeMask, FringeSize(n), OVK_FALSE)
        ReceiverMask%values = ReceiverMask%values .or. (EdgeMask%values( &
          Grid%cart%is(1):Grid%cart%ie(1),Grid%cart%is(2):Grid%cart%ie(2), &
          Grid%cart%is(3):Grid%cart%ie(3)) .and. Grid%mask%values)
      end if
      call ovkFilterGridState(Grid, ior(OVK_STATE_GRID, OVK_STATE_OCCLUDED), OVK_ALL, OcclusionMask)
      ReceiverMask%values = ReceiverMask%values .or. OcclusionMask%values
      call ovkEditGridState(Grid, State)
      do k = Grid%cart%is(3), Grid%cart%ie(3)
        do j = Grid%cart%is(2), Grid%cart%ie(2)
          do i = Grid%cart%is(1), Grid%cart%ie(1)
            if (ReceiverMask%values(i,j,k)) then
              State%values(i,j,k) = ior(State%values(i,j,k),OVK_STATE_RECEIVER)
            end if
          end do
        end do
      end do
      call ovkReleaseGridState(Grid, State)
      if (Domain%logger%verbose) then
        NumReceivers = ovkCountMask(ReceiverMask)
        if (NumReceivers > 0_lk) then
          write (*, '(5a)') "* ", trim(LargeIntToString(NumReceivers)), &
            " candidate receiver points on grid ", trim(IntToString(Grid%properties%id)), "."
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished locating candidate receiver points."
    end if

  end subroutine LocateCandidateReceivers

  subroutine ChooseDonors(Domain, ReducedDomainInfo, OverlapResolutions, DonorGridIDs)

    type(ovk_domain), intent(in) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), intent(in) :: OverlapResolutions
    type(ovk_field_int), dimension(:), intent(out) :: DonorGridIDs

    integer :: i, j, k, m, n
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:,:), pointer :: ConnectionType
    integer, dimension(:,:), pointer :: EdgePadding
    integer(lk) :: l
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap
    type(ovk_field_logical) :: ReceiverMask
    type(ovk_field_logical) :: OcclusionMask
    type(ovk_field_logical) :: OverlappedByOccludedMask
    type(ovk_field_logical) :: DonorIsOccluded
    type(ovk_field_real) :: DonorEdgeDistance
    type(ovk_field_logical) :: NoOverlapMask
    type(ovk_field_int) :: OverlapEdgeDistances
    type(ovk_field_real) :: DonorResolution
    logical :: IsOccluded
    real(rk) :: EdgeDistance
    real(rk) :: Resolution
    logical :: BetterDonor
    type(ovk_field_logical) :: OrphanMask
    type(ovk_field_int), pointer :: State
    integer :: NumWarnings
    integer, dimension(MAX_ND) :: Point
    integer(lk) :: NumOrphans

    real(rk), parameter :: TOLERANCE = 1.e-12_rk

    if (Domain%logger%verbose) then
      write (*, '(a)') "Choosing donors..."
    end if

    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    ConnectionType => ReducedDomainInfo%connection_type
    EdgePadding => ReducedDomainInfo%edge_padding

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      DonorGridIDs(n) = ovk_field_int_(Grid_n%cart, 0)
      DonorIsOccluded = ovk_field_logical_(Grid_n%cart)
      DonorEdgeDistance = ovk_field_real_(Grid_n%cart)
      DonorResolution = ovk_field_real_(Grid_n%cart)
      call ovkFilterGridState(Grid_n, OVK_STATE_RECEIVER, OVK_ANY, ReceiverMask)
      do m = 1, NumGrids
        if (ConnectionType(m,n) == OVK_CONNECTION_NONE) cycle
        Grid_m => Domain%grid(IndexToID(m))
        Overlap => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
        call ovkFilterGridState(Grid_m, OVK_STATE_OCCLUDED, OVK_ANY, OcclusionMask)
        call ovkFindOverlappedPoints(Grid_m, Grid_n, Overlap, OcclusionMask, &
          OverlappedByOccludedMask)
        NoOverlapMask = ovk_field_logical_(Grid_n%cart)
        NoOverlapMask%values = (.not. Overlap%mask%values .or. OverlappedByOccludedMask%values) &
          .and. Grid_n%mask%values
        call ovkDistanceField(NoOverlapMask, OVK_MIRROR, OverlapEdgeDistances)
        l = 1_lk
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              if (ReceiverMask%values(i,j,k) .and. Overlap%mask%values(i,j,k)) then
                if (DonorGridIDs(n)%values(i,j,k) == 0) then
                  DonorGridIDs(n)%values(i,j,k) = Grid_m%properties%id
                  DonorIsOccluded%values(i,j,k) = OverlappedByOccludedMask%values(i,j,k)
                  DonorEdgeDistance%values(i,j,k) = min(real(OverlapEdgeDistances%values(i,j,k), &
                    kind=rk)/real(max(EdgePadding(m,n),1),kind=rk),1._rk)
                  DonorResolution%values(i,j,k) = OverlapResolutions(m,n)%values(l)
                else
                  IsOccluded = OverlappedByOccludedMask%values(i,j,k)
                  EdgeDistance = min(real(OverlapEdgeDistances%values(i,j,k),kind=rk)/ &
                    real(max(EdgePadding(m,n),1),kind=rk),1._rk)
                  Resolution = OverlapResolutions(m,n)%values(l)
                  if (DonorIsOccluded%values(i,j,k) .neqv. IsOccluded) then
                    BetterDonor = .not. IsOccluded
                  else if (abs(EdgeDistance-DonorEdgeDistance%values(i,j,k)) > TOLERANCE) then
                    BetterDonor = EdgeDistance > DonorEdgeDistance%values(i,j,k)
                  else
                    BetterDonor = Resolution > DonorResolution%values(i,j,k)
                  end if
                  if (BetterDonor) then
                    DonorGridIDs(n)%values(i,j,k) = Grid_m%properties%id
                    DonorIsOccluded%values(i,j,k) = IsOccluded
                    DonorEdgeDistance%values(i,j,k) = EdgeDistance
                    DonorResolution%values(i,j,k) = Resolution
                  end if
                end if
              end if
              if (Overlap%mask%values(i,j,k)) then
                l = l + 1_lk
              end if
            end do
          end do
        end do
      end do
      OrphanMask = ovk_field_logical_(Grid_n%cart)
      OrphanMask%values = ReceiverMask%values .and. DonorGridIDs(n)%values == 0
      call ovkEditGridState(Grid_n, State)
      do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
        do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
          do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
            if (OrphanMask%values(i,j,k)) then
              State%values(i,j,k) = ior(iand(State%values(i,j,k), not(OVK_STATE_RECEIVER)), &
                OVK_STATE_ORPHAN)
            end if
          end do
        end do
      end do
      call ovkReleaseGridState(Grid_n, State)
      if (Domain%logger%verbose) then
        NumWarnings = 0
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              Point = [i,j,k]
              if (OrphanMask%values(i,j,k)) then
                if (NumWarnings <= 100) then
                  write (ERROR_UNIT, '(4a)') "WARNING: Could not find suitable donor for point ", &
                    trim(TupleToString(Point(:Grid_n%cart%nd))), " of grid ", &
                    trim(IntToString(Grid_n%properties%id))
                  if (NumWarnings == 100) then
                    write (ERROR_UNIT, '(a)') "WARNING: Further warnings suppressed."
                  end if
                  NumWarnings = NumWarnings + 1
                end if
              end if
            end do
          end do
        end do
        NumOrphans = ovkCountMask(OrphanMask)
        write (*, '(5a)') "* Done choosing donors on grid ", &
          trim(IntToString(Grid_n%properties%id)), " (", trim(LargeIntToString(NumOrphans)), &
          " orphans)."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished choosing donors."
    end if

  end subroutine ChooseDonors

  subroutine MinimizeOverlap(Domain, ReducedDomainInfo, DonorGridIDs)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_int), dimension(:), intent(in) :: DonorGridIDs

    integer :: i, j, k, m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:), pointer :: FringeSize
    logical, dimension(:,:), pointer :: OverlapMinimization
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap
    type(ovk_field_logical), dimension(:), allocatable :: OverlapMinimizationMasks
    type(ovk_field_logical) :: RemovableMask
    integer(lk) :: NumRemoved
    logical, dimension(:), allocatable :: UpdateGrid
    type(ovk_field_int), pointer :: State

    NumDims = Domain%properties%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    FringeSize => ReducedDomainInfo%fringe_size
    OverlapMinimization => ReducedDomainInfo%overlap_minimization

    if (Domain%logger%verbose) then
      write (*, '(a)') "Minimizing overlap..."
    end if

    allocate(OverlapMinimizationMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(OverlapMinimization(:,n))) then
        OverlapMinimizationMasks(n) = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          if (OverlapMinimization(m,n)) then
            OverlapMinimizationMasks(n)%values = OverlapMinimizationMasks(n)%values .or. &
              DonorGridIDs(n)%values == IndexToID(m)
          end if
        end do
        if (FringeSize(n) > 0) then
          call ovkFilterGridState(Grid_n, ior(OVK_STATE_OCCLUDED, OVK_STATE_HOLE), OVK_ANY, &
            RemovableMask)
          call ovkErode(RemovableMask, FringeSize(n), OVK_TRUE)
          OverlapMinimizationMasks(n)%values = OverlapMinimizationMasks(n)%values .and. &
            Grid_n%mask%values .and. RemovableMask%values
        else
          OverlapMinimizationMasks(n)%values = OverlapMinimizationMasks(n)%values .and. &
            Grid_n%mask%values
        end if
        if (Domain%logger%verbose) then
          NumRemoved = ovkCountMask(OverlapMinimizationMasks(n))
          if (NumRemoved > 0_lk) then
            write (*, '(5a)') "* ", trim(LargeIntToString(NumRemoved)), &
              " points removed from grid ", trim(IntToString(Grid_n%properties%id)), "."
          end if
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished minimizing overlap."
    end if

    if (Domain%logger%verbose) then
      write (*, '(a)') "Updating grids and overlap information after overlap minimization..."
    end if

    allocate(UpdateGrid(NumGrids))

    do n = 1, NumGrids
      UpdateGrid(n) = .false.
      if (any(OverlapMinimization(:,n))) then
        if (any(OverlapMinimizationMasks(n)%values)) then
          UpdateGrid(n) = .true.
        end if
      end if
    end do

    do n = 1, NumGrids
      if (UpdateGrid(n)) then
        Grid_n => Domain%grid(IndexToID(n))
        call ovkEditGridState(Grid_n, State)
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              if (OverlapMinimizationMasks(n)%values(i,j,k)) then
                State%values(i,j,k) = iand(State%values(i,j,k), not(ior(OVK_STATE_GRID, &
                  OVK_STATE_RECEIVER)))
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_HOLE, &
                  OVK_STATE_OVERLAP_MINIMIZED))
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid_n, State)
      end if
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done updating grid ", trim(IntToString(Grid_n%properties%id)), "."
      end if
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (UpdateGrid(m) .or. UpdateGrid(n)) then
          Grid_m => Domain%grid(IndexToID(m))
          Grid_n => Domain%grid(IndexToID(n))
          Overlap => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
          if (OverlapExists(Overlap)) then
            call UpdateOverlapAfterCut(Grid_m, Grid_n, Overlap)
          end if
        end if
      end do
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done updating overlap information on grid ", &
          trim(IntToString(Grid_n%properties%id)), "."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished updating grids and overlap information."
    end if

  end subroutine MinimizeOverlap

  subroutine FillConnectivity(Domain, ReducedDomainInfo, DonorGridIDs)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_int), dimension(:), intent(in) :: DonorGridIDs

    integer :: d, i, j, k, m, n
    integer(lk) :: l, p
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:,:), pointer :: ConnectionType
    integer, dimension(:,:), pointer :: InterpScheme
    integer(lk), dimension(:), allocatable :: NumConnections
    integer :: NumWarnings
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap
    type(ovk_connectivity), pointer :: Connectivity
    type(ovk_field_logical) :: ReceiverMask
    integer, dimension(MAX_ND) :: ReceiverPoint
    integer, dimension(MAX_ND,2) :: DonorExtents
    real(rk), dimension(Domain%properties%nd) :: ReceiverCoords
    real(rk), dimension(Domain%properties%nd) :: DonorCoords
    real(rk), dimension(:,:), allocatable :: DonorInterpCoefs
    logical :: Success
    character(len=STRING_LENGTH) :: DonorCellString
    character(len=STRING_LENGTH) :: DonorGridIDString
    character(len=STRING_LENGTH) :: ReceiverPointString
    character(len=STRING_LENGTH) :: ReceiverGridIDString

    if (Domain%logger%verbose) then
      write (*, '(a)') "Generating connectivity information..."
    end if

    NumDims = Domain%properties%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    ConnectionType => ReducedDomainInfo%connection_type
    InterpScheme => ReducedDomainInfo%interp_scheme

    allocate(NumConnections(NumGrids))

    do n = 1, NumGrids

      Grid_n => Domain%grid(IndexToID(n))

      call ovkFilterGridState(Grid_n, OVK_STATE_RECEIVER, OVK_ALL, ReceiverMask)

      NumConnections = 0_lk
      do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
        do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
          do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
            if (ReceiverMask%values(i,j,k)) then
              m = DonorGridIDs(n)%values(i,j,k)
              NumConnections(m) = NumConnections(m) + 1_lk
            end if
          end do
        end do
      end do

      do m = 1, NumGrids

        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE) then

          Grid_m => Domain%grid(IndexToID(m))
          Overlap => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
          Connectivity => Domain%connectivity(Grid_m%properties%id,Grid_n%properties%id)

          if (NumConnections(m) > 0_lk) then

            call ResizeConnectivity(Connectivity, NumConnections(m))

            allocate(DonorInterpCoefs(size(Connectivity%donor_interp_coefs,1), &
              size(Connectivity%donor_interp_coefs,2)))

            NumWarnings = 0

            l = 1_lk
            p = 1_lk
            do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
              do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
                do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
                  if (ReceiverMask%values(i,j,k) .and. DonorGridIDs(n)%values(i,j,k) == &
                    Grid_m%properties%id) then
                    ReceiverPoint = [i,j,k]
                    DonorExtents(:,1) = Overlap%cells(:,l)
                    DonorExtents(:NumDims,2) = Overlap%cells(:NumDims,l) + 1
                    DonorExtents(NumDims+1:,2) = 1
                    DonorCoords = Overlap%coords(:,l)
                    select case (InterpScheme(m,n))
                    case (OVK_INTERP_LINEAR)
                      do d = 1, NumDims
                        DonorInterpCoefs(:,d) = ovkInterpBasisLinear(DonorCoords(d))
                      end do
                    case (OVK_INTERP_CUBIC)
                      do d = 1, NumDims
                        ReceiverCoords(d) = Grid_n%coords(d)%values(i,j,k)
                      end do
                      call ExpandDonorCell(Grid_m, ReceiverCoords, DonorExtents, DonorCoords, &
                        Success)
                      if (Success) then
                        do d = 1, NumDims
                          DonorInterpCoefs(:,d) = ovkInterpBasisCubic(DonorCoords(d))
                        end do
                      else
                        if (Domain%logger%verbose) then
                          DonorCellString = TupleToString(DonorExtents(:Grid_m%cart%nd,1))
                          DonorGridIDString = IntToString(Grid_m%properties%id)
                          ReceiverPointString = TupleToString(ReceiverPoint(:Grid_n%cart%nd))
                          ReceiverGridIDString = IntToString(Grid_n%properties%id)
                          write (ERROR_UNIT, '(9a)') "WARNING: Could not use cubic ", &
                            "interpolation for donor cell ", trim(DonorCellString), " of grid ", &
                            trim(DonorGridIDString), " corresponding to receiver point ", &
                            trim(ReceiverPointString), " of grid ", trim(ReceiverGridIDString)
                          if (NumWarnings == 100) then
                            write (ERROR_UNIT, '(a)') "WARNING: Further warnings suppressed."
                          end if
                          NumWarnings = NumWarnings + 1
                        end if
                        DonorInterpCoefs = 0._rk
                        do d = 1, NumDims
                          DonorInterpCoefs(:2,d) = ovkInterpBasisLinear(DonorCoords(d))
                        end do
                      end if
                    end select
                    Connectivity%receiver_points(:,p) = ReceiverPoint
                    Connectivity%donor_extents(:,:,p) = DonorExtents
                    Connectivity%donor_coords(:,p) = DonorCoords
                    Connectivity%donor_interp_coefs(:,:,p) = DonorInterpCoefs
                    p = p + 1_lk
                  end if
                  if (Overlap%mask%values(i,j,k)) then
                    l = l + 1_lk
                  end if
                end do
              end do
            end do

            deallocate(DonorInterpCoefs)

            if (Domain%logger%verbose) then
              if (NumConnections(m) > 0_lk) then
                write (*, '(7a)') "* ", trim(LargeIntToString(NumConnections(m))), &
                  " donor/receiver pairs between grid ", trim(IntToString(Grid_m%properties%id)), &
                  " and grid ", trim(IntToString(Grid_n%properties%id)), "."
              end if
            end if

          end if

        end if

      end do

    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished generating connectivity information."
    end if

  end subroutine FillConnectivity

  subroutine FinalizeAssembly(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(inout) :: ReducedDomainInfo

    deallocate(ReducedDomainInfo%index_to_id)
    deallocate(ReducedDomainInfo%overlappable)
    deallocate(ReducedDomainInfo%overlap_tolerance)
    deallocate(ReducedDomainInfo%overlap_accel_quality_adjust)
    deallocate(ReducedDomainInfo%infer_boundaries)
    deallocate(ReducedDomainInfo%boundary_hole_cutting)
    deallocate(ReducedDomainInfo%occludes)
    deallocate(ReducedDomainInfo%connection_type)
    deallocate(ReducedDomainInfo%interp_scheme)
    deallocate(ReducedDomainInfo%fringe_size)
    deallocate(ReducedDomainInfo%edge_padding)
    deallocate(ReducedDomainInfo%overlap_minimization)

    call ResetDomainEdits(Domain)

  end subroutine FinalizeAssembly

  subroutine WriteGrids(Domain, IBlanks, FilePath)

    type(ovk_domain), intent(in) :: Domain
    type(ovk_field_int), dimension(:), intent(in) :: IBlanks
    character(len=*), intent(in) :: FilePath

    integer :: n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:,:), allocatable :: NumPointsAll
    type(ovk_grid), pointer :: Grid
    type(ovk_plot3d_grid_file) :: GridFile

    NumDims = Domain%properties%nd
    NumGrids = Domain%properties%ngrids

    allocate(NumPointsAll(MAX_ND,NumGrids))
    do n = 1, NumGrids
      Grid => Domain%grid(n)
      NumPointsAll(:,n) = ovkCartSize(Grid%cart)
    end do

    call ovkCreateP3D(GridFile, FilePath, NumDims=NumDims, NumGrids=NumGrids, &
      NumPointsAll=NumPointsAll, WithIBlank=.true., Verbose=.true.)

    select case (NumDims)
    case (2)
      do n = 1, NumGrids
        Grid => Domain%grid(n)
        call ovkWriteP3D(GridFile, n, Grid%coords(1), Grid%coords(2), IBlanks(n))
      end do
    case (3)
      do n = 1, NumGrids
        Grid => Domain%grid(n)
        call ovkWriteP3D(GridFile, n, Grid%coords(1), Grid%coords(2), Grid%coords(3), IBlanks(n))
      end do
    end select

    call ovkCloseP3D(GridFile)

  end subroutine WriteGrids

end module ovkAssembly
