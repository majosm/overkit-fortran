! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkAssembly

  use ovkArray
  use ovkAssemblyOptions
  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDomain
  use ovkDonorStencil
  use ovkOverlap
  use ovkOverlapAccel
  use ovkField
  use ovkFieldOps
  use ovkGlobal
  use ovkGrid
  use ovkLogger
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
    logical, dimension(:,:), pointer :: cut_boundary_holes
    integer, dimension(:,:), pointer :: occludes
    integer, dimension(:,:), pointer :: edge_padding
    integer, dimension(:), pointer :: edge_smoothing
    integer, dimension(:,:), pointer :: connection_type
    integer, dimension(:), pointer :: fringe_size
    logical, dimension(:,:), pointer :: minimize_overlap
  end type t_reduced_domain_info

contains

  subroutine ovkAssemble(Domain, AssemblyOptions)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_assembly_options), intent(in) :: AssemblyOptions

    type(t_logger) :: Logger
    integer :: ClockInitial, ClockFinal, ClockRate
    type(t_reduced_domain_info) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), allocatable :: OverlapVolumes
    type(ovk_field_logical), dimension(:,:), allocatable :: PairwiseOcclusionMasks
    type(ovk_field_int), dimension(:), allocatable :: DonorGridIDs

    Logger = Domain%logger

    if (Logger%log_status) then
      call system_clock(ClockInitial, ClockRate)
      write (Logger%status_file, '(a)') "Overset grid assembly started..."
    end if

    if (OVK_DEBUG) then
      if ( &
        any(Domain%grid_edit_ref_counts > 0) .or. &
        any(Domain%connectivity_edit_ref_counts > 0) &
      ) then
        write (ERROR_UNIT, '(a)') "ERROR: Cannot perform assembly; domain is still being edited."
        stop 1
      end if
    end if

    call InitAssembly(Domain, AssemblyOptions, ReducedDomainInfo)

    call CollideGrids(Domain, ReducedDomainInfo)

    call InferNonOverlappingBoundaries(Domain, ReducedDomainInfo)

    call CutHoles(Domain, ReducedDomainInfo)

    allocate(OverlapVolumes(ReducedDomainInfo%ngrids,ReducedDomainInfo%ngrids))

    call ComputeOverlapVolumes(Domain, ReducedDomainInfo, OverlapVolumes)

    call LocateOuterFringe(Domain, ReducedDomainInfo)

    allocate(PairwiseOcclusionMasks(ReducedDomainInfo%ngrids,ReducedDomainInfo%ngrids))

    call DetectOccludedPoints(Domain, ReducedDomainInfo, OverlapVolumes, PairwiseOcclusionMasks)

    call ApplyOverlapMinimization(Domain, ReducedDomainInfo, PairwiseOcclusionMasks)

    deallocate(PairwiseOcclusionMasks)

    call LocateReceivers(Domain, ReducedDomainInfo)

    allocate(DonorGridIDs(ReducedDomainInfo%ngrids))

    call ChooseDonors(Domain, ReducedDomainInfo, OverlapVolumes, DonorGridIDs)

    deallocate(OverlapVolumes)

    call GenerateConnectivity(Domain, ReducedDomainInfo, DonorGridIDs)

    call FinalizeAssembly(Domain, ReducedDomainInfo, AssemblyOptions)

    if (Logger%log_status) then
      call system_clock(ClockFinal, ClockRate)
      write (Logger%status_file, '(a,f0.3,a)') "Overset grid assembly finished (time: ", &
        real(ClockFinal-ClockInitial,kind=rk)/real(ClockRate,kind=rk), " seconds)."
    end if

  end subroutine ovkAssemble

  subroutine InitAssembly(Domain, AssemblyOptions, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_assembly_options), intent(in) :: AssemblyOptions
    type(t_reduced_domain_info), intent(out) :: ReducedDomainInfo

    integer :: m, n, p, q
    integer :: NumGrids
    logical :: GridExists
    logical :: GridOverlaps
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: Overlappable
    real(rk), dimension(:,:), pointer :: OverlapTolerance
    real(rk), dimension(:), pointer :: OverlapAccelQualityAdjust
    logical, dimension(:), pointer :: InferBoundaries
    logical, dimension(:,:), pointer :: CutBoundaryHoles
    integer, dimension(:,:), pointer :: Occludes
    integer, dimension(:,:), pointer :: EdgePadding
    integer, dimension(:), pointer :: EdgeSmoothing
    integer, dimension(:,:), pointer :: ConnectionType
    integer, dimension(:), pointer :: FringeSize
    logical, dimension(:,:), pointer :: MinimizeOverlap

    call PrintDomainSummary(Domain)

    ! Clear out any old assembly data
    do q = 1, Domain%ngrids
      do p = 1, Domain%ngrids
        if (ovkHasOverlap(Domain, p, q)) call ovkDestroyOverlap(Domain, p, q)
        if (ovkHasConnectivity(Domain, p, q)) call ovkDestroyConnectivity(Domain, p, q)
      end do
      call ovkResetGridState(Domain%grid(q))
    end do

    ! Ignore empty and non-overlapping grids
    NumGrids = 0
    do q = 1, Domain%ngrids
      GridExists = ovkGridExists(Domain%grid(q))
      GridOverlaps = .false.
      do p = 1, Domain%ngrids
        if (p /= q) then
          GridOverlaps = GridOverlaps .or. AssemblyOptions%overlappable(p,q) .or. &
            AssemblyOptions%overlappable(q,p)
        end if
      end do
      if (GridExists .and. GridOverlaps) then
        NumGrids = NumGrids + 1
      end if
    end do

    ReducedDomainInfo%ngrids = NumGrids

    allocate(ReducedDomainInfo%index_to_id(NumGrids))
    IndexToID => ReducedDomainInfo%index_to_id

    n = 1
    do q = 1, Domain%ngrids
      GridExists = ovkGridExists(Domain%grid(q))
      GridOverlaps = .false.
      do p = 1, Domain%ngrids
        if (p /= q) then
          GridOverlaps = GridOverlaps .or. AssemblyOptions%overlappable(p,q) .or. &
            AssemblyOptions%overlappable(q,p)
        end if
      end do
      if (GridExists .and. GridOverlaps) then
        IndexToID(n) = q
        n = n + 1
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
        if (p /= q) then
          Overlappable(m,n) = AssemblyOptions%overlappable(p,q)
        else
          Overlappable(m,n) = .false.
        end if
        ! Tolerance only used if overlappable
        if (Overlappable(m,n)) then
          OverlapTolerance(m,n) = AssemblyOptions%overlap_tolerance(p,q)
        else
          OverlapTolerance(m,n) = 0._rk
        end if
      end do
      ! Quality adjust only used if overlappable
      if (any(Overlappable(:,n))) then
        OverlapAccelQualityAdjust(n) = AssemblyOptions%overlap_accel_quality_adjust(q)
      else
        OverlapAccelQualityAdjust(n) = 0._rk
      end if
    end do

    allocate(ReducedDomainInfo%infer_boundaries(NumGrids))
    allocate(ReducedDomainInfo%cut_boundary_holes(NumGrids,NumGrids))
    InferBoundaries => ReducedDomainInfo%infer_boundaries
    CutBoundaryHoles => ReducedDomainInfo%cut_boundary_holes

    do n = 1, NumGrids
      q = IndexToID(n)
      InferBoundaries(n) = AssemblyOptions%infer_boundaries(q)
      do m = 1, NumGrids
        p = IndexToID(m)
        ! Need overlap data to cut
        if (Overlappable(m,n)) then
          CutBoundaryHoles(m,n) = AssemblyOptions%cut_boundary_holes(p,q)
        else
          CutBoundaryHoles(m,n) = .false.
        end if
      end do
    end do

    allocate(ReducedDomainInfo%occludes(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%edge_padding(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%edge_smoothing(NumGrids))
    Occludes => ReducedDomainInfo%occludes
    EdgePadding => ReducedDomainInfo%edge_padding
    EdgeSmoothing => ReducedDomainInfo%edge_smoothing

    do n = 1, NumGrids
      q = IndexToID(n)
      do m = 1, NumGrids
        p = IndexToID(m)
        ! Can't occlude if it doesn't overlap
        if (Overlappable(m,n)) then
          Occludes(m,n) = AssemblyOptions%occludes(p,q)
        else
          Occludes(m,n) = OVK_OCCLUDES_NONE
        end if
        ! Padding still needed even if occlusion is disabled (used for deciding where to switch from
        ! one donor grid to the other in 3-grid interfaces)
        EdgePadding(m,n) = AssemblyOptions%edge_padding(p,q)
      end do
      ! Smoothing only used if some occlusion exists
      if (any(Occludes(:,n) /= OVK_OCCLUDES_NONE)) then
        EdgeSmoothing(n) = AssemblyOptions%edge_smoothing(q)
      else
        EdgeSmoothing(n) = 0
      end if
    end do

    allocate(ReducedDomainInfo%connection_type(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%fringe_size(NumGrids))
    allocate(ReducedDomainInfo%minimize_overlap(NumGrids,NumGrids))
    ConnectionType => ReducedDomainInfo%connection_type
    FringeSize => ReducedDomainInfo%fringe_size
    MinimizeOverlap => ReducedDomainInfo%minimize_overlap

    do n = 1, NumGrids
      q = IndexToID(n)
      do m = 1, NumGrids
        p = IndexToID(m)
        ! Can't communicate if it doesn't overlap
        if (Overlappable(m,n)) then
          ConnectionType(m,n) = AssemblyOptions%connection_type(p,q)
        else
          ConnectionType(m,n) = OVK_CONNECTION_NONE
        end if
        ! Overlap minimization only affects occluded regions
        if (Occludes(m,n) /= OVK_OCCLUDES_NONE) then
          MinimizeOverlap(m,n) = AssemblyOptions%minimize_overlap(p,q)
        else
          MinimizeOverlap(m,n) = .false.
        end if
      end do
      ! Leave fringe size alone -- want orphans if no suitable overlap is found
      FringeSize(n) = AssemblyOptions%fringe_size(q)
    end do

  end subroutine InitAssembly

  subroutine CollideGrids(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: m, n
    type(t_logger) :: Logger
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

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Detecting overlap between grids..."
    end if

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    Overlappable => ReducedDomainInfo%overlappable
    OverlapTolerance => ReducedDomainInfo%overlap_tolerance
    OverlapAccelQualityAdjust => ReducedDomainInfo%overlap_accel_quality_adjust

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = 1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        if (Overlappable(m,n)) then
          call ovkCreateOverlap(Domain, Grid_m%id, Grid_n%id)
        end if
      end do
    end do

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
        if (Logger%log_status) then
          write (Logger%status_file, '(4a)') "* Generating overlap search accelerator on ", &
            "grid ", trim(IntToString(Grid_m%id)), "..."
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
        if (Logger%log_status) then
          write (Logger%status_file, '(4a)') "* Finished generating overlap search ", &
            "accelerator on grid ", trim(IntToString(Grid_m%id)), "."
        end if
        do n = 1, NumGrids
          Grid_n => Domain%grid(IndexToID(n))
          if (Overlappable(m,n)) then
            Overlap => Domain%overlap(Grid_m%id, Grid_n%id)
            call DetectOverlap(Overlap, OverlapAccel, Bounds(m,n), OverlapTolerance(m,n))
            if (Logger%log_status) then
              if (Overlap%noverlap > 0_lk) then
                write (Logger%status_file, '(7a)') "* Detected ", trim(LargeIntToString( &
                  Overlap%noverlap)), " points on grid ", trim(IntToString(Grid_n%id)), &
                  " overlapped by grid ", trim(IntToString(Grid_m%id)), "."
              end if
            end if
          end if
        end do
        call DestroyOverlapAccel(OverlapAccel)
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished detecting overlap between grids."
    end if

  end subroutine CollideGrids

  subroutine InferNonOverlappingBoundaries(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, m, n
    type(t_logger) :: Logger
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:), pointer :: InferBoundaries
    type(ovk_field_logical) :: InferredBoundaryMask
    type(ovk_field_logical) :: InitialBoundaryMask
    type(ovk_grid), pointer :: Grid
    type(ovk_overlap), pointer :: Overlap
    type(ovk_field_int), pointer :: State
    integer(lk) :: NumBoundaryPoints

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Inferring domain boundaries in non-overlapping regions..."
    end if

    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    InferBoundaries => ReducedDomainInfo%infer_boundaries

    do n = 1, NumGrids
      Grid => Domain%grid(IndexToID(n))
      if (InferBoundaries(n)) then
        call ovkDetectEdge(Grid%mask, OVK_INNER_EDGE, OVK_FALSE, .false., InferredBoundaryMask)
        call ovkFilterGridState(Grid, ior(OVK_STATE_GRID, OVK_STATE_DOMAIN_BOUNDARY), OVK_ALL, &
          InitialBoundaryMask)
        InferredBoundaryMask%values = InferredBoundaryMask%values .and. .not. &
          InitialBoundaryMask%values
        do m = 1, NumGrids
          Overlap => Domain%overlap(IndexToID(m),IndexToID(n))
          if (ovkOverlapExists(Overlap)) then
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
        if (Logger%log_status) then
          NumBoundaryPoints = ovkCountMask(InferredBoundaryMask)
          if (NumBoundaryPoints > 0_lk) then
            write (Logger%status_file, '(5a)') "* ", trim(LargeIntToString(NumBoundaryPoints)), &
              " points marked as boundaries on grid ", trim(IntToString(Grid%id)), "."
          end if
        end if
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(2a)') "Finished inferring domain boundaries in ", &
        "non-overlapping regions."
    end if

  end subroutine InferNonOverlappingBoundaries

  subroutine CutHoles(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, m, n
    type(t_logger) :: Logger
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: CutBoundaryHoles
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_field_logical), dimension(:), allocatable :: OwnBoundaryMasks
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

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Cutting boundary holes..."
    end if

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    CutBoundaryHoles => ReducedDomainInfo%cut_boundary_holes

    allocate(OwnBoundaryMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      call ovkFilterGridState(Grid_n, ior(OVK_STATE_GRID, OVK_STATE_DOMAIN_BOUNDARY), OVK_ALL, &
        OwnBoundaryMasks(n))
    end do

    allocate(BoundaryHoleMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(CutBoundaryHoles(:,n))) then
        BoundaryMask = ovk_field_logical_(Grid_n%cart, .false.)
        InteriorMask = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          Grid_m => Domain%grid(IndexToID(m))
          Overlap_mn => Domain%overlap(Grid_m%id,Grid_n%id)
          Overlap_nm => Domain%overlap(Grid_n%id,Grid_m%id)
          if (CutBoundaryHoles(m,n)) then
            call ovkDetectEdge(Overlap_mn%mask, OVK_OUTER_EDGE, OVK_MIRROR, .false., EdgeMask1)
            call ovkDetectEdge(Overlap_nm%mask, OVK_INNER_EDGE, OVK_FALSE, .false., EdgeMask2)
            call ovkFindOverlappingPoints(Overlap_nm, EdgeMask2, OverlappingMask)
            EdgeMask1%values = EdgeMask1%values .and. .not. OverlappingMask%values
            i = 0
            ! Explicit conversion to logical in order to work around GCC 4.7 bug
            do while (logical(any(EdgeMask1%values)))
              call ovkDilate(OverlappingMask, 1, OVK_FALSE)
              EdgeMask1%values = EdgeMask1%values .and. .not. OverlappingMask%values
              i = i + 1
            end do
            call ovkFindOverlappingPoints(Overlap_nm, OwnBoundaryMasks(m), OverlappingMask)
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
          BoundaryMask%values = BoundaryMask%values .or. OwnBoundaryMasks(n)%values
          call ovkFlood(InteriorMask, BoundaryMask)
          call ovkDetectEdge(Grid_n%mask, OVK_INNER_EDGE, OVK_FALSE, .false., EdgeMask1)
          call ovkDetectEdge(InteriorMask, OVK_OUTER_EDGE, OVK_FALSE, .false., EdgeMask2)
          SpuriousBoundaryMask = ovk_field_logical_(Grid_n%cart)
          SpuriousBoundaryMask%values = OwnBoundaryMasks(n)%values .and. EdgeMask1%values &
            .and. .not. EdgeMask2%values
          BoundaryMask%values = BoundaryMask%values .and. .not. SpuriousBoundaryMask%values
          BoundaryHoleMasks(n)%values = Grid_n%mask%values .and. .not. (InteriorMask%values .or. &
            BoundaryMask%values)
          if (Logger%log_status) then
            NumRemoved = ovkCountMask(BoundaryHoleMasks(n))
            if (NumRemoved > 0_lk) then
              write (Logger%status_file, '(5a)') "* ", trim(LargeIntToString(NumRemoved)), &
                " points removed from grid ", trim(IntToString(Grid_n%id)), "."
            end if
          end if
        end if
      end if
    end do

    deallocate(OwnBoundaryMasks)

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished cutting boundary holes."
    end if

    if (Logger%log_status) then
      write (Logger%status_file, '(2a)') "Updating grids and overlap information after ", &
        "boundary hole cutting..."
    end if

    allocate(UpdateGrid(NumGrids))

    do n = 1, NumGrids
      UpdateGrid(n) = .false.
      if (any(CutBoundaryHoles(:,n))) then
        if (any(BoundaryHoleMasks(n)%values)) then
          UpdateGrid(n) = .true.
        end if
      end if
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (UpdateGrid(n)) then
        call ovkEditGridState(Grid_n, State)
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              if (BoundaryHoleMasks(n)%values(i,j,k)) then
                State%values(i,j,k) = iand(State%values(i,j,k), not(OVK_STATE_GRID))
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_EXTERIOR, &
                  OVK_STATE_BOUNDARY_HOLE))
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid_n, State)
      end if
      if (Logger%log_status) then
        write (Logger%status_file, '(3a)') "* Done updating grid ", trim(IntToString(Grid_n%id)), "."
      end if
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (UpdateGrid(m) .or. UpdateGrid(n)) then
          Grid_m => Domain%grid(IndexToID(m))
          Grid_n => Domain%grid(IndexToID(n))
          Overlap_mn => Domain%overlap(Grid_m%id,Grid_n%id)
          if (ovkOverlapExists(Overlap_mn)) then
            call UpdateOverlapAfterCut(Overlap_mn)
          end if
        end if
      end do
      if (Logger%log_status) then
        write (Logger%status_file, '(3a)') "* Done updating overlap information on grid ", &
          trim(IntToString(Grid_n%id)), "."
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished updating grids and overlap information."
    end if

  end subroutine CutHoles

  subroutine ComputeOverlapVolumes(Domain, ReducedDomainInfo, OverlapVolumes)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), intent(out) :: OverlapVolumes

    integer :: m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID

    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = 1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        Overlap => Domain%overlap(Grid_m%id,Grid_n%id)
        if (ovkOverlapExists(Overlap)) then
          call ovkOverlapCollect(Overlap, OVK_COLLECT_INTERPOLATE, Grid_m%volumes, OverlapVolumes(m,n))
        end if
      end do
    end do

  end subroutine ComputeOverlapVolumes

  subroutine LocateOuterFringe(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, n
    type(t_logger) :: Logger
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:), pointer :: FringeSize
    type(ovk_grid), pointer :: Grid
    type(ovk_field_logical) :: OuterFringeMask
    type(ovk_field_logical) :: BoundaryMask
    type(ovk_field_logical) :: BoundaryEdgeMask
    type(ovk_field_logical) :: NonBoundaryMask
    type(ovk_field_logical) :: NonBoundaryEdgeMask
    type(ovk_field_logical) :: ExtendedOuterFringeMask
    type(ovk_field_int), pointer :: State
    integer(lk) :: NumFringe

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Locating outer fringe points..."
    end if

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    FringeSize => ReducedDomainInfo%fringe_size

    do n = 1, NumGrids
      if (FringeSize(n) > 0) then
        Grid => Domain%grid(IndexToID(n))
        OuterFringeMask = ovk_field_logical_(Grid%cart)
        call ovkFilterGridState(Grid, ior(OVK_STATE_DOMAIN_BOUNDARY, OVK_STATE_INTERNAL_BOUNDARY), &
          OVK_ANY, BoundaryMask)
        BoundaryMask%values = BoundaryMask%values .and. Grid%mask%values
        call ovkDetectEdge(BoundaryMask, OVK_OUTER_EDGE, OVK_FALSE, .true., BoundaryEdgeMask)
        NonBoundaryMask = ovk_field_logical_(Grid%cart)
        NonBoundaryMask%values = Grid%mask%values .and. .not. BoundaryMask%values
        call ovkDetectEdge(NonBoundaryMask, OVK_OUTER_EDGE, OVK_FALSE, .true., NonBoundaryEdgeMask)
        call ovkDetectEdge(Grid%mask, OVK_OUTER_EDGE, OVK_FALSE, .true., ExtendedOuterFringeMask)
        ExtendedOuterFringeMask%values = ExtendedOuterFringeMask%values .and. &
          (NonBoundaryEdgeMask%values .or. .not. BoundaryEdgeMask%values)
        call ovkDilate(ExtendedOuterFringeMask, FringeSize(n), OVK_FALSE)
        OuterFringeMask%values = Grid%mask%values .and. ExtendedOuterFringeMask%values( &
          Grid%cart%is(1):Grid%cart%ie(1),Grid%cart%is(2):Grid%cart%ie(2), &
          Grid%cart%is(3):Grid%cart%ie(3))
        call ovkEditGridState(Grid, State)
        do k = Grid%cart%is(3), Grid%cart%ie(3)
          do j = Grid%cart%is(2), Grid%cart%ie(2)
            do i = Grid%cart%is(1), Grid%cart%ie(1)
              if (OuterFringeMask%values(i,j,k)) then
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_FRINGE, &
                  OVK_STATE_OUTER_FRINGE))
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid, State)
        if (Logger%log_status) then
          NumFringe = ovkCountMask(OuterFringeMask)
          if (NumFringe > 0_lk) then
            write (Logger%status_file, '(5a)') "* ", trim(LargeIntToString(NumFringe)), &
              " outer fringe points on grid ", trim(IntToString(Grid%id)), "."
          end if
        end if
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished locating outer fringe points."
    end if

  end subroutine LocateOuterFringe

  subroutine DetectOccludedPoints(Domain, ReducedDomainInfo, OverlapVolumes, &
    PairwiseOcclusionMasks)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), intent(in) :: OverlapVolumes
    type(ovk_field_logical), dimension(:,:), intent(out) :: PairwiseOcclusionMasks

    integer :: i, j, k, m, n
    type(t_logger) :: Logger
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: Overlappable
    integer, dimension(:,:), pointer :: Occludes
    integer, dimension(:,:), pointer :: EdgePadding
    integer, dimension(:), pointer :: EdgeSmoothing
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap_mn, Overlap_nm
    type(ovk_field_logical) :: OverlappedMask_m, OverlappedMask_n
    integer(lk) :: PointCount
    type(ovk_field_logical), dimension(:,:), allocatable :: PaddingMasks
    type(ovk_field_logical) :: BaseOcclusionMask
    type(ovk_field_logical) :: OcclusionMask
    type(ovk_field_logical) :: OuterFringeMask
    type(ovk_field_int) :: OverlapEdgeDistance
    type(ovk_field_logical) :: EdgeMask
    type(ovk_field_int) :: EdgeDistance
    type(ovk_array_int) :: CellEdgeDistance
    type(ovk_field_int), pointer :: State

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Detecting pairwise occlusion..."
    end if

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    Overlappable => ReducedDomainInfo%overlappable
    Occludes => ReducedDomainInfo%occludes
    EdgePadding => ReducedDomainInfo%edge_padding
    EdgeSmoothing => ReducedDomainInfo%edge_smoothing

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = n+1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        Overlap_mn => Domain%overlap(Grid_m%id,Grid_n%id)
        Overlap_nm => Domain%overlap(Grid_n%id,Grid_m%id)
        select case (Occludes(m,n))
        case (OVK_OCCLUDES_COARSE)
          call FindCoarsePoints(Grid_n, Overlap_mn, OverlapVolumes(m,n), &
            PairwiseOcclusionMasks(m,n))
        case (OVK_OCCLUDES_ALL)
          PairwiseOcclusionMasks(m,n) = Overlap_mn%mask
        case (OVK_OCCLUDES_NONE)
          PairwiseOcclusionMasks(m,n) = ovk_field_logical_(NumDims)
        end select
        select case (Occludes(n,m))
        case (OVK_OCCLUDES_COARSE)
          call FindCoarsePoints(Grid_m, Overlap_nm, OverlapVolumes(n,m), &
            PairwiseOcclusionMasks(n,m))
        case (OVK_OCCLUDES_ALL)
          PairwiseOcclusionMasks(n,m) = Overlap_nm%mask
        case (OVK_OCCLUDES_NONE)
          PairwiseOcclusionMasks(n,m) = ovk_field_logical_(NumDims)
        end select
        if (Occludes(m,n) == OVK_OCCLUDES_COARSE .and. Occludes(n,m) == OVK_OCCLUDES_COARSE) then
          ! Exclude occluded points that are overlapped by occluded points
          call ovkFindOverlappedPoints(Overlap_mn, PairwiseOcclusionMasks(n,m), OverlappedMask_n)
          PairwiseOcclusionMasks(m,n)%values = PairwiseOcclusionMasks(m,n)%values .and. .not. &
            OverlappedMask_n%values
          call ovkFindOverlappedPoints(Overlap_nm, PairwiseOcclusionMasks(m,n), OverlappedMask_m)
          PairwiseOcclusionMasks(n,m)%values = PairwiseOcclusionMasks(n,m)%values .and. .not. &
            OverlappedMask_m%values
        end if
        if (Logger%log_status) then
          PointCount = ovkCountMask(PairwiseOcclusionMasks(m,n))
          if (PointCount > 0_lk) then
            write (Logger%status_file, '(7a)') "* ", trim(LargeIntToString(PointCount)), &
              " points on grid ", trim(IntToString(Grid_n%id)), " occluded by grid ", &
              trim(IntToString(Grid_m%id)), "."
          end if
          PointCount = ovkCountMask(PairwiseOcclusionMasks(n,m))
          if (PointCount > 0_lk) then
            write (Logger%status_file, '(7a)') "* ", trim(LargeIntToString(PointCount)), &
              " points on grid ", trim(IntToString(Grid_m%id)), " occluded by grid ", &
              trim(IntToString(Grid_n%id)), "."
          end if
        end if
      end do
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished detecting pairwise occlusion."
    end if

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Applying edge padding and smoothing..."
    end if

    allocate(PaddingMasks(NumGrids,NumGrids))

    do n = 1, NumGrids
      if (any(Occludes(:,n) /= OVK_OCCLUDES_NONE)) then
        Grid_n => Domain%grid(IndexToID(n))
        BaseOcclusionMask = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_OCCLUDES_NONE) then
            BaseOcclusionMask%values = BaseOcclusionMask%values .or. &
              PairwiseOcclusionMasks(m,n)%values
          end if
        end do
        OcclusionMask = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_OCCLUDES_NONE) then
            Grid_m => Domain%grid(IndexToID(m))
            Overlap_mn => Domain%overlap(Grid_m%id,Grid_n%id)
            call ovkFilterGridState(Grid_m, OVK_STATE_OUTER_FRINGE, OVK_ANY, OuterFringeMask)
            EdgeMask = ovk_field_logical_(Grid_m%cart)
            EdgeMask%values = OuterFringeMask%values
            if (Occludes(n,m) /= OVK_OCCLUDES_NONE) then
              EdgeMask%values = EdgeMask%values .or. PairwiseOcclusionMasks(n,m)%values .or. &
                .not. Grid_m%mask%values
            end if
            call ovkDistanceField(EdgeMask, OVK_MIRROR, EdgeDistance)
            call ovkOverlapCollect(Overlap_mn, OVK_COLLECT_MIN, EdgeDistance, CellEdgeDistance)
            OverlapEdgeDistance = ovk_field_int_(Grid_n%cart, -1)
            call ovkOverlapDisperse(Overlap_mn, OVK_DISPERSE_OVERWRITE, CellEdgeDistance, &
              OverlapEdgeDistance)
            PaddingMasks(m,n) = ovk_field_logical_(Grid_n%cart)
            PaddingMasks(m,n)%values = PairwiseOcclusionMasks(m,n)%values .and. &
              OverlapEdgeDistance%values <= EdgePadding(m,n)
            OcclusionMask%values = OcclusionMask%values .or. (PairwiseOcclusionMasks(m,n)%values &
              .and. .not. PaddingMasks(m,n)%values)
          end if
        end do
        if (EdgeSmoothing(n) > 0) then
          call ovkDilate(OcclusionMask, EdgeSmoothing(n), OVK_MIRROR)
          call ovkErode(OcclusionMask, EdgeSmoothing(n), OVK_MIRROR)
          OcclusionMask%values = OcclusionMask%values .and. BaseOcclusionMask%values
        end if
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_OCCLUDES_NONE) then
            Grid_m => Domain%grid(IndexToID(m))
            PaddingMasks(m,n)%values = PaddingMasks(m,n)%values .and. .not. OcclusionMask%values
            if (Logger%log_status) then
              PointCount = ovkCountMask(PaddingMasks(m,n))
              if (PointCount > 0_lk) then
                write (Logger%status_file, '(7a)') "* ", trim(LargeIntToString( &
                  PointCount)), " points on grid ", trim(IntToString(Grid_n%id)), &
                  " marked as not occluded by grid ", trim(IntToString(Grid_m%id)), &
                  " due to padding/smoothing."
              end if
            end if
          end if
        end do
      end if
    end do

    do n = 1, NumGrids
      if (any(Occludes(:,n) /= OVK_OCCLUDES_NONE)) then
        Grid_n => Domain%grid(IndexToID(n))
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_OCCLUDES_NONE) then
            PairwiseOcclusionMasks(m,n)%values = PairwiseOcclusionMasks(m,n)%values .and. .not. &
              PaddingMasks(m,n)%values
          end if
        end do
      end if
    end do

    deallocate(PaddingMasks)

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished applying edge padding and smoothing."
    end if

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Accumulating occlusion..."
    end if

    do n = 1, NumGrids
      if (any(Occludes(:,n) /= OVK_OCCLUDES_NONE)) then
        Grid_n => Domain%grid(IndexToID(n))
        OcclusionMask = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_OCCLUDES_NONE) then
            OcclusionMask%values = OcclusionMask%values .or. PairwiseOcclusionMasks(m,n)%values
          end if
        end do
        call ovkEditGridState(Grid_n, State)
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              if (OcclusionMask%values(i,j,k)) then
                State%values(i,j,k) = ior(State%values(i,j,k),OVK_STATE_OCCLUDED)
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid_n, State)
        if (Logger%log_status) then
          PointCount = ovkCountMask(OcclusionMask)
          if (PointCount > 0_lk) then
            write (Logger%status_file, '(5a)') "* ", trim(LargeIntToString(PointCount)), &
              " occluded points on grid ", trim(IntToString(Grid_n%id)), "."
          end if
        end if
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished accumulating occlusion."
    end if

  end subroutine DetectOccludedPoints

  subroutine FindCoarsePoints(OverlappedGrid, Overlap, OverlapVolumes, CoarseMask)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_real), intent(in) :: OverlapVolumes
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
            CoarseMask%values(i,j,k) = OverlappedGrid%volumes%values(i,j,k) > &
              (1._rk+TOLERANCE) * OverlapVolumes%values(l)
            l = l + 1_lk
          end if
        end do
      end do
    end do

  end subroutine FindCoarsePoints

  subroutine ApplyOverlapMinimization(Domain, ReducedDomainInfo, PairwiseOcclusionMasks)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_logical), dimension(:,:), intent(in) :: PairwiseOcclusionMasks

    integer :: i, j, k, m, n
    type(t_logger) :: Logger
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:,:), pointer :: Occludes
    integer, dimension(:), pointer :: FringeSize
    logical, dimension(:,:), pointer :: MinimizeOverlap
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap
    type(ovk_field_logical), dimension(:), allocatable :: OverlapMinimizationMasks
    type(ovk_field_logical), dimension(:), allocatable :: InnerFringeMasks
    type(ovk_field_logical) :: RemovableMask
    type(ovk_field_logical) :: OcclusionMask
    integer(lk) :: NumRemoved
    logical, dimension(:), allocatable :: UpdateGrid
    type(ovk_field_int), pointer :: State

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Minimizing overlap..."
    end if

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    Occludes => ReducedDomainInfo%occludes
    FringeSize => ReducedDomainInfo%fringe_size
    MinimizeOverlap => ReducedDomainInfo%minimize_overlap

    allocate(OverlapMinimizationMasks(NumGrids))
    allocate(InnerFringeMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(MinimizeOverlap(:,n))) then
        OverlapMinimizationMasks(n) = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_OCCLUDES_NONE .and. MinimizeOverlap(m,n)) then
            OverlapMinimizationMasks(n)%values = OverlapMinimizationMasks(n)%values .or. &
              PairwiseOcclusionMasks(m,n)%values
          end if
        end do
        if (FringeSize(n) > 0) then
          call ovkFilterGridState(Grid_n, OVK_STATE_OCCLUDED, OVK_ANY, OcclusionMask)
          RemovableMask = ovk_field_logical_(Grid_n%cart)
          RemovableMask%values = OcclusionMask%values .or. .not. Grid_n%mask%values
          call ovkErode(RemovableMask, FringeSize(n), OVK_TRUE)
          RemovableMask%values = RemovableMask%values .and. Grid_n%mask%values
          OverlapMinimizationMasks(n)%values = OverlapMinimizationMasks(n)%values .and. &
            RemovableMask%values
          InnerFringeMasks(n) = OverlapMinimizationMasks(n)
          call ovkDilate(InnerFringeMasks(n), FringeSize(n), OVK_FALSE)
          InnerFringeMasks(n)%values = InnerFringeMasks(n)%values .and. (Grid_n%mask%values .and. &
            .not. OverlapMinimizationMasks(n)%values)
        else
          OverlapMinimizationMasks(n)%values = OverlapMinimizationMasks(n)%values .and. &
            Grid_n%mask%values
          InnerFringeMasks(n) = ovk_field_logical_(Grid_n%cart, .false.)
        end if
        if (Logger%log_status) then
          NumRemoved = ovkCountMask(OverlapMinimizationMasks(n))
          if (NumRemoved > 0_lk) then
            write (Logger%status_file, '(5a)') "* ", trim(LargeIntToString(NumRemoved)), &
              " points removed from grid ", trim(IntToString(Grid_n%id)), "."
          end if
        end if
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished minimizing overlap."
    end if

    if (Logger%log_status) then
      write (Logger%status_file, '(2a)') "Updating grids and overlap information after ", &
        "overlap minimization..."
    end if

    allocate(UpdateGrid(NumGrids))

    do n = 1, NumGrids
      UpdateGrid(n) = .false.
      if (any(MinimizeOverlap(:,n))) then
        if (any(OverlapMinimizationMasks(n)%values)) then
          UpdateGrid(n) = .true.
        end if
      end if
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (UpdateGrid(n)) then
        call ovkEditGridState(Grid_n, State)
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              if (OverlapMinimizationMasks(n)%values(i,j,k)) then
                State%values(i,j,k) = iand(State%values(i,j,k), not(OVK_STATE_GRID))
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_EXTERIOR, &
                  OVK_STATE_OVERLAP_MINIMIZED))
              else if (InnerFringeMasks(n)%values(i,j,k)) then
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_FRINGE, &
                  OVK_STATE_INNER_FRINGE))
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid_n, State)
      end if
      if (Logger%log_status) then
        write (Logger%status_file, '(3a)') "* Done updating grid ", trim(IntToString(Grid_n%id)), "."
      end if
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (UpdateGrid(m) .or. UpdateGrid(n)) then
          Grid_m => Domain%grid(IndexToID(m))
          Grid_n => Domain%grid(IndexToID(n))
          Overlap => Domain%overlap(Grid_m%id,Grid_n%id)
          if (ovkOverlapExists(Overlap)) then
            call UpdateOverlapAfterCut(Overlap)
          end if
        end if
      end do
      if (Logger%log_status) then
        write (Logger%status_file, '(3a)') "* Done updating overlap information on grid ", &
          trim(IntToString(Grid_n%id)), "."
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished updating grids and overlap information."
    end if

  end subroutine ApplyOverlapMinimization

  subroutine LocateReceivers(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, n
    type(t_logger) :: Logger
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    type(ovk_grid), pointer :: Grid
    type(ovk_field_logical) :: FringeMask
    type(ovk_field_logical) :: OcclusionMask
    type(ovk_field_logical) :: ReceiverMask
    type(ovk_field_int), pointer :: State
    integer(lk) :: NumReceivers

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Locating receiver points..."
    end if

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id

    do n = 1, NumGrids
      Grid => Domain%grid(IndexToID(n))
      call ovkFilterGridState(Grid, ior(OVK_STATE_GRID, OVK_STATE_FRINGE), OVK_ALL, FringeMask)
      call ovkFilterGridState(Grid, ior(OVK_STATE_GRID, OVK_STATE_OCCLUDED), OVK_ALL, OcclusionMask)
      ReceiverMask = ovk_field_logical_(Grid%cart)
      ReceiverMask%values = FringeMask%values .or. OcclusionMask%values
      call ovkEditGridState(Grid, State)
      do k = Grid%cart%is(3), Grid%cart%ie(3)
        do j = Grid%cart%is(2), Grid%cart%ie(2)
          do i = Grid%cart%is(1), Grid%cart%ie(1)
            if (ReceiverMask%values(i,j,k)) then
              State%values(i,j,k) = ior(State%values(i,j,k), OVK_STATE_RECEIVER)
            end if
          end do
        end do
      end do
      call ovkReleaseGridState(Grid, State)
      if (Logger%log_status) then
        NumReceivers = ovkCountMask(ReceiverMask)
        write (Logger%status_file, '(5a)') "* ", trim(LargeIntToString(NumReceivers)), &
          " receiver points on grid ", trim(IntToString(Grid%id)), "."
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished locating receiver points."
    end if

  end subroutine LocateReceivers

  subroutine ChooseDonors(Domain, ReducedDomainInfo, OverlapVolumes, DonorGridIDs)

    type(ovk_domain), intent(in) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), intent(in) :: OverlapVolumes
    type(ovk_field_int), dimension(:), intent(out) :: DonorGridIDs

    integer :: i, j, k, m, n
    type(t_logger) :: Logger
    integer(lk) :: l
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:,:), pointer :: ConnectionType
    integer, dimension(:,:), pointer :: EdgePadding
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap
    type(ovk_field_logical), dimension(:), allocatable :: ReceiverMasks
    type(ovk_field_int), dimension(:), allocatable :: ReceiverDistances
    type(ovk_field_real) :: DonorDistance
    type(ovk_field_real) :: DonorVolume
    type(ovk_array_int) :: CellReceiverDistance
    type(ovk_field_int) :: OverlapReceiverDistance
    real(rk) :: Distance
    real(rk) :: Volume
    logical :: BetterDonor
    type(ovk_field_logical) :: OrphanMask
    type(ovk_field_int), pointer :: State
    integer, dimension(MAX_ND) :: Point
    integer :: NumWarnings
    integer(lk) :: NumOrphans

    real(rk), parameter :: TOLERANCE = 1.e-12_rk

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Choosing donors..."
    end if

    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    EdgePadding => ReducedDomainInfo%edge_padding
    ConnectionType => ReducedDomainInfo%connection_type

    allocate(ReceiverMasks(NumGrids))
    allocate(ReceiverDistances(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      call ovkFilterGridState(Grid_n, OVK_STATE_RECEIVER, OVK_ALL, ReceiverMasks(n))
      call ovkDistanceField(ReceiverMasks(n), OVK_MIRROR, ReceiverDistances(n))
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      DonorGridIDs(n) = ovk_field_int_(Grid_n%cart, 0)
      DonorDistance = ovk_field_real_(Grid_n%cart)
      DonorVolume = ovk_field_real_(Grid_n%cart)
      do m = 1, NumGrids
        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE) then
          Grid_m => Domain%grid(IndexToID(m))
          Overlap => Domain%overlap(Grid_m%id,Grid_n%id)
          call ovkOverlapCollect(Overlap, OVK_COLLECT_MIN, ReceiverDistances(m), CellReceiverDistance)
          OverlapReceiverDistance = ovk_field_int_(Grid_n%cart, -1)
          call ovkOverlapDisperse(Overlap, OVK_DISPERSE_OVERWRITE, CellReceiverDistance, &
            OverlapReceiverDistance)
          l = 1_lk
          do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
            do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
              do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
                if (ReceiverMasks(n)%values(i,j,k) .and. Overlap%mask%values(i,j,k)) then
                  if (DonorGridIDs(n)%values(i,j,k) == 0) then
                    DonorGridIDs(n)%values(i,j,k) = Grid_m%id
                    DonorDistance%values(i,j,k) = min(real(OverlapReceiverDistance%values(i,j,k), &
                      kind=rk)/real(max(EdgePadding(m,n),1),kind=rk),1._rk)
                    DonorVolume%values(i,j,k) = OverlapVolumes(m,n)%values(l)
                  else
                    Distance = min(real(OverlapReceiverDistance%values(i,j,k),kind=rk)/ &
                      real(max(EdgePadding(m,n),1),kind=rk),1._rk)
                    Volume = OverlapVolumes(m,n)%values(l)
                    if (abs(Distance-DonorDistance%values(i,j,k)) > TOLERANCE) then
                      BetterDonor = Distance > DonorDistance%values(i,j,k)
                    else
                      BetterDonor = Volume < DonorVolume%values(i,j,k)
                    end if
                    if (BetterDonor) then
                      DonorGridIDs(n)%values(i,j,k) = Grid_m%id
                      DonorDistance%values(i,j,k) = Distance
                      DonorVolume%values(i,j,k) = Volume
                    end if
                  end if
                end if
                if (Overlap%mask%values(i,j,k)) then
                  l = l + 1_lk
                end if
              end do
            end do
          end do
        end if
      end do
      OrphanMask = ovk_field_logical_(Grid_n%cart)
      OrphanMask%values = ReceiverMasks(n)%values .and. DonorGridIDs(n)%values == 0
      call ovkEditGridState(Grid_n, State)
      do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
        do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
          do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
            if (OrphanMask%values(i,j,k)) then
              State%values(i,j,k) = iand(State%values(i,j,k), not(OVK_STATE_RECEIVER))
              State%values(i,j,k) = ior(State%values(i,j,k), OVK_STATE_ORPHAN)
            end if
          end do
        end do
      end do
      call ovkReleaseGridState(Grid_n, State)
      if (Logger%log_errors) then
        NumWarnings = 0
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              Point = [i,j,k]
              if (OrphanMask%values(i,j,k)) then
                if (NumWarnings <= 100) then
                  write (Logger%error_file, '(6a)') "WARNING: Could not find suitable ", &
                    "donor for point ", trim(TupleToString(Point(:Grid_n%nd))), " of grid ", &
                    trim(IntToString(Grid_n%id)), "."
                  if (NumWarnings == 100) then
                    write (Logger%error_file, '(a)') "WARNING: Further warnings suppressed."
                  end if
                  NumWarnings = NumWarnings + 1
                end if
              end if
            end do
          end do
        end do
        NumOrphans = ovkCountMask(OrphanMask)
        write (Logger%status_file, '(5a)') "* Done choosing donors on grid ", &
          trim(IntToString(Grid_n%id)), " (", trim(LargeIntToString(NumOrphans)), " orphans found)."
      end if
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished choosing donors."
    end if

  end subroutine ChooseDonors

  subroutine GenerateConnectivity(Domain, ReducedDomainInfo, DonorGridIDs)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_int), dimension(:), intent(in) :: DonorGridIDs

    type t_donor_grid_info
      type(t_donor_stencil), dimension(:), allocatable :: stencil
      integer, dimension(:), allocatable :: stencil_index
    end type t_donor_grid_info

    integer :: i, j, m, n
    type(t_logger) :: Logger
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:,:), pointer :: ConnectionType
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap
    type(ovk_connectivity), pointer :: Connectivity
    type(t_donor_grid_info), dimension(:), allocatable, target :: DonorGridInfo
    integer :: MinConnectionType, MaxConnectionType
    logical, dimension(:), allocatable :: HasConnectionType
    integer :: NumConnectionTypes
    type(ovk_field_logical) :: ReceiverMask
    type(t_donor_stencil), pointer :: DonorStencil
    integer(lk) :: NumConnections

    Logger = Domain%logger

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Generating connectivity information..."
    end if

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    ConnectionType => ReducedDomainInfo%connection_type

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = 1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE) then
          call ovkCreateConnectivity(Domain, Grid_m%id, Grid_n%id)
        end if
      end do
    end do

    allocate(DonorGridInfo(NumGrids))

    do m = 1, NumGrids
      Grid_m => Domain%grid(IndexToID(m))
      MinConnectionType = huge(0)
      MaxConnectionType = -huge(0)
      do n = 1, NumGrids
        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE) then
          MinConnectionType = min(MinConnectionType, ConnectionType(m,n))
          MaxConnectionType = max(MaxConnectionType, ConnectionType(m,n))
        end if
      end do
      if (MaxConnectionType < MinConnectionType) cycle
      allocate(HasConnectionType(MinConnectionType:MaxConnectionType))
      HasConnectionType = .false.
      do n = 1, NumGrids
        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE) then
          HasConnectionType(ConnectionType(m,n)) = .true.
        end if
      end do
      NumConnectionTypes = count(HasConnectionType)
      allocate(DonorGridInfo(m)%stencil(NumConnectionTypes))
      allocate(DonorGridInfo(m)%stencil_index(MinConnectionType:MaxConnectionType))
      j = 1
      DonorGridInfo(m)%stencil_index = 0
      do i = MinConnectionType, MaxConnectionType
        if (HasConnectionType(i)) then
          call CreateDonorStencil(DonorGridInfo(m)%stencil(j), Grid_m, i)
          DonorGridInfo(m)%stencil_index(i) = j
          j = j + 1
        end if
      end do
      deallocate(HasConnectionType)
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = 1, NumGrids
        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE) then
          Grid_m => Domain%grid(IndexToID(m))
          Overlap => Domain%overlap(Grid_m%id,Grid_n%id)
          Connectivity => Domain%connectivity(Grid_m%id,Grid_n%id)
          ReceiverMask = ovk_field_logical_(Grid_n%cart)
          ReceiverMask%values = DonorGridIDs(n)%values == Grid_m%id
          DonorStencil => DonorGridInfo(m)%stencil(DonorGridInfo(m)%stencil_index(ConnectionType(m,n)))
          call FillConnectivity(Connectivity, Overlap, DonorStencil, ReceiverMask)
          if (Logger%log_status) then
            call ovkGetConnectivityCount(Connectivity, NumConnections)
            if (NumConnections > 0_lk) then
              write (Logger%status_file, '(7a)') "* ", trim(LargeIntToString( &
                NumConnections)), " donor/receiver pairs between grid ", &
                trim(IntToString(Grid_m%id)), " and grid ", trim(IntToString(Grid_n%id)), "."
            end if
          end if
        end if
      end do
    end do

    if (Logger%log_status) then
      write (Logger%status_file, '(a)') "Finished generating connectivity information."
    end if

  end subroutine GenerateConnectivity

  subroutine FinalizeAssembly(Domain, ReducedDomainInfo, AssemblyOptions)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(inout) :: ReducedDomainInfo
    type(ovk_assembly_options), intent(in) :: AssemblyOptions

    deallocate(ReducedDomainInfo%index_to_id)
    deallocate(ReducedDomainInfo%overlappable)
    deallocate(ReducedDomainInfo%overlap_tolerance)
    deallocate(ReducedDomainInfo%overlap_accel_quality_adjust)
    deallocate(ReducedDomainInfo%infer_boundaries)
    deallocate(ReducedDomainInfo%cut_boundary_holes)
    deallocate(ReducedDomainInfo%occludes)
    deallocate(ReducedDomainInfo%edge_padding)
    deallocate(ReducedDomainInfo%edge_smoothing)
    deallocate(ReducedDomainInfo%connection_type)
    deallocate(ReducedDomainInfo%fringe_size)
    deallocate(ReducedDomainInfo%minimize_overlap)

    call ResetDomainEdits(Domain)

    Domain%cached_assembly_options = AssemblyOptions

  end subroutine FinalizeAssembly

#if false
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

    NumDims = Domain%nd
    NumGrids = Domain%ngrids

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
#endif

end module ovkAssembly
