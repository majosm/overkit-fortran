! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkAssembly

  use ovkArray
  use ovkAssemblyOptions
  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDomain
  use ovkOverlap
  use ovkOverlapAccel
  use ovkField
  use ovkFieldOps
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

    allocate(OverlapResolutions(ReducedDomainInfo%ngrids,ReducedDomainInfo%ngrids))

    call ComputeOverlapResolutions(Domain, ReducedDomainInfo, OverlapResolutions)

    call LocateOuterFringe(Domain, ReducedDomainInfo)

    call DetectOccludedPoints(Domain, ReducedDomainInfo, OverlapResolutions)

    allocate(DonorGridIDs(ReducedDomainInfo%ngrids))

    call ChooseDonors(Domain, ReducedDomainInfo, OverlapResolutions, DonorGridIDs)

    deallocate(OverlapResolutions)

    call ApplyOverlapMinimization(Domain, ReducedDomainInfo, DonorGridIDs)

    call LocateReceivers(Domain, ReducedDomainInfo, DonorGridIDs)

    call FillConnectivity(Domain, ReducedDomainInfo, DonorGridIDs)

    call FinalizeAssembly(Domain, ReducedDomainInfo, AssemblyOptions)

    if (Domain%logger%verbose) then
      call system_clock(ClockFinal, ClockRate)
      write (*, '(a,f0.3,a)') "Overset grid assembly finished (time: ", &
        real(ClockFinal-ClockInitial,kind=rk)/real(ClockRate,kind=rk), " seconds)."
    end if

  end subroutine ovkAssemble

  subroutine InitAssembly(Domain, AssemblyOptions, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_assembly_options), intent(in) :: AssemblyOptions
    type(t_reduced_domain_info), intent(out) :: ReducedDomainInfo

    integer :: m, n, p, q
    type(ovk_grid), pointer :: Grid_m, Grid_n
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
    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = 1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        if (ovkHasOverlap(Domain, Grid_m%id, Grid_n%id)) then
          call DestroyOverlap(Domain%overlap(Grid_m%id,Grid_n%id))
        end if
        if (ovkHasConnectivity(Domain, Grid_m%id, Grid_n%id)) then
          call DestroyConnectivity(Domain%connectivity(Grid_m%id,Grid_n%id))
        end if
      end do
      call ovkResetGridState(Grid_n)
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
          call CreateOverlap(Domain%overlap(Grid_m%id,Grid_n%id), &
            Grid_m%id, Grid_n%id, Domain%logger, Grid_n%cart)
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
        if (Domain%logger%verbose) then
          write (*, '(3a)') "* Generating overlap search accelerator on grid ", &
            trim(IntToString(Grid_m%id)), "..."
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
            trim(IntToString(Grid_m%id)), "."
        end if
        do n = 1, NumGrids
          Grid_n => Domain%grid(IndexToID(n))
          if (Overlappable(m,n)) then
            Overlap => Domain%overlap(Grid_m%id, Grid_n%id)
            call DetectOverlap(Grid_m, Grid_n, OverlapAccel, Bounds(m,n), OverlapTolerance(m,n), &
              Overlap)
            if (Domain%logger%verbose) then
              if (Overlap%noverlap > 0_lk) then
                write (*, '(7a)') "* Detected ", trim(LargeIntToString( &
                  Overlap%noverlap)), " points on grid ", &
                  trim(IntToString(Grid_n%id)), " overlapped by grid ", &
                  trim(IntToString(Grid_m%id)), "."
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
        if (Domain%logger%verbose) then
          NumBoundaryPoints = ovkCountMask(InferredBoundaryMask)
          if (NumBoundaryPoints > 0_lk) then
            write (*, '(5a)') "* ", trim(LargeIntToString(NumBoundaryPoints)), &
              " points marked as boundaries on grid ", trim(IntToString(Grid%id)), "."
          end if
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished inferring domain boundaries in non-overlapping regions."
    end if

  end subroutine InferNonOverlappingBoundaries

  subroutine CutHoles(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: CutBoundaryHoles
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

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    CutBoundaryHoles => ReducedDomainInfo%cut_boundary_holes

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
                " points removed from grid ", trim(IntToString(Grid_n%id)), "."
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
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_HOLE, &
                  OVK_STATE_BOUNDARY_HOLE))
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid_n, State)
      end if
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done updating grid ", trim(IntToString(Grid_n%id)), "."
      end if
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (UpdateGrid(m) .or. UpdateGrid(n)) then
          Grid_m => Domain%grid(IndexToID(m))
          Grid_n => Domain%grid(IndexToID(n))
          Overlap_mn => Domain%overlap(Grid_m%id,Grid_n%id)
          if (ovkOverlapExists(Overlap_mn)) then
            call UpdateOverlapAfterCut(Grid_m, Grid_n, Overlap_mn)
          end if
        end if
      end do
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done updating overlap information on grid ", &
          trim(IntToString(Grid_n%id)), "."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished updating grids and overlap information."
    end if

  end subroutine CutHoles

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

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = 1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        Overlap => Domain%overlap(Grid_m%id,Grid_n%id)
        if (ovkOverlapExists(Overlap)) then
          call ovkOverlapCollect(Grid_m, Overlap, OVK_COLLECT_INTERPOLATE, Grid_m%resolution, &
            OverlapResolutions(m,n))
        end if
      end do
    end do

  end subroutine ComputeOverlapResolutions

  subroutine LocateOuterFringe(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, n
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

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    FringeSize => ReducedDomainInfo%fringe_size

    if (Domain%logger%verbose) then
      write (*, '(a)') "Locating outer fringe points..."
    end if

    do n = 1, NumGrids
      if (FringeSize(n) > 0) then
        Grid => Domain%grid(IndexToID(n))
        OuterFringeMask = ovk_field_logical_(Grid%cart)
        BoundaryMask = ovk_field_logical_(Grid%cart)
        BoundaryMask%values = Grid%boundary_mask%values .or. Grid%internal_boundary_mask%values
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
        if (Domain%logger%verbose) then
          NumFringe = ovkCountMask(OuterFringeMask)
          if (NumFringe > 0_lk) then
            write (*, '(5a)') "* ", trim(LargeIntToString(NumFringe)), &
              " outer fringe points on grid ", trim(IntToString(Grid%id)), "."
          end if
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished locating outer fringe points."
    end if

  end subroutine LocateOuterFringe

  subroutine DetectOccludedPoints(Domain, ReducedDomainInfo, OverlapResolutions)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), intent(in) :: OverlapResolutions

    integer :: i, j, k, m, n
    integer(lk) :: l
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: Overlappable
    integer, dimension(:,:), pointer :: Occludes
    integer, dimension(:,:), pointer :: EdgePadding
    integer, dimension(:), pointer :: EdgeSmoothing
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap_mn, Overlap_nm
    type(ovk_field_logical), dimension(:,:), allocatable :: PairwiseOcclusionMasks
    type(ovk_field_logical) :: OverlappedMask_m, OverlappedMask_n
    integer(lk) :: PointCount
    type(ovk_field_logical), dimension(:), allocatable :: OcclusionMasks
    type(ovk_field_int), dimension(:), allocatable :: FinestOverlappingGrid
    type(ovk_field_logical) :: OuterFringeMask
    type(ovk_field_real) :: BestResolutions
    type(ovk_field_int) :: OverlapEdgeDistance
    type(ovk_field_logical) :: EdgeMask
    type(ovk_field_int) :: EdgeDistance
    type(ovk_array_int) :: CellEdgeDistance
    type(ovk_field_logical) :: PaddingMask
    type(ovk_field_logical) :: OverlapMask
    type(ovk_field_int), pointer :: State

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    Overlappable => ReducedDomainInfo%overlappable
    Occludes => ReducedDomainInfo%occludes
    EdgePadding => ReducedDomainInfo%edge_padding
    EdgeSmoothing => ReducedDomainInfo%edge_smoothing

    if (Domain%logger%verbose) then
      write (*, '(a)') "Detecting pairwise occlusion..."
    end if

    allocate(PairwiseOcclusionMasks(NumGrids,NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = n+1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        Overlap_mn => Domain%overlap(Grid_m%id,Grid_n%id)
        Overlap_nm => Domain%overlap(Grid_n%id,Grid_m%id)
        select case (Occludes(m,n))
        case (OVK_OCCLUDES_COARSE)
          call FindCoarsePoints(Grid_n, Overlap_mn, OverlapResolutions(m,n), &
            PairwiseOcclusionMasks(m,n))
        case (OVK_OCCLUDES_ALL)
          PairwiseOcclusionMasks(m,n) = Overlap_mn%mask
        case (OVK_OCCLUDES_NONE)
          PairwiseOcclusionMasks(m,n) = ovk_field_logical_(NumDims)
        end select
        select case (Occludes(n,m))
        case (OVK_OCCLUDES_COARSE)
          call FindCoarsePoints(Grid_m, Overlap_nm, OverlapResolutions(n,m), &
            PairwiseOcclusionMasks(n,m))
        case (OVK_OCCLUDES_ALL)
          PairwiseOcclusionMasks(n,m) = Overlap_nm%mask
        case (OVK_OCCLUDES_NONE)
          PairwiseOcclusionMasks(n,m) = ovk_field_logical_(NumDims)
        end select
        if (Occludes(m,n) == OVK_OCCLUDES_COARSE .and. Occludes(n,m) == OVK_OCCLUDES_COARSE) then
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
              trim(IntToString(Grid_n%id)), " occluded by grid ", &
              trim(IntToString(Grid_m%id)), "."
          end if
          PointCount = ovkCountMask(PairwiseOcclusionMasks(n,m))
          if (PointCount > 0_lk) then
            write (*, '(7a)') "* ", trim(LargeIntToString(PointCount)), " points on grid ", &
              trim(IntToString(Grid_m%id)), " occluded by grid ", &
              trim(IntToString(Grid_n%id)), "."
          end if
        end if
      end do
    end do

    allocate(OcclusionMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      OcclusionMasks(n) = ovk_field_logical_(Grid_n%cart)
      if (any(Occludes(:,n) /= OVK_OCCLUDES_NONE)) then
        OcclusionMasks(n)%values = .not. Grid_n%mask%values
        do m = 1, NumGrids
          if (Occludes(m,n) /= OVK_OCCLUDES_NONE) then
            OcclusionMasks(n)%values = OcclusionMasks(n)%values .or. &
              PairwiseOcclusionMasks(m,n)%values
          end if
        end do
      else
        OcclusionMasks(n)%values = .false.
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished detecting pairwise occlusion."
    end if

    if (Domain%logger%verbose) then
      write (*, '(a)') "Applying edge padding..."
    end if

    allocate(FinestOverlappingGrid(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      call ovkFilterGridState(Grid_n, OVK_STATE_OUTER_FRINGE, OVK_ANY, OuterFringeMask)
      FinestOverlappingGrid(n) = ovk_field_int_(Grid_n%cart, -1)
      BestResolutions = ovk_field_real_(Grid_n%cart, 0._rk)
      do m = 1, NumGrids
        if (Overlappable(m,n)) then
          Grid_m => Domain%grid(IndexToID(m))
          Overlap_mn => Domain%overlap(Grid_m%id,Grid_n%id)
          l = 1_lk
          do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
            do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
              do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
                if (Overlap_mn%mask%values(i,j,k)) then
                  if ((OcclusionMasks(n)%values(i,j,k) .or. OuterFringeMask%values(i,j,k)) .and. &
                    OverlapResolutions(m,n)%values(l) > BestResolutions%values(i,j,k)) then
                    BestResolutions%values(i,j,k) = OverlapResolutions(m,n)%values(l)
                    FinestOverlappingGrid(n)%values(i,j,k) = Grid_m%id
                  end if
                  l = l + 1_lk
                end if
              end do
            end do
          end do
        end if
      end do
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(Occludes(:,n) /= OVK_OCCLUDES_NONE)) then
        OcclusionMasks(n)%values = .false.
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
            call ovkOverlapCollect(Grid_m, Overlap_mn, OVK_COLLECT_MIN, EdgeDistance, &
              CellEdgeDistance)
            OverlapEdgeDistance = ovk_field_int_(Grid_n%cart, -1)
            call ovkOverlapDisperse(Grid_n, Overlap_mn, OVK_DISPERSE_OVERWRITE, CellEdgeDistance, &
              OverlapEdgeDistance)
            PaddingMask = ovk_field_logical_(Grid_n%cart, .false.)
            do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
              do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
                do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
                  if (PairwiseOcclusionMasks(m,n)%values(i,j,k) .and. &
                    OverlapEdgeDistance%values(i,j,k) <= EdgePadding(m,n)) then
                    PaddingMask%values(i,j,k) = .true.
                  end if
                end do
              end do
            end do
            OcclusionMasks(n)%values = OcclusionMasks(n)%values .or. &
              (PairwiseOcclusionMasks(m,n)%values .and. .not. PaddingMask%values)
            if (Domain%logger%verbose) then
              PointCount = ovkCountMask(PaddingMask)
              if (PointCount > 0_lk) then
                write (*, '(7a)') "* ", trim(LargeIntToString(PointCount)), " points on grid ", &
                  trim(IntToString(Grid_n%id)), " marked as not occluded by grid ", &
                  trim(IntToString(Grid_m%id)), " due to padding."
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
      write (*, '(a)') "Applying edge smoothing..."
    end if

    do n = 1, NumGrids
      if (EdgeSmoothing(n) > 0) then
        Grid_n => Domain%grid(IndexToID(n))
        if (any(Occludes(:,n) /= OVK_OCCLUDES_NONE)) then
          OverlapMask = ovk_field_logical_(Grid_n%cart, .false.)
          do m = 1, NumGrids
            if (Overlappable(m,n)) then
              Grid_m => Domain%grid(IndexToID(m))
              Overlap_mn => Domain%overlap(Grid_m%id,Grid_n%id)
              OverlapMask%values = OverlapMask%values .or. Overlap_mn%mask%values
            end if
          end do
          call ovkDilate(OcclusionMasks(n), EdgeSmoothing(n), OVK_MIRROR)
          call ovkErode(OcclusionMasks(n), EdgeSmoothing(n), OVK_MIRROR)
          OcclusionMasks(n)%values = OcclusionMasks(n)%values .and. OverlapMask%values
        end if
        if (Domain%logger%verbose) then
          write (*, '(3a)') "* Done applying edge smoothing to grid ", &
            trim(IntToString(Grid_n%id)), "."
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished applying edge smoothing."
    end if

    if (Domain%logger%verbose) then
      write (*, '(a)') "Accumulating occlusion..."
    end if

    do n = 1, NumGrids
      if (any(Occludes(:,n) /= OVK_OCCLUDES_NONE)) then
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
              " occluded points on grid ", trim(IntToString(Grid_n%id)), "."
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

  subroutine ChooseDonors(Domain, ReducedDomainInfo, OverlapResolutions, DonorGridIDs)

    type(ovk_domain), intent(in) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_array_real), dimension(:,:), intent(in) :: OverlapResolutions
    type(ovk_field_int), dimension(:), intent(out) :: DonorGridIDs

    integer :: i, j, k, m, n
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
    type(ovk_field_real) :: DonorResolution
    type(ovk_array_int) :: CellReceiverDistance
    type(ovk_field_int) :: OverlapReceiverDistance
    real(rk) :: Distance
    real(rk) :: Resolution
    logical :: BetterDonor

    real(rk), parameter :: TOLERANCE = 1.e-12_rk

    if (Domain%logger%verbose) then
      write (*, '(a)') "Choosing donors..."
    end if

    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    EdgePadding => ReducedDomainInfo%edge_padding
    ConnectionType => ReducedDomainInfo%connection_type

    allocate(ReceiverMasks(NumGrids))
    allocate(ReceiverDistances(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      call ovkFilterGridState(Grid_n, ior(OVK_STATE_OUTER_FRINGE, OVK_STATE_OCCLUDED), OVK_ANY, &
        ReceiverMasks(n))
      call ovkDistanceField(ReceiverMasks(n), OVK_MIRROR, ReceiverDistances(n))
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      DonorGridIDs(n) = ovk_field_int_(Grid_n%cart, 0)
      DonorDistance = ovk_field_real_(Grid_n%cart)
      DonorResolution = ovk_field_real_(Grid_n%cart)
      do m = 1, NumGrids
        if (ConnectionType(m,n) /= OVK_CONNECTION_NONE) then
          Grid_m => Domain%grid(IndexToID(m))
          Overlap => Domain%overlap(Grid_m%id,Grid_n%id)
          call ovkOverlapCollect(Grid_m, Overlap, OVK_COLLECT_MIN, ReceiverDistances(m), &
            CellReceiverDistance)
          OverlapReceiverDistance = ovk_field_int_(Grid_n%cart, -1)
          call ovkOverlapDisperse(Grid_n, Overlap, OVK_DISPERSE_OVERWRITE, CellReceiverDistance, &
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
                    DonorResolution%values(i,j,k) = OverlapResolutions(m,n)%values(l)
                  else
                    Distance = min(real(OverlapReceiverDistance%values(i,j,k),kind=rk)/ &
                      real(max(EdgePadding(m,n),1),kind=rk),1._rk)
                    Resolution = OverlapResolutions(m,n)%values(l)
                    if (abs(Distance-DonorDistance%values(i,j,k)) > TOLERANCE) then
                      BetterDonor = Distance > DonorDistance%values(i,j,k)
                    else
                      BetterDonor = Resolution > DonorResolution%values(i,j,k)
                    end if
                    if (BetterDonor) then
                      DonorGridIDs(n)%values(i,j,k) = Grid_m%id
                      DonorDistance%values(i,j,k) = Distance
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
        end if
      end do
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done choosing donors on grid ", &
          trim(IntToString(Grid_n%id)), "."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished choosing donors."
    end if

  end subroutine ChooseDonors

  subroutine ApplyOverlapMinimization(Domain, ReducedDomainInfo, DonorGridIDs)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_int), dimension(:), intent(in) :: DonorGridIDs

    integer :: i, j, k, m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
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

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    FringeSize => ReducedDomainInfo%fringe_size
    MinimizeOverlap => ReducedDomainInfo%minimize_overlap

    if (Domain%logger%verbose) then
      write (*, '(a)') "Minimizing overlap..."
    end if

    allocate(OverlapMinimizationMasks(NumGrids))
    allocate(InnerFringeMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(MinimizeOverlap(:,n))) then
        OverlapMinimizationMasks(n) = ovk_field_logical_(Grid_n%cart, .false.)
        do m = 1, NumGrids
          if (MinimizeOverlap(m,n)) then
            OverlapMinimizationMasks(n)%values = OverlapMinimizationMasks(n)%values .or. &
              DonorGridIDs(n)%values == IndexToID(m)
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
        if (Domain%logger%verbose) then
          NumRemoved = ovkCountMask(OverlapMinimizationMasks(n))
          if (NumRemoved > 0_lk) then
            write (*, '(5a)') "* ", trim(LargeIntToString(NumRemoved)), &
              " points removed from grid ", trim(IntToString(Grid_n%id)), "."
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
                State%values(i,j,k) = ior(State%values(i,j,k), ior(OVK_STATE_HOLE, &
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
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done updating grid ", trim(IntToString(Grid_n%id)), "."
      end if
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (UpdateGrid(m) .or. UpdateGrid(n)) then
          Grid_m => Domain%grid(IndexToID(m))
          Grid_n => Domain%grid(IndexToID(n))
          Overlap => Domain%overlap(Grid_m%id,Grid_n%id)
          if (ovkOverlapExists(Overlap)) then
            call UpdateOverlapAfterCut(Grid_m, Grid_n, Overlap)
          end if
        end if
      end do
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done updating overlap information on grid ", &
          trim(IntToString(Grid_n%id)), "."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished updating grids and overlap information."
    end if

  end subroutine ApplyOverlapMinimization

  subroutine LocateReceivers(Domain, ReducedDomainInfo, DonorGridIDs)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_int), dimension(:), intent(in) :: DonorGridIDs

    integer :: i, j, k, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    type(ovk_grid), pointer :: Grid
    type(ovk_field_logical) :: FringeMask
    type(ovk_field_logical) :: OcclusionMask
    type(ovk_field_logical) :: ReceiverMask
    type(ovk_field_logical) :: OrphanMask
    type(ovk_field_int), pointer :: State
    integer :: NumWarnings
    integer, dimension(MAX_ND) :: Point
    integer(lk) :: NumReceivers
    integer(lk) :: NumOrphans

    NumDims = Domain%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id

    if (Domain%logger%verbose) then
      write (*, '(a)') "Locating receiver points..."
    end if

    do n = 1, NumGrids
      Grid => Domain%grid(IndexToID(n))
      call ovkFilterGridState(Grid, ior(OVK_STATE_GRID, OVK_STATE_FRINGE), OVK_ALL, FringeMask)
      call ovkFilterGridState(Grid, ior(OVK_STATE_GRID, OVK_STATE_OCCLUDED), OVK_ALL, OcclusionMask)
      ReceiverMask = ovk_field_logical_(Grid%cart)
      ReceiverMask%values = FringeMask%values .or. OcclusionMask%values
      OrphanMask = ovk_field_logical_(Grid%cart)
      OrphanMask%values = ReceiverMask%values .and. DonorGridIDs(n)%values == 0
      ReceiverMask%values = ReceiverMask%values .and. .not. OrphanMask%values
      call ovkEditGridState(Grid, State)
      do k = Grid%cart%is(3), Grid%cart%ie(3)
        do j = Grid%cart%is(2), Grid%cart%ie(2)
          do i = Grid%cart%is(1), Grid%cart%ie(1)
            if (ReceiverMask%values(i,j,k)) then
              State%values(i,j,k) = ior(State%values(i,j,k), OVK_STATE_RECEIVER)
            else if (OrphanMask%values(i,j,k)) then
              State%values(i,j,k) = ior(State%values(i,j,k), OVK_STATE_ORPHAN)
            end if
          end do
        end do
      end do
      call ovkReleaseGridState(Grid, State)
      if (Domain%logger%verbose) then
        NumWarnings = 0
        do k = Grid%cart%is(3), Grid%cart%ie(3)
          do j = Grid%cart%is(2), Grid%cart%ie(2)
            do i = Grid%cart%is(1), Grid%cart%ie(1)
              Point = [i,j,k]
              if (OrphanMask%values(i,j,k)) then
                if (NumWarnings <= 100) then
                  write (ERROR_UNIT, '(4a)') "WARNING: Could not find suitable donor for point ", &
                    trim(TupleToString(Point(:Grid%nd))), " of grid ", &
                    trim(IntToString(Grid%id))
                  if (NumWarnings == 100) then
                    write (ERROR_UNIT, '(a)') "WARNING: Further warnings suppressed."
                  end if
                  NumWarnings = NumWarnings + 1
                end if
              end if
            end do
          end do
        end do
        NumReceivers = ovkCountMask(ReceiverMask)
        NumOrphans = ovkCountMask(OrphanMask)
        write (*, '(7a)') "* ", trim(LargeIntToString(NumReceivers)), &
          " receiver points on grid ", trim(IntToString(Grid%id)), " (", &
          trim(LargeIntToString(NumOrphans)), " orphans)."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished locating receiver points."
    end if

  end subroutine LocateReceivers

  subroutine FillConnectivity(Domain, ReducedDomainInfo, DonorGridIDs)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_int), dimension(:), intent(in) :: DonorGridIDs

    integer :: i, j, k, m, n
    integer(lk) :: l, p
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:,:), pointer :: ConnectionType
    integer(lk), dimension(:), allocatable :: NumConnections
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap
    type(ovk_connectivity), pointer :: Connectivity
    type(ovk_field_logical) :: ReceiverMask
    integer, dimension(MAX_ND) :: DonorSize
    integer :: MaxDonorSize
    integer, dimension(MAX_ND,2) :: DonorExtents
    real(rk), dimension(Domain%nd) :: DonorCoords
    real(rk), dimension(:,:), allocatable :: DonorInterpCoefs
    integer, dimension(MAX_ND) :: ReceiverPoint

    if (Domain%logger%verbose) then
      write (*, '(a)') "Generating connectivity information..."
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
          call CreateConnectivity(Domain%connectivity(Grid_m%id, Grid_n%id), &
            Grid_m%id, Grid_n%id, Domain%logger, NumDims)
        end if
      end do
    end do

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
          Overlap => Domain%overlap(Grid_m%id,Grid_n%id)
          Connectivity => Domain%connectivity(Grid_m%id,Grid_n%id)

          if (NumConnections(m) > 0_lk) then

            DonorSize = 1
            DonorSize(:NumDims) = ovkDonorSize(NumDims, ConnectionType(m,n))
            MaxDonorSize = maxval(DonorSize)

            call ResizeConnectivity(Connectivity, NumConnections(m), MaxDonorSize)

            allocate(DonorInterpCoefs(MaxDonorSize,NumDims))

            l = 1_lk
            p = 1_lk
            DonorExtents = 1
            do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
              do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
                do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
                  if (ReceiverMask%values(i,j,k) .and. DonorGridIDs(n)%values(i,j,k) == &
                    Grid_m%id) then
                    ReceiverPoint = [i,j,k]
                    call ovkFindDonor(Grid_m, Grid_n, Overlap, ReceiverPoint, l, &
                      ConnectionType(m,n), DonorExtents, DonorCoords, DonorInterpCoefs)
                    Connectivity%donor_extents(:,:,p) = DonorExtents
                    Connectivity%donor_coords(:,p) = DonorCoords
                    Connectivity%donor_interp_coefs(:,:,p) = DonorInterpCoefs
                    Connectivity%receiver_points(:,p) = ReceiverPoint
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
                  " donor/receiver pairs between grid ", trim(IntToString(Grid_m%id)), &
                  " and grid ", trim(IntToString(Grid_n%id)), "."
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
