! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkAssembly

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
    logical, dimension(:,:), pointer :: overlap_hole_cutting
    integer, dimension(:,:), pointer :: connection_type
    integer, dimension(:,:), pointer :: interp_scheme
    integer, dimension(:), pointer :: fringe_size
    integer, dimension(:,:), pointer :: interface_padding
  end type t_reduced_domain_info

contains

  subroutine ovkAssemble(Domain)

    type(ovk_domain), intent(inout) :: Domain

    integer :: ClockInitial, ClockFinal, ClockRate
    type(t_reduced_domain_info) :: ReducedDomainInfo

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

    call CutHoles(Domain, ReducedDomainInfo)

    call GenerateConnectivity(Domain, ReducedDomainInfo)

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
    logical, dimension(:,:), pointer :: OverlapHoleCutting
    integer, dimension(:,:), pointer :: ConnectionType
    integer, dimension(:,:), pointer :: InterpScheme
    integer, dimension(:), pointer :: FringeSize
    integer, dimension(:,:), pointer :: InterfacePadding

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
    allocate(ReducedDomainInfo%overlap_hole_cutting(NumGrids,NumGrids))
    InferBoundaries => ReducedDomainInfo%infer_boundaries
    BoundaryHoleCutting => ReducedDomainInfo%boundary_hole_cutting
    OverlapHoleCutting => ReducedDomainInfo%overlap_hole_cutting

    do n = 1, NumGrids
      q = IndexToID(n)
      InferBoundaries(n) = Domain%properties%infer_boundaries(q)
      do m = 1, NumGrids
        p = IndexToID(m)
        BoundaryHoleCutting(m,n) = Domain%properties%boundary_hole_cutting(p,q)
        OverlapHoleCutting(m,n) = Domain%properties%overlap_hole_cutting(p,q)
      end do
    end do

    allocate(ReducedDomainInfo%connection_type(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%interp_scheme(NumGrids,NumGrids))
    allocate(ReducedDomainInfo%fringe_size(NumGrids))
    allocate(ReducedDomainInfo%interface_padding(NumGrids,NumGrids))
    ConnectionType => ReducedDomainInfo%connection_type
    InterpScheme => ReducedDomainInfo%interp_scheme
    FringeSize => ReducedDomainInfo%fringe_size
    InterfacePadding => ReducedDomainInfo%interface_padding

    do n = 1, NumGrids
      q = IndexToID(n)
      FringeSize(n) = Domain%properties%fringe_size(q)
      do m = 1, NumGrids
        p = IndexToID(m)
        ConnectionType(m,n) = Domain%properties%connection_type(p,q)
        InterpScheme(m,n) = Domain%properties%interp_scheme(p,q)
        InterfacePadding(m,n) = Domain%properties%interface_padding(p,q)
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
        do m = 1, NumGrids
          Overlap => Domain%overlap(IndexToID(m),IndexToID(n))
          if (OverlapExists(Overlap)) then
            InferredBoundaryMask%values = InferredBoundaryMask%values .and. .not. &
              Overlap%mask%values
          end if
        end do
        call ovkFilterGridState(Grid, OVK_STATE_DOMAIN_BOUNDARY, OVK_ALL, InitialBoundaryMask)
        InferredBoundaryMask%values = InferredBoundaryMask%values .and. .not. &
          InitialBoundaryMask%values
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

  subroutine CutHoles(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, m, n
    integer :: NumDims
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    logical, dimension(:,:), pointer :: BoundaryHoleCutting
    logical, dimension(:,:), pointer :: OverlapHoleCutting
    integer, dimension(:,:), pointer :: ConnectionType
    integer, dimension(:), pointer :: FringeSize
    integer, dimension(:,:), pointer :: InterfacePadding
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap_mn, Overlap_nm
    type(ovk_field_logical), dimension(:), allocatable :: BoundaryHoleMasks
    type(ovk_field_logical), dimension(:), allocatable :: OverlapHoleMasks
    type(ovk_field_logical) :: EdgeMask1, EdgeMask2
    type(ovk_field_logical) :: BoundaryMask
    type(ovk_field_logical) :: InteriorMask
    type(ovk_field_logical) :: SpuriousBoundaryMask
    integer(lk) :: NumRemoved
    type(ovk_field_logical) :: OverlappingMask
    type(ovk_field_logical) :: CoarseMask_m, CoarseMask_n
    type(ovk_field_logical) :: OverlappedMask_m, OverlappedMask_n
    type(ovk_field_logical), dimension(:), allocatable :: FineMasks
    type(ovk_field_logical), dimension(:,:), allocatable :: PairwiseFineMasks
    type(ovk_field_logical) :: CoarseMask
    type(ovk_field_logical) :: CoarseEdgeMask
    type(ovk_field_logical) :: PairwiseCoarseMask
    type(ovk_field_logical) :: PaddingEdgeMask
    type(ovk_field_logical) :: CutMask
    type(ovk_field_int), pointer :: State

    if (Domain%logger%verbose) then
      write (*, '(a)') "Cutting boundary holes..."
    end if

    NumDims = Domain%properties%nd
    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    BoundaryHoleCutting => ReducedDomainInfo%boundary_hole_cutting
    OverlapHoleCutting => ReducedDomainInfo%overlap_hole_cutting
    ConnectionType => ReducedDomainInfo%connection_type
    FringeSize => ReducedDomainInfo%fringe_size
    InterfacePadding => ReducedDomainInfo%interface_padding

    allocate(BoundaryHoleMasks(NumGrids))
    allocate(OverlapHoleMasks(NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      BoundaryHoleMasks(n) = ovk_field_logical_(Grid_n%cart)
      OverlapHoleMasks(n) = ovk_field_logical_(Grid_n%cart)
    end do

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
        if (any(InteriorMask%values)) then
          BoundaryMask%values = BoundaryMask%values .or. Grid_n%boundary_mask%values
          call ovkFloodFill(InteriorMask, BoundaryMask)
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
      write (*, '(a)') "Cutting overlap holes..."
    end if

    allocate(FineMasks(NumGrids))
    allocate(PairwiseFineMasks(NumGrids,NumGrids))

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(OverlapHoleCutting(:,n))) then
        FineMasks(n) = Grid_n%mask
        do m = 1, NumGrids
          if (OverlapHoleCutting(m,n)) then
            PairwiseFineMasks(m,n) = Grid_n%mask
          else
            PairwiseFineMasks(m,n) = ovk_field_logical_(NumDims)
          end if
        end do
      else
        FineMasks(n) = ovk_field_logical_(NumDims)
        PairwiseFineMasks(:,n) = ovk_field_logical_(NumDims)
      end if
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      do m = n+1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        Overlap_mn => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
        Overlap_nm => Domain%overlap(Grid_n%properties%id,Grid_m%properties%id)
        if (OverlapHoleCutting(m,n) .and. OverlapHoleCutting(n,m)) then
          call FindCoarsePoints(Grid_m, Overlap_nm, CoarseMask_m)
          call FindCoarsePoints(Grid_n, Overlap_mn, CoarseMask_n)
          ! If any coarse points on grid m are overlapped by coarse points on grid n, change
          ! them to fine


          ! This seems weird -- check it



          if (any(CoarseMask_m%values)) then
            call ovkFindOverlappedPoints(Grid_n, Grid_m, Overlap_nm, CoarseMask_n, &
              OverlappedMask_m)
            OverlappedMask_m%values = OverlappedMask_m%values .and. CoarseMask_m%values
            FineMasks(m)%values = FineMasks(m)%values .or. OverlappedMask_m%values
          end if
          if (any(CoarseMask_n%values)) then
            call ovkFindOverlappedPoints(Grid_m, Grid_n, Overlap_mn, CoarseMask_m, &
              OverlappedMask_n)
            OverlappedMask_n%values = OverlappedMask_n%values .and. CoarseMask_n%values
            FineMasks(n)%values = FineMasks(n)%values .or. OverlappedMask_n%values
          end if
          FineMasks(m)%values = FineMasks(m)%values .and. .not. CoarseMask_m%values
          PairwiseFineMasks(n,m)%values = PairwiseFineMasks(n,m)%values .and. .not. &
            CoarseMask_m%values
          FineMasks(n)%values = FineMasks(n)%values .and. .not. CoarseMask_n%values
          PairwiseFineMasks(m,n)%values = PairwiseFineMasks(m,n)%values .and. .not. &
            CoarseMask_n%values




        else if (OverlapHoleCutting(m,n)) then
          ! Behave as if grid m is finer everywhere
          FineMasks(n)%values = FineMasks(n)%values .and. .not. Overlap_mn%mask%values
          PairwiseFineMasks(m,n)%values = PairwiseFineMasks(m,n)%values .and. .not. &
            Overlap_mn%mask%values
        else if (OverlapHoleCutting(n,m)) then
          ! Behave as if grid n is finer everywhere
          FineMasks(m)%values = FineMasks(m)%values .and. .not. Overlap_mn%mask%values
          PairwiseFineMasks(n,m)%values = PairwiseFineMasks(n,m)%values .and. .not. &
            Overlap_nm%mask%values
        end if
      end do
    end do

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      if (any(OverlapHoleCutting(:,n))) then
        CoarseMask = ovk_field_logical_(Grid_n%cart)
        CoarseMask%values = .not. FineMasks(n)%values .or. BoundaryHoleMasks(n)%values
        CutMask = ovk_field_logical_(Grid_n%cart)
        CutMask%values = Grid_n%mask%values .and. CoarseMask%values
        call ovkDetectEdge(CoarseMask, OVK_OUTER_EDGE, OVK_MIRROR, .false., CoarseEdgeMask)
        do m = 1, NumGrids
          Grid_m => Domain%grid(IndexToID(m))
          Overlap_mn => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
          if (OverlapHoleCutting(m,n) .and. ConnectionType(m,n) == OVK_CONNECTION_FRINGE) then
            PairwiseCoarseMask = ovk_field_logical_(Grid_n%cart)
            PairwiseCoarseMask%values = .not. PairwiseFineMasks(m,n)%values
            call ovkDetectEdge(PairwiseCoarseMask, OVK_OUTER_EDGE, OVK_MIRROR, .false., &
              PaddingEdgeMask)
            PaddingEdgeMask%values = PaddingEdgeMask%values .and. CoarseEdgeMask%values
            call ovkDilate(PaddingEdgeMask, FringeSize(n)+InterfacePadding(m,n), OVK_FALSE)
            CutMask%values = CutMask%values .and. .not. &
              (Overlap_mn%mask%values .and. PaddingEdgeMask%values)
          end if
        end do
        CutMask%values = CutMask%values .and. .not. BoundaryHoleMasks(n)%values
        OverlapHoleMasks(n)%values = CutMask%values
        if (Domain%logger%verbose) then
          NumRemoved = ovkCountMask(OverlapHoleMasks(n))
          if (NumRemoved > 0_lk) then
            write (*, '(5a)') "* ", trim(LargeIntToString(NumRemoved)), &
              " points removed from grid ", trim(IntToString(Grid_n%properties%id)), "."
          end if
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished cutting overlap holes."
    end if

    if (Domain%logger%verbose) then
      write (*, '(a)') "Updating grids and overlap information after hole cutting..."
    end if

    do n = 1, NumGrids
      if (any(OverlapHoleCutting(:,n)) .or. any(BoundaryHoleCutting(:,n))) then
        Grid_n => Domain%grid(IndexToID(n))
        call ovkEditGridState(Grid_n, State)
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              if (BoundaryHoleMasks(n)%values(i,j,k) .or. OverlapHoleMasks(n)%values(i,j,k)) then
                State%values(i,j,k) = ior(iand(State%values(i,j,k),not(OVK_STATE_GRID)), &
                  OVK_STATE_HOLE)
                if (BoundaryHoleMasks(n)%values(i,j,k)) then
                  State%values(i,j,k) = ior(State%values(i,j,k),OVK_STATE_BOUNDARY_HOLE)
                end if
                if (OverlapHoleMasks(n)%values(i,j,k)) then
                  State%values(i,j,k) = ior(State%values(i,j,k),OVK_STATE_OVERLAP_HOLE)
                end if
              end if
            end do
          end do
        end do
        call ovkReleaseGridState(Grid_n, State)
        if (Domain%logger%verbose) then
          write (*, '(3a)') "* Done updating grid ", trim(IntToString(Grid_n%properties%id)), "."
        end if
      end if
    end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        Grid_m => Domain%grid(IndexToID(m))
        Grid_n => Domain%grid(IndexToID(n))
        Overlap_mn => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
        if (OverlapExists(Overlap_mn)) then
          call UpdateOverlapAfterCut(Grid_m, Grid_n, Overlap_mn)
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

  end subroutine CutHoles

  subroutine GenerateConnectivity(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(inout) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    type(ovk_field_int), dimension(:), allocatable :: DonorGridIDs

    allocate(DonorGridIDs(ReducedDomainInfo%ngrids))

    call LocateReceivers(Domain, ReducedDomainInfo)
    call ChooseDonors(Domain, ReducedDomainInfo, DonorGridIDs)
    call DetectOrphans(Domain, ReducedDomainInfo, DonorGridIDs)
    call FillConnectivity(Domain, ReducedDomainInfo, DonorGridIDs)

  end subroutine GenerateConnectivity

  subroutine LocateReceivers(Domain, ReducedDomainInfo)

    type(ovk_domain), intent(in) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo

    integer :: i, j, k, n
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:), pointer :: FringeSize
    type(ovk_grid), pointer :: Grid
    type(ovk_field_logical) :: EdgeMask
    type(ovk_field_logical) :: NonSubsetMask
    type(ovk_field_logical) :: NonSubsetEdgeMask
    type(ovk_field_logical) :: SubsetMask
    type(ovk_field_logical) :: SubsetEdgeMask
    type(ovk_field_logical) :: ReceiverMask
    type(ovk_field_int), pointer :: State
    integer(lk) :: NumReceivers

    if (Domain%logger%verbose) then
      write (*, '(a)') "Locating receiver points..."
    end if

    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    FringeSize => ReducedDomainInfo%fringe_size

    do n = 1, NumGrids
      Grid => Domain%grid(IndexToID(n))
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
      ReceiverMask = ovk_field_logical_(Grid%cart)
      ReceiverMask%values = Grid%mask%values .and. EdgeMask%values(&
        Grid%cart%is(1):Grid%cart%ie(1),Grid%cart%is(2):Grid%cart%ie(2), &
        Grid%cart%is(3):Grid%cart%ie(3))
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
            " receiver points on grid ", trim(IntToString(Grid%properties%id)), "."
        end if
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished locating receiver points."
    end if

  end subroutine LocateReceivers

  subroutine ChooseDonors(Domain, ReducedDomainInfo, DonorGridIDs)

    type(ovk_domain), intent(in) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_int), dimension(:), intent(out) :: DonorGridIDs

    integer :: i, j, k, m, n
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    integer, dimension(:,:), pointer :: ConnectionType
    integer, dimension(:), pointer :: FringeSize
    integer, dimension(:,:), pointer :: InterfacePadding
    integer(lk) :: l
    type(ovk_grid), pointer :: Grid_m, Grid_n
    type(ovk_overlap), pointer :: Overlap
    logical :: IsReceiver
    type(ovk_field_real) :: DonorEdgeDistance
    type(ovk_field_real) :: DonorResolution
    real(rk) :: EdgeDistance
    real(rk) :: Resolution
    logical :: BetterDonor

    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    if (Domain%logger%verbose) then
      write (*, '(a)') "Choosing donors..."
    end if

    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id
    ConnectionType => ReducedDomainInfo%connection_type
    FringeSize => ReducedDomainInfo%fringe_size
    InterfacePadding => ReducedDomainInfo%interface_padding

    do n = 1, NumGrids
      Grid_n => Domain%grid(IndexToID(n))
      DonorGridIDs(n) = ovk_field_int_(Grid_n%cart, 0)
      DonorEdgeDistance = ovk_field_real_(Grid_n%cart)
      DonorResolution = ovk_field_real_(Grid_n%cart)
      do m = 1, NumGrids
        if (ConnectionType(m,n) == OVK_CONNECTION_NONE) cycle
        Grid_m => Domain%grid(IndexToID(m))
        Overlap => Domain%overlap(Grid_m%properties%id,Grid_n%properties%id)
        l = 1_lk
        do k = Grid_n%cart%is(3), Grid_n%cart%ie(3)
          do j = Grid_n%cart%is(2), Grid_n%cart%ie(2)
            do i = Grid_n%cart%is(1), Grid_n%cart%ie(1)
              IsReceiver = iand(Grid_n%state%values(i,j,k),OVK_STATE_RECEIVER) /= 0
              if (IsReceiver .and. Overlap%mask%values(i,j,k)) then
                if (DonorGridIDs(n)%values(i,j,k) == 0) then
                  DonorGridIDs(n)%values(i,j,k) = Grid_m%properties%id
                  DonorEdgeDistance%values(i,j,k) = real(Overlap%edge_dists(l),kind=rk)/ &
                    real(max(FringeSize(m)+InterfacePadding(n,m),1),kind=rk)
                  DonorResolution%values(i,j,k) = Overlap%resolutions(l)
                else
                  EdgeDistance = real(Overlap%edge_dists(l),kind=rk)/ &
                    real(max(FringeSize(m)+InterfacePadding(n,m),1),kind=rk)
                  Resolution = Overlap%resolutions(l)
                  if (EdgeDistance > 1._rk-TOLERANCE .and. DonorEdgeDistance%values(i,j,k) > &
                    1._rk-TOLERANCE) then
                    BetterDonor = Resolution > DonorResolution%values(i,j,k)
                  else
                    BetterDonor = EdgeDistance > DonorEdgeDistance%values(i,j,k)
                  end if
                  if (BetterDonor) then
                    DonorGridIDs(n)%values(i,j,k) = Grid_m%properties%id
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
      if (Domain%logger%verbose) then
        write (*, '(3a)') "* Done choosing donors on grid ", &
          trim(IntToString(Grid_n%properties%id)), "."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished choosing donors."
    end if

  end subroutine ChooseDonors

  subroutine DetectOrphans(Domain, ReducedDomainInfo, DonorGridIDs)

    type(ovk_domain), intent(in) :: Domain
    type(t_reduced_domain_info), intent(in) :: ReducedDomainInfo
    type(ovk_field_int), dimension(:), intent(in) :: DonorGridIDs

    integer :: i, j, k, n
    integer :: NumGrids
    integer, dimension(:), pointer :: IndexToID
    type(ovk_grid), pointer :: Grid
    type(ovk_field_logical) :: ReceiverMask
    type(ovk_field_logical) :: OrphanMask
    type(ovk_field_int), pointer :: State
    integer :: NumWarnings
    integer, dimension(MAX_ND) :: ReceiverPoint
    integer(lk) :: NumOrphans

    if (Domain%logger%verbose) then
      write (*, '(a)') "Checking for orphans..."
    end if

    NumGrids = ReducedDomainInfo%ngrids
    IndexToID => ReducedDomainInfo%index_to_id

    do n = 1, NumGrids
      Grid => Domain%grid(IndexToID(n))
      call ovkFilterGridState(Grid, OVK_STATE_RECEIVER, OVK_ALL, ReceiverMask)
      OrphanMask = ovk_field_logical_(Grid%cart)
      OrphanMask%values = ReceiverMask%values .and. DonorGridIDs(n)%values == 0
      call ovkEditGridState(Grid, State)
      do k = Grid%cart%is(3), Grid%cart%ie(3)
        do j = Grid%cart%is(2), Grid%cart%ie(2)
          do i = Grid%cart%is(1), Grid%cart%ie(1)
            if (OrphanMask%values(i,j,k)) then
              State%values(i,j,k) = ior(iand(State%values(i,j,k),not(OVK_STATE_RECEIVER)), &
                OVK_STATE_ORPHAN)
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
              ReceiverPoint = [i,j,k]
              if (OrphanMask%values(i,j,k)) then
                if (NumWarnings <= 100) then
                  write (ERROR_UNIT, '(4a)') "WARNING: Orphan detected at point ", &
                    trim(TupleToString(ReceiverPoint(:Grid%cart%nd))), " of grid ", &
                    trim(IntToString(Grid%properties%id))
                  if (NumWarnings == 100) then
                    write (ERROR_UNIT, '(a)') "Further warnings suppressed."
                  end if
                  NumWarnings = NumWarnings + 1
                end if
              end if
            end do
          end do
        end do
        NumOrphans = ovkCountMask(OrphanMask)
        write (*, '(5a)') "* ", trim(LargeIntToString(NumOrphans)), &
          " orphan points found on grid ", trim(IntToString(Grid%properties%id)), "."
      end if
    end do

    if (Domain%logger%verbose) then
      write (*, '(a)') "Finished checking for orphans."
    end if

  end subroutine DetectOrphans

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
                            write (ERROR_UNIT, '(a)') "Further warnings suppressed."
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
    deallocate(ReducedDomainInfo%overlap_hole_cutting)
    deallocate(ReducedDomainInfo%connection_type)
    deallocate(ReducedDomainInfo%interp_scheme)
    deallocate(ReducedDomainInfo%fringe_size)
    deallocate(ReducedDomainInfo%interface_padding)

    call ResetDomainEdits(Domain)

  end subroutine FinalizeAssembly

end module ovkAssembly
