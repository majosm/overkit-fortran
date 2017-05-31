! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkOverset

  use ovkBoundingBox
  use ovkCart
  use ovkDonorAccel
  use ovkDonors
  use ovkField
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkInterp
  use ovkMask
  implicit none

  private

  ! API
  public :: ovkAssembleOverset
  public :: ovkPartitionReceivers
  public :: ovkGenerateOverlapOptimizationMask

contains

  subroutine ovkAssembleOverset(Grids, InterpData, FringeSize, InterpScheme, AllowOverlap, &
    AllowBoundaryHoleCutting, AllowOverlapHoleCutting, AllowInterpolation, DisjointFringes, &
    FringePadding, OverlapTolerance, HoleMasks, OrphanMasks)

    type(ovk_grid), dimension(:), intent(inout) :: Grids
    type(ovk_interp), dimension(size(Grids)), intent(out) :: InterpData
    integer, dimension(size(Grids)), intent(in), optional :: FringeSize
    integer, dimension(size(Grids)), intent(in), optional :: InterpScheme
    logical, dimension(size(Grids),size(Grids)), intent(in), optional :: AllowOverlap
    logical, dimension(size(Grids),size(Grids)), intent(in), optional :: AllowBoundaryHoleCutting
    logical, dimension(size(Grids),size(Grids)), intent(in), optional :: AllowOverlapHoleCutting
    logical, dimension(size(Grids),size(Grids)), intent(in), optional :: AllowInterpolation
    logical, dimension(size(Grids),size(Grids)), intent(in), optional :: DisjointFringes
    integer, dimension(size(Grids),size(Grids)), intent(in), optional :: FringePadding
    real(rk), dimension(size(Grids),size(Grids)), intent(in), optional :: OverlapTolerance
    type(ovk_field_logical), dimension(size(Grids)), intent(out), optional :: HoleMasks
    type(ovk_field_logical), dimension(size(Grids)), intent(out), optional :: OrphanMasks

    integer, dimension(:), allocatable :: FringeSize_
    integer, dimension(:), allocatable :: InterpScheme_
    logical, dimension(:,:), allocatable :: AllowOverlap_
    logical, dimension(:,:), allocatable :: AllowBoundaryHoleCutting_
    logical, dimension(:,:), allocatable :: AllowOverlapHoleCutting_
    logical, dimension(:,:), allocatable :: AllowInterpolation_
    logical, dimension(:,:), allocatable :: DisjointFringes_
    integer, dimension(:,:), allocatable :: FringePadding_
    real(rk), dimension(:,:), allocatable :: OverlapTolerance_
    integer :: i, j, k, l, m, n, p, q, r
    integer :: ClockInitial, ClockFinal, ClockRate
    character(len=STRING_LENGTH) :: NumPointsTotalString
    character(len=STRING_LENGTH) :: iSString, iEString, jSString, jEString, kSString, kEString
    type(ovk_bbox) :: Bounds
    type(ovk_bbox), dimension(:), allocatable :: OverlapBounds
    integer :: PaddingMax, PaddingSum
    real(rk) :: PaddingFrac
    integer :: Padding1, Padding2
    integer(lk) :: NumReceivers
    integer(lk) :: NumOuterReceivers, NumInnerReceivers
    integer(lk) :: NumInvalidatedDonors, NumInvalidatedDonors1, NumInvalidatedDonors2
    integer(lk) :: NumOrphans
    integer(lk) :: NumRemovedPoints
    type(ovk_donor_accel) :: DonorAccel
    type(ovk_donors), dimension(:,:), allocatable :: PairwiseDonors
    type(ovk_donors), dimension(:), allocatable :: Donors
    type(ovk_field_logical), dimension(:), allocatable :: OuterReceiverMasks
    type(ovk_field_logical), dimension(:), allocatable :: InnerReceiverMasks
    type(ovk_field_logical) :: DonorMask
    type(ovk_field_logical) :: EdgeMask1
    type(ovk_field_logical) :: EdgeMask2
    type(ovk_field_logical) :: BoundaryMask
    type(ovk_field_logical) :: InteriorMask
    type(ovk_field_logical) :: ExteriorMask
    type(ovk_field_logical) :: CoarseToFineMask
    type(ovk_field_logical) :: NearCrossoverMask1, NearCrossoverMask2
    type(ovk_field_logical) :: OverlapOptimizationMask
    type(ovk_field_logical) :: ReceiverMask
    type(ovk_field_logical) :: OrphanMask

    type(ovk_field_logical), dimension(:), allocatable :: ValidCellMasks
    type(ovk_field_logical) :: ValidCellInnerEdgeMask
    type(ovk_field_logical), dimension(:), allocatable :: ExpandableCellMasks
    type(ovk_field_int), dimension(:), allocatable :: CellQualities
    logical :: ExpandableCell
    integer, dimension(MAX_ND) :: ReceiverPoint
    integer, dimension(MAX_ND) :: DonorCell
    integer, dimension(MAX_ND) :: NeighborCellLower, NeighborCellUpper
    logical :: AwayFromBoundary
    integer, dimension(MAX_ND) :: NeighborCell
    integer, dimension(MAX_ND) :: DonorCellShift
    integer :: BestCellQuality
    integer :: CellQuality
    integer, dimension(MAX_ND) :: ExpandedDonorCell
    real(rk), dimension(Grids(1)%cart%nd) :: ReceiverCoords
    real(rk), dimension(Grids(1)%cart%nd) :: DonorCellCoords
    real(rk), dimension(Grids(1)%cart%nd) :: ExpandedDonorCellCoords

    allocate(FringeSize_(size(Grids)))
    if (present(FringeSize)) then
      FringeSize_ = FringeSize
    else
      FringeSize_ = 2
    end if

    allocate(InterpScheme_(size(Grids)))
    if (present(InterpScheme)) then
      InterpScheme_ = InterpScheme
    else
      InterpScheme_ = OVK_INTERP_LINEAR
    end if

    allocate(AllowOverlap_(size(Grids),size(Grids)))
    if (present(AllowOverlap)) then
      AllowOverlap_ = AllowOverlap
    else
      AllowOverlap_ = .true.
    end if

    allocate(AllowBoundaryHoleCutting_(size(Grids),size(Grids)))
    if (present(AllowBoundaryHoleCutting)) then
      AllowBoundaryHoleCutting_ = AllowBoundaryHoleCutting
    else
      AllowBoundaryHoleCutting_ = .true.
    end if

    allocate(AllowOverlapHoleCutting_(size(Grids),size(Grids)))
    if (present(AllowOverlapHoleCutting)) then
      AllowOverlapHoleCutting_ = AllowOverlapHoleCutting
    else
      AllowOverlapHoleCutting_ = .true.
    end if

    allocate(AllowInterpolation_(size(Grids),size(Grids)))
    if (present(AllowInterpolation)) then
      AllowInterpolation_ = AllowInterpolation
    else
      AllowInterpolation_ = .true.
    end if

    allocate(DisjointFringes_(size(Grids),size(Grids)))
    if (present(DisjointFringes)) then
      DisjointFringes_ = DisjointFringes
    else
      DisjointFringes_ = .true.
    end if

    allocate(FringePadding_(size(Grids),size(Grids)))
    if (present(FringePadding)) then
      FringePadding_ = FringePadding
    else
      FringePadding_ = 0
    end if

    allocate(OverlapTolerance_(size(Grids),size(Grids)))
    if (present(OverlapTolerance)) then
      OverlapTolerance_ = OverlapTolerance
    else
      OverlapTolerance_ = 0._rk
    end if

    if (OVK_VERBOSE) then
      call system_clock(ClockInitial, ClockRate)
      write (*, '(a)') "Overset grid assembly started..."
    end if

    do m = 1, size(Grids)
      Grids(m)%id = m
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Grid info:"
      write (*, '(3a)') "* Dimension: ", trim(IntToString(Grids(1)%cart%nd)), "D"
      write (*, '(2a)') "* Number of grids: ", trim(IntToString(size(Grids)))
      NumPointsTotalString = LargeIntToString(sum([(ovkCartCount(Grids(m)%cart),m=1,size(Grids))]))
      write (*, '(2a)') "* Total number of grid points: ", trim(NumPointsTotalString)
      do m = 1, size(Grids)
        NumPointsTotalString = LargeIntToString(ovkCartCount(Grids(m)%cart))
        iSString = IntToString(Grids(m)%cart%is(1))
        iEString = IntToString(Grids(m)%cart%ie(1))
        jSString = IntToString(Grids(m)%cart%is(2))
        jEString = IntToString(Grids(m)%cart%ie(2))
        kSString = IntToString(Grids(m)%cart%is(3))
        kEString = IntToString(Grids(m)%cart%ie(3))
        write (*, '(3a)', advance="no") "* Grid ", trim(IntToString(Grids(m)%id)), ": "
        write (*, '(2a)', advance="no") trim(NumPointsTotalString), " points "
        select case (Grids(m)%cart%nd)
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

    allocate(OverlapBounds(size(Grids)))
    do n = 1, size(Grids)
      OverlapBounds(n) = ovk_bbox_(Grids(n)%cart%nd)
      AllowOverlap_(n,n) = .false.
      do m = 1, size(Grids)
        if (AllowOverlap_(m,n)) then
          Bounds = ovkBBIntersect(Grids(m)%bounds, Grids(n)%bounds)
          if (.not. ovkBBIsEmpty(Bounds)) then
            OverlapBounds(n) = ovkBBUnion(OverlapBounds(n), Bounds)
          else
            AllowOverlap_(m,n) = .false.
          end if
        end if
      end do
    end do

    do n = 1, size(Grids)
      do m = 1, size(Grids)
        if (.not. AllowOverlap_(m,n)) then
          AllowBoundaryHoleCutting_(m,n) = .false.
          AllowOverlapHoleCutting_(m,n) = .false.
          AllowInterpolation_(m,n) = .false.
        end if
      end do
    end do

    allocate(PairwiseDonors(size(Grids),size(Grids)))

    do m = 1, size(Grids)
      if (OVK_VERBOSE) then
        write (*, '(3a)') "* Generating donor search accelerator on grid ", trim(IntToString(m)), &
          "..."
      end if
      call ovkGenerateDonorAccel(Grids(m), DonorAccel, Bounds=OverlapBounds(m), &
        OverlapTolerance=maxval(OverlapTolerance_(m,:)))
      do n = 1, size(Grids)
        if (AllowOverlap_(m,n)) then
          call ovkFindDonors(Grids(m), Grids(n), DonorAccel, PairwiseDonors(m,n), &
            OverlapTolerance=OverlapTolerance_(m,n))
          if (OVK_VERBOSE) then
            NumReceivers = ovkCountMask(PairwiseDonors(m,n)%valid_mask)
            write (*, '(7a)') "* ", trim(LargeIntToString(NumReceivers)), &
              " candidate donors from grid ", trim(IntToString(m)), " to grid ", &
              trim(IntToString(n)), " found."
          end if
        else
          call ovkMakeDonors(PairwiseDonors(m,n), Grids(n)%cart)
          PairwiseDonors(m,n)%valid_mask%values = .false.
        end if
      end do
      call ovkDestroyDonorAccel(DonorAccel)
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished searching for candidate donor/receiver pairs."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Cutting boundary holes..."
    end if

    if (present(HoleMasks)) then
      do m = 1, size(Grids)
        HoleMasks(m) = ovk_field_logical_(Grids(m)%cart, .false.)
      end do
    end if

    do n = 1, size(Grids)
      if (any(AllowBoundaryHoleCutting_(:,n))) then
        BoundaryMask = Grids(n)%boundary_mask
        InteriorMask = ovk_field_logical_(Grids(n)%cart, .false.)
        do m = 1, size(Grids)
          if (AllowBoundaryHoleCutting_(m,n)) then
            call ovkFindMaskEdge(PairwiseDonors(m,n)%valid_mask, OVK_EDGE_TYPE_INNER, EdgeMask1)
            call ovkGenerateDonorMask(Grids(m), Grids(n), PairwiseDonors(m,n), DonorMask, &
              ReceiverSubset=EdgeMask1)
            call ovkFindMaskEdge(PairwiseDonors(n,m)%valid_mask, OVK_EDGE_TYPE_INNER, EdgeMask2)
            DonorMask%values = DonorMask%values .and. .not. EdgeMask2%values
            i = 0
            do while (any(DonorMask%values))
              call ovkGrowMask(EdgeMask2, 1)
              DonorMask%values = DonorMask%values .and. .not. EdgeMask2%values
              i = i + 1
            end do
            EdgeMask2 = Grids(m)%boundary_mask
            do j = 1, i
              call ovkGrowMask(EdgeMask2, 1)
            end do
            call ovkGenerateReceiverMask(Grids(n), Grids(m), PairwiseDonors(m,n), ReceiverMask, &
              DonorSubset=EdgeMask2)
            call ovkGrowMask(ReceiverMask, 1)
            ReceiverMask%values = ReceiverMask%values .and. EdgeMask1%values
            BoundaryMask%values = BoundaryMask%values .or. ReceiverMask%values
            call ovkFindMaskEdge(ReceiverMask, OVK_EDGE_TYPE_OUTER, EdgeMask1)
            InteriorMask%values = InteriorMask%values .or. (EdgeMask1%values( &
              Grids(n)%cart%is(1):Grids(n)%cart%ie(1),Grids(n)%cart%is(2):Grids(n)%cart%ie(2), &
              Grids(n)%cart%is(3):Grids(n)%cart%ie(3)) .and. PairwiseDonors(m,n)%valid_mask%values)
          end if
        end do
        call ovkFillMask(InteriorMask, BoundaryMask)
        InteriorMask%values = InteriorMask%values .or. BoundaryMask%values
        call ovkFindMaskEdge(InteriorMask, OVK_EDGE_TYPE_INNER, EdgeMask1)
        BoundaryMask%values = BoundaryMask%values .and. EdgeMask1%values
        call ovkFindMaskEdge(BoundaryMask, OVK_EDGE_TYPE_OUTER, EdgeMask1)
        ExteriorMask = ovk_field_logical_(Grids(n)%cart)
        ExteriorMask%values = EdgeMask1%values(Grids(n)%cart%is(1):Grids(n)%cart%ie(1), &
          Grids(n)%cart%is(2):Grids(n)%cart%ie(2),Grids(n)%cart%is(3):Grids(n)%cart%ie(3)) .and. &
          .not. InteriorMask%values
        call ovkFillMask(ExteriorMask, BoundaryMask)
        Grids(n)%grid_mask%values = Grids(n)%grid_mask%values .and. .not. ExteriorMask%values
        Grids(n)%boundary_mask%values = Grids(n)%boundary_mask%values .and. .not. &
          ExteriorMask%values
        do m = 1, size(Grids)
          if (AllowOverlap_(m,n)) then
            PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
              .not. ExteriorMask%values
          end if
        end do
        do m = 1, size(Grids)
          if (AllowOverlap_(n,m)) then
            call ovkGenerateReceiverMask(Grids(m), Grids(n), PairwiseDonors(n,m), ReceiverMask, &
              DonorSubset=ExteriorMask)
            PairwiseDonors(n,m)%valid_mask%values = PairwiseDonors(n,m)%valid_mask%values .and. &
              .not. ReceiverMask%values
          end if
        end do
        if (present(HoleMasks)) then
          HoleMasks(n)%values = HoleMasks(n)%values .or. ExteriorMask%values
        end if
        if (OVK_VERBOSE) then
          NumRemovedPoints = ovkCountMask(ExteriorMask)
        end if
      else
        if (OVK_VERBOSE) then
          NumRemovedPoints = 0_lk
        end if
      end if
      if (OVK_VERBOSE) then
        write (*, '(5a)') "* ", trim(LargeIntToString(NumRemovedPoints)), &
          " points removed from grid ", trim(IntToString(n)), "."
      end if
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished cutting boundary holes."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Partitioning receivers..."
    end if

    allocate(OuterReceiverMasks(size(Grids)))
    allocate(InnerReceiverMasks(size(Grids)))

    do n = 1, size(Grids)
      call ovkPartitionReceivers(Grids(n), PairwiseDonors(:,n), OuterReceiverMasks(n), &
        InnerReceiverMasks(n), FringeSize_(n))
      if (OVK_VERBOSE) then
        NumOuterReceivers = ovkCountMask(OuterReceiverMasks(n))
        NumInnerReceivers = ovkCountMask(InnerReceiverMasks(n))
        write (*, '(7a)') "* Partitioned receivers on grid ", trim(IntToString(n)), " into ", &
          trim(LargeIntToString(NumOuterReceivers)), " near-edge points and ", &
          trim(LargeIntToString(NumInnerReceivers)), " interior points."
      end if
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished partitioning receivers."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Invalidating coarse-to-fine donors..."
    end if

    do n = 1, size(Grids)
      do m = 1, size(Grids)
        if (AllowOverlapHoleCutting_(m,n) .and. AllowOverlapHoleCutting_(n,m)) then
          call ovkGenerateCoarseToFineMask(Grids(n), PairwiseDonors(m,n), CoarseToFineMask, &
            Subset=InnerReceiverMasks(n))
          PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
            .not. CoarseToFineMask%values
          if (OVK_VERBOSE) then
            NumInvalidatedDonors = ovkCountMask(CoarseToFineMask)
            write (*, '(7a)') "* ", trim(LargeIntToString(NumInvalidatedDonors)), &
              " donors from grid ", trim(IntToString(m)), " to grid ", trim(IntToString(n)), &
              " invalidated."
          end if
        else if (.not. AllowOverlapHoleCutting_(m,n)) then
          ! Behave as if grid m is coarser everywhere
          PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
              .not. InnerReceiverMasks(n)%values
        end if
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished invalidating coarse-to-fine donors."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Choosing best donor/receiver pairs near grid edges..."
    end if

    do n = 1, size(Grids)
      do m = 1, size(Grids)
        if (AllowInterpolation_(m,n) .and. DisjointFringes_(m,n)) then
          call ovkGenerateReceiverMask(Grids(n), Grids(m), PairwiseDonors(m,n), &
            ReceiverMask, DonorSubset=OuterReceiverMasks(m))
          PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
            .not. ReceiverMask%values
        end if
      end do
    end do

    do n = 1, size(Grids)
      call ovkChooseDonors(PairwiseDonors(:,n), Subset=OuterReceiverMasks(n))
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished choosing best donor/receiver pairs near grid edges."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Invalidating cyclic donors and applying fringe padding..."
    end if

    do n = 1, size(Grids)
      do m = n+1, size(Grids)
        if (AllowInterpolation_(m,n) .and. AllowInterpolation_(n,m)) then
          PaddingMax = max(FringePadding_(m,n), FringePadding_(n,m))
          PaddingSum = FringePadding_(m,n) + FringePadding_(n,m)
          PaddingFrac = real(FringePadding_(n,m),kind=rk)/real(max(PaddingSum,1),kind=rk)
          Padding1 = int(real(PaddingMax,kind=rk) * PaddingFrac + 0.5_rk)
          Padding2 = PaddingMax - Padding1
          if (Padding1 > 0 .or. Padding2 > 0 .or. DisjointFringes_(m,n) .or. &
            DisjointFringes_(n,m)) then
            call ovkGenerateNearCrossoverMask(Grids(m), Grids(n), PairwiseDonors(n,m), &
              PairwiseDonors(m,n), Padding1, NearCrossoverMask1, Subset1=InnerReceiverMasks(m))
            call ovkGenerateNearCrossoverMask(Grids(n), Grids(m), PairwiseDonors(m,n), &
              PairwiseDonors(n,m), Padding2, NearCrossoverMask2, Subset1=InnerReceiverMasks(n))
            PairwiseDonors(n,m)%valid_mask%values = PairwiseDonors(n,m)%valid_mask%values .and. &
              .not. NearCrossoverMask1%values
            PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
              .not. NearCrossoverMask2%values
            if (OVK_VERBOSE) then
              NumInvalidatedDonors1 = ovkCountMask(NearCrossoverMask1)
              NumInvalidatedDonors2 = ovkCountMask(NearCrossoverMask2)
            end if
            call ovkGenerateNearCrossoverMask(Grids(m), Grids(n), PairwiseDonors(n,m), &
              PairwiseDonors(m,n), FringePadding_(n,m), NearCrossoverMask1, &
              Subset1=InnerReceiverMasks(m))
            call ovkGenerateNearCrossoverMask(Grids(n), Grids(m), PairwiseDonors(m,n), &
              PairwiseDonors(n,m), FringePadding_(m,n), NearCrossoverMask2, &
              Subset1=InnerReceiverMasks(n))
            PairwiseDonors(n,m)%valid_mask%values = PairwiseDonors(n,m)%valid_mask%values .and. &
              .not. NearCrossoverMask1%values
            PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
              .not. NearCrossoverMask2%values
            if (OVK_VERBOSE) then
              NumInvalidatedDonors1 = NumInvalidatedDonors1 + ovkCountMask(NearCrossoverMask1)
              NumInvalidatedDonors2 = NumInvalidatedDonors2 + ovkCountMask(NearCrossoverMask2)
            end if
          else
            if (OVK_VERBOSE) then
              NumInvalidatedDonors1 = 0
              NumInvalidatedDonors2 = 0
            end if
          end if
          if (OVK_VERBOSE) then
            write (*, '(7a)') "* ", trim(LargeIntToString(NumInvalidatedDonors1)), &
              " donor/receiver pairs from grid ", trim(IntToString(n)), " to ", &
              trim(IntToString(m)), " invalidated."
            write (*, '(7a)') "* ", trim(LargeIntToString(NumInvalidatedDonors2)), &
              " donor/receiver pairs from grid ", trim(IntToString(m)), " to ", &
              trim(IntToString(n)), " invalidated."
          end if
        end if
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished invalidating cyclic donors and applying fringe padding."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Choosing best donor/receiver pairs on grid interior..."
    end if

    do n = 1, size(Grids)
      call ovkChooseDonors(PairwiseDonors(:,n), Subset=InnerReceiverMasks(n))
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished choosing best donor/receiver pairs on grid interior."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Merging remaining candidate donor/receiver pairs into final set..."
    end if

    allocate(Donors(size(Grids)))

    do m = 1, size(Grids)
      call ovkMergeDonors(PairwiseDonors(:,m), Donors(m))
      if (OVK_VERBOSE) then
        NumReceivers = ovkCountMask(Donors(m)%valid_mask)
        write (*, '(5a)') "* ", trim(LargeIntToString(NumReceivers)), &
          " receiver points on grid ", trim(IntToString(m)), "."
      end if
    end do

    deallocate(PairwiseDonors)

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished merging candidate donor/receiver pairs."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Cutting overlap holes..."
    end if

    do n = 1, size(Grids)
      if (any(AllowOverlapHoleCutting_(:,n))) then
        call ovkGenerateOverlapOptimizationMask(Grids(n), Donors(n), OverlapOptimizationMask, &
          FringeSize=FringeSize_(n))
        do k = Grids(n)%cart%is(3), Grids(n)%cart%ie(3)
          do j = Grids(n)%cart%is(2), Grids(n)%cart%ie(2)
            do i = Grids(n)%cart%is(1), Grids(n)%cart%ie(1)
              if (OverlapOptimizationMask%values(i,j,k)) then
                Grids(n)%grid_mask%values(i,j,k) = .false.
                Grids(n)%boundary_mask%values(i,j,k) = .false.
                Donors(n)%valid_mask%values(i,j,k) = .false.
                OuterReceiverMasks(n)%values(i,j,k) = .false.
                InnerReceiverMasks(n)%values(i,j,k) = .false.
              end if
            end do
          end do
        end do
        if (present(HoleMasks)) then
          HoleMasks(n)%values = HoleMasks(n)%values .or. OverlapOptimizationMask%values
        end if
        if (OVK_VERBOSE) then
          NumRemovedPoints = ovkCountMask(OverlapOptimizationMask)
          write (*, '(5a)') "* ", trim(LargeIntToString(NumRemovedPoints)), &
            " points removed from grid ", trim(IntToString(n)), "."
        end if
      end if
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished cutting overlap holes."
    end if

    ! TODO: Decide whether it makes sense for InterpScheme_ to be non-scalar (and if so, should
    ! it be per-grid or per-grid-pair?)
    if (any(InterpScheme_ == OVK_INTERP_CUBIC)) then

      if (OVK_VERBOSE) then
        write (*, '(a)') "Expanding donor cells for cubic interpolation..."
      end if

      allocate(ValidCellMasks(size(Grids)))
      allocate(ExpandableCellMasks(size(Grids)))
      allocate(CellQualities(size(Grids)))

      do m = 1, size(Grids)
        call GenerateValidCellMask(Grids(m), ValidCellMasks(m))
        call ovkFindMaskEdge(ValidCellMasks(m), OVK_EDGE_TYPE_INNER, ValidCellInnerEdgeMask)
        ExpandableCellMasks(m) = ovk_field_logical_(Grids(m)%cell_cart)
        ExpandableCellMasks(m)%values = ValidCellMasks(m)%values .and. .not. &
          ValidCellInnerEdgeMask%values
        call CountNeighbors(ExpandableCellMasks(m), CellQualities(m))
      end do

      do n = 1, size(Grids)
        if (InterpScheme_(n) /= OVK_INTERP_CUBIC) cycle
        do k = Grids(n)%cart%is(3), Grids(n)%cart%ie(3)
          do j = Grids(n)%cart%is(2), Grids(n)%cart%ie(2)
            do i = Grids(n)%cart%is(1), Grids(n)%cart%ie(1)
              ReceiverPoint = [i,j,k]
              if (Donors(n)%valid_mask%values(i,j,k)) then
                m = Donors(n)%grid_ids%values(i,j,k)
                do l = 1, Grids(m)%cart%nd
                  DonorCell(l) = Donors(n)%cells(l)%values(i,j,k)
                  DonorCellCoords(l) = Donors(n)%cell_coords(l)%values(i,j,k)
                end do
                DonorCell(Grids(m)%cell_cart%nd+1:) = 1
                if (ExpandableCellMasks(m)%values(DonorCell(1),DonorCell(2),DonorCell(3))) then
                  ExpandableCell = .true.
                else
                  BestCellQuality = 0
                  DonorCellShift = 0
                  NeighborCellLower(:Grids(m)%cart%nd) = DonorCell(:Grids(m)%cart%nd)-1
                  NeighborCellLower(Grids(m)%cart%nd+1:) = DonorCell(Grids(m)%cart%nd+1:)
                  NeighborCellUpper(:Grids(m)%cart%nd) = DonorCell(:Grids(m)%cart%nd)+1
                  NeighborCellUpper(Grids(m)%cart%nd+1:) = DonorCell(Grids(m)%cart%nd+1:)
                  AwayFromBoundary = ovkCartContains(Grids(m)%cell_cart, NeighborCellLower) .and. &
                    ovkCartContains(Grids(m)%cell_cart, NeighborCellUpper)
                  if (AwayFromBoundary) then
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
                          NeighborCell(:Grids(m)%cell_cart%nd) = ovkCartPeriodicAdjust( &
                            Grids(m)%cell_cart, NeighborCell)
                          if (ovkCartContains(Grids(m)%cell_cart, NeighborCell)) then
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
                    DonorCellCoords = DonorCellCoords - real(DonorCellShift(:Grids(m)%cell_cart%nd),kind=rk)
                  else
                    ExpandableCell = .false.
                    if (OVK_VERBOSE) then
                      write (ERROR_UNIT, '(8a)') "WARNING: Could not use cubic interpolation for donor cell ", &
                        trim(TupleToString(DonorCell(:Grids(m)%cart%nd))), " of grid ", trim(IntToString(m)), &
                        " corresponding to receiver point ", trim(TupleToString(ReceiverPoint(:Grids(n)%cart%nd))), &
                        " of grid ", trim(IntToString(n))
                    end if
                  end if
                end if
                if (ExpandableCell) then
                  do l = 1, Grids(n)%cart%nd
                    ReceiverCoords(l) = Grids(n)%xyz(l)%values(ReceiverPoint(1),ReceiverPoint(2), &
                      ReceiverPoint(3))
                  end do
                  call ExpandDonorCell(Grids(m), DonorCell, DonorCellCoords, &
                    ReceiverCoords, ExpandedDonorCell, ExpandedDonorCellCoords)
                  ExpandedDonorCell(Grids(m)%cart%nd+1:) = 1
                  do l = 1, Grids(m)%cart%nd
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

      do n = 1, size(Grids)
        do m = 1, size(Grids)
          if (AllowInterpolation_(m,n) .and. DisjointFringes_(m,n)) then
            call ovkGenerateReceiverMask(Grids(n), Grids(m), Donors(n), ReceiverMask, &
              DonorSubset=Donors(m)%valid_mask)
            Donors(n)%valid_mask%values = Donors(n)%valid_mask%values .and. .not. &
              ReceiverMask%values
          end if
        end do
      end do

      if (OVK_VERBOSE) then
        write (*, '(a)') "Finished expanding donor cells."
      end if

    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Searching for orphan points..."
    end if

    do m = 1, size(Grids)
      call ovkGenerateOrphanMask(Donors(m), OuterReceiverMasks(m), OrphanMask)
      if (present(OrphanMasks)) then
        OrphanMasks(m) = OrphanMask
      end if
      if (OVK_VERBOSE) then
        do k = Grids(m)%cart%is(3), Grids(m)%cart%ie(3)
          do j = Grids(m)%cart%is(2), Grids(m)%cart%ie(2)
            do i = Grids(m)%cart%is(1), Grids(m)%cart%ie(1)
              ReceiverPoint = [i,j,k]
              if (OrphanMask%values(i,j,k)) then
                write (ERROR_UNIT, '(4a)') "WARNING: Orphan detected at point ", &
                  trim(TupleToString(ReceiverPoint(:Grids(m)%cart%nd))), " of grid ", &
                  trim(IntToString(m))
              end if
            end do
          end do
        end do
        NumOrphans = ovkCountMask(OrphanMask)
        write (*, '(5a)') "* ", trim(LargeIntToString(NumOrphans)), &
          " orphan points found on grid ", trim(IntToString(m)), "."
      end if
    end do

    deallocate(OuterReceiverMasks)
    deallocate(InnerReceiverMasks)

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished searching for orphan points."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Generating interpolation data..."
    end if

    do m = 1, size(Grids)
      call ovkGenerateInterpData(Donors(m), InterpData(m), InterpScheme=InterpScheme_(m))
      if (OVK_VERBOSE) then
        write (*, '(3a)') "* Generated interpolation data for grid ", trim(IntToString(m)), "."
      end if
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished generating interpolation data."
    end if

    if (OVK_VERBOSE) then
      call system_clock(ClockFinal, ClockRate)
      write (*, '(a,f0.3,a)') "Overset grid assembly finished (time: ", &
        real(ClockFinal-ClockInitial,kind=rk)/real(ClockRate,kind=rk), " seconds)."
    end if

  end subroutine ovkAssembleOverset

  subroutine ovkPartitionReceivers(Grid, CandidateDonors, OuterReceiverMask, InnerReceiverMask, &
    OuterWidth)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_donors), dimension(:), intent(in) :: CandidateDonors
    type(ovk_field_logical), intent(out) :: OuterReceiverMask
    type(ovk_field_logical), intent(out) :: InnerReceiverMask
    integer, intent(in) :: OuterWidth

    type(ovk_field_logical) :: GridOuterEdgeMask
    type(ovk_field_logical) :: OverlapMask
    type(ovk_field_logical) :: OverlapOuterEdgeMask
    type(ovk_field_logical) :: NonBoundaryMask
    type(ovk_field_logical) :: NonBoundaryOuterEdgeMask

    call ovkFindMaskEdge(Grid%grid_mask, OVK_EDGE_TYPE_OUTER, GridOuterEdgeMask)

    NonBoundaryMask = ovk_field_logical_(Grid%cart)
    NonBoundaryMask%values = .not. Grid%boundary_mask%values

    call ovkFindMaskEdge(NonBoundaryMask, OVK_EDGE_TYPE_OUTER, NonBoundaryOuterEdgeMask)

    GridOuterEdgeMask%values = GridOuterEdgeMask%values .and. NonBoundaryOuterEdgeMask%values

    call ovkGenerateOverlapMask(CandidateDonors, OverlapMask)
    call ovkFindMaskEdge(OverlapMask, OVK_EDGE_TYPE_OUTER, OverlapOuterEdgeMask)

    OverlapOuterEdgeMask%values = OverlapOuterEdgeMask%values .and. GridOuterEdgeMask%values

    call ovkGrowMask(OverlapOuterEdgeMask, OuterWidth)

    OuterReceiverMask = ovk_field_logical_(Grid%cart)
    OuterReceiverMask%values = OverlapMask%values .and. OverlapOuterEdgeMask%values( &
      Grid%cart%is(1):Grid%cart%ie(1),Grid%cart%is(2):Grid%cart%ie(2), &
      Grid%cart%is(3):Grid%cart%ie(3))

    InnerReceiverMask = ovk_field_logical_(Grid%cart)
    InnerReceiverMask%values = OverlapMask%values .and. .not. OverlapOuterEdgeMask%values( &
      Grid%cart%is(1):Grid%cart%ie(1),Grid%cart%is(2):Grid%cart%ie(2), &
      Grid%cart%is(3):Grid%cart%ie(3))

  end subroutine ovkPartitionReceivers

  subroutine ovkGenerateOverlapOptimizationMask(Grid, Donors, OverlapOptimizationMask, FringeSize)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_donors), intent(in) :: Donors
    type(ovk_field_logical), intent(out) :: OverlapOptimizationMask
    integer, intent(in), optional :: FringeSize

    integer :: FringeSize_

    if (present(FringeSize)) then
      FringeSize_ = FringeSize
    else
      FringeSize_ = 2
    end if

    OverlapOptimizationMask = ovk_field_logical_(Grid%cart)
    OverlapOptimizationMask%values = .not. Grid%grid_mask%values .or. Donors%valid_mask%values

    call ovkGrowMask(OverlapOptimizationMask, -FringeSize_)

    OverlapOptimizationMask%values = OverlapOptimizationMask%values .and. Grid%grid_mask%values

  end subroutine ovkGenerateOverlapOptimizationMask

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
    logical :: AwayFromBoundary
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
          AwayFromBoundary = ovkCartContains(Mask%cart, NeighborLower) .and. &
            ovkCartContains(Mask%cart, NeighborUpper)
          if (AwayFromBoundary) then
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
    logical :: AwayFromBoundary
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

    AwayFromBoundary = ovkCartContains(DonorGrid%cart, ExpandedDonorLower) .and. &
      ovkCartContains(DonorGrid%cart, ExpandedDonorUpper)

    ExpandedDonorCell = ovkCartPeriodicAdjust(DonorGrid%cell_cart, ExpandedDonorLower)

    if (AwayFromBoundary) then
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
        VertexCoords(:,i) = ovkPeriodicExtend(DonorGrid%cart, DonorGrid%periodic_length, Vertex, &
          PrincipalCoords)
      end do
    end if

    Success = .true.

    select case (DonorGrid%cart%nd)
    case (2)
      select case (DonorGrid%grid_type)
      case (OVK_GRID_TYPE_CARTESIAN)
        ExpandedDonorCellCoords = ovkRectangleIsoInverseCubic(VertexCoords, ReceiverCoords)
      case default
        ExpandedDonorCellCoords = ovkQuadIsoInverseCubic(VertexCoords, ReceiverCoords, &
          Guess=DonorCellCoords, Success=Success)
      end select
    case (3)
      select case (DonorGrid%grid_type)
      case (OVK_GRID_TYPE_CARTESIAN)
        ExpandedDonorCellCoords = ovkCuboidIsoInverseCubic(VertexCoords, ReceiverCoords)
      case default
        ExpandedDonorCellCoords = ovkHexahedronIsoInverseCubic(VertexCoords, ReceiverCoords, &
          Guess=DonorCellCoords, Success=Success)
      end select
    end select

    if (OVK_VERBOSE) then
      if (.not. Success) then
        write (ERROR_UNIT, '(6a)') "WARNING: Failed to compute isoparametric coordinates of ", &
          trim(CoordsToString(ReceiverCoords)), " in expanded donor cell surrounding cell ", &
          trim(TupleToString(DonorCell)), " of grid ", trim(IntToString(DonorGrid%id))
      end if
    end if

  end subroutine ExpandDonorCell

end module ovkOverset
