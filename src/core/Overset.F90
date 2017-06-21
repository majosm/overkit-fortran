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
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkInterp
  use ovkMask
  implicit none

  private

  ! API
  public :: ovkAssemble
  public :: ovkPartitionReceivers
  public :: ovkGenerateOverlapOptimizationMask

contains

  subroutine ovkAssemble(Assembler)

    type(ovk_assembler), intent(inout) :: Assembler

    integer :: i, j, k, l, m, n, p, q, r
    integer :: ClockInitial, ClockFinal, ClockRate
    integer :: NumDims
    integer :: NumGrids
    type(ovk_grid), dimension(:), pointer :: Grids
    integer, dimension(:), allocatable :: InterpScheme
    integer, dimension(:), allocatable :: FringeSize
    integer, dimension(:), allocatable :: FringePadding
    character(len=STRING_LENGTH) :: NumPointsTotalString
    character(len=STRING_LENGTH) :: iSString, iEString, jSString, jEString, kSString, kEString
    type(ovk_bbox) :: Bounds
    type(ovk_bbox), dimension(:), allocatable :: OverlapBounds
    real(rk), dimension(:), allocatable :: MaxOverlapTolerance
    integer :: Connection1, Connection2
    logical :: Disjoint1, Disjoint2
    integer :: Padding1, Padding2
    integer :: PaddingMax, PaddingSum
    real(rk) :: PaddingFrac
    integer :: WeightedPadding1, WeightedPadding2
    integer(lk) :: NumReceivers
    integer(lk) :: NumOuterReceivers, NumInnerReceivers
    integer(lk) :: NumInvalidatedDonors, NumInvalidatedDonors1, NumInvalidatedDonors2
    integer(lk) :: NumOrphans
    integer(lk) :: NumRemovedPoints
    type(ovk_donor_accel) :: DonorAccel
    type(ovk_donors), dimension(:,:), allocatable :: PairwiseDonors
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
    type(ovk_donors), dimension(:), pointer :: Donors
    type(ovk_connectivity), pointer :: Connectivity
    integer :: StencilSize
    type(ovk_interp), dimension(:), pointer :: InterpData

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

    if (OVK_VERBOSE) then
      call system_clock(ClockInitial, ClockRate)
      write (*, '(a)') "Overset grid assembly started..."
    end if

    NumDims = Assembler%properties%nd
    NumGrids = Assembler%properties%ngrids
    Grids => Assembler%domain%grids

!     ! Fringe size is currently assumed to be constant for each receiver grid
!     allocate(FringeSize(NumGrids))
!     do n_ = 1, NumActiveGrids
!       n = Assembler%grid_id(n_)
!       FringeSize(n) = 0
!       do m_ = 1, NumActiveGrids
!         m = Assembler%grid_id(m_)
!         if (Assembler%graph%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
!           FringeSize(n) = Assembler%graph%fringe_size(m,n)
!           exit
!         end if
!       end do
!     end do

    ! Interp scheme is currently assumed to be constant for each receiver grid
    allocate(InterpScheme(NumGrids))
    do n = 1, NumGrids
      InterpScheme(n) = OVK_INTERP_LINEAR
      do m = 1, NumGrids
        if (Assembler%graph%interp_scheme(m,n) == OVK_INTERP_CUBIC) then
          InterpScheme(n) = OVK_INTERP_CUBIC
          exit
        end if
      end do
    end do

    ! Fringe size is currently assumed to be constant for each receiver grid
    allocate(FringeSize(NumGrids))
    do n = 1, NumGrids
      FringeSize(n) = 0
      do m = 1, NumGrids
        if (Assembler%graph%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
          FringeSize(n) = Assembler%graph%fringe_size(m,n)
          exit
        end if
      end do
    end do

    ! Fringe padding is currently assumed to be constant for each receiver grid
    allocate(FringePadding(NumGrids))
    do n = 1, NumGrids
      FringePadding(n) = 0
      do m = 1, NumGrids
        if (Assembler%graph%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
          FringePadding(n) = Assembler%graph%fringe_padding(m,n)
          exit
        end if
      end do
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Grid info:"
      write (*, '(3a)') "* Dimension: ", trim(IntToString(NumDims)), "D"
      write (*, '(2a)') "* Number of grids: ", trim(IntToString(NumGrids))
      NumPointsTotalString = LargeIntToString(sum([(ovkCartCount(Grids(m)%cart),m=1,NumGrids)]))
      write (*, '(2a)') "* Total number of grid points: ", trim(NumPointsTotalString)
      do m = 1, NumGrids
        NumPointsTotalString = LargeIntToString(ovkCartCount(Grids(m)%cart))
        iSString = IntToString(Grids(m)%cart%is(1))
        iEString = IntToString(Grids(m)%cart%ie(1))
        jSString = IntToString(Grids(m)%cart%is(2))
        jEString = IntToString(Grids(m)%cart%ie(2))
        kSString = IntToString(Grids(m)%cart%is(3))
        kEString = IntToString(Grids(m)%cart%ie(3))
        write (*, '(3a)', advance="no") "* Grid ", trim(IntToString(Grids(m)%properties%id)), ": "
        write (*, '(2a)', advance="no") trim(NumPointsTotalString), " points "
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

!     if (present(DebugMasks)) then
!       do n = 1, NumGrids
!         do m = 1, size(DebugMasks,1)
!           DebugMasks(m,n) = ovk_field_logical_(Grids(n)%cart, .false.)
!         end do
!       end do
!     end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Searching for candidate donor/receiver pairs..."
    end if

    allocate(OverlapBounds(NumGrids))
    allocate(MaxOverlapTolerance(NumGrids))
    do m = 1, NumGrids
      OverlapBounds(m) = ovk_bbox_(NumDims)
      MaxOverlapTolerance(m) = 0._rk
      do n = 1, NumGrids
        if (Assembler%graph%overlap(m,n)) then
          Bounds = ovkBBIntersect(Grids(m)%bounds, Grids(n)%bounds)
          OverlapBounds(m) = ovkBBUnion(OverlapBounds(m), Bounds)
          MaxOverlapTolerance(m) = max(MaxOverlapTolerance(m), &
            Assembler%graph%overlap_tolerance(m,n))
        end if
      end do
    end do

    allocate(PairwiseDonors(NumGrids,NumGrids))

    do m = 1, NumGrids
      if (OVK_VERBOSE) then
        write (*, '(3a)') "* Generating donor search accelerator on grid ", trim(IntToString(m)), &
          "..."
      end if
      call ovkGenerateDonorAccel(Grids(m), DonorAccel, Bounds=OverlapBounds(m), &
        OverlapTolerance=MaxOverlapTolerance(m))
      do n = 1, NumGrids
        if (Assembler%graph%overlap(m,n)) then
          call ovkFindDonors(Grids(m), Grids(n), DonorAccel, PairwiseDonors(m,n), &
            OverlapTolerance=Assembler%graph%overlap_tolerance(m,n))
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

    deallocate(OverlapBounds)
    deallocate(MaxOverlapTolerance)

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished searching for candidate donor/receiver pairs."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Cutting boundary holes..."
    end if

    do n = 1, NumGrids
      if (any(Assembler%graph%boundary_hole_cutting(:,n))) then
        BoundaryMask = Grids(n)%boundary_mask
        InteriorMask = ovk_field_logical_(Grids(n)%cart, .false.)
        do m = 1, NumGrids
          if (Assembler%graph%boundary_hole_cutting(m,n)) then
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
        do m = 1, NumGrids
          if (Assembler%graph%overlap(m,n)) then
            PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
              .not. ExteriorMask%values
          end if
        end do
        do m = 1, NumGrids
          if (Assembler%graph%overlap(n,m)) then
            call ovkGenerateReceiverMask(Grids(m), Grids(n), PairwiseDonors(n,m), ReceiverMask, &
              DonorSubset=ExteriorMask)
            PairwiseDonors(n,m)%valid_mask%values = PairwiseDonors(n,m)%valid_mask%values .and. &
              .not. ReceiverMask%values
          end if
        end do
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

    allocate(OuterReceiverMasks(NumGrids))
    allocate(InnerReceiverMasks(NumGrids))

    do n = 1, NumGrids
      call ovkPartitionReceivers(Grids(n), PairwiseDonors(:,n), OuterReceiverMasks(n), &
        InnerReceiverMasks(n), FringeSize(n))
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

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (Assembler%graph%overlap_hole_cutting(m,n) .and. &
          Assembler%graph%overlap_hole_cutting(n,m)) then
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
        else if (.not. Assembler%graph%overlap_hole_cutting(m,n)) then
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

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (Assembler%graph%connection_type(m,n) == OVK_CONNECTION_FRINGE .and. &
          Assembler%graph%disjoint_connection(m,n)) then
          call ovkGenerateReceiverMask(Grids(n), Grids(m), PairwiseDonors(m,n), &
            ReceiverMask, DonorSubset=OuterReceiverMasks(m))
          PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
            .not. ReceiverMask%values
        end if
      end do
    end do

    do n = 1, NumGrids
      call ovkChooseDonors(PairwiseDonors(:,n), Subset=OuterReceiverMasks(n))
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished choosing best donor/receiver pairs near grid edges."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Invalidating cyclic donors and applying fringe padding..."
    end if

    do n = 1, NumGrids
      do m = n+1, NumGrids
        Connection1 = Assembler%graph%connection_type(m,n)
        Connection2 = Assembler%graph%connection_type(n,m)
        Disjoint1 = Assembler%graph%disjoint_connection(m,n)
        Disjoint2 = Assembler%graph%disjoint_connection(n,m)
        if (Connection1 == OVK_CONNECTION_FRINGE .and. Connection2 == OVK_CONNECTION_FRINGE) then
          Padding1 = Assembler%graph%fringe_padding(m,n)
          Padding2 = Assembler%graph%fringe_padding(n,m)
          PaddingMax = max(Padding1, Padding2)
          PaddingSum = Padding1 + Padding2
          PaddingFrac = real(Padding2,kind=rk)/real(max(PaddingSum,1),kind=rk)
          WeightedPadding1 = int(real(PaddingMax,kind=rk) * PaddingFrac + 0.5_rk)
          WeightedPadding2 = PaddingMax - WeightedPadding1
          if (WeightedPadding1 > 0 .or. WeightedPadding2 > 0 .or. Disjoint1 .or. Disjoint2) then
            call ovkGenerateNearCrossoverMask(Grids(m), Grids(n), PairwiseDonors(n,m), &
              PairwiseDonors(m,n), WeightedPadding1, NearCrossoverMask1, Subset1=InnerReceiverMasks(m))
            call ovkGenerateNearCrossoverMask(Grids(n), Grids(m), PairwiseDonors(m,n), &
              PairwiseDonors(n,m), WeightedPadding2, NearCrossoverMask2, Subset1=InnerReceiverMasks(n))
            PairwiseDonors(n,m)%valid_mask%values = PairwiseDonors(n,m)%valid_mask%values .and. &
              .not. NearCrossoverMask1%values
            PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
              .not. NearCrossoverMask2%values
            if (OVK_VERBOSE) then
              NumInvalidatedDonors1 = ovkCountMask(NearCrossoverMask1)
              NumInvalidatedDonors2 = ovkCountMask(NearCrossoverMask2)
            end if
            call ovkGenerateNearCrossoverMask(Grids(m), Grids(n), PairwiseDonors(n,m), &
              PairwiseDonors(m,n), Padding1, NearCrossoverMask1, &
              Subset1=InnerReceiverMasks(m))
            call ovkGenerateNearCrossoverMask(Grids(n), Grids(m), PairwiseDonors(m,n), &
              PairwiseDonors(n,m), Padding2, NearCrossoverMask2, &
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

    do n = 1, NumGrids
      call ovkChooseDonors(PairwiseDonors(:,n), Subset=InnerReceiverMasks(n))
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished choosing best donor/receiver pairs on grid interior."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Merging remaining candidate donor/receiver pairs into final set..."
    end if

    Donors => Assembler%donors

    do m = 1, NumGrids
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

    do n = 1, NumGrids
      if (any(Assembler%graph%overlap_hole_cutting(:,n))) then
        call ovkGenerateOverlapOptimizationMask(Grids(n), Donors(n), OverlapOptimizationMask, &
          FringeSize=FringeSize(n))
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

    ! TODO: Need to update this to work with per-pair interp scheme parameter
    if (any(Assembler%graph%interp_scheme == OVK_INTERP_CUBIC)) then

      if (OVK_VERBOSE) then
        write (*, '(a)') "Expanding donor cells for cubic interpolation..."
      end if

      allocate(ValidCellMasks(NumGrids))
      allocate(ExpandableCellMasks(NumGrids))
      allocate(CellQualities(NumGrids))

      do m = 1, NumGrids
        call GenerateValidCellMask(Grids(m), ValidCellMasks(m))
        call ovkFindMaskEdge(ValidCellMasks(m), OVK_EDGE_TYPE_INNER, ValidCellInnerEdgeMask)
        ExpandableCellMasks(m) = ovk_field_logical_(Grids(m)%cell_cart)
        ExpandableCellMasks(m)%values = ValidCellMasks(m)%values .and. .not. &
          ValidCellInnerEdgeMask%values
        call CountNeighbors(ExpandableCellMasks(m), CellQualities(m))
      end do

      do n = 1, NumGrids
        if (all(Assembler%graph%interp_scheme(:,n) /= OVK_INTERP_CUBIC)) cycle
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
                  AwayFromEdge = ovkCartContains(Grids(m)%cell_cart, NeighborCellLower) .and. &
                    ovkCartContains(Grids(m)%cell_cart, NeighborCellUpper)
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

      do n = 1, NumGrids
        do m = 1, NumGrids
          if (Assembler%graph%connection_type(m,n) /= OVK_CONNECTION_NONE .and. &
            Assembler%graph%disjoint_connection(m,n)) then
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

    do n = 1, NumGrids
      call ovkGenerateOrphanMask(Donors(n), OuterReceiverMasks(n), OrphanMask)
!       if (present(OrphanMasks)) then
!         OrphanMasks(n) = OrphanMask
!       end if
      if (OVK_VERBOSE) then
        do k = Grids(n)%cart%is(3), Grids(n)%cart%ie(3)
          do j = Grids(n)%cart%is(2), Grids(n)%cart%ie(2)
            do i = Grids(n)%cart%is(1), Grids(n)%cart%ie(1)
              ReceiverPoint = [i,j,k]
              if (OrphanMask%values(i,j,k)) then
                write (ERROR_UNIT, '(4a)') "WARNING: Orphan detected at point ", &
                  trim(TupleToString(ReceiverPoint(:Grids(n)%cart%nd))), " of grid ", &
                  trim(IntToString(n))
              end if
            end do
          end do
        end do
        NumOrphans = ovkCountMask(OrphanMask)
        write (*, '(5a)') "* ", trim(LargeIntToString(NumOrphans)), &
          " orphan points found on grid ", trim(IntToString(n)), "."
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

    call ovkEditAssemblerConnectivity(Assembler, Connectivity)

    do n = 1, NumGrids
      StencilSize = 2
      do m = 1, NumGrids
        if (Assembler%graph%connection_type(m,n) /= OVK_CONNECTION_NONE .and. &
          Assembler%graph%interp_scheme(m,n) == OVK_INTERP_CUBIC) then
          StencilSize = 4
          exit
        end if
      end do
      call ovkCreateConnectivityInterpData(Connectivity, n, Grids(n), StencilSize)
    end do

    InterpData => Connectivity%interp_data

    do n = 1, NumGrids
      call ovkFillInterpData(InterpData(n), Donors(n), InterpScheme=InterpScheme(n))
      if (OVK_VERBOSE) then
        write (*, '(3a)') "* Generated interpolation data for grid ", trim(IntToString(n)), "."
      end if
    end do

    call ovkReleaseAssemblerConnectivity(Assembler, Connectivity)

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished generating interpolation data."
    end if

    if (OVK_VERBOSE) then
      call system_clock(ClockFinal, ClockRate)
      write (*, '(a,f0.3,a)') "Overset grid assembly finished (time: ", &
        real(ClockFinal-ClockInitial,kind=rk)/real(ClockRate,kind=rk), " seconds)."
    end if

  end subroutine ovkAssemble

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

    select case (DonorGrid%cart%nd)
    case (2)
      select case (DonorGrid%properties%geometry_type)
      case (OVK_GRID_GEOMETRY_CARTESIAN)
        ExpandedDonorCellCoords = ovkRectangleIsoInverseCubic(VertexCoords, ReceiverCoords)
      case default
        ExpandedDonorCellCoords = ovkQuadIsoInverseCubic(VertexCoords, ReceiverCoords, &
          Guess=DonorCellCoords, Success=Success)
      end select
    case (3)
      select case (DonorGrid%properties%geometry_type)
      case (OVK_GRID_GEOMETRY_CARTESIAN)
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
          trim(TupleToString(DonorCell)), " of grid ", trim(IntToString(DonorGrid%properties%id))
      end if
    end if

  end subroutine ExpandDonorCell

end module ovkOverset
