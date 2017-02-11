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

  subroutine ovkAssembleOverset(Grids, InterpData, FringeSize, DisjointFringes, FringePadding, &
    InterpScheme, AllowInterpolation, OptimizeOverlap, OverlapTolerance, HoleMasks, OrphanMasks)

    type(ovk_grid), dimension(:), intent(inout) :: Grids
    type(ovk_interp), dimension(size(Grids)), intent(out) :: InterpData
    integer, dimension(size(Grids)), intent(in), optional :: FringeSize
    logical, dimension(size(Grids),size(Grids)), intent(in), optional :: DisjointFringes
    integer, dimension(size(Grids),size(Grids)), intent(in), optional :: FringePadding
    integer, dimension(size(Grids)), intent(in), optional :: InterpScheme
    logical, dimension(size(Grids),size(Grids)), intent(in), optional :: AllowInterpolation
    logical, dimension(size(Grids)), intent(in), optional :: OptimizeOverlap
    real(rk), dimension(size(Grids)), intent(in), optional :: OverlapTolerance
    type(ovk_field_logical), dimension(size(Grids)), intent(out), optional :: HoleMasks
    type(ovk_field_logical), dimension(size(Grids)), intent(out), optional :: OrphanMasks

    integer, dimension(:), allocatable :: FringeSize_
    logical, dimension(:,:), allocatable :: DisjointFringes_
    integer, dimension(:,:), allocatable :: FringePadding_
    integer, dimension(:), allocatable :: InterpScheme_
    logical, dimension(:,:), allocatable :: AllowInterpolation_
    logical, dimension(:), allocatable :: OptimizeOverlap_
    real(rk), dimension(:), allocatable :: OverlapTolerance_
    integer :: i, j, k, l, m, n, p, q, r
    integer :: ClockInitial, ClockFinal, ClockRate
    character(len=STRING_LENGTH) :: nPointsTotalString
    character(len=STRING_LENGTH) :: iSString, iEString, jSString, jEString, kSString, kEString
    type(ovk_bbox) :: Bounds
    type(ovk_bbox), dimension(:), allocatable :: OverlapBounds
    integer :: PaddingMax, PaddingSum
    real(rk) :: PaddingFrac
    integer :: Padding1, Padding2
    integer(lk) :: nReceivers
    integer(lk) :: nOuterReceivers, nInnerReceivers
    integer(lk) :: nInvalidatedDonors, nInvalidatedDonors1, nInvalidatedDonors2
    integer(lk) :: nOrphans
    integer(lk) :: nRemovedPoints
    type(ovk_donor_accel) :: DonorAccel
    type(ovk_donors), dimension(:,:), allocatable :: PairwiseDonors
    type(ovk_donors), dimension(:), allocatable :: Donors
    type(ovk_field_logical), dimension(:), allocatable :: OuterReceiverMasks
    type(ovk_field_logical), dimension(:), allocatable :: InnerReceiverMasks
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

    allocate(InterpScheme_(size(Grids)))
    if (present(InterpScheme)) then
      InterpScheme_ = InterpScheme
    else
      InterpScheme_ = OVK_INTERP_LINEAR
    end if

    allocate(AllowInterpolation_(size(Grids),size(Grids)))
    if (present(AllowInterpolation)) then
      AllowInterpolation_ = AllowInterpolation
    else
      AllowInterpolation_ = .true.
    end if

    allocate(OptimizeOverlap_(size(Grids)))
    if (present(OptimizeOverlap)) then
      OptimizeOverlap_ = OptimizeOverlap
    else
      OptimizeOverlap_ = .true.
    end if

    allocate(OverlapTolerance_(size(Grids)))
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
      nPointsTotalString = LargeIntToString(sum([(ovkCartCount(Grids(m)%cart),m=1,size(Grids))]))
      write (*, '(2a)') "* Total number of grid points: ", trim(nPointsTotalString)
      do m = 1, size(Grids)
        nPointsTotalString = LargeIntToString(ovkCartCount(Grids(m)%cart))
        iSString = IntToString(Grids(m)%cart%is(1))
        iEString = IntToString(Grids(m)%cart%ie(1))
        jSString = IntToString(Grids(m)%cart%is(2))
        jEString = IntToString(Grids(m)%cart%ie(2))
        kSString = IntToString(Grids(m)%cart%is(3))
        kEString = IntToString(Grids(m)%cart%ie(3))
        write (*, '(3a)', advance="no") "* Grid ", trim(IntToString(Grids(m)%id)), ": "
        write (*, '(2a)', advance="no") trim(nPointsTotalString), " points "
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

    do m = 1, size(Grids)
      AllowInterpolation_(m,m) = .false.
    end do

    allocate(OverlapBounds(size(Grids)))
    OverlapBounds = ovk_bbox_(Grids(1)%cart%nd)
    do n = 1, size(Grids)
      do m = 1, size(Grids)
        if (AllowInterpolation_(m,n)) then
          Bounds = ovkBBIntersect(Grids(m)%bounds, Grids(n)%bounds)
          OverlapBounds(n) = ovkBBUnion(OverlapBounds(n), Bounds)
          if (ovkBBIsEmpty(Bounds)) then
            AllowInterpolation_(m,n) = .false.
          end if
        end if
      end do
    end do

    allocate(PairwiseDonors(size(Grids),size(Grids)))

    do n = 1, size(Grids)
      if (OVK_VERBOSE) then
        write (*, '(3a)') "* Generating donor search accelerator on grid ", trim(IntToString(n)), &
          "..."
      end if
      call ovkGenerateDonorAccel(Grids(n), DonorAccel, Bounds=OverlapBounds(n), &
        OverlapTolerance=OverlapTolerance_(n))
      do m = 1, size(Grids)
        if (AllowInterpolation_(n,m)) then
          call FindDonors(Grids(n), Grids(m), DonorAccel, PairwiseDonors(n,m))
          if (OVK_VERBOSE) then
            nReceivers = ovkCountMask(PairwiseDonors(n,m)%valid_mask)
            write (*, '(7a)') "* ", trim(LargeIntToString(nReceivers)), &
              " candidate donors from grid ", trim(IntToString(n)), " to grid ", &
              trim(IntToString(m)), " found."
          end if
        else
          call ovkMakeDonors(PairwiseDonors(n,m), Grids(m)%cart)
          PairwiseDonors(n,m)%valid_mask%values = .false.
        end if
      end do
      call ovkDestroyDonorAccel(DonorAccel)
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished searching for candidate donor/receiver pairs."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Partitioning receivers..."
    end if

    allocate(OuterReceiverMasks(size(Grids)))
    allocate(InnerReceiverMasks(size(Grids)))

    do m = 1, size(Grids)
      call ovkPartitionReceivers(Grids(m), PairwiseDonors(:,m), OuterReceiverMasks(m), &
        InnerReceiverMasks(m), FringeSize_(m))
      if (OVK_VERBOSE) then
        nOuterReceivers = ovkCountMask(OuterReceiverMasks(m))
        nInnerReceivers = ovkCountMask(InnerReceiverMasks(m))
        write (*, '(7a)') "* Partitioned receivers on grid ", trim(IntToString(m)), " into ", &
          trim(LargeIntToString(nOuterReceivers)), " near-edge points and ", &
          trim(LargeIntToString(nInnerReceivers)), " interior points."
      end if
    end do

    do n = 1, size(Grids)
      if (.not. OptimizeOverlap_(n)) then
        do m = 1, size(Grids)
          if (AllowInterpolation_(m,n)) then
            PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
              .not. InnerReceiverMasks(n)%values
          end if
        end do
        InnerReceiverMasks(n)%values = .false.
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
        if (AllowInterpolation_(m,n) .and. OptimizeOverlap_(m)) then
          call ovkGenerateCoarseToFineMask(Grids(n), PairwiseDonors(m,n), CoarseToFineMask, &
            Subset=InnerReceiverMasks(n))
          PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
            .not. CoarseToFineMask%values
          if (OVK_VERBOSE) then
            nInvalidatedDonors = ovkCountMask(CoarseToFineMask)
            write (*, '(7a)') "* ", trim(LargeIntToString(nInvalidatedDonors)), &
              " donors from grid ", trim(IntToString(m)), " to grid ", trim(IntToString(n)), &
              " invalidated."
          end if
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
        if (AllowInterpolation_(m,n) .and. AllowInterpolation_(n,m)) then
          if (DisjointFringes_(m,n) .or. DisjointFringes_(n,m)) then
            call ovkGenerateReceiverMask(Grids(n), Grids(m), PairwiseDonors(m,n), &
              ReceiverMask, DonorSubset=OuterReceiverMasks(m))
            PairwiseDonors(m,n)%valid_mask%values = PairwiseDonors(m,n)%valid_mask%values .and. &
              .not. ReceiverMask%values
          end if
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
              nInvalidatedDonors1 = ovkCountMask(NearCrossoverMask1)
              nInvalidatedDonors2 = ovkCountMask(NearCrossoverMask2)
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
              nInvalidatedDonors1 = nInvalidatedDonors1 + ovkCountMask(NearCrossoverMask1)
              nInvalidatedDonors2 = nInvalidatedDonors2 + ovkCountMask(NearCrossoverMask2)
            end if
          else
            if (OVK_VERBOSE) then
              nInvalidatedDonors1 = 0
              nInvalidatedDonors2 = 0
            end if
          end if
          if (OVK_VERBOSE) then
            write (*, '(7a)') "* ", trim(LargeIntToString(nInvalidatedDonors1)), &
              " donor/receiver pairs from grid ", trim(IntToString(n)), " to ", &
              trim(IntToString(m)), " invalidated."
            write (*, '(7a)') "* ", trim(LargeIntToString(nInvalidatedDonors2)), &
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
        nReceivers = ovkCountMask(Donors(m)%valid_mask)
        write (*, '(5a)') "* ", trim(LargeIntToString(nReceivers)), &
          " receiver points on grid ", trim(IntToString(m)), "."
      end if
    end do

    deallocate(PairwiseDonors)

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished merging candidate donor/receiver pairs."
    end if

    if (OVK_VERBOSE) then
      write (*, '(a)') "Cutting holes and optimizing overlap..."
    end if

    do m = 1, size(Grids)
      if (OptimizeOverlap_(m)) then
        call ovkGenerateOverlapOptimizationMask(Grids(m), Donors(m), OverlapOptimizationMask, &
          FringeSize=FringeSize_(m))
        do k = Grids(m)%cart%is(3), Grids(m)%cart%ie(3)
          do j = Grids(m)%cart%is(2), Grids(m)%cart%ie(2)
            do i = Grids(m)%cart%is(1), Grids(m)%cart%ie(1)
              if (OverlapOptimizationMask%values(i,j,k)) then
                Grids(m)%grid_mask%values(i,j,k) = .false.
                Grids(m)%boundary_mask%values(i,j,k) = .false.
                Donors(m)%valid_mask%values(i,j,k) = .false.
                OuterReceiverMasks(m)%values(i,j,k) = .false.
                InnerReceiverMasks(m)%values(i,j,k) = .false.
              end if
            end do
          end do
        end do
        if (present(HoleMasks)) then
          HoleMasks(m) = OverlapOptimizationMask
        end if
        if (OVK_VERBOSE) then
          nRemovedPoints = ovkCountMask(OverlapOptimizationMask)
          write (*, '(5a)') "* ", trim(LargeIntToString(nRemovedPoints)), &
            " points removed from grid ", trim(IntToString(m)), "."
        end if
      else
        if (present(HoleMasks)) then
          HoleMasks(m) = ovk_field_logical_(Grids(m)%cart, .false.)
        end if
      end if
    end do

    if (OVK_VERBOSE) then
      write (*, '(a)') "Finished cutting holes and optimizing overlap."
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

      do m = 1, size(Grids)
        if (InterpScheme_(m) /= OVK_INTERP_CUBIC) cycle
        do k = Grids(m)%cart%is(3), Grids(m)%cart%ie(3)
          do j = Grids(m)%cart%is(2), Grids(m)%cart%ie(2)
            do i = Grids(m)%cart%is(1), Grids(m)%cart%ie(1)
              ReceiverPoint = [i,j,k]
              if (Donors(m)%valid_mask%values(i,j,k)) then
                n = Donors(m)%grid_ids%values(i,j,k)
                do l = 1, Grids(n)%cart%nd
                  DonorCell(l) = Donors(m)%cells(l)%values(i,j,k)
                  DonorCellCoords(l) = Donors(m)%cell_coords(l)%values(i,j,k)
                end do
                DonorCell(Grids(n)%cell_cart%nd+1:) = 1
                if (ExpandableCellMasks(n)%values(DonorCell(1),DonorCell(2),DonorCell(3))) then
                  ExpandableCell = .true.
                else
                  BestCellQuality = 0
                  DonorCellShift = 0
                  NeighborCellLower(:Grids(n)%cart%nd) = DonorCell(:Grids(n)%cart%nd)-1
                  NeighborCellLower(Grids(n)%cart%nd+1:) = DonorCell(Grids(n)%cart%nd+1:)
                  NeighborCellUpper(:Grids(n)%cart%nd) = DonorCell(:Grids(n)%cart%nd)+1
                  NeighborCellUpper(Grids(n)%cart%nd+1:) = DonorCell(Grids(n)%cart%nd+1:)
                  AwayFromBoundary = ovkCartContains(Grids(n)%cell_cart, NeighborCellLower) .and. &
                    ovkCartContains(Grids(n)%cell_cart, NeighborCellUpper)
                  if (AwayFromBoundary) then
                    do r = NeighborCellLower(3), NeighborCellUpper(3)
                      do q = NeighborCellLower(2), NeighborCellUpper(2)
                        do p = NeighborCellLower(1), NeighborCellUpper(1)
                          NeighborCell = [p,q,r]
                          if (ExpandableCellMasks(n)%values(NeighborCell(1),NeighborCell(2), &
                            NeighborCell(3))) then
                            CellQuality = CellQualities(n)%values(NeighborCell(1),NeighborCell(2),NeighborCell(3))
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
                          NeighborCell(:Grids(n)%cell_cart%nd) = ovkCartPeriodicAdjust( &
                            Grids(n)%cell_cart, NeighborCell)
                          if (ovkCartContains(Grids(n)%cell_cart, NeighborCell)) then
                            if (ExpandableCellMasks(n)%values(NeighborCell(1),NeighborCell(2), &
                              NeighborCell(3))) then
                              CellQuality = CellQualities(n)%values(NeighborCell(1),NeighborCell(2),NeighborCell(3))
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
                    DonorCellCoords = DonorCellCoords - real(DonorCellShift(:Grids(n)%cell_cart%nd),kind=rk)
                  else
                    ExpandableCell = .false.
                    if (OVK_VERBOSE) then
                      write (ERROR_UNIT, '(8a)') "WARNING: Could not use cubic interpolation for donor cell ", &
                        trim(TupleToString(DonorCell(:Grids(n)%cart%nd))), " of grid ", trim(IntToString(n)), &
                        " corresponding to receiver point ", trim(TupleToString(ReceiverPoint(:Grids(m)%cart%nd))), &
                        " of grid ", trim(IntToString(m))
                    end if
                  end if
                end if
                if (ExpandableCell) then
                  do l = 1, Grids(m)%cart%nd
                    ReceiverCoords(l) = Grids(m)%xyz(l)%values(ReceiverPoint(1),ReceiverPoint(2), &
                      ReceiverPoint(3))
                  end do
                  call ExpandDonorCell(Grids(n), DonorCell, DonorCellCoords, &
                    ReceiverCoords, ExpandedDonorCell, ExpandedDonorCellCoords)
                  ExpandedDonorCell(Grids(n)%cart%nd+1:) = 1
                  do l = 1, Grids(n)%cart%nd
                    Donors(m)%cells(l)%values(i,j,k) = ExpandedDonorCell(l)
                    Donors(m)%cell_coords(l)%values(i,j,k) = ExpandedDonorCellCoords(l)
                  end do
                  Donors(m)%cell_extents%values(i,j,k) = 4
                end if
              end if
            end do
          end do
        end do
      end do

      do n = 1, size(Grids)
        do m = 1, size(Grids)
          if (AllowInterpolation_(m,n) .and. AllowInterpolation_(n,m)) then
            if (DisjointFringes_(m,n) .or. DisjointFringes_(n,m)) then
              call ovkGenerateReceiverMask(Grids(n), Grids(m), Donors(n), ReceiverMask, &
                DonorSubset=Donors(m)%valid_mask)
              Donors(n)%valid_mask%values = Donors(n)%valid_mask%values .and. .not. &
                ReceiverMask%values
            end if
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
        nOrphans = ovkCountMask(OrphanMask)
        write (*, '(5a)') "* ", trim(LargeIntToString(nOrphans)), &
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
    type(ovk_field_logical) :: BoundaryOuterEdgeMask
    type(ovk_field_logical) :: NonBoundaryMask
    type(ovk_field_logical) :: NonBoundaryOuterEdgeMask

    call ovkFindMaskEdge(Grid%grid_mask, OVK_EDGE_TYPE_OUTER, GridOuterEdgeMask)

    NonBoundaryMask = ovk_field_logical_(Grid%cart)
    NonBoundaryMask%values = .not. Grid%boundary_mask%values

    call ovkFindMaskEdge(Grid%boundary_mask, OVK_EDGE_TYPE_OUTER, BoundaryOuterEdgeMask)
    call ovkFindMaskEdge(NonBoundaryMask, OVK_EDGE_TYPE_OUTER, NonBoundaryOuterEdgeMask)

    GridOuterEdgeMask%values = GridOuterEdgeMask%values .and. .not. &
      (BoundaryOuterEdgeMask%values .and. .not. NonBoundaryOuterEdgeMask%values)

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
    integer :: i, j, k

    if (present(FringeSize)) then
      FringeSize_ = FringeSize
    else
      FringeSize_ = 2
    end if

    OverlapOptimizationMask = ovk_field_logical_(Grid%cart)

    do k = Grid%cart%is(3), Grid%cart%ie(3)
      do j = Grid%cart%is(2), Grid%cart%ie(2)
        do i = Grid%cart%is(1), Grid%cart%ie(1)
          OverlapOptimizationMask%values(i,j,k) = .not. Grid%grid_mask%values(i,j,k) .or. &
            Donors%valid_mask%values(i,j,k)
        end do
      end do
    end do

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
