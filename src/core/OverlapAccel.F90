! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkOverlapAccel

  use ovkBoundingBox
  use ovkCart
  use ovkField
  use ovkFieldOps
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkHashGrid
  use ovkLogger
  implicit none

  private

  ! Internal API
  public :: t_overlap_accel
  public :: t_overlap_accel_
  public :: CreateOverlapAccel
  public :: DestroyOverlapAccel
  public :: PopulateOverlapAccel
  public :: FindOverlappingCell

  type t_node
    real(rk) :: split
    integer :: split_dir
    type(t_node), pointer :: left_child
    type(t_node), pointer :: right_child
    type(t_hash_grid), pointer :: hash_grid
  end type t_node

  type t_overlap_accel
    type(t_noconstruct) :: noconstruct
    type(ovk_grid), pointer :: grid
    integer :: nd
    type(ovk_bbox) :: bounds
    type(t_node), pointer :: root
  end type t_overlap_accel

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_overlap_accel_
    module procedure t_overlap_accel_Default
  end interface t_overlap_accel_

contains

  pure function t_overlap_accel_Default() result(Accel)

    type(t_overlap_accel) :: Accel

    nullify(Accel%grid)
    Accel%nd = 2
    Accel%bounds = ovk_bbox_()
    nullify(Accel%root)

  end function t_overlap_accel_Default

  subroutine CreateOverlapAccel(Accel, Grid)

    type(t_overlap_accel), intent(out) :: Accel
    type(ovk_grid), pointer, intent(in) :: Grid

    Accel%grid => Grid
    Accel%nd = Grid%nd
    Accel%bounds = ovk_bbox_(Grid%nd)
    nullify(Accel%root)

  end subroutine CreateOverlapAccel

  subroutine DestroyOverlapAccel(Accel)

    type(t_overlap_accel), intent(inout) :: Accel

    if (associated(Accel%root)) then
      call DestroyOverlapAccelNode(Accel%root)
      deallocate(Accel%root)
    end if

  end subroutine DestroyOverlapAccel

  subroutine PopulateOverlapAccel(Accel, Bounds, MaxOverlapTolerance, NumCellsLeaf, &
    MaxNodeUnoccupiedVolume, MaxNodeCellVolumeVariation, BinScale)

    type(t_overlap_accel), intent(inout) :: Accel
    type(ovk_bbox), intent(in) :: Bounds
    real(rk), intent(in) :: MaxOverlapTolerance
    integer(lk), intent(in) :: NumCellsLeaf
    real(rk), intent(in) :: MaxNodeUnoccupiedVolume
    real(rk), intent(in) :: MaxNodeCellVolumeVariation
    real(rk), intent(in) :: BinScale

    integer :: d, i, j, k
    integer :: NumDims
    type(ovk_grid), pointer :: Grid
    type(t_logger) :: Logger
    type(ovk_field_logical) :: GridCellOverlapMask
    type(ovk_field_real), dimension(:), allocatable :: GridCellLower, GridCellUpper
    real(rk), dimension(MAX_ND) :: AccelLower, AccelUpper
    integer, dimension(MAX_ND) :: Cell
    type(ovk_bbox) :: GridCellBounds
    integer(lk) :: NumOverlappingCells, iNextOverlappingCell
    integer, dimension(:,:), allocatable :: OverlappingCells
    integer(lk), dimension(:), allocatable :: CellIndices
    integer :: MaxDepth
    real(rk), dimension(Accel%nd) :: BoundsSize
    real(rk) :: BoundsVolume

    real(rk), parameter :: TOLERANCE = 1.e-12_rk

    NumDims = Accel%nd
    Grid => Accel%grid
    Logger = Grid%logger

    if (associated(Accel%root)) then
      call DestroyOverlapAccelNode(Accel%root)
      deallocate(Accel%root)
    end if

    GridCellOverlapMask = ovk_field_logical_(Grid%cell_cart, .false.)

    allocate(GridCellLower(NumDims))
    allocate(GridCellUpper(NumDims))
    do d = 1, NumDims
      GridCellLower(d) = ovk_field_real_(Grid%cell_cart, 0._rk)
      GridCellUpper(d) = ovk_field_real_(Grid%cell_cart, 0._rk)
    end do

    AccelLower = huge(0._rk)
    AccelUpper = -huge(0._rk)

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  FIRSTPRIVATE(Bounds, MaxOverlapTolerance) &
!$OMP&  SHARED(Grid, GridCellOverlapMask, GridCellLower, GridCellUpper) &
!$OMP&  REDUCTION(min:AccelLower) &
!$OMP&  REDUCTION(max:AccelUpper)
    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          Cell = [i,j,k]
          if (Grid%cell_mask%values(Cell(1),Cell(2),Cell(3))) then
            GridCellBounds = ovkGridCellBounds(Grid, Cell)
            GridCellBounds = ovkBBScale(GridCellBounds, 1._rk + 2._rk*MaxOverlapTolerance)
            GridCellOverlapMask%values(i,j,k) = ovkBBOverlaps(Bounds, GridCellBounds)
            if (GridCellOverlapMask%values(i,j,k)) then
              do d = 1, NumDims
                GridCellLower(d)%values(i,j,k) = GridCellBounds%b(d)
                GridCellUpper(d)%values(i,j,k) = GridCellBounds%e(d)
              end do
              AccelLower = min(AccelLower, GridCellBounds%b)
              AccelUpper = max(AccelUpper, GridCellBounds%e)
            end if
          end if
        end do
      end do
    end do
!$OMP END PARALLEL DO

    Accel%bounds = ovk_bbox_(NumDims, AccelLower, AccelUpper)

    NumOverlappingCells = ovkCountMask(GridCellOverlapMask)

    if (NumOverlappingCells > 0_lk) then

      allocate(OverlappingCells(NumOverlappingCells,MAX_ND))
      allocate(CellIndices(NumOverlappingCells))

      iNextOverlappingCell = 1_lk
      do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
        do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
          do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
            Cell = [i,j,k]
            if (GridCellOverlapMask%values(i,j,k)) then
              OverlappingCells(iNextOverlappingCell,:) = Cell
              CellIndices(iNextOverlappingCell) = iNextOverlappingCell
              iNextOverlappingCell = iNextOverlappingCell + 1_lk
            end if
          end do
        end do
      end do

      MaxDepth = 0
      do while (ishft(NumOverlappingCells/NumCellsLeaf, -MaxDepth) /= 0)
        MaxDepth = MaxDepth + 1
      end do

      BoundsSize = ovkBBSize(Accel%bounds)
      BoundsVolume = product(BoundsSize)

      allocate(Accel%root)
      call GenerateOverlapAccelNode(Accel%root, Grid%cell_cart, Accel%bounds, BoundsVolume, &
        GridCellLower, GridCellUpper, Grid%cell_volumes, OverlappingCells, CellIndices, 0, &
        MaxDepth, NumCellsLeaf, MaxNodeUnoccupiedVolume, MaxNodeCellVolumeVariation, BinScale)

    end if

    if (Logger%log_status) then
      call PrintStats(Accel)
    end if

  end subroutine PopulateOverlapAccel

  recursive subroutine GenerateOverlapAccelNode(Node, CellCart, Bounds, BoundsVolume, &
    AllCellLowers, AllCellUppers, AllCellVolumes, OverlappingCells, CellIndices, Depth, MaxDepth, &
    NumCellsLeaf, MaxUnoccupiedVolume, MaxCellVolumeVariation, BinScale)

    type(t_node), intent(out) :: Node
    type(ovk_cart), intent(in) :: CellCart
    type(ovk_bbox), intent(in) :: Bounds
    real(rk), intent(in) :: BoundsVolume
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellLowers, AllCellUppers
    type(ovk_field_real), intent(in) :: AllCellVolumes
    integer, dimension(:,:), intent(in) :: OverlappingCells
    integer(lk), dimension(:), intent(in) :: CellIndices
    integer, intent(in) :: Depth
    integer, intent(in) :: MaxDepth
    integer(lk), intent(in) :: NumCellsLeaf
    real(rk), intent(in) :: MaxUnoccupiedVolume
    real(rk), intent(in) :: MaxCellVolumeVariation
    real(rk), intent(in) :: BinScale

    integer :: i, j, k, d
    integer(lk) :: l
    integer(lk) :: NumCells
    integer, dimension(MAX_ND) :: Cell
    real(rk), dimension(MAX_ND) :: CellLower, CellUpper
    real(rk), dimension(MAX_ND) :: NodeBoundsLower, NodeBoundsUpper
    type(ovk_bbox) :: NodeBounds
    real(rk), dimension(CellCart%nd) :: NodeBoundsSize
    real(rk) :: NodeBoundsVolume
    real(rk) :: CellVolume, MeanCellVolume, CellVolumeVariation
    real(rk), dimension(CellCart%nd) :: CellBoundsSize, MeanCellBoundsSize
    real(rk) :: UnoccupiedVolume
    logical :: OccupiedEnough, UniformEnough, NotTooBig
    logical :: LeafNode
    integer :: SplitDir
    real(rk) :: Split
    integer(lk), dimension(:), allocatable :: LeftCellIndices, RightCellIndices
    real(rk), dimension(CellCart%nd) :: BinSize
    integer, dimension(CellCart%nd) :: NumBins
    type(ovk_cart) :: BinCart
    type(ovk_field_large_int) :: BinCellsStart, BinCellsEnd
    integer(lk), dimension(:), allocatable :: BinCells
    type(ovk_field_int) :: NumBinElements
    integer(lk) :: BinStart, BinEnd

    integer(lk), parameter :: MaxHashGridSize = 2_lk**26

    nullify(Node%hash_grid)
    nullify(Node%left_child)
    nullify(Node%right_child)

    NumCells = size(CellIndices,kind=lk)

    NodeBoundsLower = Bounds%e
    NodeBoundsUpper = Bounds%b

! This isn't working with OpenMP for some reason
! !$OMP PARALLEL DO &
! !$OMP&  DEFAULT(PRIVATE) &
! !$OMP&  SHARED(CellIndices, OverlappingCells, AllCellLowers, AllCellUppers) &
! !$OMP&  REDUCTION(min:NodeBoundsLower) &
! !$OMP&  REDUCTION(max:NodeBoundsUpper)
    do l = 1_lk, NumCells
      Cell = OverlappingCells(CellIndices(l),:)
      do d = 1, CellCart%nd
        CellLower(d) = AllCellLowers(d)%values(Cell(1),Cell(2),Cell(3))
        CellUpper(d) = AllCellUppers(d)%values(Cell(1),Cell(2),Cell(3))
      end do
      CellLower(CellCart%nd+1:) = 0._rk
      CellUpper(CellCart%nd+1:) = 0._rk
      NodeBoundsLower = min(NodeBoundsLower, CellLower)
      NodeBoundsUpper = max(NodeBoundsUpper, CellUpper)
    end do
! !$OMP END PARALLEL DO

    NodeBounds = ovk_bbox_(CellCart%nd, NodeBoundsLower, NodeBoundsUpper)
    NodeBoundsSize = ovkBBSize(NodeBounds)
    NodeBoundsVolume = product(NodeBoundsSize)

    if (Depth == MaxDepth .or. NumCells <= NumCellsLeaf) then

      LeafNode = .true.

    else

      MeanCellVolume = 0._rk
      do l = 1_lk, NumCells
        Cell = OverlappingCells(CellIndices(l),:)
        CellVolume = AllCellVolumes%values(Cell(1),Cell(2),Cell(3))
        MeanCellVolume = MeanCellVolume + CellVolume
      end do
      MeanCellVolume = MeanCellVolume/real(NumCells,kind=rk)

      CellVolumeVariation = 0._rk
      UnoccupiedVolume = NodeBoundsVolume
      do l = 1_lk, NumCells
        Cell = OverlappingCells(CellIndices(l),:)
        CellVolume = AllCellVolumes%values(Cell(1),Cell(2),Cell(3))
        CellVolumeVariation = CellVolumeVariation + (CellVolume - MeanCellVolume)**2
        UnoccupiedVolume = UnoccupiedVolume - CellVolume
      end do
      UnoccupiedVolume = max(UnoccupiedVolume,0._rk)
      CellVolumeVariation = sqrt(CellVolumeVariation/real(NumCells-1_lk,kind=rk))/MeanCellVolume

      OccupiedEnough = UnoccupiedVolume < MaxUnoccupiedVolume*BoundsVolume
      UniformEnough = CellVolumeVariation <= MaxCellVolumeVariation
      NotTooBig = NumCells <= MaxHashGridSize

      LeafNode = OccupiedEnough .and. UniformEnough .and. NotTooBig

    end if

    if (.not. LeafNode) then

      call FindSplit(CellCart, NodeBounds, AllCellLowers, AllCellUppers, OverlappingCells, &
        CellIndices, SplitDir, Split)

      Node%split_dir = SplitDir
      Node%split = Split

      allocate(Node%left_child)
      allocate(Node%right_child)

      call DivideCells(CellCart, AllCellLowers, AllCellUppers, OverlappingCells, CellIndices, &
        SplitDir, Split, LeftCellIndices, RightCellIndices)

      call GenerateOverlapAccelNode(Node%left_child, CellCart, Bounds, BoundsVolume, &
        AllCellLowers, AllCellUppers, AllCellVolumes, OverlappingCells, LeftCellIndices, Depth+1, &
        MaxDepth, NumCellsLeaf, MaxUnoccupiedVolume, MaxCellVolumeVariation, BinScale)
      call GenerateOverlapAccelNode(Node%right_child, CellCart, Bounds, BoundsVolume, &
        AllCellLowers, AllCellUppers, AllCellVolumes, OverlappingCells, RightCellIndices, Depth+1, &
        MaxDepth, NumCellsLeaf, MaxUnoccupiedVolume, MaxCellVolumeVariation, BinScale)

    else

      MeanCellBoundsSize = 0._rk
! This isn't working with OpenMP for some reason
! !$OMP PARALLEL DO &
! !$OMP&  DEFAULT(PRIVATE) &
! !$OMP&  SHARED(CellIndices, OverlappingCells, AllCellLowers, AllCellUppers) &
! !$OMP&  REDUCTION(+:MeanCellBoundsSize)
      do l = 1_lk, NumCells
        Cell = OverlappingCells(CellIndices(l),:)
        do d = 1, CellCart%nd
          CellBoundsSize(d) = max(AllCellUppers(d)%values(Cell(1),Cell(2),Cell(3)) - &
            AllCellLowers(d)%values(Cell(1),Cell(2),Cell(3)), 0._rk)
        end do
        MeanCellBoundsSize = MeanCellBoundsSize + CellBoundsSize
      end do
! !$OMP END PARALLEL DO
      MeanCellBoundsSize = MeanCellBoundsSize/real(max(NumCells,1_lk),kind=rk)

      BinSize = BinScale * MeanCellBoundsSize

      NumBins = merge(max(int(ceiling(NodeBoundsSize/BinSize)),1), 1, NumCells > 1_lk)

      BinCart = ovk_cart_(CellCart%nd, NumBins)

      call DistributeCellsToBins(BinCart, NodeBounds, CellCart, AllCellLowers, AllCellUppers, &
        OverlappingCells, CellIndices, BinCellsStart, BinCellsEnd, BinCells)

      NumBinElements = ovk_field_int_(BinCart)
      NumBinElements%values = int(BinCellsEnd%values - BinCellsStart%values + 1_lk)

      allocate(Node%hash_grid)
      Node%hash_grid = t_hash_grid_(BinCart, NodeBounds, NumBinElements)

      l = 1_lk
      do k = BinCart%is(3), BinCart%ie(3)
        do j = BinCart%is(2), BinCart%ie(2)
          do i = BinCart%is(1), BinCart%ie(1)
            BinStart = Node%hash_grid%bin_start(l)
            BinEnd = Node%hash_grid%bin_start(l+1_lk)-1_lk
            if (BinEnd >= BinStart) then
              Node%hash_grid%bin_contents(BinStart:BinEnd) = BinCells(BinCellsStart%values(i,j,k): &
                BinCellsEnd%values(i,j,k))
            end if
            l = l + 1_lk
          end do
        end do
      end do

    end if

  end subroutine GenerateOverlapAccelNode

  subroutine FindSplit(CellCart, Bounds, AllCellLowers, AllCellUppers, OverlappingCells, &
    CellIndices, SplitDir, Split)

    type(ovk_cart), intent(in) :: CellCart
    type(ovk_bbox), intent(in) :: Bounds
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellLowers, AllCellUppers
    integer, dimension(:,:), intent(in) :: OverlappingCells
    integer(lk), dimension(:), intent(in) :: CellIndices
    integer, intent(out) :: SplitDir
    real(rk), intent(out) :: Split

    integer(lk) :: l
    integer(lk) :: NumCells
    real(rk), dimension(CellCart%nd) :: BoundsSize
    integer, dimension(MAX_ND) :: Cell
    real(rk) :: MeanCenter

    NumCells = size(CellIndices,kind=lk)

    BoundsSize = ovkBBSize(Bounds)

    SplitDir = maxloc(BoundsSize,dim=1)

    MeanCenter = 0._rk
    do l = 1_lk, NumCells
      Cell = OverlappingCells(CellIndices(l),:)
      MeanCenter = MeanCenter + 0.5_rk*(AllCellLowers(SplitDir)%values(Cell(1),Cell(2),Cell(3)) + &
        AllCellUppers(SplitDir)%values(Cell(1),Cell(2),Cell(3)))
    end do
    MeanCenter = MeanCenter/max(real(NumCells,kind=rk), 1._rk)

    Split = MeanCenter

  end subroutine FindSplit

  subroutine DivideCells(CellCart, AllCellLowers, AllCellUppers, OverlappingCells, CellIndices, &
    SplitDir, Split, LeftCellIndices, RightCellIndices)

    type(ovk_cart), intent(in) :: CellCart
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellLowers, AllCellUppers
    integer, dimension(:,:), intent(in) :: OverlappingCells
    integer(lk), dimension(:), intent(in) :: CellIndices
    integer, intent(in) :: SplitDir
    real(rk), intent(in) :: Split
    integer(lk), dimension(:), allocatable :: LeftCellIndices
    integer(lk), dimension(:), allocatable :: RightCellIndices

    integer(lk) :: l
    integer(lk) :: NumCells
    integer, dimension(MAX_ND) :: Cell
    integer(lk) :: NumLeftCells, NumRightCells
    integer(lk) :: NextLeftCell, NextRightCell

    NumCells = size(CellIndices,kind=lk)

    NumLeftCells = 0_lk
    NumRightCells = 0_lk
    do l = 1_lk, NumCells
      Cell = OverlappingCells(CellIndices(l),:)
      if (AllCellLowers(SplitDir)%values(Cell(1),Cell(2),Cell(3)) <= Split) then
        NumLeftCells = NumLeftCells + 1_lk
      end if
      if (AllCellUppers(SplitDir)%values(Cell(1),Cell(2),Cell(3)) >= Split) then
        NumRightCells = NumRightCells + 1_lk
      end if
    end do

    allocate(LeftCellIndices(NumLeftCells))
    allocate(RightCellIndices(NumRightCells))

    NextLeftCell = 1_lk
    NextRightCell = 1_lk
    do l = 1_lk, NumCells
      Cell = OverlappingCells(CellIndices(l),:)
      if (AllCellLowers(SplitDir)%values(Cell(1),Cell(2),Cell(3)) <= Split) then
        LeftCellIndices(NextLeftCell) = CellIndices(l)
        NextLeftCell = NextLeftCell + 1_lk
      end if
      if (AllCellUppers(SplitDir)%values(Cell(1),Cell(2),Cell(3)) >= Split) then
        RightCellIndices(NextRightCell) = CellIndices(l)
        NextRightCell = NextRightCell + 1_lk
      end if
    end do

  end subroutine DivideCells

  subroutine DistributeCellsToBins(BinCart, Bounds, CellCart, AllCellLowers, AllCellUppers, &
    OverlappingCells, CellIndices, BinCellsStart, BinCellsEnd, BinCells)

    type(ovk_cart), intent(in) :: BinCart
    type(ovk_bbox), intent(in) :: Bounds
    type(ovk_cart), intent(in) :: CellCart
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellLowers, AllCellUppers
    integer, dimension(:,:), intent(in) :: OverlappingCells
    integer(lk), dimension(:), intent(in) :: CellIndices
    type(ovk_field_large_int), intent(out) :: BinCellsStart, BinCellsEnd
    integer(lk), dimension(:), intent(out), allocatable :: BinCells

    integer :: i, j, k, d
    integer(lk) :: l, m
    integer(lk) :: NumCells
    real(rk), dimension(CellCart%nd) :: Origin
    real(rk), dimension(CellCart%nd) :: BinSize
    type(ovk_field_large_int) :: NumCellsInBin
    integer(lk) :: NumBinCells
    integer, dimension(MAX_ND) :: Cell
    real(rk), dimension(CellCart%nd) :: CellLower, CellUpper
    integer, dimension(MAX_ND) :: LowerBin, UpperBin
    type(ovk_field_large_int) :: NextCellInBin
    integer(lk) :: BinCellsIndex

    NumCells = size(CellIndices,kind=lk)

    Origin = Bounds%b(:BinCart%nd)
    BinSize = ovkBBSize(Bounds)/real(ovkCartSize(BinCart),kind=rk)

    NumCellsInBin = ovk_field_large_int_(BinCart, 0_lk)
    NumBinCells = 0

    do l = 1_lk, NumCells
      Cell = OverlappingCells(CellIndices(l),:)
      do d = 1, CellCart%nd
        CellLower(d) = AllCellLowers(d)%values(Cell(1),Cell(2),Cell(3))
        CellUpper(d) = AllCellUppers(d)%values(Cell(1),Cell(2),Cell(3))
      end do
      LowerBin(:BinCart%nd) = int(ovkCartesianGridCell(Origin, BinSize, CellLower))
      LowerBin(:BinCart%nd) = ovkCartClamp(BinCart, LowerBin)
      LowerBin(BinCart%nd+1:) = 1
      UpperBin(:BinCart%nd) = int(ovkCartesianGridCell(Origin, BinSize, CellUpper))
      UpperBin(:BinCart%nd) = ovkCartClamp(BinCart, UpperBin)
      UpperBin(BinCart%nd+1:) = 1
      do k = LowerBin(3), UpperBin(3)
        do j = LowerBin(2), UpperBin(2)
          do i = LowerBin(1), UpperBin(1)
            NumCellsInBin%values(i,j,k) = NumCellsInBin%values(i,j,k) + 1_lk
            NumBinCells = NumBinCells + 1_lk
          end do
        end do
      end do
    end do

    BinCellsStart = ovk_field_large_int_(BinCart)
    BinCellsEnd = ovk_field_large_int_(BinCart)

    m = 1_lk
    do k = BinCart%is(3), BinCart%ie(3)
      do j = BinCart%is(2), BinCart%ie(2)
        do i = BinCart%is(1), BinCart%ie(1)
          BinCellsStart%values(i,j,k) = m
          BinCellsEnd%values(i,j,k) = m + NumCellsInBin%values(i,j,k)-1
          m = m + NumCellsInBin%values(i,j,k)
        end do
      end do
    end do

    NextCellInBin = ovk_field_large_int_(BinCart, 1_lk)

    allocate(BinCells(NumBinCells))

    do l = 1_lk, NumCells
      Cell = OverlappingCells(CellIndices(l),:)
      do d = 1, CellCart%nd
        CellLower(d) = AllCellLowers(d)%values(Cell(1),Cell(2),Cell(3))
        CellUpper(d) = AllCellUppers(d)%values(Cell(1),Cell(2),Cell(3))
      end do
      LowerBin(:BinCart%nd) = int(ovkCartesianGridCell(Origin, BinSize, CellLower))
      LowerBin(:BinCart%nd) = ovkCartClamp(BinCart, LowerBin)
      LowerBin(BinCart%nd+1:) = 1
      UpperBin(:BinCart%nd) = int(ovkCartesianGridCell(Origin, BinSize, CellUpper))
      UpperBin(:BinCart%nd) = ovkCartClamp(BinCart, UpperBin)
      UpperBin(BinCart%nd+1:) = 1
      do k = LowerBin(3), UpperBin(3)
        do j = LowerBin(2), UpperBin(2)
          do i = LowerBin(1), UpperBin(1)
            BinCellsIndex = BinCellsStart%values(i,j,k) + NextCellInBin%values(i,j,k)-1
            BinCells(BinCellsIndex) = ovkCartTupleToIndex(CellCart, Cell)
            NextCellInBin%values(i,j,k) = NextCellInBin%values(i,j,k) + 1
          end do
        end do
      end do
    end do

  end subroutine DistributeCellsToBins

  recursive subroutine DestroyOverlapAccelNode(Node)

    type(t_node), intent(inout) :: Node

    if (associated(Node%hash_grid)) then
      deallocate(Node%hash_grid)
    end if

    if (associated(Node%left_child)) then
      call DestroyOverlapAccelNode(Node%left_child)
      deallocate(Node%left_child)
    end if

    if (associated(Node%right_child)) then
      call DestroyOverlapAccelNode(Node%right_child)
      deallocate(Node%right_child)
    end if

  end subroutine DestroyOverlapAccelNode

  function FindOverlappingCell(Accel, Coords, OverlapTolerance) result(Cell)

    type(t_overlap_accel), intent(in) :: Accel
    real(rk), dimension(Accel%nd), intent(in) :: Coords
    real(rk), intent(in) :: OverlapTolerance
    integer, dimension(Accel%nd) :: Cell

    type(ovk_grid), pointer :: Grid

    Grid => Accel%grid

    Cell = Grid%cart%is(:Grid%nd)-1

    if (ovkBBContainsPoint(Accel%bounds, Coords)) then
      Cell = FindOverlappingCellInNode(Accel%root, Grid, Coords, OverlapTolerance)
    end if

  end function FindOverlappingCell

  recursive function FindOverlappingCellInNode(Node, Grid, Coords, OverlapTolerance) result(Cell)

    type(t_node), intent(in) :: Node
    type(ovk_grid), intent(in) :: Grid
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    real(rk), intent(in) :: OverlapTolerance
    integer, dimension(Grid%nd) :: Cell

    integer(lk) :: l, m
    integer, dimension(Grid%nd) :: Bin
    integer(lk) :: BinStart, BinEnd
    logical :: LeafNode
    integer, dimension(Grid%nd) :: CandidateCell

    Cell = Grid%cell_cart%is(:Grid%nd) - 1

    LeafNode = associated(Node%hash_grid)

    if (.not. LeafNode) then

      if (Coords(Node%split_dir) <= Node%split) then
        Cell = FindOverlappingCellInNode(Node%left_child, Grid, Coords, OverlapTolerance)
      else
        Cell = FindOverlappingCellInNode(Node%right_child, Grid, Coords, OverlapTolerance)
      end if

    else

      Bin = HashGridBin(Node%hash_grid, Coords)

      if (ovkCartContains(Node%hash_grid%cart, Bin)) then

        l = ovkCartTupleToIndex(Node%hash_grid%cart, Bin)
        BinStart = Node%hash_grid%bin_start(l)
        BinEnd = Node%hash_grid%bin_start(l+1_lk)-1_lk

        do m = BinStart, BinEnd
          CandidateCell = ovkCartIndexToTuple(Grid%cell_cart, Node%hash_grid%bin_contents(m))
          if (ovkOverlapsGridCell(Grid, CandidateCell, Coords, OverlapTolerance)) then
            Cell = CandidateCell
            return
          end if
        end do

      end if

    end if

  end function FindOverlappingCellInNode

  subroutine PrintStats(Accel)

    type(t_overlap_accel), intent(in) :: Accel

    type(t_logger) :: Logger
    integer(lk) :: NumLeaves
    integer(lk) :: TotalLeafDepth
    integer(lk) :: NumBins
    integer(lk) :: NumNonEmptyBins
    integer :: MinBinEntries, MaxBinEntries
    integer(lk) :: TotalBinEntries
    real(rk) :: PercentFilled
    real(rk) :: AvgEntriesPerBin
    real(rk) :: AvgLeafDepth
    integer(lk), dimension(10) :: Histogram

    Logger = Accel%grid%logger

    if (associated(Accel%root)) then
      call LeafStats(Accel%root, NumLeaves, TotalLeafDepth)
      call BinStats(Accel%root, NumBins, NumNonEmptyBins, MinBinEntries, MaxBinEntries, &
        TotalBinEntries)
    else
      NumLeaves = 0_lk
      NumBins = 0_lk
    end if

    write (Logger%status_file, '(2a)') "* Number of leaf nodes: ", trim(LargeIntToString(NumLeaves))
    if (NumLeaves > 0_lk) then
      AvgLeafDepth = real(TotalLeafDepth,kind=rk)/real(NumLeaves,kind=rk)
      write (Logger%status_file, '(a,f10.4)') "* Average leaf node depth: ", AvgLeafDepth
    end if
    write (Logger%status_file, '(2a)') "* Number of bins: ", trim(LargeIntToString(NumBins))
    if (NumBins > 0_lk) then
      PercentFilled = 100._rk * real(NumNonEmptyBins,kind=rk)/real(NumBins,kind=rk)
      AvgEntriesPerBin = real(TotalBinEntries,kind=rk)/real(NumNonEmptyBins,kind=rk)
      call BinEntryHistogram(Accel%root, MinBinEntries, MaxBinEntries, 10, Histogram)
      write (Logger%status_file, '(3a,f8.4,a)') "* Number of non-empty bins: ", &
        trim(LargeIntToString(NumNonEmptyBins)), " (", PercentFilled, "%)"
      write (Logger%status_file, '(a,f10.4)') "* Average cells per non-empty bin: ", AvgEntriesPerBin
      write (Logger%status_file, '(2a)') "* Smallest number of cells per bin: ", &
        trim(IntToString(MinBinEntries))
      write (Logger%status_file, '(2a)') "* Largest number of cells per bin: ", &
        trim(IntToString(MaxBinEntries))
      write (Logger%status_file, '(a)') "* Bin cell count histogram:"
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(1))), "  "
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(2))), "  "
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(3))), "  "
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(4))), "  "
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(5))), "  "
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(6))), "  "
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(7))), "  "
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(8))), "  "
      write (Logger%status_file, '(2a)', advance='no') trim(LargeIntToString(Histogram(9))), "  "
      write (Logger%status_file, '(a)') trim(LargeIntToString(Histogram(10)))
    end if

  end subroutine PrintStats

  recursive subroutine LeafStats(Node, NumLeaves, TotalLeafDepth)

    type(t_node), intent(in) :: Node
    integer(lk), intent(out) :: NumLeaves
    integer(lk), intent(out) :: TotalLeafDepth

    logical :: LeafNode
    integer(lk) :: NumLeavesLeft, NumLeavesRight
    integer(lk) :: TotalLeafDepthLeft, TotalLeafDepthRight

    LeafNode = associated(Node%hash_grid)

    if (LeafNode) then
      TotalLeafDepth = 0_lk
      NumLeaves = 1_lk
    else
      TotalLeafDepth = 0_lk
      NumLeaves = 0_lk
      call LeafStats(Node%left_child, NumLeavesLeft, TotalLeafDepthLeft)
      call LeafStats(Node%right_child, NumLeavesRight, TotalLeafDepthRight)
      NumLeaves = NumLeaves + NumLeavesLeft + NumLeavesRight
      TotalLeafDepth = TotalLeafDepth + TotalLeafDepthLeft + TotalLeafDepthRight + &
        NumLeavesLeft + NumLeavesRight
    end if

  end subroutine LeafStats

  recursive subroutine BinStats(Node, NumBins, NumNonEmptyBins, MinBinEntries, MaxBinEntries, &
    TotalBinEntries)

    type(t_node), intent(in) :: Node
    integer(lk), intent(out) :: NumBins
    integer(lk), intent(out) :: NumNonEmptyBins
    integer, intent(out) :: MinBinEntries
    integer, intent(out) :: MaxBinEntries
    integer(lk), intent(out) :: TotalBinEntries

    logical :: LeafNode
    integer(lk) :: NumBinsLeft, NumBinsRight
    integer(lk) :: NumNonEmptyBinsLeft, NumNonEmptyBinsRight
    integer :: MinBinEntriesLeft, MinBinEntriesRight
    integer :: MaxBinEntriesLeft, MaxBinEntriesRight
    integer(lk) :: TotalBinEntriesLeft, TotalBinEntriesRight

    LeafNode = associated(Node%hash_grid)

    if (LeafNode) then
      call HashGridStats(Node%hash_grid, NumBins, NumNonEmptyBins, MinBinEntries, &
        MaxBinEntries, TotalBinEntries)
    else
      NumBins = 0_lk
      NumNonEmptyBins = 0_lk
      MinBinEntries = huge(0)
      MaxBinEntries = 0
      TotalBinEntries = 0_lk
      call BinStats(Node%left_child, NumBinsLeft, NumNonEmptyBinsLeft, &
        MinBinEntriesLeft, MaxBinEntriesLeft, TotalBinEntriesLeft)
      call BinStats(Node%right_child, NumBinsRight, NumNonEmptyBinsRight, &
        MinBinEntriesRight, MaxBinEntriesRight, TotalBinEntriesRight)
      NumBins = NumBins + NumBinsLeft + NumBinsRight
      NumNonEmptyBins = NumNonEmptyBins + NumNonEmptyBinsLeft + NumNonEmptyBinsRight
      MinBinEntries = min(MinBinEntries, min(MinBinEntriesLeft, MinBinEntriesRight))
      MaxBinEntries = max(MaxBinEntries, max(MaxBinEntriesLeft, MaxBinEntriesRight))
      TotalBinEntries = TotalBinEntries + TotalBinEntriesLeft + TotalBinEntriesRight
    end if

  end subroutine BinStats

  recursive subroutine BinEntryHistogram(Node, Lower, Upper, N, Histogram)

    type(t_node), intent(in) :: Node
    integer, intent(in) :: Lower
    integer, intent(in) :: Upper
    integer, intent(in) :: N
    integer(lk), dimension(N), intent(out) :: Histogram

    logical :: LeafNode
    integer(lk), dimension(N) :: HistogramChild

    LeafNode = associated(Node%hash_grid)

    if (LeafNode) then
      call HashGridHistogram(Node%hash_grid, Lower, Upper, N, Histogram)
    else
      Histogram = 0_lk
      call BinEntryHistogram(Node%left_child, Lower, Upper, N, HistogramChild)
      Histogram = Histogram + HistogramChild
      call BinEntryHistogram(Node%right_child, Lower, Upper, N, HistogramChild)
      Histogram = Histogram + HistogramChild
    end if

  end subroutine BinEntryHistogram

end module ovkOverlapAccel
