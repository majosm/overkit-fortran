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
  implicit none

  private

  ! Internal
  public :: t_overlap_accel
  public :: t_overlap_accel_
  public :: GenerateOverlapAccel
  public :: DestroyOverlapAccel
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
    integer :: nd
    type(ovk_bbox) :: bounds
    real(rk) :: max_cell_size_deviation
    real(rk) :: bin_scale
    integer :: max_depth
    type(t_node), pointer :: root
  end type t_overlap_accel

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_overlap_accel_
    module procedure t_overlap_accel_Default
    module procedure t_overlap_accel_Empty
  end interface t_overlap_accel_

contains

  pure function t_overlap_accel_Default() result(Accel)

    type(t_overlap_accel) :: Accel

    Accel = t_overlap_accel_Empty(2)

  end function t_overlap_accel_Default

  pure function t_overlap_accel_Empty(NumDims) result(Accel)

    integer, intent(in) :: NumDims
    type(t_overlap_accel) :: Accel

    Accel%nd = NumDims
    Accel%bounds = ovk_bbox_(NumDims)
    Accel%max_cell_size_deviation = 0._rk
    Accel%bin_scale = 0._rk
    Accel%max_depth = 0
    nullify(Accel%root)

  end function t_overlap_accel_Empty

  subroutine GenerateOverlapAccel(Grid, Accel, Bounds, MaxOverlapTolerance, QualityAdjust)

    type(ovk_grid), intent(in) :: Grid
    type(t_overlap_accel), intent(out) :: Accel
    type(ovk_bbox), intent(in) :: Bounds
    real(rk), intent(in) :: MaxOverlapTolerance
    real(rk), intent(in) :: QualityAdjust

    real(rk) :: MaxCellSizeDeviation
    real(rk) :: BinScale
    integer :: MaxDepth
    integer :: i, j, k, d
    type(ovk_field_logical) :: GridCellOverlapMask
    type(ovk_field_real), dimension(:), allocatable :: GridCellLower, GridCellUpper
    integer, dimension(MAX_ND) :: Cell
    type(ovk_bbox) :: GridCellBounds
    integer(lk) :: NumOverlappingCells, iNextOverlappingCell
    integer(lk), dimension(:), allocatable :: OverlappingCells

    MaxCellSizeDeviation = 0.5_rk
    BinScale = 0.5_rk**QualityAdjust

    if (.not. ovkBBOverlaps(Bounds, Grid%bounds)) then
      Accel = t_overlap_accel_(Grid%cart%nd)
      return
    end if

    GridCellOverlapMask = ovk_field_logical_(Grid%cell_cart, .false.)

    allocate(GridCellLower(Grid%cell_cart%nd))
    allocate(GridCellUpper(Grid%cell_cart%nd))
    do d = 1, Grid%cell_cart%nd
      GridCellLower(d) = ovk_field_real_(Grid%cell_cart, 0._rk)
      GridCellUpper(d) = ovk_field_real_(Grid%cell_cart, 0._rk)
    end do

    Accel%nd = Grid%cart%nd
    Accel%bounds = ovk_bbox_(Grid%cart%nd)

    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          Cell = [i,j,k]
          if (ovkGridCellExists(Grid, Cell)) then
            GridCellBounds = ovkGridCellBounds(Grid, Cell)
            GridCellBounds = ovkBBScale(GridCellBounds, 1._rk + 2._rk * MaxOverlapTolerance)
            GridCellOverlapMask%values(i,j,k) = ovkBBOverlaps(Bounds, GridCellBounds)
            if (GridCellOverlapMask%values(i,j,k)) then
              do d = 1, Grid%cell_cart%nd
                GridCellLower(d)%values(i,j,k) = GridCellBounds%b(d)
                GridCellUpper(d)%values(i,j,k) = GridCellBounds%e(d)
              end do
              Accel%bounds = ovkBBUnion(Accel%bounds, GridCellBounds)
            end if
          end if
        end do
      end do
    end do

    NumOverlappingCells = ovkCountMask(GridCellOverlapMask)

    if (NumOverlappingCells == 0) then
      Accel = t_overlap_accel_(Grid%cart%nd)
      return
    end if

    allocate(OverlappingCells(NumOverlappingCells))

    iNextOverlappingCell = 1_lk
    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          Cell = [i,j,k]
          if (GridCellOverlapMask%values(i,j,k)) then
            OverlappingCells(iNextOverlappingCell) = ovkCartTupleToIndex(Grid%cell_cart, Cell)
            iNextOverlappingCell = iNextOverlappingCell + 1_lk
          end if
        end do
      end do
    end do

    ! Stop when the number of cells is likely to be small (say 1000)
    MaxDepth = int(ceiling(log(max(real(NumOverlappingCells,kind=rk)/1000._rk,1._rk))/log(2._rk)))

    Accel%max_cell_size_deviation = MaxCellSizeDeviation
    Accel%bin_scale = BinScale
    Accel%max_depth = MaxDepth

    allocate(Accel%root)
    call GenerateOverlapAccelNode(Accel%root, Grid%cell_cart, Accel%bounds, GridCellLower, &
      GridCellUpper, OverlappingCells, Accel%max_cell_size_deviation, Accel%bin_scale, &
      Accel%max_depth, 0)

    if (Grid%logger%verbose) then
      call PrintStats(Accel)
    end if

  end subroutine GenerateOverlapAccel

  subroutine DestroyOverlapAccel(Accel)

    type(t_overlap_accel), intent(inout) :: Accel

    if (associated(Accel%root)) then
      call DestroyOverlapAccelNode(Accel%root)
      deallocate(Accel%root)
    end if

  end subroutine DestroyOverlapAccel

  recursive subroutine GenerateOverlapAccelNode(Node, CellCart, Bounds, AllCellLower, &
    AllCellUpper, OverlappingCells, MaxCellSizeDeviation, BinScale, MaxDepth, Depth)

    type(t_node), intent(out) :: Node
    type(ovk_cart), intent(in) :: CellCart
    type(ovk_bbox), intent(in) :: Bounds
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellLower
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellUpper
    integer(lk), dimension(:), intent(in) :: OverlappingCells
    real(rk), intent(in) :: MaxCellSizeDeviation
    real(rk), intent(in) :: BinScale
    integer, intent(in) :: MaxDepth
    integer, intent(in) :: Depth

    integer :: i, j, k, d
    integer(lk) :: l
    integer(lk) :: NumOverlappingCells
    integer, dimension(MAX_ND) :: Cell
    real(rk), dimension(CellCart%nd) :: CellSize, MeanCellSize, CellSizeDeviation
    logical :: LeafNode
    real(rk), dimension(CellCart%nd) :: BoundsSize
    integer :: SplitDir
    real(rk) :: Split
    type(ovk_bbox) :: LeftBounds, RightBounds
    integer(lk), dimension(:), allocatable :: LeftCells, RightCells
    integer, dimension(CellCart%nd) :: NumBins
    type(ovk_cart) :: BinCart
    type(ovk_field_large_int) :: BinCellsStart, BinCellsEnd
    integer(lk), dimension(:), allocatable :: BinCells
    type(ovk_field_int) :: NumBinElements
    type(ovk_bbox) :: HashGridBounds
    type(ovk_bbox) :: CellBounds
    integer(lk) :: BinStart, BinEnd

    nullify(Node%hash_grid)
    nullify(Node%left_child)
    nullify(Node%right_child)

    NumOverlappingCells = size(OverlappingCells,kind=lk)

    MeanCellSize = 0._rk
    do l = 1_lk, NumOverlappingCells
      Cell(:CellCart%nd) = ovkCartIndexToTuple(CellCart, OverlappingCells(l))
      Cell(CellCart%nd+1:) = 1
      do d = 1, CellCart%nd
        CellSize(d) = max(AllCellUpper(d)%values(Cell(1),Cell(2),Cell(3)) - &
          AllCellLower(d)%values(Cell(1),Cell(2),Cell(3)), 0._rk)
      end do
      MeanCellSize = MeanCellSize + CellSize
    end do
    MeanCellSize = MeanCellSize/real(NumOverlappingCells,kind=rk)

    if (NumOverlappingCells > 1_lk) then
      CellSizeDeviation = 0._rk
      do l = 1_lk, NumOverlappingCells
        Cell(:CellCart%nd) = ovkCartIndexToTuple(CellCart, OverlappingCells(l))
        Cell(CellCart%nd+1:) = 1
        do d = 1, CellCart%nd
          CellSize(d) = max(AllCellUpper(d)%values(Cell(1),Cell(2),Cell(3)) - &
            AllCellLower(d)%values(Cell(1),Cell(2),Cell(3)), 0._rk)
        end do
        CellSizeDeviation = CellSizeDeviation + (CellSize - MeanCellSize)**2
      end do
      CellSizeDeviation = sqrt(CellSizeDeviation/real(NumOverlappingCells-1_lk,kind=rk))
    else
      CellSizeDeviation = 0._rk
    end if

    CellSizeDeviation = CellSizeDeviation/MeanCellSize

    ! Don't want to generate a hash grid if the set of cells is highly uniform or very large
    LeafNode = (all(CellSizeDeviation <= MaxCellSizeDeviation) .and. NumOverlappingCells <= 2**26) &
      .or. Depth == MaxDepth

    if (.not. LeafNode) then

      BoundsSize = ovkBBSize(Bounds)
      SplitDir = maxloc(BoundsSize,dim=1)

      Split = FindSplit(CellCart, Bounds, AllCellLower, AllCellUpper, OverlappingCells, SplitDir)

      Node%split_dir = SplitDir
      Node%split = Split

      allocate(Node%left_child)
      allocate(Node%right_child)

      call DivideCells(CellCart, AllCellLower, AllCellUpper, OverlappingCells, SplitDir, Split, &
        LeftCells, RightCells)

      LeftBounds = Bounds
      LeftBounds%e(SplitDir) = Split
      RightBounds = Bounds
      RightBounds%b(SplitDir) = Split

      call GenerateOverlapAccelNode(Node%left_child, CellCart, LeftBounds, AllCellLower, &
        AllCellUpper, LeftCells, MaxCellSizeDeviation, BinScale, MaxDepth, Depth+1)
      call GenerateOverlapAccelNode(Node%right_child, CellCart, RightBounds, AllCellLower, &
        AllCellUpper, RightCells, MaxCellSizeDeviation, BinScale, MaxDepth, Depth+1)

    else

      NumBins = max(int(ceiling(ovkBBSize(Bounds)/(BinScale * MeanCellSize))),1)
      BinCart = ovk_cart_(CellCart%nd, NumBins)

      HashGridBounds = ovk_bbox_(CellCart%nd)
      do l = 1, NumOverlappingCells
        Cell(:CellCart%nd) = ovkCartIndexToTuple(CellCart, OverlappingCells(l))
        Cell(CellCart%nd+1:) = 1
        CellBounds = ovk_bbox_(CellCart%nd)
        do d = 1, CellCart%nd
          CellBounds%b(d) = AllCellLower(d)%values(Cell(1),Cell(2),Cell(3))
          CellBounds%e(d) = AllCellUpper(d)%values(Cell(1),Cell(2),Cell(3))
        end do
        HashGridBounds = ovkBBUnion(HashGridBounds, CellBounds)
      end do

      call DistributeCellsToBins(BinCart, HashGridBounds, CellCart, AllCellLower, AllCellUpper, &
        OverlappingCells, BinCellsStart, BinCellsEnd, BinCells)

      NumBinElements = ovk_field_int_(BinCart)
      NumBinElements%values = int(BinCellsEnd%values - BinCellsStart%values + 1_lk)

      allocate(Node%hash_grid)
      Node%hash_grid = t_hash_grid_(BinCart, HashGridBounds, NumBinElements)

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

  function FindSplit(CellCart, Bounds, AllCellLower, AllCellUpper, OverlappingCells, SplitDir) &
    result(Split)

    type(ovk_cart), intent(in) :: CellCart
    type(ovk_bbox), intent(in) :: Bounds
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellLower, AllCellUpper
    integer(lk), dimension(:), intent(in) :: OverlappingCells
    integer, intent(in) :: SplitDir
    real(rk) :: Split

    integer(lk) :: l
    integer(lk) :: NumOverlappingCells
    integer(lk) :: NumBins
    real(rk) :: Origin
    real(rk) :: BinSize
    integer(lk), dimension(:), allocatable :: NumCellsPerBin
    integer, dimension(MAX_ND) :: Cell
    real(rk) :: Coord
    integer(lk) :: Bin
    integer(lk) :: iSplit
    integer(lk) :: NumCellsBelowSplit

    NumOverlappingCells = size(OverlappingCells,kind=lk)

    NumBins = NumOverlappingCells/10_lk

    Origin = Bounds%b(SplitDir)
    BinSize = (Bounds%e(SplitDir)-Bounds%b(SplitDir))/real(NumBins,kind=rk)

    allocate(NumCellsPerBin(NumBins))
    NumCellsPerBin = 0_lk

    do l = 1_lk, NumOverlappingCells
      Cell(:CellCart%nd) = ovkCartIndexToTuple(CellCart, OverlappingCells(l))
      Cell(CellCart%nd+1:) = 1
      Coord = 0.5_rk * (AllCellLower(SplitDir)%values(Cell(1),Cell(2),Cell(3)) + &
        AllCellUpper(SplitDir)%values(Cell(1),Cell(2),Cell(3)))
      Bin = min(max(ovkCartesianGridCell(Origin, BinSize, Coord), 1_lk), NumBins)
      NumCellsPerBin(Bin) = NumCellsPerBin(Bin) + 1_lk
    end do

    iSplit = 0_lk
    NumCellsBelowSplit = 0_lk
    do while (NumCellsBelowSplit < NumOverlappingCells/2_lk)
      iSplit = iSplit + 1_lk
      NumCellsBelowSplit = NumCellsBelowSplit + NumCellsPerBin(iSplit)
    end do

    Split = Origin + BinSize * real(iSplit-1_lk,kind=rk)

  end function FindSplit

  subroutine DivideCells(CellCart, AllCellLower, AllCellUpper, OverlappingCells, SplitDir, Split, &
    LeftCells, RightCells)

    type(ovk_cart), intent(in) :: CellCart
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellLower, AllCellUpper
    integer(lk), dimension(:), intent(in) :: OverlappingCells
    integer, intent(in) :: SplitDir
    real(rk), intent(in) :: Split
    integer(lk), dimension(:), allocatable :: LeftCells
    integer(lk), dimension(:), allocatable :: RightCells

    integer(lk) :: l
    integer(lk) :: NumOverlappingCells
    integer, dimension(MAX_ND) :: Cell
    integer(lk) :: NumLeftCells, NumRightCells
    integer(lk) :: NextLeftCell, NextRightCell

    NumOverlappingCells = size(OverlappingCells,kind=lk)

    NumLeftCells = 0_lk
    NumRightCells = 0_lk
    do l = 1_lk, NumOverlappingCells
      Cell(:CellCart%nd) = ovkCartIndexToTuple(CellCart, OverlappingCells(l))
      Cell(CellCart%nd+1:) = 1
      if (AllCellLower(SplitDir)%values(Cell(1),Cell(2),Cell(3)) <= Split) then
        NumLeftCells = NumLeftCells + 1_lk
      end if
      if (AllCellUpper(SplitDir)%values(Cell(1),Cell(2),Cell(3)) >= Split) then
        NumRightCells = NumRightCells + 1_lk
      end if
    end do

    allocate(LeftCells(NumLeftCells))
    allocate(RightCells(NumRightCells))

    NextLeftCell = 1_lk
    NextRightCell = 1_lk
    do l = 1_lk, NumOverlappingCells
      Cell(:CellCart%nd) = ovkCartIndexToTuple(CellCart, OverlappingCells(l))
      Cell(CellCart%nd+1:) = 1
      if (AllCellLower(SplitDir)%values(Cell(1),Cell(2),Cell(3)) <= Split) then
        LeftCells(NextLeftCell) = OverlappingCells(l)
        NextLeftCell = NextLeftCell + 1_lk
      end if
      if (AllCellUpper(SplitDir)%values(Cell(1),Cell(2),Cell(3)) >= Split) then
        RightCells(NextRightCell) = OverlappingCells(l)
        NextRightCell = NextRightCell + 1_lk
      end if
    end do

  end subroutine DivideCells

  subroutine DistributeCellsToBins(BinCart, Bounds, CellCart, AllCellLower, AllCellUpper, &
    OverlappingCells, BinCellsStart, BinCellsEnd, BinCells)

    type(ovk_cart), intent(in) :: BinCart
    type(ovk_bbox), intent(in) :: Bounds
    type(ovk_cart), intent(in) :: CellCart
    type(ovk_field_real), dimension(CellCart%nd), intent(in) :: AllCellLower, AllCellUpper
    integer(lk), dimension(:), intent(in) :: OverlappingCells
    type(ovk_field_large_int), intent(out) :: BinCellsStart, BinCellsEnd
    integer(lk), dimension(:), intent(out), allocatable :: BinCells

    integer :: i, j, k, d
    integer(lk) :: l, m
    integer(lk) :: NumOverlappingCells
    real(rk), dimension(CellCart%nd) :: Origin
    real(rk), dimension(CellCart%nd) :: BinSize
    type(ovk_field_large_int) :: NumCellsInBin
    integer(lk) :: NumBinCells
    integer, dimension(MAX_ND) :: Cell
    real(rk), dimension(CellCart%nd) :: CellLower, CellUpper
    integer, dimension(MAX_ND) :: LowerBin, UpperBin
    type(ovk_field_large_int) :: NextCellInBin
    integer(lk) :: BinCellsIndex

    NumOverlappingCells = size(OverlappingCells,kind=lk)

    Origin = Bounds%b(:BinCart%nd)
    BinSize = ovkBBSize(Bounds)/real(ovkCartSize(BinCart),kind=rk)

    NumCellsInBin = ovk_field_large_int_(BinCart, 0_lk)
    NumBinCells = 0

    do l = 1_lk, NumOverlappingCells
      Cell(:CellCart%nd) = ovkCartIndexToTuple(CellCart, OverlappingCells(l))
      Cell(CellCart%nd+1:) = 1
      do d = 1, CellCart%nd
        CellLower(d) = AllCellLower(d)%values(Cell(1),Cell(2),Cell(3))
        CellUpper(d) = AllCellUpper(d)%values(Cell(1),Cell(2),Cell(3))
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

    do l = 1_lk, NumOverlappingCells
      Cell(:CellCart%nd) = ovkCartIndexToTuple(CellCart, OverlappingCells(l))
      Cell(CellCart%nd+1:) = 1
      do d = 1, CellCart%nd
        CellLower(d) = AllCellLower(d)%values(Cell(1),Cell(2),Cell(3))
        CellUpper(d) = AllCellUpper(d)%values(Cell(1),Cell(2),Cell(3))
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
            BinCells(BinCellsIndex) = OverlappingCells(l)
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

  function FindOverlappingCell(Grid, Accel, Coords, OverlapTolerance) result(Cell)

    type(ovk_grid), intent(in) :: Grid
    type(t_overlap_accel), intent(in) :: Accel
    real(rk), dimension(Grid%cart%nd), intent(in) :: Coords
    real(rk), intent(in) :: OverlapTolerance
    integer, dimension(Grid%cart%nd) :: Cell

    Cell = Grid%cart%is(:Grid%cart%nd)-1

    if (ovkBBContainsPoint(Accel%bounds, Coords)) then
      Cell = FindOverlappingCellInNode(Grid, Accel%root, OverlapTolerance, Coords)
    end if

  end function FindOverlappingCell

  recursive function FindOverlappingCellInNode(Grid, Node, OverlapTolerance, Coords) result(Cell)

    type(ovk_grid), intent(in) :: Grid
    type(t_node), intent(in) :: Node
    real(rk), intent(in) :: OverlapTolerance
    real(rk), dimension(Grid%cart%nd), intent(in) :: Coords
    integer, dimension(Grid%cart%nd) :: Cell

    integer(lk) :: l, m
    integer, dimension(Grid%cart%nd) :: Bin
    integer(lk) :: BinStart, BinEnd
    logical :: LeafNode
    integer, dimension(Grid%cart%nd) :: CandidateCell

    Cell = Grid%cell_cart%is(:Grid%cell_cart%nd) - 1

    LeafNode = associated(Node%hash_grid)

    if (.not. LeafNode) then

      if (Coords(Node%split_dir) <= Node%split) then
        Cell = FindOverlappingCellInNode(Grid, Node%left_child, OverlapTolerance, Coords)
      else
        Cell = FindOverlappingCellInNode(Grid, Node%right_child, OverlapTolerance, Coords)
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

    call LeafStats(Accel%root, NumLeaves, TotalLeafDepth)

    AvgLeafDepth = real(TotalLeafDepth,kind=rk)/real(NumLeaves,kind=rk)

    call BinStats(Accel%root, NumBins, NumNonEmptyBins, MinBinEntries, MaxBinEntries, &
      TotalBinEntries)

    PercentFilled = 100._rk * real(NumNonEmptyBins,kind=rk)/real(NumBins,kind=rk)
    AvgEntriesPerBin = real(TotalBinEntries,kind=rk)/real(NumNonEmptyBins,kind=rk)

    write (*, '(2a)') "* Number of leaf nodes: ", trim(LargeIntToString(NumLeaves))
    write (*, '(a,f10.4)') "* Average leaf node depth: ", AvgLeafDepth
    write (*, '(2a)') "* Number of bins: ", trim(LargeIntToString(NumBins))
    write (*, '(3a,f8.4,a)') "* Number of non-empty bins: ", trim(LargeIntToString(NumNonEmptyBins)), &
      " (", PercentFilled, "%)"
    write (*, '(a,f10.4)') "* Average cells per non-empty bin: ", AvgEntriesPerBin
    write (*, '(2a)') "* Smallest number of cells per bin: ", trim(IntToString(MinBinEntries))
    write (*, '(2a)') "* Largest number of cells per bin: ", trim(IntToString(MaxBinEntries))

    call BinEntryHistogram(Accel%root, MinBinEntries, MaxBinEntries, 10, Histogram)

    write (*, '(a)') "* Bin cell count histogram:"
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(1))), "  "
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(2))), "  "
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(3))), "  "
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(4))), "  "
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(5))), "  "
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(6))), "  "
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(7))), "  "
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(8))), "  "
    write (*, '(2a)', advance='no') trim(LargeIntToString(Histogram(9))), "  "
    write (*, '(a)') trim(LargeIntToString(Histogram(10)))

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
