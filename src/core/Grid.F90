! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkGrid

  use ovkBoundingBox
  use ovkCart
  use ovkField
  use ovkGeometry
  use ovkGlobal
  use ovkMask
  implicit none

  private

  ! API
  public :: ovk_grid
  public :: ovk_grid_
  public :: ovkMakeGrid
  public :: ovkDestroyGrid
  public :: ovkGetCellVertexData
  public :: ovkOverlapsCell
  public :: ovkCoordsInCell
  public :: ovkCellSize
  public :: ovkAvgCellSizeAroundPoint
  public :: ovkGenerateBBOverlapMask
  public :: ovkPeriodicExtend
  public :: OVK_GRID_TYPE_CARTESIAN
  public :: OVK_GRID_TYPE_CARTESIAN_ROTATED
  public :: OVK_GRID_TYPE_RECTILINEAR
  public :: OVK_GRID_TYPE_RECTILINEAR_ROTATED
  public :: OVK_GRID_TYPE_CURVILINEAR

  type ovk_grid
    type(ovk_cart) :: cart
    type(ovk_cart) :: cell_cart
    real(rk), dimension(MAX_ND) :: periodic_length
    integer :: grid_type
    integer :: id
    type(ovk_bbox) :: bounds
    type(ovk_field_real), dimension(:), allocatable :: xyz
    type(ovk_field_logical) :: grid_mask
    type(ovk_field_logical) :: boundary_mask
    type(ovk_field_real) :: cell_sizes
  end type ovk_grid

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_grid_
    module procedure ovk_grid_Default
  end interface ovk_grid_

  integer, parameter :: OVK_GRID_TYPE_CARTESIAN = 1
  integer, parameter :: OVK_GRID_TYPE_CARTESIAN_ROTATED = 2
  integer, parameter :: OVK_GRID_TYPE_RECTILINEAR = 3
  integer, parameter :: OVK_GRID_TYPE_RECTILINEAR_ROTATED = 4
  integer, parameter :: OVK_GRID_TYPE_CURVILINEAR = 5

contains

  function ovk_grid_Default(NumDims) result(Grid)

    integer, intent(in) :: NumDims
    type(ovk_grid) :: Grid

    Grid%cart = ovk_cart_(NumDims)
    Grid%cell_cart = ovk_cart_(NumDims)
    Grid%periodic_length = 0._rk
    Grid%grid_type = OVK_GRID_TYPE_CURVILINEAR
    Grid%id = 0
    Grid%bounds = ovk_bbox_(NumDims)
    Grid%grid_mask = ovk_field_logical_(NumDims)
    Grid%boundary_mask = ovk_field_logical_(NumDims)
    Grid%cell_sizes = ovk_field_real_(NumDims)

  end function ovk_grid_Default

  subroutine ovkMakeGrid(Grid, Cart, Coords, PeriodicLength, GridMask, BoundaryMask, GridType)

    type(ovk_grid), intent(out) :: Grid
    type(ovk_cart), intent(in) :: Cart
    type(ovk_field_real), dimension(Cart%nd), intent(in) :: Coords
    real(rk), dimension(Cart%nd), intent(in), optional :: PeriodicLength
    type(ovk_field_logical), intent(in), optional :: GridMask
    type(ovk_field_logical), intent(in), optional :: BoundaryMask
    integer, intent(in), optional :: GridType

    real(rk), dimension(MAX_ND) :: PeriodicLength_
    integer :: GridType_
    integer :: dir, idir, jdir
    integer :: i, j, k, l
    integer, dimension(MAX_ND) :: Point, AdjustedPoint
    real(rk), dimension(Cart%nd) :: PeriodicCoords
    integer, dimension(MAX_ND) :: Cell
    real(rk), dimension(Cart%nd,2**Cart%nd) :: VertexCoords
    logical, dimension(2**Cart%nd) :: VertexGridMaskValues

    if (present(PeriodicLength)) then
      PeriodicLength_(:Cart%nd) = PeriodicLength
      PeriodicLength_(Cart%nd+1:) = 0._rk
    else
      PeriodicLength_ = 0._rk
    end if

    if (present(GridType)) then
      GridType_ = GridType
    else
      GridType_ = OVK_GRID_TYPE_CURVILINEAR
    end if

    Grid%cart = ovkCartConvertPeriodicStorage(Cart, OVK_NO_OVERLAP_PERIODIC)
    Grid%cell_cart = ovkCartConvertPointToCell(Grid%cart)
    Grid%periodic_length = PeriodicLength_
    Grid%grid_type = GridType_
    Grid%id = 0

    allocate(Grid%xyz(Grid%cart%nd))
    do dir = 1, Grid%cart%nd
      if (Coords(dir)%cart == Grid%cart) then
        Grid%xyz(dir) = Coords(dir)
      else
        Grid%xyz(dir) = ovk_field_real_(Grid%cart)
        Grid%xyz(dir)%values = Coords(dir)%values(Grid%cart%is(1):Grid%cart%ie(1), &
          Grid%cart%is(2):Grid%cart%ie(2),Grid%cart%is(3):Grid%cart%ie(3))
      end if
    end do

    Grid%bounds = ovk_bbox_(Grid%cart%nd)
    do dir = 1, Grid%cart%nd
      Grid%bounds%b(dir) = minval(Grid%xyz(dir)%values)
      Grid%bounds%e(dir) = maxval(Grid%xyz(dir)%values)
    end do

    ! Add contribution to bounds from periodic points
    do dir = 1, Grid%cart%nd
      if (Grid%cart%periodic(dir)) then
        idir = modulo((dir+1)-1,MAX_ND) + 1
        jdir = modulo((dir+2)-1,MAX_ND) + 1
        do j = Grid%cart%is(jdir), Grid%cart%ie(jdir)
          do i = Grid%cart%is(idir), Grid%cart%ie(idir)
            Point(dir) = Grid%cart%ie(dir)+1
            Point(idir) = i
            Point(jdir) = j
            AdjustedPoint(dir) = Grid%cart%is(dir)
            AdjustedPoint(idir) = i
            AdjustedPoint(jdir) = j
            do l = 1, Grid%cart%nd
              PeriodicCoords(l) = Grid%xyz(l)%values(AdjustedPoint(1),AdjustedPoint(2), &
                AdjustedPoint(3))
            end do
            PeriodicCoords = ovkPeriodicExtend(Grid%cart, Grid%periodic_length, Point, &
              PeriodicCoords)
            Grid%bounds = ovkBBExtend(Grid%bounds, PeriodicCoords)
          end do
        end do
      end if
    end do

    if (present(GridMask)) then
      if (GridMask%cart == Grid%cart) then
        Grid%grid_mask = GridMask
      else
        Grid%grid_mask = ovk_field_logical_(Grid%cart)
        Grid%grid_mask%values = GridMask%values(Grid%cart%is(1):Grid%cart%ie(1), &
          Grid%cart%is(2):Grid%cart%ie(2),Grid%cart%is(3):Grid%cart%ie(3))
      end if
    else
      Grid%grid_mask = ovk_field_logical_(Grid%cart, .true.)
    end if

    if (present(BoundaryMask)) then
      if (BoundaryMask%cart == Grid%cart) then
        Grid%boundary_mask = BoundaryMask
      else
        Grid%boundary_mask = ovk_field_logical_(Grid%cart)
        Grid%boundary_mask%values = BoundaryMask%values(Grid%cart%is(1):Grid%cart%ie(1), &
          Grid%cart%is(2):Grid%cart%ie(2),Grid%cart%is(3):Grid%cart%ie(3))
      end if
    else
      Grid%boundary_mask = ovk_field_logical_(Grid%cart, .false.)
    end if

    Grid%cell_sizes = ovk_field_real_(Grid%cell_cart)

    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          Cell = [i,j,k]
          call ovkGetCellVertexData(Grid, Cell, VertexCoords=VertexCoords, &
            VertexGridMaskValues=VertexGridMaskValues)
          if (all(VertexGridMaskValues)) then
            select case (Grid%cart%nd)
            case (2)
              select case (Grid%grid_type)
              case (OVK_GRID_TYPE_CARTESIAN,OVK_GRID_TYPE_RECTILINEAR)
                Grid%cell_sizes%values(i,j,k) = ovkRectangleSize(VertexCoords)
              case (OVK_GRID_TYPE_CURVILINEAR)
                Grid%cell_sizes%values(i,j,k) = ovkQuadSize(VertexCoords)
              end select
            case (3)
              select case (Grid%grid_type)
              case (OVK_GRID_TYPE_CARTESIAN,OVK_GRID_TYPE_RECTILINEAR)
                Grid%cell_sizes%values(i,j,k) = ovkCuboidSize(VertexCoords)
              case (OVK_GRID_TYPE_CURVILINEAR)
                Grid%cell_sizes%values(i,j,k) = ovkHexahedronSize(VertexCoords)
              end select
            end select
          else
            Grid%cell_sizes%values(i,j,k) = 0._rk
          end if
        end do
      end do
    end do

  end subroutine ovkMakeGrid

  subroutine ovkDestroyGrid(Grid)

    type(ovk_grid), intent(inout) :: Grid

    deallocate(Grid%xyz)

    Grid%grid_mask = ovk_field_logical_(Grid%cart%nd)
    Grid%boundary_mask = ovk_field_logical_(Grid%cart%nd)
    Grid%cell_sizes = ovk_field_real_(Grid%cart%nd)

  end subroutine ovkDestroyGrid

  subroutine ovkGetCellVertexData(Grid, Cell, VertexCoords, VertexGridMaskValues)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%cart%nd), intent(in) :: Cell
    real(rk), dimension(Grid%cart%nd,2**Grid%cart%nd), intent(out), optional :: VertexCoords
    logical, dimension(2**Grid%cart%nd), intent(out), optional :: VertexGridMaskValues

    integer :: i, j, k, l, m
    integer, dimension(MAX_ND) :: VertexLower, VertexUpper
    logical :: AwayFromBoundary
    integer, dimension(MAX_ND) :: Vertex
    integer, dimension(MAX_ND) :: AdjustedVertex
    real(rk), dimension(Grid%cart%nd) :: PrincipalCoords

    VertexLower(:Grid%cart%nd) = Cell
    VertexLower(Grid%cart%nd+1:) = 1
    VertexUpper(:Grid%cart%nd) = Cell+1
    VertexUpper(Grid%cart%nd+1:) = 1

    AwayFromBoundary = ovkCartContains(Grid%cart, VertexLower) .and. &
      ovkCartContains(Grid%cart, VertexUpper)

    if (AwayFromBoundary) then
      l = 1
      do k = VertexLower(3), VertexUpper(3)
        do j = VertexLower(2), VertexUpper(2)
          do i = VertexLower(1), VertexUpper(1)
            Vertex = [i,j,k]
            if (present(VertexCoords)) then
              do m = 1, Grid%cart%nd
                VertexCoords(m,l) = Grid%xyz(m)%values(Vertex(1),Vertex(2),Vertex(3))
              end do
            end if
            if (present(VertexGridMaskValues)) then
              VertexGridMaskValues(l) = Grid%grid_mask%values(Vertex(1),Vertex(2),Vertex(3))
            end if
            l = l + 1
          end do
        end do
      end do
    else
      l = 1
      do k = VertexLower(3), VertexUpper(3)
        do j = VertexLower(2), VertexUpper(2)
          do i = VertexLower(1), VertexUpper(1)
            Vertex = [i,j,k]
            AdjustedVertex(:Grid%cart%nd) = ovkCartPeriodicAdjust(Grid%cart, Vertex)
            AdjustedVertex(Grid%cart%nd+1:) = 1
            if (present(VertexCoords)) then
              do m = 1, Grid%cart%nd
                PrincipalCoords(m) = Grid%xyz(m)%values(AdjustedVertex(1),AdjustedVertex(2), &
                  AdjustedVertex(3))
              end do
              VertexCoords(:,l) = ovkPeriodicExtend(Grid%cart, Grid%periodic_length, Vertex, &
                PrincipalCoords)
            end if
            if (present(VertexGridMaskValues)) then
              VertexGridMaskValues(l) = Grid%grid_mask%values(AdjustedVertex(1),AdjustedVertex(2), &
                AdjustedVertex(3))
            end if
            l = l + 1
          end do
        end do
      end do
    end if

  end subroutine ovkGetCellVertexData

  function ovkOverlapsCell(Grid, Cell, Coords) result(Overlaps)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%cart%nd), intent(in) :: Cell
    real(rk), dimension(Grid%cart%nd), intent(in) :: Coords
    logical :: Overlaps

    real(rk), dimension(Grid%cart%nd,2**Grid%cart%nd) :: VertexCoords
    logical, dimension(2**Grid%cart%nd) :: VertexGridMaskValues

    call ovkGetCellVertexData(Grid, Cell, VertexCoords=VertexCoords, &
      VertexGridMaskValues=VertexGridMaskValues)

    if (.not. all(VertexGridMaskValues)) then
      Overlaps = .false.
      return
    end if

    select case (Grid%cart%nd)
    case (2)
      select case (Grid%grid_type)
      case (OVK_GRID_TYPE_CARTESIAN,OVK_GRID_TYPE_RECTILINEAR)
        Overlaps = ovkOverlapsRectangle(VertexCoords, Coords)
      case (OVK_GRID_TYPE_CURVILINEAR)
        Overlaps = ovkOverlapsQuad(VertexCoords, Coords)
      end select
    case (3)
      select case (Grid%grid_type)
      case (OVK_GRID_TYPE_CARTESIAN,OVK_GRID_TYPE_RECTILINEAR)
        Overlaps = ovkOverlapsCuboid(VertexCoords, Coords)
      case (OVK_GRID_TYPE_CURVILINEAR)
        Overlaps = ovkOverlapsHexahedron(VertexCoords, Coords)
      end select
    end select

  end function ovkOverlapsCell

  function ovkCoordsInCell(Grid, Cell, Coords) result(CoordsInCell)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%cart%nd), intent(in) :: Cell
    real(rk), dimension(Grid%cart%nd), intent(in) :: Coords
    real(rk), dimension(Grid%cart%nd) :: CoordsInCell

    real(rk), dimension(Grid%cart%nd,2**Grid%cart%nd) :: VertexCoords
    logical, dimension(2**Grid%cart%nd) :: VertexGridMaskValues
    logical :: Success

    call ovkGetCellVertexData(Grid, Cell, VertexCoords=VertexCoords, &
      VertexGridMaskValues=VertexGridMaskValues)

    if (.not. all(VertexGridMaskValues)) then
      ! Pick a number large enough to not represent any valid cell coordinates, but small enough
      ! that it can be subsequently manipulated without causing floating point errors
      CoordsInCell = 1000._rk
      return
    end if

    Success = .true.

    select case (Grid%cart%nd)
    case (2)
      select case (Grid%grid_type)
      case (OVK_GRID_TYPE_CARTESIAN,OVK_GRID_TYPE_RECTILINEAR)
        CoordsInCell = ovkRectangleIsoInverseLinear(VertexCoords, Coords)
      case (OVK_GRID_TYPE_CURVILINEAR)
        CoordsInCell = ovkQuadIsoInverseLinear(VertexCoords, Coords, Success=Success)
      end select
    case (3)
      select case (Grid%grid_type)
      case (OVK_GRID_TYPE_CARTESIAN,OVK_GRID_TYPE_RECTILINEAR)
        CoordsInCell = ovkCuboidIsoInverseLinear(VertexCoords, Coords)
      case (OVK_GRID_TYPE_CURVILINEAR)
        CoordsInCell = ovkHexahedronIsoInverseLinear(VertexCoords, Coords, Success=Success)
      end select
    end select

    if (OVK_VERBOSE) then
      if (.not. Success) then
        write (ERROR_UNIT, '(6a)') "WARNING: Failed to compute isoparametric coordinates of ", &
          trim(CoordsToString(Coords)), " in cell ", trim(TupleToString(Cell)), &
          " of grid ", trim(IntToString(Grid%id))
      end if
    end if

  end function ovkCoordsInCell

  function ovkCellSize(Grid, Cell) result(CellSize)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%cart%nd), intent(in) :: Cell
    real(rk) :: CellSize

    integer, dimension(MAX_ND) :: PaddedCell

    PaddedCell(:Grid%cart%nd) = Cell
    PaddedCell(Grid%cart%nd+1:) = 1

    CellSize = Grid%cell_sizes%values(PaddedCell(1),PaddedCell(2),PaddedCell(3))

  end function ovkCellSize

  function ovkAvgCellSizeAroundPoint(Grid, Point) result(AvgCellSize)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%cart%nd), intent(in) :: Point
    real(rk) :: AvgCellSize

    integer :: i, j, k
    integer, dimension(MAX_ND) :: NeighborCellLower, NeighborCellUpper
    integer :: NumCells
    logical :: AwayFromBoundary
    integer, dimension(MAX_ND) :: NeighborCell

    NeighborCellLower(:Grid%cart%nd) = Point-1
    NeighborCellLower(Grid%cart%nd+1:) = 1
    NeighborCellUpper(:Grid%cart%nd) = Point
    NeighborCellUpper(Grid%cart%nd+1:) = 1

    AvgCellSize = 0._rk
    NumCells = 0

    AwayFromBoundary = ovkCartContains(Grid%cell_cart, NeighborCellLower) .and. &
      ovkCartContains(Grid%cell_cart, NeighborCellUpper)

    if (AwayFromBoundary) then
      do k = NeighborCellLower(3), NeighborCellUpper(3)
        do j = NeighborCellLower(2), NeighborCellUpper(2)
          do i = NeighborCellLower(1), NeighborCellUpper(1)
            NeighborCell = [i,j,k]
            AvgCellSize = AvgCellSize + Grid%cell_sizes%values(NeighborCell(1),NeighborCell(2), &
              NeighborCell(3))
            NumCells = NumCells + 1
          end do
        end do
      end do
    else
      do k = NeighborCellLower(3), NeighborCellUpper(3)
        do j = NeighborCellLower(2), NeighborCellUpper(2)
          do i = NeighborCellLower(1), NeighborCellUpper(1)
            NeighborCell = [i,j,k]
            NeighborCell(:Grid%cart%nd) = ovkCartPeriodicAdjust(Grid%cell_cart, NeighborCell)
            if (ovkCartContains(Grid%cell_cart, NeighborCell)) then
              AvgCellSize = AvgCellSize + Grid%cell_sizes%values(NeighborCell(1),NeighborCell(2), &
                NeighborCell(3))
              NumCells = NumCells + 1
            end if
          end do
        end do
      end do
    end if

    AvgCellSize = AvgCellSize/real(NumCells, kind=rk)

  end function ovkAvgCellSizeAroundPoint

  subroutine ovkGenerateBBOverlapMask(Grid, Bounds, BBOverlapMask)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_bbox), intent(in) :: Bounds
    type(ovk_field_logical), intent(out) :: BBOverlapMask

    integer :: i, j, k, l
    real(rk), dimension(Grid%cart%nd) :: Coords

    BBOverlapMask = ovk_field_logical_(Grid%cart)

    do k = Grid%cart%is(3), Grid%cart%ie(3)
      do j = Grid%cart%is(2), Grid%cart%ie(2)
        do i = Grid%cart%is(1), Grid%cart%ie(1)
          do l = 1, Grid%cart%nd
            Coords(l) = Grid%xyz(l)%values(i,j,k)
          end do
          BBOverlapMask%values(i,j,k) = ovkBBContainsPoint(Bounds, Coords)
        end do
      end do
    end do

  end subroutine ovkGenerateBBOverlapMask

  pure function ovkPeriodicExtend(Cart, PeriodicLength, Point, Coords) result(ExtendedCoords)

    type(ovk_cart), intent(in) :: Cart
    real(rk), dimension(Cart%nd), intent(in) :: PeriodicLength
    integer, dimension(Cart%nd), intent(in) :: Point
    real(rk), dimension(Cart%nd), intent(in) :: Coords
    real(rk), dimension(Cart%nd) :: ExtendedCoords

    real(rk), dimension(Cart%nd) :: PositiveAdjustment, NegativeAdjustment

    PositiveAdjustment = merge(PeriodicLength, 0._rk, Cart%periodic(:Cart%nd) .and. &
      Point > Cart%ie(:Cart%nd))
    NegativeAdjustment = merge(PeriodicLength, 0._rk, Cart%periodic(:Cart%nd) .and. &
      Point < Cart%is(:Cart%nd))

    ExtendedCoords = Coords + PositiveAdjustment - NegativeAdjustment

  end function ovkPeriodicExtend

end module ovkGrid
