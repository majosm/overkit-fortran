! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkOverlap

  use ovkBoundingBox
  use ovkCart
  use ovkOverlapAccel
  use ovkField
  use ovkFieldOps
  use ovkGlobal
  use ovkGrid
  use ovkLogger
  implicit none

  private

  ! API
  public :: ovk_overlap
  public :: ovk_overlap_properties
  public :: ovkGetOverlapProperties
  public :: ovkGetOverlapCart
  public :: ovkGetOverlapBounds
  public :: ovkGetOverlapMask
  public :: ovkGetOverlapCells
  public :: ovkGetOverlapCoords
  public :: ovkGetOverlapResolutions
  public :: ovkGetOverlapEdgeDistances
  public :: ovkFindOverlappingPoints
  public :: ovkFindOverlappedPoints
  public :: ovkGetOverlapPropertyOverlappingGridID
  public :: ovkGetOverlapPropertyOverlappedGridID
  public :: ovkGetOverlapPropertyDimension
  public :: ovkGetOverlapPropertySize
  public :: ovkGetOverlapPropertyPeriodicity
  public :: ovkGetOverlapPropertyPeriodicStorage
  public :: ovkGetOverlapPropertyNumOverlapped

  ! Internal
  public :: ovk_overlap_
  public :: CreateOverlap
  public :: DestroyOverlap
  public :: OverlapExists
  public :: DetectOverlap
  public :: ResetOverlap
  public :: UpdateOverlapAfterCut
  public :: FindCoarsePoints

  type ovk_overlap_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: overlapping_grid_id
    integer :: overlapped_grid_id
    integer :: nd
    integer, dimension(MAX_ND) :: npoints
    logical, dimension(MAX_ND) :: periodic
    integer :: periodic_storage
    integer(lk) :: noverlap
  end type ovk_overlap_properties

  type ovk_overlap
    type(t_noconstruct) :: noconstruct
    type(ovk_overlap_properties), pointer :: properties
    type(t_logger), pointer :: logger
    type(ovk_cart) :: cart
    type(ovk_bbox) :: bounds
    type(ovk_field_logical), pointer :: mask
    integer, dimension(:,:), pointer :: cells
    real(rk), dimension(:,:), pointer :: coords
    real(rk), dimension(:), pointer :: resolutions
    integer, dimension(:), pointer :: edge_dists
  end type ovk_overlap

contains

  pure function ovk_overlap_() result(Overlap)

    type(ovk_overlap) :: Overlap

    nullify(Overlap%properties)
    nullify(Overlap%logger)
    Overlap%cart = ovk_cart_()
    Overlap%bounds = ovk_bbox_()
    nullify(Overlap%mask)
    nullify(Overlap%cells)
    nullify(Overlap%coords)
    nullify(Overlap%resolutions)
    nullify(Overlap%edge_dists)

  end function ovk_overlap_

  subroutine CreateOverlap(Overlap, OverlappingGridID, OverlappedGridID, Logger, Cart)

    type(ovk_overlap), intent(out) :: Overlap
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    type(t_logger), pointer, intent(in) :: Logger
    type(ovk_cart), intent(in) :: Cart

    allocate(Overlap%properties)
    Overlap%properties = ovk_overlap_properties_(Cart%nd)
    Overlap%properties%overlapping_grid_id = OverlappingGridID
    Overlap%properties%overlapped_grid_id = OverlappedGridID
    Overlap%properties%npoints(:Cart%nd) = ovkCartSize(Cart)
    Overlap%properties%periodic = Cart%periodic
    Overlap%properties%periodic_storage = Cart%periodic_storage

    Overlap%logger => Logger

    Overlap%cart = Cart

    Overlap%bounds = ovk_bbox_(Cart%nd)

    allocate(Overlap%mask)
    Overlap%mask = ovk_field_logical_(Cart, .false.)

    allocate(Overlap%cells(MAX_ND,0))
    allocate(Overlap%coords(Cart%nd,0))
    allocate(Overlap%resolutions(0))
    allocate(Overlap%edge_dists(0))

  end subroutine CreateOverlap

  subroutine DestroyOverlap(Overlap)

    type(ovk_overlap), intent(inout) :: Overlap

    if (.not. OverlapExists(Overlap)) return

    deallocate(Overlap%mask)

    deallocate(Overlap%cells)
    deallocate(Overlap%coords)
    deallocate(Overlap%resolutions)
    deallocate(Overlap%edge_dists)

    deallocate(Overlap%properties)

  end subroutine DestroyOverlap

  function OverlapExists(Overlap) result(Exists)

    type(ovk_overlap), intent(in) :: Overlap
    logical :: Exists

    Exists = associated(Overlap%properties)

  end function OverlapExists

  subroutine ovkGetOverlapProperties(Overlap, Properties)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_overlap_properties), pointer, intent(out) :: Properties

    Properties => Overlap%properties

  end subroutine ovkGetOverlapProperties

  subroutine ovkGetOverlapCart(Overlap, Cart)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_cart), intent(out) :: Cart

    Cart = Overlap%cart

  end subroutine ovkGetOverlapCart

  subroutine ovkGetOverlapBounds(Overlap, Bounds)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_bbox), intent(out) :: Bounds

    Bounds = Overlap%bounds

  end subroutine ovkGetOverlapBounds

  subroutine ovkGetOverlapMask(Overlap, OverlapMask)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), pointer, intent(out) :: OverlapMask

    OverlapMask => Overlap%mask

  end subroutine ovkGetOverlapMask

  subroutine ovkGetOverlapCells(Overlap, Cells)

    type(ovk_overlap), intent(in) :: Overlap
    integer, dimension(:,:), pointer, intent(out) :: Cells

    Cells => Overlap%cells

  end subroutine ovkGetOverlapCells

  subroutine ovkGetOverlapCoords(Overlap, Coords)

    type(ovk_overlap), intent(in) :: Overlap
    real(rk), dimension(:,:), pointer, intent(out) :: Coords

    Coords => Overlap%coords

  end subroutine ovkGetOverlapCoords

  subroutine ovkGetOverlapResolutions(Overlap, Resolutions)

    type(ovk_overlap), intent(in) :: Overlap
    real(rk), dimension(:), pointer, intent(out) :: Resolutions

    Resolutions => Overlap%resolutions

  end subroutine ovkGetOverlapResolutions

  subroutine ovkGetOverlapEdgeDistances(Overlap, EdgeDistances)

    type(ovk_overlap), intent(in) :: Overlap
    integer, dimension(:), pointer, intent(out) :: EdgeDistances

    EdgeDistances => Overlap%edge_dists

  end subroutine ovkGetOverlapEdgeDistances

  subroutine DetectOverlap(OverlappingGrid, OverlappedGrid, OverlapAccel, OverlapBounds, &
    OverlapTolerance, Overlap)

    type(ovk_grid), intent(in) :: OverlappingGrid, OverlappedGrid
    type(t_overlap_accel), intent(in) :: OverlapAccel
    type(ovk_bbox), intent(in) :: OverlapBounds
    real(rk), intent(in) :: OverlapTolerance
    type(ovk_overlap), intent(inout) :: Overlap

    integer :: d, i, j, k
    integer(lk) :: l
    type(ovk_bbox) :: Bounds
    type(ovk_field_large_int) :: OverlappingCells
    integer(lk) :: NumOverlappedPoints
    type(ovk_field_large_int) :: OverlapIndices
    real(rk), dimension(OverlappedGrid%cart%nd) :: OverlappedCoords
    integer, dimension(MAX_ND) :: Cell
    real(rk), dimension(OverlappingGrid%cart%nd) :: CoordsInCell

    Overlap%mask%values = .false.

    deallocate(Overlap%cells)
    deallocate(Overlap%coords)
    deallocate(Overlap%resolutions)
    deallocate(Overlap%edge_dists)

    Bounds = ovkBBScale(OverlapBounds, 1._rk + OverlapTolerance)

    NumOverlappedPoints = 0_lk

    if (.not. ovkBBIsEmpty(Bounds)) then

      OverlappingCells = ovk_field_large_int_(OverlappedGrid%cart)

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  FIRSTPRIVATE(Bounds, OverlapTolerance) &
!$OMP&  SHARED(OverlappingGrid, OverlappedGrid, OverlapAccel, Overlap, OverlappingCells)
      do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
        do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
          do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
            if (OverlappedGrid%mask%values(i,j,k)) then
              do d = 1, OverlappedGrid%cart%nd
                OverlappedCoords(d) = OverlappedGrid%coords(d)%values(i,j,k)
              end do
              if (ovkBBContainsPoint(Bounds, OverlappedCoords)) then
                Cell(:OverlappingGrid%cart%nd) = FindOverlappingCell(OverlappingGrid, &
                  OverlapAccel, OverlappedCoords, OverlapTolerance)
                Cell(OverlappingGrid%cart%nd+1:) = 1
                if (ovkCartContains(OverlappingGrid%cell_cart, Cell)) then
                  Overlap%mask%values(i,j,k) = .true.
                  OverlappingCells%values(i,j,k) = ovkCartTupleToIndex(OverlappingGrid%cart, Cell)
                end if
              end if
            end if
          end do
        end do
      end do
!$OMP END PARALLEL DO

      NumOverlappedPoints = ovkCountMask(Overlap%mask)

    end if

    Overlap%properties%noverlap = NumOverlappedPoints

    allocate(Overlap%cells(MAX_ND,NumOverlappedPoints))
    allocate(Overlap%coords(Overlap%cart%nd,NumOverlappedPoints))
    allocate(Overlap%resolutions(NumOverlappedPoints))
    allocate(Overlap%edge_dists(NumOverlappedPoints))

    Overlap%bounds = ovk_bbox_(Overlap%cart%nd)

    if (NumOverlappedPoints > 0_lk) then

      OverlapIndices = ovk_field_large_int_(Overlap%cart)

      l = 1_lk
      do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
        do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
          do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
            if (Overlap%mask%values(i,j,k)) then
              OverlapIndices%values(i,j,k) = l
              l = l + 1_lk
            end if
          end do
        end do
      end do

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  SHARED(OverlappingGrid, OverlappedGrid, OverlapAccel, Overlap, OverlappingCells, OverlapIndices)
      do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
        do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
          do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
            if (Overlap%mask%values(i,j,k)) then
              l = OverlapIndices%values(i,j,k)
              do d = 1, OverlappedGrid%cart%nd
                OverlappedCoords(d) = OverlappedGrid%coords(d)%values(i,j,k)
              end do
              Overlap%bounds = ovkBBExtend(Overlap%bounds, OverlappedCoords)
              Cell(:OverlappingGrid%cart%nd) = ovkCartIndexToTuple(OverlappingGrid%cart, &
                OverlappingCells%values(i,j,k))
              Cell(OverlappingGrid%cart%nd+1:) = 1
              CoordsInCell = ovkCoordsInGridCell(OverlappingGrid, Cell, OverlappedCoords)
              Overlap%cells(:,l) = Cell
              Overlap%coords(:,l) = CoordsInCell
              Overlap%resolutions(l) = ovkGridResolution(OverlappingGrid, Cell, CoordsInCell)
              Overlap%edge_dists(l) = OverlappingGrid%cell_edge_dist%values(Cell(1),Cell(2),Cell(3))
            end if
          end do
        end do
      end do
!$OMP END PARALLEL DO

    end if

  end subroutine DetectOverlap

  subroutine ResetOverlap(Overlap)

    type(ovk_overlap), intent(inout) :: Overlap

    Overlap%mask%values = .false.

    deallocate(Overlap%cells)
    deallocate(Overlap%coords)
    deallocate(Overlap%resolutions)
    deallocate(Overlap%edge_dists)

    Overlap%properties%noverlap = 0_lk

    allocate(Overlap%cells(MAX_ND,0))
    allocate(Overlap%coords(Overlap%cart%nd,0))
    allocate(Overlap%resolutions(0))
    allocate(Overlap%edge_dists(0))

    Overlap%bounds = ovk_bbox_(Overlap%cart%nd)

  end subroutine ResetOverlap

  subroutine UpdateOverlapAfterCut(OverlappingGrid, OverlappedGrid, Overlap)

    type(ovk_grid), intent(in) :: OverlappingGrid, OverlappedGrid
    type(ovk_overlap), intent(inout) :: Overlap

    integer :: i, j, k
    integer(lk) :: l_old, l
    type(ovk_field_logical) :: OldOverlapMask
    integer, dimension(:,:), pointer :: OldCells
    real(rk), dimension(:,:), pointer :: OldCoords
    real(rk), dimension(:), pointer :: OldResolutions
    integer, dimension(:), pointer :: OldEdgeDistances
    type(ovk_field_logical) :: HoleMask
    type(ovk_field_logical) :: OverlappedByHoleMask
    integer(lk) :: NumOverlappedPoints
    integer, dimension(MAX_ND) :: Cell

    if (Overlap%properties%noverlap == 0_lk) return

    OldOverlapMask = Overlap%mask
    OldCells => Overlap%cells
    OldCoords => Overlap%coords
    OldResolutions => Overlap%resolutions
    OldEdgeDistances => Overlap%edge_dists

    HoleMask = ovk_field_logical_(OverlappingGrid%cart)
    HoleMask%values = .not. OverlappingGrid%mask%values

    call ovkFindOverlappedPoints(OverlappingGrid, OverlappedGrid, Overlap, HoleMask, &
      OverlappedByHoleMask)

    Overlap%mask%values = Overlap%mask%values .and. OverlappedGrid%mask%values .and. &
      .not. OverlappedByHoleMask%values

    NumOverlappedPoints = ovkCountMask(Overlap%mask)

    allocate(Overlap%cells(MAX_ND,NumOverlappedPoints))
    allocate(Overlap%coords(OverlappedGrid%cart%nd,NumOverlappedPoints))
    allocate(Overlap%resolutions(NumOverlappedPoints))
    allocate(Overlap%edge_dists(NumOverlappedPoints))

    Overlap%properties%noverlap = NumOverlappedPoints

    l_old = 1_lk
    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (OldOverlapMask%values(i,j,k)) then
            if (Overlap%mask%values(i,j,k)) then
              Overlap%cells(:,l) = OldCells(:,l_old)
              Overlap%coords(:,l) = OldCoords(:,l_old)
              Overlap%resolutions(l) = OldResolutions(l_old)
              Overlap%edge_dists(l) = OldEdgeDistances(l_old)
              l = l + 1_lk
            end if
            l_old = l_old + 1_lk
          end if
        end do
      end do
    end do

    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (Overlap%mask%values(i,j,k)) then
            Cell = Overlap%cells(:,l)
            Overlap%edge_dists(l) = OverlappingGrid%cell_edge_dist%values(Cell(1),Cell(2),Cell(3))
            l = l + 1_lk
          end if
        end do
      end do
    end do

  end subroutine UpdateOverlapAfterCut

  subroutine ovkFindOverlappingPoints(OverlappingGrid, OverlappedGrid, Overlap, OverlappedSubset, &
    OverlappingMask)

    type(ovk_grid), intent(in) :: OverlappingGrid, OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappedSubset
    type(ovk_field_logical), intent(out) :: OverlappingMask

    integer :: i, j, k, m, n, o
    integer(lk) :: l
    type(ovk_field_large_int) :: OverlapIndices
    integer, dimension(MAX_ND) :: CellLower
    integer, dimension(MAX_ND) :: CellUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Vertex

    OverlappingMask = ovk_field_logical_(OverlappingGrid%cart, .false.)

    OverlapIndices = ovk_field_large_int_(OverlappedGrid%cart)

    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (Overlap%mask%values(i,j,k)) then
            OverlapIndices%values(i,j,k) = l
            l = l + 1_lk
          end if
        end do
      end do
    end do

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  SHARED(OverlappingGrid, OverlappedGrid, Overlap, OverlappedSubset, OverlapIndices, OverlappingMask)
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)

          if (Overlap%mask%values(i,j,k)) then

            if (OverlappedSubset%values(i,j,k)) then

              l = OverlapIndices%values(i,j,k)

              CellLower = Overlap%cells(:,l)
              CellUpper(:OverlappingGrid%cart%nd) = CellLower(:OverlappingGrid%cart%nd)+1
              CellUpper(OverlappingGrid%cart%nd+1:) = 1

              AwayFromEdge = ovkCartContains(OverlappingGrid%cart, CellLower) .and. &
                ovkCartContains(OverlappingGrid%cart, CellUpper)

              if (AwayFromEdge) then

                do o = CellLower(3), CellUpper(3)
                  do n = CellLower(2), CellUpper(2)
                    do m = CellLower(1), CellUpper(1)
                      Vertex = [m,n,o]
                      OverlappingMask%values(Vertex(1),Vertex(2),Vertex(3)) = .true.
                    end do
                  end do
                end do

              else

                do o = CellLower(3), CellUpper(3)
                  do n = CellLower(2), CellUpper(2)
                    do m = CellLower(1), CellUpper(1)
                      Vertex = [m,n,o]
                      Vertex(:OverlappingGrid%cart%nd) = ovkCartPeriodicAdjust(OverlappingGrid%cart, &
                        Vertex)
                      OverlappingMask%values(Vertex(1),Vertex(2),Vertex(3)) = .true.
                    end do
                  end do
                end do

              end if

            end if

          end if

        end do
      end do
    end do
!$OMP END PARALLEL DO

    call ovkFieldPeriodicFill(OverlappingMask)

  end subroutine ovkFindOverlappingPoints

  subroutine ovkFindOverlappedPoints(OverlappingGrid, OverlappedGrid, Overlap, OverlappingSubset, &
    OverlappedMask)

    type(ovk_grid), intent(in) :: OverlappingGrid, OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingSubset
    type(ovk_field_logical), intent(out) :: OverlappedMask

    integer :: i, j, k, m, n, o
    integer(lk) :: l
    type(ovk_field_large_int) :: OverlapIndices
    integer, dimension(MAX_ND) :: CellLower
    integer, dimension(MAX_ND) :: CellUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Vertex

    OverlappedMask = ovk_field_logical_(OverlappedGrid%cart, .false.)

    OverlapIndices = ovk_field_large_int_(OverlappedGrid%cart)

    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (Overlap%mask%values(i,j,k)) then
            OverlapIndices%values(i,j,k) = l
            l = l + 1_lk
          end if
        end do
      end do
    end do

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  SHARED(OverlappingGrid, OverlappedGrid, Overlap, OverlappingSubset, OverlapIndices, OverlappedMask)
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)

          if (Overlap%mask%values(i,j,k)) then

            l = OverlapIndices%values(i,j,k)

            CellLower = Overlap%cells(:,l)
            CellUpper(:OverlappingGrid%cart%nd) = CellLower(:OverlappingGrid%cart%nd)+1
            CellUpper(OverlappingGrid%cart%nd+1:) = 1

            AwayFromEdge = ovkCartContains(OverlappingGrid%cart, CellLower) .and. &
              ovkCartContains(OverlappingGrid%cart, CellUpper)

            if (AwayFromEdge) then

              L1: &
              do o = CellLower(3), CellUpper(3)
                do n = CellLower(2), CellUpper(2)
                  do m = CellLower(1), CellUpper(1)
                    Vertex = [m,n,o]
                    if (OverlappingSubset%values(Vertex(1),Vertex(2),Vertex(3))) then
                      OverlappedMask%values(i,j,k) = .true.
                      exit L1
                    end if
                  end do
                end do
              end do L1

            else

              L2: &
              do o = CellLower(3), CellUpper(3)
                do n = CellLower(2), CellUpper(2)
                  do m = CellLower(1), CellUpper(1)
                    Vertex = [m,n,o]
                    Vertex(:OverlappingGrid%cart%nd) = ovkCartPeriodicAdjust(OverlappingGrid%cart, &
                      Vertex)
                    if (OverlappingSubset%values(Vertex(1),Vertex(2),Vertex(3))) then
                      OverlappedMask%values(i,j,k) = .true.
                      exit L2
                    end if
                  end do
                end do
              end do L2

            end if

          end if

        end do
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine ovkFindOverlappedPoints

  subroutine FindCoarsePoints(OverlappedGrid, Overlap, CoarseMask)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
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
              (1._rk-TOLERANCE) * Overlap%resolutions(l)
            l = l + 1_lk
          end if
        end do
      end do
    end do

  end subroutine FindCoarsePoints

  function ovk_overlap_properties_(NumDims) result(Properties)

    integer, intent(in) :: NumDims
    type(ovk_overlap_properties) :: Properties

    Properties%overlapping_grid_id = 0
    Properties%overlapped_grid_id = 0
    Properties%nd = NumDims
    Properties%npoints(:NumDims) = 0
    Properties%npoints(NumDims+1:) = 1
    Properties%periodic = .false.
    Properties%periodic_storage = OVK_NO_OVERLAP_PERIODIC
    Properties%noverlap = 0

  end function ovk_overlap_properties_

  subroutine ovkGetOverlapPropertyOverlappingGridID(Properties, OverlappingGridID)

    type(ovk_overlap_properties), intent(in) :: Properties
    integer, intent(out) :: OverlappingGridID

    OverlappingGridID = Properties%overlapping_grid_id

  end subroutine ovkGetOverlapPropertyOverlappingGridID

  subroutine ovkGetOverlapPropertyOverlappedGridID(Properties, OverlappedGridID)

    type(ovk_overlap_properties), intent(in) :: Properties
    integer, intent(out) :: OverlappedGridID

    OverlappedGridID = Properties%overlapped_grid_id

  end subroutine ovkGetOverlapPropertyOverlappedGridID

  subroutine ovkGetOverlapPropertyDimension(Properties, NumDims)

    type(ovk_overlap_properties), intent(in) :: Properties
    integer, intent(out) :: NumDims

    NumDims = Properties%nd

  end subroutine ovkGetOverlapPropertyDimension

  subroutine ovkGetOverlapPropertySize(Properties, NumPoints)

    type(ovk_overlap_properties), intent(in) :: Properties
    integer, dimension(Properties%nd), intent(out) :: NumPoints

    NumPoints = Properties%npoints(:Properties%nd)

  end subroutine ovkGetOverlapPropertySize

  subroutine ovkGetOverlapPropertyPeriodicity(Properties, Periodic)

    type(ovk_overlap_properties), intent(in) :: Properties
    logical, dimension(Properties%nd), intent(out) :: Periodic

    Periodic = Properties%periodic(:Properties%nd)

  end subroutine ovkGetOverlapPropertyPeriodicity

  subroutine ovkGetOverlapPropertyPeriodicStorage(Properties, PeriodicStorage)

    type(ovk_overlap_properties), intent(in) :: Properties
    integer, intent(out) :: PeriodicStorage

    PeriodicStorage = Properties%periodic_storage

  end subroutine ovkGetOverlapPropertyPeriodicStorage

  subroutine ovkGetOverlapPropertyNumOverlapped(Properties, NumOverlapped)

    type(ovk_overlap_properties), intent(in) :: Properties
    integer(lk), intent(out) :: NumOverlapped

    NumOverlapped = Properties%noverlap

  end subroutine ovkGetOverlapPropertyNumOverlapped

end module ovkOverlap
