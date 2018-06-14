! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkOverlap

  use ovkArray
  use ovkBoundingBox
  use ovkCart
  use ovkOverlapAccel
  use ovkField
  use ovkFieldOps
  use ovkGeometry
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
  public :: ovkFindOverlappingPoints
  public :: ovkFindOverlappedPoints
  public :: ovkOverlapCollect
  public :: ovkOverlapDisperse
  public :: ovkGetOverlapPropertyOverlappingGridID
  public :: ovkGetOverlapPropertyOverlappedGridID
  public :: ovkGetOverlapPropertyDimension
  public :: ovkGetOverlapPropertySize
  public :: ovkGetOverlapPropertyPeriodicity
  public :: ovkGetOverlapPropertyPeriodicStorage
  public :: ovkGetOverlapPropertyNumOverlapped
  public :: OVK_COLLECT_SIMPLE
  public :: OVK_COLLECT_MIN
  public :: OVK_COLLECT_MAX
  public :: OVK_COLLECT_NONE
  public :: OVK_COLLECT_ANY
  public :: OVK_COLLECT_NOT_ALL
  public :: OVK_COLLECT_ALL
  public :: OVK_COLLECT_INTERPOLATE
  public :: OVK_DISPERSE_OVERWRITE

  ! Internal
  public :: ovk_overlap_
  public :: CreateOverlap
  public :: DestroyOverlap
  public :: OverlapExists
  public :: DetectOverlap
  public :: ResetOverlap
  public :: UpdateOverlapAfterCut

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
  end type ovk_overlap

  interface ovkOverlapCollect
    module procedure ovkOverlapCollect_Integer
    module procedure ovkOverlapCollect_LargeInteger
    module procedure ovkOverlapCollect_Real
    module procedure ovkOverlapCollect_Logical
  end interface ovkOverlapCollect

  interface ovkOverlapDisperse
    module procedure ovkOverlapDisperse_Integer
    module procedure ovkOverlapDisperse_LargeInteger
    module procedure ovkOverlapDisperse_Real
    module procedure ovkOverlapDisperse_Logical
  end interface ovkOverlapDisperse

  integer, parameter :: OVK_COLLECT_SIMPLE = 1
  integer, parameter :: OVK_COLLECT_MIN = 2
  integer, parameter :: OVK_COLLECT_MAX = 3
  integer, parameter :: OVK_COLLECT_NONE = 4
  integer, parameter :: OVK_COLLECT_ANY = 5
  integer, parameter :: OVK_COLLECT_NOT_ALL = 6
  integer, parameter :: OVK_COLLECT_ALL = 7
  integer, parameter :: OVK_COLLECT_INTERPOLATE = 8

  integer, parameter :: OVK_DISPERSE_OVERWRITE = 1

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

  end subroutine CreateOverlap

  subroutine DestroyOverlap(Overlap)

    type(ovk_overlap), intent(inout) :: Overlap

    if (.not. OverlapExists(Overlap)) return

    deallocate(Overlap%mask)

    deallocate(Overlap%cells)
    deallocate(Overlap%coords)

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

    Overlap%properties%noverlap = 0_lk

    allocate(Overlap%cells(MAX_ND,0))
    allocate(Overlap%coords(Overlap%cart%nd,0))

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
    type(ovk_field_logical) :: HoleMask
    type(ovk_field_logical) :: OverlappedByHoleMask
    integer(lk) :: NumOverlappedPoints

    if (Overlap%properties%noverlap == 0_lk) return

    OldOverlapMask = Overlap%mask
    OldCells => Overlap%cells
    OldCoords => Overlap%coords

    HoleMask = ovk_field_logical_(OverlappingGrid%cart)
    HoleMask%values = .not. OverlappingGrid%mask%values

    call ovkFindOverlappedPoints(OverlappingGrid, OverlappedGrid, Overlap, HoleMask, &
      OverlappedByHoleMask)

    Overlap%mask%values = Overlap%mask%values .and. OverlappedGrid%mask%values .and. &
      .not. OverlappedByHoleMask%values

    NumOverlappedPoints = ovkCountMask(Overlap%mask)

    allocate(Overlap%cells(MAX_ND,NumOverlappedPoints))
    allocate(Overlap%coords(OverlappedGrid%cart%nd,NumOverlappedPoints))

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
              l = l + 1_lk
            end if
            l_old = l_old + 1_lk
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
    type(ovk_cart) :: PrincipalCart
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

    PrincipalCart = ovkCartConvertPeriodicStorage(OverlappingGrid%cart, OVK_NO_OVERLAP_PERIODIC)

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  FIRSTPRIVATE(PrincipalCart) &
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
              AwayFromEdge = ovkCartContains(PrincipalCart, CellLower) .and. &
                ovkCartContains(PrincipalCart, CellUpper)
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
                      Vertex(:OverlappingGrid%cart%nd) = ovkCartPeriodicAdjust(PrincipalCart, &
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

    call ovkFieldPeriodicFill(OverlappingMask, PrincipalCart)

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

  subroutine ovkOverlapCollect_Integer(OverlappingGrid, Overlap, CollectOp, OverlappingGridData, &
    CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: CollectOp
    type(ovk_field_int), intent(in) :: OverlappingGridData
    type(ovk_array_int), intent(out) :: CollectedData

    select case (CollectOp)
    case (OVK_COLLECT_SIMPLE)
      call CollectSimple_Integer(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MIN)
      call CollectMin_Integer(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MAX)
      call CollectMax_Integer(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid collect operation."
      stop 1
    end select

  end subroutine ovkOverlapCollect_Integer

  subroutine ovkOverlapCollect_LargeInteger(OverlappingGrid, Overlap, CollectOp, &
    OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: CollectOp
    type(ovk_field_large_int), intent(in) :: OverlappingGridData
    type(ovk_array_large_int), intent(out) :: CollectedData

    select case (CollectOp)
    case (OVK_COLLECT_SIMPLE)
      call CollectSimple_LargeInteger(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MIN)
      call CollectMin_LargeInteger(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MAX)
      call CollectMax_LargeInteger(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid collect operation."
      stop 1
    end select

  end subroutine ovkOverlapCollect_LargeInteger

  subroutine ovkOverlapCollect_Real(OverlappingGrid, Overlap, CollectOp, OverlappingGridData, &
    CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: CollectOp
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    select case (CollectOp)
    case (OVK_COLLECT_SIMPLE)
      call CollectSimple_Real(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MIN)
      call CollectMin_Real(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MAX)
      call CollectMax_Real(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_INTERPOLATE)
      call CollectInterpolate(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid collect operation."
      stop 1
    end select

  end subroutine ovkOverlapCollect_Real

  subroutine ovkOverlapCollect_Logical(OverlappingGrid, Overlap, CollectOp, OverlappingGridData, &
    CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: CollectOp
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    select case (CollectOp)
    case (OVK_COLLECT_SIMPLE)
      call CollectSimple_Logical(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_NONE)
      call CollectNone(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_ANY)
      call CollectAny(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_NOT_ALL)
      call CollectNotAll(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_ALL)
      call CollectAll(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid collect operation."
      stop 1
    end select

  end subroutine ovkOverlapCollect_Logical

  subroutine CollectSimple_Integer(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_int), intent(in) :: OverlappingGridData
    type(ovk_array_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Cell

    CollectedData = ovk_array_int_(Overlap%properties%noverlap)

    do l = 1_lk, Overlap%properties%noverlap
      Cell = Overlap%cells(:,l)
      CollectedData%values(l) = OverlappingGridData%values(Cell(1),Cell(2),Cell(3))
    end do

  end subroutine CollectSimple_Integer

  subroutine CollectSimple_LargeInteger(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_large_int), intent(in) :: OverlappingGridData
    type(ovk_array_large_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Cell

    CollectedData = ovk_array_large_int_(Overlap%properties%noverlap)

    do l = 1_lk, Overlap%properties%noverlap
      Cell = Overlap%cells(:,l)
      CollectedData%values(l) = OverlappingGridData%values(Cell(1),Cell(2),Cell(3))
    end do

  end subroutine CollectSimple_LargeInteger

  subroutine CollectSimple_Real(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Cell

    CollectedData = ovk_array_real_(Overlap%properties%noverlap)

    do l = 1_lk, Overlap%properties%noverlap
      Cell = Overlap%cells(:,l)
      CollectedData%values(l) = OverlappingGridData%values(Cell(1),Cell(2),Cell(3))
    end do

  end subroutine CollectSimple_Real

  subroutine CollectSimple_Logical(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Cell

    CollectedData = ovk_array_logical_(Overlap%properties%noverlap)

    do l = 1_lk, Overlap%properties%noverlap
      Cell = Overlap%cells(:,l)
      CollectedData%values(l) = OverlappingGridData%values(Cell(1),Cell(2),Cell(3))
    end do

  end subroutine CollectSimple_Logical

  subroutine CollectMin_Integer(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_int), intent(in) :: OverlappingGridData
    type(ovk_array_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    integer, dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_int_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = minval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMin_Integer

  subroutine CollectMin_LargeInteger(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_large_int), intent(in) :: OverlappingGridData
    type(ovk_array_large_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    integer(lk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_large_int_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = minval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMin_LargeInteger

  subroutine CollectMin_Real(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    real(rk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_real_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = minval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMin_Real

  subroutine CollectMax_Integer(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_int), intent(in) :: OverlappingGridData
    type(ovk_array_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    integer, dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_int_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = maxval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMax_Integer

  subroutine CollectMax_LargeInteger(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_large_int), intent(in) :: OverlappingGridData
    type(ovk_array_large_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    integer(lk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_large_int_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = maxval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMax_LargeInteger

  subroutine CollectMax_Real(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    real(rk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_real_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = maxval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMax_Real

  subroutine CollectInterpolate(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    integer :: i, j, k
    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    real(rk), dimension(MAX_ND,0:1) :: InterpBasis
    real(rk), dimension(0:1,0:1,0:1) :: VertexData
    real(rk) :: Weight

    CollectedData = ovk_array_real_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      InterpBasis(1,:) = ovkInterpBasisLinear(Overlap%coords(1,l))
      InterpBasis(2,:) = ovkInterpBasisLinear(Overlap%coords(2,l))
      if (OverlappingGrid%cart%nd == 3) then
        InterpBasis(3,:) = ovkInterpBasisLinear(Overlap%coords(3,l))
      else
        InterpBasis(3,:) = ovkInterpBasisLinear(0._rk)
      end if
      CollectedData%values(l) = 0._rk
      do k = Lower(3), Upper(3)
        do j = Lower(2), Upper(2)
          do i = Lower(1), Upper(1)
            Weight = InterpBasis(1,i) * InterpBasis(2,j) * InterpBasis(3,k)
            CollectedData%values(l) = CollectedData%values(l) + VertexData(i,j,k) * Weight
          end do
        end do
      end do
    end do

  end subroutine CollectInterpolate

  subroutine CollectNone(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    logical(bk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_logical_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = .not. any(VertexData)
    end do

  end subroutine CollectNone

  subroutine CollectAny(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    logical(bk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_logical_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = any(VertexData)
    end do

  end subroutine CollectAny

  subroutine CollectNotAll(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    logical(bk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_logical_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = .not. all(VertexData)
    end do

  end subroutine CollectNotAll

  subroutine CollectAll(OverlappingGrid, Overlap, OverlappingGridData, CollectedData)

    type(ovk_grid), intent(in) :: OverlappingGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    logical(bk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_logical_(Overlap%properties%noverlap)

    Lower = 0
    Upper(:OverlappingGrid%cart%nd) = 1
    Upper(OverlappingGrid%cart%nd+1:) = 0

    do l = 1_lk, Overlap%properties%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = all(VertexData)
    end do

  end subroutine CollectAll

  subroutine ovkOverlapDisperse_Integer(OverlappedGrid, Overlap, DisperseOp, CollectedData, &
    OverlappedGridData)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: DisperseOp
    type(ovk_array_int), intent(in) :: CollectedData
    type(ovk_field_int), intent(inout) :: OverlappedGridData

    select case (DisperseOp)
    case (OVK_DISPERSE_OVERWRITE)
      call DisperseOverwrite_Integer(OverlappedGrid, Overlap, CollectedData, OverlappedGridData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid disperse operation."
      stop 1
    end select

  end subroutine ovkOverlapDisperse_Integer

  subroutine ovkOverlapDisperse_LargeInteger(OverlappedGrid, Overlap, DisperseOp, &
    CollectedData, OverlappedGridData)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: DisperseOp
    type(ovk_array_large_int), intent(in) :: CollectedData
    type(ovk_field_large_int), intent(inout) :: OverlappedGridData

    select case (DisperseOp)
    case (OVK_DISPERSE_OVERWRITE)
      call DisperseOverwrite_LargeInteger(OverlappedGrid, Overlap, CollectedData, OverlappedGridData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid disperse operation."
      stop 1
    end select

  end subroutine ovkOverlapDisperse_LargeInteger

  subroutine ovkOverlapDisperse_Real(OverlappedGrid, Overlap, DisperseOp, CollectedData, &
    OverlappedGridData)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: DisperseOp
    type(ovk_array_real), intent(in) :: CollectedData
    type(ovk_field_real), intent(inout) :: OverlappedGridData

    select case (DisperseOp)
    case (OVK_DISPERSE_OVERWRITE)
      call DisperseOverwrite_Real(OverlappedGrid, Overlap, CollectedData, OverlappedGridData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid disperse operation."
      stop 1
    end select

  end subroutine ovkOverlapDisperse_Real

  subroutine ovkOverlapDisperse_Logical(OverlappedGrid, Overlap, DisperseOp, CollectedData, &
    OverlappedGridData)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: DisperseOp
    type(ovk_array_logical), intent(in) :: CollectedData
    type(ovk_field_logical), intent(inout) :: OverlappedGridData

    select case (DisperseOp)
    case (OVK_DISPERSE_OVERWRITE)
      call DisperseOverwrite_Logical(OverlappedGrid, Overlap, CollectedData, OverlappedGridData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid disperse operation."
      stop 1
    end select

  end subroutine ovkOverlapDisperse_Logical

  subroutine DisperseOverwrite_Integer(OverlappedGrid, Overlap, CollectedData, OverlappedGridData)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_int), intent(in) :: CollectedData
    type(ovk_field_int), intent(inout) :: OverlappedGridData

    integer :: i, j, k
    integer(lk) :: l

    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (Overlap%mask%values(i,j,k)) then
            OverlappedGridData%values(i,j,k) = CollectedData%values(l)
            l = l + 1_lk
          end if
        end do
      end do
    end do

  end subroutine DisperseOverwrite_Integer

  subroutine DisperseOverwrite_LargeInteger(OverlappedGrid, Overlap, CollectedData, OverlappedGridData)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_large_int), intent(in) :: CollectedData
    type(ovk_field_large_int), intent(inout) :: OverlappedGridData

    integer :: i, j, k
    integer(lk) :: l

    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (Overlap%mask%values(i,j,k)) then
            OverlappedGridData%values(i,j,k) = CollectedData%values(l)
            l = l + 1_lk
          end if
        end do
      end do
    end do

  end subroutine DisperseOverwrite_LargeInteger

  subroutine DisperseOverwrite_Real(OverlappedGrid, Overlap, CollectedData, OverlappedGridData)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_real), intent(in) :: CollectedData
    type(ovk_field_real), intent(inout) :: OverlappedGridData

    integer :: i, j, k
    integer(lk) :: l

    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (Overlap%mask%values(i,j,k)) then
            OverlappedGridData%values(i,j,k) = CollectedData%values(l)
            l = l + 1_lk
          end if
        end do
      end do
    end do

  end subroutine DisperseOverwrite_Real

  subroutine DisperseOverwrite_Logical(OverlappedGrid, Overlap, CollectedData, OverlappedGridData)

    type(ovk_grid), intent(in) :: OverlappedGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_logical), intent(in) :: CollectedData
    type(ovk_field_logical), intent(inout) :: OverlappedGridData

    integer :: i, j, k
    integer(lk) :: l

    l = 1_lk
    do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
      do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
        do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
          if (Overlap%mask%values(i,j,k)) then
            OverlappedGridData%values(i,j,k) = CollectedData%values(l)
            l = l + 1_lk
          end if
        end do
      end do
    end do

  end subroutine DisperseOverwrite_Logical

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
