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
  public :: ovkOverlapExists
  public :: ovkResetOverlap
  public :: ovkGetOverlapOverlappingGrid
  public :: ovkGetOverlapOverlappedGrid
  public :: ovkGetOverlapDimension
  public :: ovkGetOverlapCount
  public :: ovkGetOverlapMask
  public :: ovkGetOverlapCells
  public :: ovkGetOverlapCoords
  public :: ovkFindOverlappingPoints
  public :: ovkFindOverlappedPoints
  public :: ovkOverlapCollect
  public :: ovkOverlapDisperse
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
  public :: DetectOverlap
  public :: UpdateOverlapAfterCut

  type ovk_overlap
    type(t_noconstruct) :: noconstruct
    type(t_existence_flag) :: existence_flag
    type(t_logger), pointer :: logger
    type(ovk_grid), pointer :: overlapping_grid
    type(ovk_grid), pointer :: overlapped_grid
    integer :: nd
    integer(lk) :: noverlap
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

    nullify(Overlap%logger)
    nullify(Overlap%overlapping_grid)
    nullify(Overlap%overlapped_grid)
    Overlap%nd = 2
    Overlap%noverlap = 0_lk
    nullify(Overlap%mask)
    nullify(Overlap%cells)
    nullify(Overlap%coords)

    call SetExists(Overlap%existence_flag, .false.)

  end function ovk_overlap_

  subroutine CreateOverlap(Overlap, Logger, OverlappingGrid, OverlappedGrid)

    type(ovk_overlap), intent(out) :: Overlap
    type(t_logger), pointer, intent(in) :: Logger
    type(ovk_grid), pointer, intent(in) :: OverlappingGrid, OverlappedGrid

    integer :: NumDims

    NumDims = OverlappingGrid%nd

    Overlap%logger => Logger
    Overlap%overlapping_grid => OverlappingGrid
    Overlap%overlapped_grid => OverlappedGrid
    Overlap%nd = NumDims
    Overlap%noverlap = 0_lk

    allocate(Overlap%mask)
    Overlap%mask = ovk_field_logical_(OverlappedGrid%cart, .false.)

    allocate(Overlap%cells(MAX_ND,0))
    allocate(Overlap%coords(NumDims,0))

    call SetExists(Overlap%existence_flag, .true.)

  end subroutine CreateOverlap

  subroutine DestroyOverlap(Overlap)

    type(ovk_overlap), intent(inout) :: Overlap

    if (.not. ovkOverlapExists(Overlap)) return

    call SetExists(Overlap%existence_flag, .false.)

    deallocate(Overlap%mask)
    deallocate(Overlap%cells)
    deallocate(Overlap%coords)

  end subroutine DestroyOverlap

  function ovkOverlapExists(Overlap) result(Exists)

    type(ovk_overlap), intent(in) :: Overlap
    logical :: Exists

    Exists = CheckExists(Overlap%existence_flag)

  end function ovkOverlapExists

  subroutine ovkResetOverlap(Overlap, OverlappedMask)

    type(ovk_overlap), intent(inout) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappedMask

    Overlap%mask%values = OverlappedMask%values

    deallocate(Overlap%cells)
    deallocate(Overlap%coords)

    Overlap%noverlap = ovkCountMask(Overlap%mask)

    allocate(Overlap%cells(MAX_ND,Overlap%noverlap))

    Overlap%cells(Overlap%nd+1:,:) = 1

    allocate(Overlap%coords(Overlap%nd,Overlap%noverlap))

  end subroutine ovkResetOverlap

  subroutine ovkGetOverlapOverlappingGrid(Overlap, OverlappingGrid)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_grid), pointer, intent(out) :: OverlappingGrid

    OverlappingGrid => Overlap%overlapping_grid

  end subroutine ovkGetOverlapOverlappingGrid

  subroutine ovkGetOverlapOverlappedGrid(Overlap, OverlappedGrid)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_grid), pointer, intent(out) :: OverlappedGrid

    OverlappedGrid => Overlap%overlapped_grid

  end subroutine ovkGetOverlapOverlappedGrid

  subroutine ovkGetOverlapDimension(Overlap, NumDims)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(out) :: NumDims

    NumDims = Overlap%nd

  end subroutine ovkGetOverlapDimension

  subroutine ovkGetOverlapCount(Overlap, NumOverlapped)

    type(ovk_overlap), intent(in) :: Overlap
    integer(lk), intent(out) :: NumOverlapped

    NumOverlapped = Overlap%noverlap

  end subroutine ovkGetOverlapCount

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

  subroutine DetectOverlap(Overlap, OverlapAccel, OverlapBounds, OverlapTolerance)

    type(ovk_overlap), intent(inout) :: Overlap
    type(t_overlap_accel), intent(in) :: OverlapAccel
    type(ovk_bbox), intent(in) :: OverlapBounds
    real(rk), intent(in) :: OverlapTolerance

    integer :: d, i, j, k
    integer(lk) :: l
    integer :: NumDims
    type(ovk_grid), pointer :: OverlappingGrid, OverlappedGrid
    type(ovk_field_logical) :: OverlappedMask
    type(ovk_bbox) :: Bounds
    type(ovk_field_large_int) :: OverlappingCells
    real(rk), dimension(Overlap%nd) :: OverlappedCoords
    integer, dimension(MAX_ND) :: Cell
    integer, dimension(MAX_ND) :: Point
    real(rk), dimension(Overlap%nd) :: CoordsInCell
    logical :: Success
    integer :: NumWarnings
    character(len=STRING_LENGTH) :: PointString
    character(len=STRING_LENGTH) :: CellString
    character(len=STRING_LENGTH) :: OverlappingGridIDString
    character(len=STRING_LENGTH) :: OverlappedGridIDString

    NumDims = Overlap%nd
    OverlappingGrid => Overlap%overlapping_grid
    OverlappedGrid => Overlap%overlapped_grid

    OverlappedMask = ovk_field_logical_(OverlappedGrid%cart, .false.)

    Bounds = ovkBBScale(OverlapBounds, 1._rk + OverlapTolerance)

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
              do d = 1, NumDims
                OverlappedCoords(d) = OverlappedGrid%coords(d)%values(i,j,k)
              end do
              if (ovkBBContainsPoint(Bounds, OverlappedCoords)) then
                Cell(:NumDims) = FindOverlappingCell(OverlappingGrid, OverlapAccel, &
                  OverlappedCoords, OverlapTolerance)
                Cell(NumDims+1:) = 1
                if (ovkCartContains(OverlappingGrid%cell_cart, Cell)) then
                  OverlappedMask%values(i,j,k) = .true.
                  OverlappingCells%values(i,j,k) = ovkCartTupleToIndex(OverlappingGrid%cart, Cell)
                end if
              end if
            end if
          end do
        end do
      end do
!$OMP END PARALLEL DO

    end if

    call ovkResetOverlap(Overlap, OverlappedMask)

    if (Overlap%noverlap > 0_lk) then

      NumWarnings = 0
      l = 1_lk
      do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
        do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
          do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
            if (Overlap%mask%values(i,j,k)) then
              Point = [i,j,k]
              do d = 1, NumDims
                OverlappedCoords(d) = OverlappedGrid%coords(d)%values(i,j,k)
              end do
              Overlap%cells(:NumDims,l) = ovkCartIndexToTuple(OverlappingGrid%cart, &
                OverlappingCells%values(i,j,k))
              Overlap%cells(NumDims+1:,l) = 1
              CoordsInCell = ovkCoordsInGridCell(OverlappingGrid, Overlap%cells(:,l), &
                OverlappedCoords, Success=Success)
              if (Success) then
                Overlap%coords(:,l) = CoordsInCell
              else
                Overlap%coords(:,l) = 0.5_rk
                if (Overlap%logger%verbose) then
                  if (NumWarnings <= 100) then
                    PointString = TupleToString(Point(:NumDims))
                    CellString = TupleToString(Overlap%cells(:NumDims,l))
                    OverlappingGridIDString = IntToString(OverlappingGrid%id)
                    OverlappedGridIDString = IntToString(OverlappedGrid%id)
                    write (ERROR_UNIT, '(9a)') "WARNING: Failed to compute local coordinates of point ", &
                      trim(PointString), " of grid ", trim(OverlappedGridIDString), " inside cell ", &
                      trim(CellString), " of grid ", trim(OverlappingGridIDString), "."
                    if (NumWarnings == 100) then
                      write (ERROR_UNIT, '(a)') "WARNING: Further warnings suppressed."
                    end if
                    NumWarnings = NumWarnings + 1
                  end if
                end if
              end if
              l = l + 1_lk
            end if
          end do
        end do
      end do

    end if

  end subroutine DetectOverlap

  subroutine UpdateOverlapAfterCut(Overlap)

    type(ovk_overlap), intent(inout) :: Overlap

    integer :: i, j, k
    integer(lk) :: l_old, l
    type(ovk_grid), pointer :: OverlappingGrid, OverlappedGrid
    type(ovk_field_logical) :: OldOverlappedMask
    integer, dimension(:,:), allocatable :: OldCells
    real(rk), dimension(:,:), allocatable :: OldCoords
    type(ovk_field_logical) :: HoleMask
    type(ovk_field_logical) :: OverlappedByHoleMask
    type(ovk_field_logical) :: OverlappedMask

    if (Overlap%noverlap == 0_lk) return

    OverlappingGrid => Overlap%overlapping_grid
    OverlappedGrid => Overlap%overlapped_grid

    OldOverlappedMask = Overlap%mask

    allocate(OldCells(MAX_ND,Overlap%noverlap))
    OldCells = Overlap%cells

    allocate(OldCoords(Overlap%nd,Overlap%noverlap))
    OldCoords = Overlap%coords

    HoleMask = ovk_field_logical_(OverlappingGrid%cart)
    HoleMask%values = .not. OverlappingGrid%mask%values

    call ovkFindOverlappedPoints(Overlap, HoleMask, OverlappedByHoleMask)

    OverlappedMask = ovk_field_logical_(OverlappedGrid%cart)
    OverlappedMask%values = OldOverlappedMask%values .and. OverlappedGrid%mask%values .and. &
      .not. OverlappedByHoleMask%values

    call ovkResetOverlap(Overlap, OverlappedMask)

    if (Overlap%noverlap > 0_lk) then

      l_old = 1_lk
      l = 1_lk
      do k = OverlappedGrid%cart%is(3), OverlappedGrid%cart%ie(3)
        do j = OverlappedGrid%cart%is(2), OverlappedGrid%cart%ie(2)
          do i = OverlappedGrid%cart%is(1), OverlappedGrid%cart%ie(1)
            if (OldOverlappedMask%values(i,j,k)) then
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

    end if

  end subroutine UpdateOverlapAfterCut

  subroutine ovkFindOverlappingPoints(Overlap, OverlappedSubset, OverlappingMask)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappedSubset
    type(ovk_field_logical), intent(out) :: OverlappingMask

    integer :: i, j, k, m, n, o
    integer(lk) :: l
    type(ovk_grid), pointer :: OverlappingGrid, OverlappedGrid
    type(ovk_cart) :: PrincipalCart
    type(ovk_field_large_int) :: OverlapIndices
    integer, dimension(MAX_ND) :: CellLower
    integer, dimension(MAX_ND) :: CellUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Vertex

    OverlappingGrid => Overlap%overlapping_grid
    OverlappedGrid => Overlap%overlapped_grid

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
              CellUpper(:Overlap%nd) = CellLower(:Overlap%nd)+1
              CellUpper(Overlap%nd+1:) = 1
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
                      Vertex(:Overlap%nd) = ovkCartPeriodicAdjust(PrincipalCart, Vertex)
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

  subroutine ovkFindOverlappedPoints(Overlap, OverlappingSubset, OverlappedMask)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingSubset
    type(ovk_field_logical), intent(out) :: OverlappedMask

    integer :: i, j, k, m, n, o
    integer(lk) :: l
    type(ovk_grid), pointer :: OverlappingGrid, OverlappedGrid
    type(ovk_field_large_int) :: OverlapIndices
    integer, dimension(MAX_ND) :: CellLower
    integer, dimension(MAX_ND) :: CellUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Vertex

    OverlappingGrid => Overlap%overlapping_grid
    OverlappedGrid => Overlap%overlapped_grid

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
            CellUpper(:Overlap%nd) = CellLower(:Overlap%nd)+1
            CellUpper(Overlap%nd+1:) = 1
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
                    Vertex(:Overlap%nd) = ovkCartPeriodicAdjust(OverlappingGrid%cart, Vertex)
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

  subroutine ovkOverlapCollect_Integer(Overlap, CollectOp, OverlappingGridData, &
    CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: CollectOp
    type(ovk_field_int), intent(in) :: OverlappingGridData
    type(ovk_array_int), intent(out) :: CollectedData

    select case (CollectOp)
    case (OVK_COLLECT_SIMPLE)
      call CollectSimple_Integer(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MIN)
      call CollectMin_Integer(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MAX)
      call CollectMax_Integer(Overlap, OverlappingGridData, CollectedData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid collect operation."
      stop 1
    end select

  end subroutine ovkOverlapCollect_Integer

  subroutine ovkOverlapCollect_LargeInteger(Overlap, CollectOp, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: CollectOp
    type(ovk_field_large_int), intent(in) :: OverlappingGridData
    type(ovk_array_large_int), intent(out) :: CollectedData

    select case (CollectOp)
    case (OVK_COLLECT_SIMPLE)
      call CollectSimple_LargeInteger(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MIN)
      call CollectMin_LargeInteger(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MAX)
      call CollectMax_LargeInteger(Overlap, OverlappingGridData, CollectedData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid collect operation."
      stop 1
    end select

  end subroutine ovkOverlapCollect_LargeInteger

  subroutine ovkOverlapCollect_Real(Overlap, CollectOp, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: CollectOp
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    select case (CollectOp)
    case (OVK_COLLECT_SIMPLE)
      call CollectSimple_Real(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MIN)
      call CollectMin_Real(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_MAX)
      call CollectMax_Real(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_INTERPOLATE)
      call CollectInterpolate(Overlap, OverlappingGridData, CollectedData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid collect operation."
      stop 1
    end select

  end subroutine ovkOverlapCollect_Real

  subroutine ovkOverlapCollect_Logical(Overlap, CollectOp, OverlappingGridData, &
    CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: CollectOp
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    select case (CollectOp)
    case (OVK_COLLECT_SIMPLE)
      call CollectSimple_Logical(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_NONE)
      call CollectNone(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_ANY)
      call CollectAny(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_NOT_ALL)
      call CollectNotAll(Overlap, OverlappingGridData, CollectedData)
    case (OVK_COLLECT_ALL)
      call CollectAll(Overlap, OverlappingGridData, CollectedData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid collect operation."
      stop 1
    end select

  end subroutine ovkOverlapCollect_Logical

  subroutine CollectSimple_Integer(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_int), intent(in) :: OverlappingGridData
    type(ovk_array_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Cell

    CollectedData = ovk_array_int_(Overlap%noverlap)

    do l = 1_lk, Overlap%noverlap
      Cell = Overlap%cells(:,l)
      CollectedData%values(l) = OverlappingGridData%values(Cell(1),Cell(2),Cell(3))
    end do

  end subroutine CollectSimple_Integer

  subroutine CollectSimple_LargeInteger(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_large_int), intent(in) :: OverlappingGridData
    type(ovk_array_large_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Cell

    CollectedData = ovk_array_large_int_(Overlap%noverlap)

    do l = 1_lk, Overlap%noverlap
      Cell = Overlap%cells(:,l)
      CollectedData%values(l) = OverlappingGridData%values(Cell(1),Cell(2),Cell(3))
    end do

  end subroutine CollectSimple_LargeInteger

  subroutine CollectSimple_Real(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Cell

    CollectedData = ovk_array_real_(Overlap%noverlap)

    do l = 1_lk, Overlap%noverlap
      Cell = Overlap%cells(:,l)
      CollectedData%values(l) = OverlappingGridData%values(Cell(1),Cell(2),Cell(3))
    end do

  end subroutine CollectSimple_Real

  subroutine CollectSimple_Logical(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Cell

    CollectedData = ovk_array_logical_(Overlap%noverlap)

    do l = 1_lk, Overlap%noverlap
      Cell = Overlap%cells(:,l)
      CollectedData%values(l) = OverlappingGridData%values(Cell(1),Cell(2),Cell(3))
    end do

  end subroutine CollectSimple_Logical

  subroutine CollectMin_Integer(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_int), intent(in) :: OverlappingGridData
    type(ovk_array_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    integer, dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_int_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = minval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMin_Integer

  subroutine CollectMin_LargeInteger(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_large_int), intent(in) :: OverlappingGridData
    type(ovk_array_large_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    integer(lk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_large_int_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = minval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMin_LargeInteger

  subroutine CollectMin_Real(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    real(rk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_real_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = minval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMin_Real

  subroutine CollectMax_Integer(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_int), intent(in) :: OverlappingGridData
    type(ovk_array_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    integer, dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_int_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = maxval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMax_Integer

  subroutine CollectMax_LargeInteger(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_large_int), intent(in) :: OverlappingGridData
    type(ovk_array_large_int), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    integer(lk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_large_int_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = maxval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMax_LargeInteger

  subroutine CollectMax_Real(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_real), intent(in) :: OverlappingGridData
    type(ovk_array_real), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    real(rk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_real_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = maxval(VertexData(Lower(1):Upper(1),Lower(2):Upper(2), &
        Lower(3):Upper(3)))
    end do

  end subroutine CollectMax_Real

  subroutine CollectInterpolate(Overlap, OverlappingGridData, CollectedData)

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

    CollectedData = ovk_array_real_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      InterpBasis(1,:) = ovkInterpBasisLinear(Overlap%coords(1,l))
      InterpBasis(2,:) = ovkInterpBasisLinear(Overlap%coords(2,l))
      if (Overlap%nd == 3) then
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

  subroutine CollectNone(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    logical(bk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_logical_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = .not. any(VertexData)
    end do

  end subroutine CollectNone

  subroutine CollectAny(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    logical(bk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_logical_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = any(VertexData)
    end do

  end subroutine CollectAny

  subroutine CollectNotAll(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    logical(bk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_logical_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = .not. all(VertexData)
    end do

  end subroutine CollectNotAll

  subroutine CollectAll(Overlap, OverlappingGridData, CollectedData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_field_logical), intent(in) :: OverlappingGridData
    type(ovk_array_logical), intent(out) :: CollectedData

    integer(lk) :: l
    integer, dimension(MAX_ND) :: Lower, Upper
    integer, dimension(MAX_ND) :: VertexStart, VertexEnd
    logical(bk), dimension(0:1,0:1,0:1) :: VertexData

    CollectedData = ovk_array_logical_(Overlap%noverlap)

    Lower = 0
    Upper(:Overlap%nd) = 1
    Upper(Overlap%nd+1:) = 0

    do l = 1_lk, Overlap%noverlap
      VertexStart = Overlap%cells(:,l) + Lower
      VertexEnd = Overlap%cells(:,l) + Upper
      call ovkGetFieldPatch(OverlappingGridData, VertexStart, VertexEnd, VertexData)
      CollectedData%values(l) = all(VertexData)
    end do

  end subroutine CollectAll

  subroutine ovkOverlapDisperse_Integer(Overlap, DisperseOp, CollectedData, &
    OverlappedGridData)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: DisperseOp
    type(ovk_array_int), intent(in) :: CollectedData
    type(ovk_field_int), intent(inout) :: OverlappedGridData

    select case (DisperseOp)
    case (OVK_DISPERSE_OVERWRITE)
      call DisperseOverwrite_Integer(Overlap, CollectedData, OverlappedGridData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid disperse operation."
      stop 1
    end select

  end subroutine ovkOverlapDisperse_Integer

  subroutine ovkOverlapDisperse_LargeInteger(Overlap, DisperseOp, &
    CollectedData, OverlappedGridData)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: DisperseOp
    type(ovk_array_large_int), intent(in) :: CollectedData
    type(ovk_field_large_int), intent(inout) :: OverlappedGridData

    select case (DisperseOp)
    case (OVK_DISPERSE_OVERWRITE)
      call DisperseOverwrite_LargeInteger(Overlap, CollectedData, OverlappedGridData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid disperse operation."
      stop 1
    end select

  end subroutine ovkOverlapDisperse_LargeInteger

  subroutine ovkOverlapDisperse_Real(Overlap, DisperseOp, CollectedData, &
    OverlappedGridData)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: DisperseOp
    type(ovk_array_real), intent(in) :: CollectedData
    type(ovk_field_real), intent(inout) :: OverlappedGridData

    select case (DisperseOp)
    case (OVK_DISPERSE_OVERWRITE)
      call DisperseOverwrite_Real(Overlap, CollectedData, OverlappedGridData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid disperse operation."
      stop 1
    end select

  end subroutine ovkOverlapDisperse_Real

  subroutine ovkOverlapDisperse_Logical(Overlap, DisperseOp, CollectedData, &
    OverlappedGridData)

    type(ovk_overlap), intent(in) :: Overlap
    integer, intent(in) :: DisperseOp
    type(ovk_array_logical), intent(in) :: CollectedData
    type(ovk_field_logical), intent(inout) :: OverlappedGridData

    select case (DisperseOp)
    case (OVK_DISPERSE_OVERWRITE)
      call DisperseOverwrite_Logical(Overlap, CollectedData, OverlappedGridData)
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Invalid disperse operation."
      stop 1
    end select

  end subroutine ovkOverlapDisperse_Logical

  subroutine DisperseOverwrite_Integer(Overlap, CollectedData, OverlappedGridData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_int), intent(in) :: CollectedData
    type(ovk_field_int), intent(inout) :: OverlappedGridData

    integer :: i, j, k
    integer(lk) :: l
    type(ovk_grid), pointer :: OverlappedGrid

    OverlappedGrid => Overlap%overlapped_grid

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

  subroutine DisperseOverwrite_LargeInteger(Overlap, CollectedData, OverlappedGridData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_large_int), intent(in) :: CollectedData
    type(ovk_field_large_int), intent(inout) :: OverlappedGridData

    integer :: i, j, k
    integer(lk) :: l
    type(ovk_grid), pointer :: OverlappedGrid

    OverlappedGrid => Overlap%overlapped_grid

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

  subroutine DisperseOverwrite_Real(Overlap, CollectedData, OverlappedGridData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_real), intent(in) :: CollectedData
    type(ovk_field_real), intent(inout) :: OverlappedGridData

    integer :: i, j, k
    integer(lk) :: l
    type(ovk_grid), pointer :: OverlappedGrid

    OverlappedGrid => Overlap%overlapped_grid

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

  subroutine DisperseOverwrite_Logical(Overlap, CollectedData, OverlappedGridData)

    type(ovk_overlap), intent(in) :: Overlap
    type(ovk_array_logical), intent(in) :: CollectedData
    type(ovk_field_logical), intent(inout) :: OverlappedGridData

    integer :: i, j, k
    integer(lk) :: l
    type(ovk_grid), pointer :: OverlappedGrid

    OverlappedGrid => Overlap%overlapped_grid

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

end module ovkOverlap
