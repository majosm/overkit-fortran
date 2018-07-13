! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkGrid

  use ovkBoundingBox
  use ovkCart
  use ovkField
  use ovkFieldOps
  use ovkGeometry
  use ovkGlobal
  use ovkLogger
  implicit none

  private

  ! API
  public :: ovk_grid
  public :: ovkGridExists
  public :: ovkGetGridID
  public :: ovkGetGridDimension
  public :: ovkGetGridSize
  public :: ovkGetGridPeriodicity
  public :: ovkGetGridPeriodicStorage
  public :: ovkGetGridPeriodicLength
  public :: ovkGetGridGeometryType
  public :: ovkGetGridCart
  public :: ovkGetGridCoords
  public :: ovkEditGridCoords
  public :: ovkReleaseGridCoords
  public :: ovkGetGridState
  public :: ovkEditGridState
  public :: ovkReleaseGridState
  public :: ovkResetGridState
  public :: ovkFilterGridState
  public :: ovkGridCellExists
  public :: ovkGridCellBounds
  public :: ovkOverlapsGridCell
  public :: ovkCoordsInGridCell
  public :: ovkCoordsInCubicGridCell
  public :: ovkGridResolution
  public :: ovkGenerateBBOverlapMask
  public :: ovkPeriodicExtend
  public :: ovkExportGridCoords
  public :: OVK_GRID_GEOMETRY_CARTESIAN
  public :: OVK_GRID_GEOMETRY_RECTILINEAR
  public :: OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN
  public :: OVK_GRID_GEOMETRY_ORIENTED_RECTILINEAR
  public :: OVK_GRID_GEOMETRY_CURVILINEAR
  public :: OVK_STATE_GRID
  public :: OVK_STATE_INTERIOR
  public :: OVK_STATE_BOUNDARY
  public :: OVK_STATE_DOMAIN_BOUNDARY
  public :: OVK_STATE_INFERRED_DOMAIN_BOUNDARY
  public :: OVK_STATE_INTERNAL_BOUNDARY
  public :: OVK_STATE_HOLE
  public :: OVK_STATE_BOUNDARY_HOLE
  public :: OVK_STATE_FRINGE
  public :: OVK_STATE_OUTER_FRINGE
  public :: OVK_STATE_INNER_FRINGE
  public :: OVK_STATE_OCCLUDED
  public :: OVK_STATE_OVERLAP_MINIMIZED
  public :: OVK_STATE_RECEIVER
  public :: OVK_STATE_ORPHAN
  public :: OVK_STATE_DEBUG1
  public :: OVK_STATE_DEBUG2
  public :: OVK_STATE_DEBUG3
  public :: OVK_STATE_DEBUG4
  public :: OVK_STATE_DEBUG5
  public :: OVK_INTERIOR_POINT
  public :: OVK_DOMAIN_BOUNDARY_POINT
  public :: OVK_INTERNAL_BOUNDARY_POINT
  public :: OVK_HOLE_POINT

  ! Internal
  public :: ovk_grid_
  public :: t_grid_edits
  public :: CreateGrid
  public :: DestroyGrid
  public :: GetGridEdits
  public :: ResetGridEdits

  type t_grid_edits
    logical :: coords
    logical :: mask
    logical :: boundary_mask
    logical :: internal_boundary_mask
  end type t_grid_edits

  type ovk_grid
    type(t_noconstruct) :: noconstruct
    type(t_existence_flag) :: existence_flag
    type(t_logger), pointer :: logger
    integer :: id
    integer :: nd
    integer, dimension(MAX_ND) :: npoints
    logical, dimension(MAX_ND) :: periodic
    integer :: periodic_storage
    real(rk), dimension(MAX_ND) :: periodic_length
    integer :: geometry_type
    type(ovk_cart) :: cart
    type(ovk_cart) :: cell_cart
    type(ovk_field_real), dimension(:), pointer :: coords
    integer :: coords_edit_ref_count
    type(ovk_field_int), pointer :: state
    type(ovk_field_int), pointer :: prev_state
    integer :: state_edit_ref_count
    type(t_grid_edits), pointer :: edits
    type(ovk_bbox) :: bounds
    type(ovk_field_logical) :: mask
    type(ovk_field_logical) :: boundary_mask
    type(ovk_field_logical) :: internal_boundary_mask
    type(ovk_field_logical) :: cell_mask
    type(ovk_field_real) :: cell_volumes
    type(ovk_field_real) :: resolution
    type(ovk_field_int) :: cell_edge_dist
  end type ovk_grid

  integer, parameter :: OVK_GRID_GEOMETRY_CARTESIAN = 1
  integer, parameter :: OVK_GRID_GEOMETRY_RECTILINEAR = 2
  integer, parameter :: OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN = 3
  integer, parameter :: OVK_GRID_GEOMETRY_ORIENTED_RECTILINEAR = 4
  integer, parameter :: OVK_GRID_GEOMETRY_CURVILINEAR = 5

  integer, parameter :: OVK_STATE_GRID = ishft(1,0)
  integer, parameter :: OVK_STATE_INTERIOR = ishft(1,1)
  integer, parameter :: OVK_STATE_BOUNDARY = ishft(1,2)
  integer, parameter :: OVK_STATE_DOMAIN_BOUNDARY = ishft(1,3)
  integer, parameter :: OVK_STATE_INFERRED_DOMAIN_BOUNDARY = ishft(1,4)
  integer, parameter :: OVK_STATE_INTERNAL_BOUNDARY = ishft(1,5)
  integer, parameter :: OVK_STATE_HOLE = ishft(1,6)
  integer, parameter :: OVK_STATE_BOUNDARY_HOLE = ishft(1,7)
  integer, parameter :: OVK_STATE_FRINGE = ishft(1,8)
  integer, parameter :: OVK_STATE_OUTER_FRINGE = ishft(1,9)
  integer, parameter :: OVK_STATE_INNER_FRINGE = ishft(1,10)
  integer, parameter :: OVK_STATE_OCCLUDED = ishft(1,11)
  integer, parameter :: OVK_STATE_OVERLAP_MINIMIZED = ishft(1,12)
  integer, parameter :: OVK_STATE_RECEIVER = ishft(1,13)
  integer, parameter :: OVK_STATE_ORPHAN = ishft(1,14)
  integer, parameter :: OVK_STATE_DEBUG1 = ishft(1,15)
  integer, parameter :: OVK_STATE_DEBUG2 = ishft(1,16)
  integer, parameter :: OVK_STATE_DEBUG3 = ishft(1,17)
  integer, parameter :: OVK_STATE_DEBUG4 = ishft(1,18)
  integer, parameter :: OVK_STATE_DEBUG5 = ishft(1,19)

  integer, parameter :: OVK_INTERIOR_POINT = ior(OVK_STATE_GRID,OVK_STATE_INTERIOR)
  integer, parameter :: OVK_DOMAIN_BOUNDARY_POINT = ior(OVK_STATE_GRID,ior(OVK_STATE_BOUNDARY, &
    OVK_STATE_DOMAIN_BOUNDARY))
  integer, parameter :: OVK_INTERNAL_BOUNDARY_POINT = ior(OVK_STATE_GRID,ior(OVK_STATE_BOUNDARY, &
    OVK_STATE_INTERNAL_BOUNDARY))
  integer, parameter :: OVK_HOLE_POINT = OVK_STATE_HOLE

contains

  function ovk_grid_() result(Grid)

    type(ovk_grid) :: Grid

    nullify(Grid%logger)
    Grid%id = 0
    Grid%nd = 2
    Grid%npoints = [0,0,1]
    Grid%periodic = .false.
    Grid%periodic_storage = OVK_NO_OVERLAP_PERIODIC
    Grid%periodic_length = 0._rk
    Grid%geometry_type = OVK_GRID_GEOMETRY_CURVILINEAR
    Grid%cart = ovk_cart_()
    Grid%cell_cart = ovk_cart_()
    nullify(Grid%coords)
    Grid%coords_edit_ref_count = 0
    nullify(Grid%state)
    Grid%state_edit_ref_count = 0
    nullify(Grid%prev_state)
    nullify(Grid%edits)
    Grid%bounds = ovk_bbox_()
    Grid%mask = ovk_field_logical_()
    Grid%boundary_mask = ovk_field_logical_()
    Grid%internal_boundary_mask = ovk_field_logical_()
    Grid%cell_mask = ovk_field_logical_()
    Grid%cell_volumes = ovk_field_real_()
    Grid%resolution = ovk_field_real_()
    Grid%cell_edge_dist = ovk_field_int_()

    call SetExists(Grid%existence_flag, .false.)

  end function ovk_grid_

  subroutine CreateGrid(Grid, ID, Logger, Cart, PeriodicLength, GeometryType)

    type(ovk_grid), intent(out) :: Grid
    integer, intent(in) :: ID
    type(t_logger), pointer, intent(in) :: Logger
    type(ovk_cart), intent(in) :: Cart
    real(rk), dimension(Cart%nd), intent(in) :: PeriodicLength
    integer, intent(in) :: GeometryType

    integer :: d
    type(ovk_cart) :: CellEdgeDistCart

    Grid%logger => Logger

    Grid%id = ID
    Grid%nd = Cart%nd
    Grid%npoints(:Cart%nd) = ovkCartSize(Cart)
    Grid%npoints(Cart%nd+1:) = 1
    Grid%periodic = Cart%periodic
    Grid%periodic_storage = Cart%periodic_storage
    Grid%periodic_length(:Cart%nd) = PeriodicLength
    Grid%periodic_length(Cart%nd+1:) = 0._rk
    Grid%geometry_type = GeometryType

    Grid%cart = Cart

    Grid%cell_cart = ovkCartPointToCell(Grid%cart)

    allocate(Grid%coords(Grid%nd))
    do d = 1, Grid%nd
      Grid%coords(d) = ovk_field_real_(Grid%cart, 0._rk)
    end do

    Grid%coords_edit_ref_count = 0

    allocate(Grid%state)
    Grid%state = ovk_field_int_(Grid%cart, OVK_INTERIOR_POINT)

    nullify(Grid%prev_state)

    Grid%state_edit_ref_count = 0

    allocate(Grid%edits)
    Grid%edits = t_grid_edits_(Cart%nd)

    Grid%bounds = ovk_bbox_(Grid%nd)
    if (ovkCartCount(Grid%cart) > 0) then
      Grid%bounds%b(:Grid%nd) = 0._rk
      Grid%bounds%e(:Grid%nd) = 0._rk
    end if

    Grid%mask = ovk_field_logical_(Grid%cart, .true.)
    Grid%boundary_mask = ovk_field_logical_(Grid%cart, .false.)
    Grid%internal_boundary_mask = ovk_field_logical_(Grid%cart, .false.)
    Grid%cell_mask = ovk_field_logical_(Grid%cell_cart, .true.)
    Grid%cell_volumes = ovk_field_real_(Grid%cell_cart, 0._rk)
    Grid%resolution = ovk_field_real_(Grid%cart, 0._rk)

    CellEdgeDistCart = Grid%cell_cart
    CellEdgeDistCart%is(:Grid%nd) = CellEdgeDistCart%is(:Grid%nd) - merge(0, 1, &
      Grid%cart%periodic(:Grid%nd))
    CellEdgeDistCart%ie(:Grid%nd) = CellEdgeDistCart%ie(:Grid%nd) + merge(0, 1, &
      Grid%cart%periodic(:Grid%nd))

    Grid%cell_edge_dist = ovk_field_int_(CellEdgeDistCart)

    call UpdateCellEdgeDistance(Grid)

    call SetExists(Grid%existence_flag, .true.)

  end subroutine CreateGrid

  subroutine DestroyGrid(Grid)

    type(ovk_grid), intent(inout) :: Grid

    if (.not. ovkGridExists(Grid)) return

    call SetExists(Grid%existence_flag, .false.)

    Grid%mask = ovk_field_logical_()
    Grid%boundary_mask = ovk_field_logical_()
    Grid%internal_boundary_mask = ovk_field_logical_()
    Grid%cell_mask = ovk_field_logical_()
    Grid%cell_volumes = ovk_field_real_()
    Grid%resolution = ovk_field_real_()
    Grid%cell_edge_dist = ovk_field_int_()

    deallocate(Grid%edits)

    deallocate(Grid%coords)
    deallocate(Grid%state)
    if (associated(Grid%prev_state)) deallocate(Grid%prev_state)

  end subroutine DestroyGrid

  function ovkGridExists(Grid) result(Exists)

    type(ovk_grid), intent(in) :: Grid
    logical :: Exists

    Exists = CheckExists(Grid%existence_flag)

  end function ovkGridExists

  subroutine ovkGetGridID(Grid, ID)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(out) :: ID

    ID = Grid%id

  end subroutine ovkGetGridID

  subroutine ovkGetGridDimension(Grid, NumDims)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(out) :: NumDims

    NumDims = Grid%nd

  end subroutine ovkGetGridDimension

  subroutine ovkGetGridSize(Grid, NumPoints)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(out) :: NumPoints

    NumPoints = Grid%npoints(:Grid%nd)

  end subroutine ovkGetGridSize

  subroutine ovkGetGridPeriodicity(Grid, Periodic)

    type(ovk_grid), intent(in) :: Grid
    logical, dimension(Grid%nd), intent(out) :: Periodic

    Periodic = Grid%periodic(:Grid%nd)

  end subroutine ovkGetGridPeriodicity

  subroutine ovkGetGridPeriodicStorage(Grid, PeriodicStorage)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(out) :: PeriodicStorage

    PeriodicStorage = Grid%periodic_storage

  end subroutine ovkGetGridPeriodicStorage

  subroutine ovkGetGridPeriodicLength(Grid, PeriodicLength)

    type(ovk_grid), intent(in) :: Grid
    real(rk), dimension(Grid%nd), intent(out) :: PeriodicLength

    PeriodicLength = Grid%periodic_length(:Grid%nd)

  end subroutine ovkGetGridPeriodicLength

  subroutine ovkGetGridGeometryType(Grid, GeometryType)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(out) :: GeometryType

    GeometryType = Grid%geometry_type

  end subroutine ovkGetGridGeometryType

  subroutine ovkGetGridCart(Grid, Cart)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_cart), intent(out) :: Cart

    Cart = Grid%cart

  end subroutine ovkGetGridCart

  subroutine ovkGetGridCoords(Grid, DimIndex, Coords)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(in) :: DimIndex
    type(ovk_field_real), pointer, intent(out) :: Coords

    Coords => Grid%coords(DimIndex)

  end subroutine ovkGetGridCoords

  subroutine ovkEditGridCoords(Grid, DimIndex, Coords)

    type(ovk_grid), intent(inout) :: Grid
    integer, intent(in) :: DimIndex
    type(ovk_field_real), pointer, intent(out) :: Coords

    logical :: Success, StartEdit

    if (DimIndex >= 1 .and. DimIndex <= Grid%nd) then

      call TryEditCoords(Grid, Success, StartEdit)

      if (Success) then
        Coords => Grid%coords(DimIndex)
      else
        ! This shouldn't happen, since there are currently no mutually exclusive edit
        ! operations; leaving this here just in case that changes later on
        if (OVK_DEBUG) then
          stop 1
        end if
        nullify(Coords)
      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to edit coordinates; invalid dimension."
        stop 1
      end if

      nullify(Coords)

    end if

  end subroutine ovkEditGridCoords

  subroutine ovkReleaseGridCoords(Grid, Coords)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_real), pointer, intent(inout) :: Coords

    integer :: d
    integer :: DimIndex
    logical :: Success, EndEdit

    DimIndex = 0
    do d = 1, Grid%nd
      if (associated(Coords, Grid%coords(d))) then
        DimIndex = d
        exit
      end if
    end do

    if (DimIndex /= 0) then

      call TryReleaseCoords(Grid, Success, EndEdit)

      if (Success) then
        if (EndEdit) then
          Grid%edits%coords = .true.
          if (.not. EditingState(Grid)) then
            call UpdateBounds(Grid)
            call UpdateResolution(Grid)
          end if
        end if
      else
        if (OVK_DEBUG) then
          write (ERROR_UNIT, '(2a)') "ERROR: Unable to release coordinates; not currently ", &
            "being edited."
          stop 1
        end if
      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release coordinates; invalid pointer."
        stop 1
      end if

    end if

    nullify(Coords)

  end subroutine ovkReleaseGridCoords

  subroutine ovkGetGridState(Grid, State)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_field_int), pointer, intent(out) :: State

    State => Grid%state

  end subroutine ovkGetGridState

  subroutine ovkEditGridState(Grid, State)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_int), pointer, intent(out) :: State

    logical :: Success, StartEdit

    call TryEditState(Grid, Success, StartEdit)

    if (Success) then
      allocate(Grid%prev_state)
      Grid%prev_state = Grid%state
      State => Grid%state
    else
      ! This shouldn't happen, since there are currently no mutually exclusive edit
      ! operations; leaving this here just in case that changes later on
      if (OVK_DEBUG) then
        stop 1
      end if
      nullify(State)
    end if

  end subroutine ovkEditGridState

  subroutine ovkReleaseGridState(Grid, State)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_int), pointer, intent(inout) :: State

    integer :: i, j, k
    logical :: Success, EndEdit
    type(ovk_field_int), pointer :: PrevState
    logical :: ModifiedMask
    logical :: ModifiedBoundaryMask
    logical :: ModifiedInternalBoundaryMask
    integer :: StateValue, PrevStateValue

    if (associated(State, Grid%state)) then

      call TryReleaseState(Grid, Success, EndEdit)

      if (Success) then

        if (EndEdit) then

          PrevState => Grid%prev_state
          nullify(Grid%prev_state)

          ModifiedMask = .false.
          L1: &
          do k = Grid%cart%is(3), Grid%cart%ie(3)
            do j = Grid%cart%is(2), Grid%cart%ie(2)
              do i = Grid%cart%is(1), Grid%cart%ie(1)
                if (State%values(i,j,k) /= PrevState%values(i,j,k)) then
                  StateValue = iand(State%values(i,j,k),OVK_STATE_GRID)
                  PrevStateValue = iand(PrevState%values(i,j,k),OVK_STATE_GRID)
                  if (StateValue /= PrevStateValue) then
                    ModifiedMask = .true.
                    exit L1
                  end if
                end if
              end do
            end do
          end do L1

          ModifiedBoundaryMask = .false.
          L2: &
          do k = Grid%cart%is(3), Grid%cart%ie(3)
            do j = Grid%cart%is(2), Grid%cart%ie(2)
              do i = Grid%cart%is(1), Grid%cart%ie(1)
                if (State%values(i,j,k) /= PrevState%values(i,j,k)) then
                  StateValue = iand(State%values(i,j,k),OVK_STATE_DOMAIN_BOUNDARY)
                  PrevStateValue = iand(PrevState%values(i,j,k),OVK_STATE_DOMAIN_BOUNDARY)
                  if (StateValue /= PrevStateValue) then
                    ModifiedBoundaryMask = .true.
                    exit L2
                  end if
                end if
              end do
            end do
          end do L2

          ModifiedInternalBoundaryMask = .false.
          L3: &
          do k = Grid%cart%is(3), Grid%cart%ie(3)
            do j = Grid%cart%is(2), Grid%cart%ie(2)
              do i = Grid%cart%is(1), Grid%cart%ie(1)
                if (State%values(i,j,k) /= PrevState%values(i,j,k)) then
                  StateValue = iand(State%values(i,j,k),OVK_STATE_INTERNAL_BOUNDARY)
                  PrevStateValue = iand(PrevState%values(i,j,k),OVK_STATE_INTERNAL_BOUNDARY)
                  if (StateValue /= PrevStateValue) then
                    ModifiedInternalBoundaryMask = .true.
                    exit L3
                  end if
                end if
              end do
            end do
          end do L3

          if (ModifiedMask) then
            call UpdateMask(Grid)
            call UpdateCellMask(Grid)
            call UpdateCellEdgeDistance(Grid)
            if (.not. EditingCoords(Grid)) then
              call UpdateBounds(Grid)
              call UpdateResolution(Grid)
            end if
            Grid%edits%mask = .true.
          end if

          if (ModifiedBoundaryMask) then
            call UpdateBoundaryMask(Grid)
            Grid%edits%boundary_mask = .true.
          end if

          if (ModifiedInternalBoundaryMask) then
            call UpdateInternalBoundaryMask(Grid)
            Grid%edits%internal_boundary_mask = .true.
          end if

          deallocate(PrevState)

        end if

      else

        if (OVK_DEBUG) then
          write (ERROR_UNIT, '(2a)') "ERROR: Unable to release state; not currently ", &
            "being edited."
          stop 1
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release state; invalid pointer."
        stop 1
      end if

    end if

    nullify(State)

  end subroutine ovkReleaseGridState

  subroutine ovkResetGridState(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: i, j, k
    type(ovk_field_int), pointer :: State
    integer :: AssemblyStates

    call ovkEditGridState(Grid, State)

    AssemblyStates = 0
    AssemblyStates = ior(AssemblyStates, OVK_STATE_INFERRED_DOMAIN_BOUNDARY)
    AssemblyStates = ior(AssemblyStates, OVK_STATE_BOUNDARY_HOLE)
    AssemblyStates = ior(AssemblyStates, OVK_STATE_FRINGE)
    AssemblyStates = ior(AssemblyStates, OVK_STATE_OUTER_FRINGE)
    AssemblyStates = ior(AssemblyStates, OVK_STATE_INNER_FRINGE)
    AssemblyStates = ior(AssemblyStates, OVK_STATE_OCCLUDED)
    AssemblyStates = ior(AssemblyStates, OVK_STATE_OVERLAP_MINIMIZED)
    AssemblyStates = ior(AssemblyStates, OVK_STATE_RECEIVER)
    AssemblyStates = ior(AssemblyStates, OVK_STATE_ORPHAN)

    do k = Grid%cart%is(3), Grid%cart%ie(3)
      do j = Grid%cart%is(2), Grid%cart%ie(2)
        do i = Grid%cart%is(1), Grid%cart%ie(1)
          ! Reset inferred boundaries
          if (iand(State%values(i,j,k),OVK_STATE_INFERRED_DOMAIN_BOUNDARY) /= 0) then
            State%values(i,j,k) = iand(State%values(i,j,k), not(OVK_STATE_DOMAIN_BOUNDARY))
          end if
          ! Reset holes
          if (iand(State%values(i,j,k),ior(OVK_STATE_BOUNDARY_HOLE, OVK_STATE_OVERLAP_MINIMIZED)) &
            /= 0) then
            State%values(i,j,k) = iand(State%values(i,j,k), not(OVK_STATE_HOLE))
            State%values(i,j,k) = ior(State%values(i,j,k), OVK_STATE_GRID)
          end if
          ! Remove other assembly states
          State%values(i,j,k) = iand(State%values(i,j,k),not(AssemblyStates))
        end do
      end do
    end do

    call ovkReleaseGridState(Grid, State)

  end subroutine ovkResetGridState

  subroutine ovkFilterGridState(Grid, StateBits, FilterMode, MatchingMask)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(in) :: StateBits
    integer, intent(in) :: FilterMode
    type(ovk_field_logical), intent(out) :: MatchingMask

    integer :: i, j, k
    type(ovk_field_int), pointer :: State

    State => Grid%state

    MatchingMask = ovk_field_logical_(State%cart)

    select case (FilterMode)
    case (OVK_NONE)
      do k = State%cart%is(3), State%cart%ie(3)
        do j = State%cart%is(2), State%cart%ie(2)
          do i = State%cart%is(1), State%cart%ie(1)
            MatchingMask%values(i,j,k) = iand(State%values(i,j,k), StateBits) == 0
          end do
        end do
      end do
    case (OVK_ANY)
      do k = State%cart%is(3), State%cart%ie(3)
        do j = State%cart%is(2), State%cart%ie(2)
          do i = State%cart%is(1), State%cart%ie(1)
            MatchingMask%values(i,j,k) = iand(State%values(i,j,k), StateBits) /= 0
          end do
        end do
      end do
    case (OVK_ALL)
      do k = State%cart%is(3), State%cart%ie(3)
        do j = State%cart%is(2), State%cart%ie(2)
          do i = State%cart%is(1), State%cart%ie(1)
            MatchingMask%values(i,j,k) = iand(State%values(i,j,k), StateBits) == StateBits
          end do
        end do
      end do
    end select

  end subroutine ovkFilterGridState

  function EditingCoords(Grid) result(Editing)

    type(ovk_grid), intent(in) :: Grid
    logical :: Editing

    Editing = Grid%coords_edit_ref_count > 0

  end function EditingCoords

  function EditingState(Grid) result(Editing)

    type(ovk_grid), intent(in) :: Grid
    logical :: Editing

    Editing = Grid%state_edit_ref_count > 0

  end function EditingState

  subroutine TryEditCoords(Grid, Success, StartEdit)

    type(ovk_grid), intent(inout) :: Grid
    logical, intent(out) :: Success
    logical, intent(out) :: StartEdit

    ! Currently no mutually exclusive edit operations; leaving this here just in case that changes
    ! later on
    Success = .true.

    if (Success) then
      StartEdit = Grid%coords_edit_ref_count == 0
      Grid%coords_edit_ref_count = Grid%coords_edit_ref_count + 1
    else
      StartEdit = .false.
    end if

  end subroutine TryEditCoords

  subroutine TryReleaseCoords(Grid, Success, EndEdit)

    type(ovk_grid), intent(inout) :: Grid
    logical, intent(out) :: Success
    logical, intent(out) :: EndEdit

    Success = EditingCoords(Grid)

    if (Success) then
      Grid%coords_edit_ref_count = Grid%coords_edit_ref_count - 1
      EndEdit = Grid%coords_edit_ref_count == 0
    else
      EndEdit = .false.
    end if

  end subroutine TryReleaseCoords

  subroutine TryEditState(Grid, Success, StartEdit)

    type(ovk_grid), intent(inout) :: Grid
    logical, intent(out) :: Success
    logical, intent(out) :: StartEdit

    ! Currently no mutually exclusive edit operations; leaving this here just in case that changes
    ! later on
    Success = .true.

    if (Success) then
      StartEdit = Grid%state_edit_ref_count == 0
      Grid%state_edit_ref_count = Grid%state_edit_ref_count + 1
    else
      StartEdit = .false.
    end if

  end subroutine TryEditState

  subroutine TryReleaseState(Grid, Success, EndEdit)

    type(ovk_grid), intent(inout) :: Grid
    logical, intent(out) :: Success
    logical, intent(out) :: EndEdit

    Success = EditingState(Grid)

    if (Success) then
      Grid%state_edit_ref_count = Grid%state_edit_ref_count - 1
      EndEdit = Grid%state_edit_ref_count == 0
    else
      EndEdit = .false.
    end if

  end subroutine TryReleaseState

  subroutine UpdateBounds(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: i, j, d, l
    integer :: idir, jdir
    integer, dimension(MAX_ND) :: Point, AdjustedPoint
    real(rk), dimension(Grid%nd) :: PeriodicCoords

    Grid%bounds = ovk_bbox_(Grid%nd)
    do d = 1, Grid%nd
      Grid%bounds%b(d) = minval(Grid%coords(d)%values)
      Grid%bounds%e(d) = maxval(Grid%coords(d)%values)
    end do

    ! Add contribution from periodic points
    if (Grid%cart%periodic_storage == OVK_NO_OVERLAP_PERIODIC .and. &
      any(Grid%cart%periodic .and. Grid%periodic_length > 0._rk)) then
      do d = 1, Grid%nd
        if (Grid%cart%periodic(d)) then
          idir = modulo((d+1)-1,MAX_ND) + 1
          jdir = modulo((d+2)-1,MAX_ND) + 1
          do j = Grid%cart%is(jdir), Grid%cart%ie(jdir)
            do i = Grid%cart%is(idir), Grid%cart%ie(idir)
              Point(d) = Grid%cart%ie(d)+1
              Point(idir) = i
              Point(jdir) = j
              AdjustedPoint(d) = Grid%cart%is(d)
              AdjustedPoint(idir) = i
              AdjustedPoint(jdir) = j
              do l = 1, Grid%nd
                PeriodicCoords(l) = Grid%coords(l)%values(AdjustedPoint(1),AdjustedPoint(2), &
                  AdjustedPoint(3))
              end do
              PeriodicCoords = ovkPeriodicExtend(Grid%cart, Grid%periodic_length, &
                Point, PeriodicCoords)
              Grid%bounds = ovkBBExtend(Grid%bounds, PeriodicCoords)
            end do
          end do
        end if
      end do
    end if

  end subroutine UpdateBounds

  subroutine UpdateMask(Grid)

    type(ovk_grid), intent(inout) :: Grid

    call ovkFilterGridState(Grid, OVK_STATE_GRID, OVK_ALL, Grid%mask)

  end subroutine UpdateMask

  subroutine UpdateBoundaryMask(Grid)

    type(ovk_grid), intent(inout) :: Grid

    call ovkFilterGridState(Grid, OVK_STATE_DOMAIN_BOUNDARY, OVK_ALL, Grid%boundary_mask)

  end subroutine UpdateBoundaryMask

  subroutine UpdateInternalBoundaryMask(Grid)

    type(ovk_grid), intent(inout) :: Grid

    call ovkFilterGridState(Grid, OVK_STATE_INTERNAL_BOUNDARY, OVK_ALL, &
      Grid%internal_boundary_mask)

  end subroutine UpdateInternalBoundaryMask

  subroutine UpdateCellMask(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: i, j, k, m, n, o
    integer, dimension(MAX_ND) :: VertexLower, VertexUpper
    integer, dimension(MAX_ND) :: Vertex
    logical :: AwayFromEdge

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  SHARED(Grid)
    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          VertexLower = [i,j,k]
          VertexUpper(:Grid%nd) = VertexLower(:Grid%nd) + 1
          VertexUpper(Grid%nd+1:) = 1
          AwayFromEdge = ovkCartContains(Grid%cart, VertexUpper)
          if (AwayFromEdge) then
            Grid%cell_mask%values(i,j,k) = .true.
            L1: &
            do o = VertexLower(3), VertexUpper(3)
              do n = VertexLower(2), VertexUpper(2)
                do m = VertexLower(1), VertexUpper(1)
                  if (.not. Grid%mask%values(m,n,o)) then
                    Grid%cell_mask%values(i,j,k) = .false.
                    exit L1
                  end if
                end do
              end do
            end do L1
          else
            Grid%cell_mask%values(i,j,k) = .true.
            L2: &
            do o = VertexLower(3), VertexUpper(3)
              do n = VertexLower(2), VertexUpper(2)
                do m = VertexLower(1), VertexUpper(1)
                  Vertex = [m,n,o]
                  Vertex(:Grid%nd) = ovkCartPeriodicAdjust(Grid%cart, Vertex)
                  if (.not. Grid%mask%values(Vertex(1),Vertex(2),Vertex(3))) then
                    Grid%cell_mask%values(i,j,k) = .false.
                    exit L2
                  end if
                end do
              end do
            end do L2
          end if
        end do
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine UpdateCellMask

  subroutine UpdateResolution(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: i, j, k, m, n, o
    integer, dimension(MAX_ND) :: Cell
    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: NeighborCellLower, NeighborCellUpper
    integer, dimension(MAX_ND) :: Neighbor
    real(rk) :: AvgCellVolume
    integer :: NumCells
    logical :: AwayFromEdge

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  SHARED(Grid)
    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          Cell = [i,j,k]
          if (Grid%cell_mask%values(i,j,k)) then
            call GetCellVertexCoords(Grid, Cell, VertexCoords)
            select case (Grid%geometry_type)
            case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
              select case (Grid%nd)
              case (2)
                Grid%cell_volumes%values(i,j,k) = ovkRectangleSize(VertexCoords)
              case (3)
                Grid%cell_volumes%values(i,j,k) = ovkCuboidSize(VertexCoords)
              end select
            case (OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN,OVK_GRID_GEOMETRY_ORIENTED_RECTILINEAR)
              select case (Grid%nd)
              case (2)
                Grid%cell_volumes%values(i,j,k) = ovkOrientedRectangleSize(VertexCoords)
              case (3)
                Grid%cell_volumes%values(i,j,k) = ovkOrientedCuboidSize(VertexCoords)
              end select
            case default
              select case (Grid%nd)
              case (2)
                Grid%cell_volumes%values(i,j,k) = ovkQuadSize(VertexCoords)
              case (3)
                Grid%cell_volumes%values(i,j,k) = ovkHexahedronSize(VertexCoords)
              end select
            end select
          else
            Grid%cell_volumes%values(i,j,k) = 0._rk
          end if
        end do
      end do
    end do
!$OMP END PARALLEL DO

    ! Compute the grid resolution at each point by averaging the sizes of neighboring cells
!$OMP PARALLEL DO &
!$OMP&  DEFAULT(PRIVATE) &
!$OMP&  SHARED(Grid)
    do k = Grid%cart%is(3), Grid%cart%ie(3)
      do j = Grid%cart%is(2), Grid%cart%ie(2)
        do i = Grid%cart%is(1), Grid%cart%ie(1)
          Point = [i,j,k]
          NeighborCellLower(:Grid%nd) = Point(:Grid%nd)-1
          NeighborCellLower(Grid%nd+1:) = 1
          NeighborCellUpper(:Grid%nd) = Point(:Grid%nd)
          NeighborCellUpper(Grid%nd+1:) = 1
          AvgCellVolume = 0._rk
          NumCells = 0
          AwayFromEdge = ovkCartContains(Grid%cell_cart, NeighborCellLower) .and. &
            ovkCartContains(Grid%cell_cart, NeighborCellUpper)
          if (AwayFromEdge) then
            do o = NeighborCellLower(3), NeighborCellUpper(3)
              do n = NeighborCellLower(2), NeighborCellUpper(2)
                do m = NeighborCellLower(1), NeighborCellUpper(1)
                  if (Grid%cell_mask%values(m,n,o)) then
                    AvgCellVolume = AvgCellVolume + Grid%cell_volumes%values(m,n,o)
                    NumCells = NumCells + 1
                  end if
                end do
              end do
            end do
          else
            do o = NeighborCellLower(3), NeighborCellUpper(3)
              do n = NeighborCellLower(2), NeighborCellUpper(2)
                do m = NeighborCellLower(1), NeighborCellUpper(1)
                  Neighbor = [m,n,o]
                  Neighbor(:Grid%nd) = ovkCartPeriodicAdjust(Grid%cell_cart, Neighbor)
                  if (ovkCartContains(Grid%cell_cart, Neighbor)) then
                    if (Grid%cell_mask%values(Neighbor(1),Neighbor(2),Neighbor(3))) then
                      AvgCellVolume = AvgCellVolume + Grid%cell_volumes%values(Neighbor(1), &
                        Neighbor(2),Neighbor(3))
                      NumCells = NumCells + 1
                    end if
                  end if
                end do
              end do
            end do
          end if
          AvgCellVolume = AvgCellVolume/real(max(NumCells,1), kind=rk)
          if (AvgCellVolume > 0._rk) then
            Grid%resolution%values(i,j,k) = 1._rk/AvgCellVolume
          else
            Grid%resolution%values(i,j,k) = 0._rk
          end if
        end do
      end do
    end do
!$OMP END PARALLEL DO

  end subroutine UpdateResolution

  subroutine UpdateCellEdgeDistance(Grid)

    type(ovk_grid), intent(inout) :: Grid

    type(ovk_field_logical) :: NotMask

    NotMask = ovk_field_logical_(Grid%cell_edge_dist%cart, .true.)
    NotMask%values(Grid%cell_cart%is(1):Grid%cell_cart%ie(1), &
      Grid%cell_cart%is(2):Grid%cell_cart%ie(2),Grid%cell_cart%is(3):Grid%cell_cart%ie(3)) = &
      .not. Grid%cell_mask%values

    call ovkDistanceField(NotMask, OVK_TRUE, Grid%cell_edge_dist)

  end subroutine UpdateCellEdgeDistance

  subroutine GetGridEdits(Grid, Edits)

    type(ovk_grid), intent(in) :: Grid
    type(t_grid_edits), pointer, intent(out) :: Edits

    Edits => Grid%edits

  end subroutine GetGridEdits

  subroutine ResetGridEdits(Grid)

    type(ovk_grid), intent(inout) :: Grid

    Grid%edits%coords = .false.
    Grid%edits%mask = .false.
    Grid%edits%boundary_mask = .false.
    Grid%edits%internal_boundary_mask = .false.

  end subroutine ResetGridEdits

  subroutine GetCellVertexCoords(Grid, Cell, VertexCoords)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    real(rk), dimension(Grid%nd,2**Grid%nd), intent(out) :: VertexCoords

    integer :: i, j, k, l, m
    integer, dimension(MAX_ND) :: VertexLower, VertexUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Vertex
    integer, dimension(MAX_ND) :: AdjustedVertex
    real(rk), dimension(Grid%nd) :: PrincipalCoords

    VertexLower(:Grid%nd) = Cell
    VertexLower(Grid%nd+1:) = 1
    VertexUpper(:Grid%nd) = Cell+1
    VertexUpper(Grid%nd+1:) = 1

    AwayFromEdge = ovkCartContains(Grid%cart, VertexLower) .and. &
      ovkCartContains(Grid%cart, VertexUpper)

    if (AwayFromEdge) then
      l = 1
      do k = VertexLower(3), VertexUpper(3)
        do j = VertexLower(2), VertexUpper(2)
          do i = VertexLower(1), VertexUpper(1)
            do m = 1, Grid%nd
              VertexCoords(m,l) = Grid%coords(m)%values(i,j,k)
            end do
            l = l + 1
          end do
        end do
      end do
    else
      AdjustedVertex = 1
      l = 1
      do k = VertexLower(3), VertexUpper(3)
        do j = VertexLower(2), VertexUpper(2)
          do i = VertexLower(1), VertexUpper(1)
            Vertex = [i,j,k]
            AdjustedVertex(:Grid%nd) = ovkCartPeriodicAdjust(Grid%cart, Vertex)
            do m = 1, Grid%nd
              PrincipalCoords(m) = Grid%coords(m)%values(AdjustedVertex(1),AdjustedVertex(2), &
                AdjustedVertex(3))
            end do
            VertexCoords(:,l) = ovkPeriodicExtend(Grid%cart, Grid%periodic_length, &
              Vertex, PrincipalCoords)
            l = l + 1
          end do
        end do
      end do
    end if

  end subroutine GetCellVertexCoords

  subroutine GetCubicCellVertexCoords(Grid, Cell, VertexCoords)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    real(rk), dimension(Grid%nd,4**Grid%nd), intent(out) :: VertexCoords

    integer :: i, j, k, l, m
    integer, dimension(MAX_ND) :: VertexLower, VertexUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Vertex
    integer, dimension(MAX_ND) :: AdjustedVertex
    real(rk), dimension(Grid%nd) :: PrincipalCoords

    VertexLower(:Grid%nd) = Cell
    VertexLower(Grid%nd+1:) = 1
    VertexUpper(:Grid%nd) = Cell+3
    VertexUpper(Grid%nd+1:) = 1

    AwayFromEdge = ovkCartContains(Grid%cart, VertexLower) .and. &
      ovkCartContains(Grid%cart, VertexUpper)

    if (AwayFromEdge) then
      l = 1
      do k = VertexLower(3), VertexUpper(3)
        do j = VertexLower(2), VertexUpper(2)
          do i = VertexLower(1), VertexUpper(1)
            do m = 1, Grid%nd
              VertexCoords(m,l) = Grid%coords(m)%values(i,j,k)
            end do
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
            AdjustedVertex(:Grid%nd) = ovkCartPeriodicAdjust(Grid%cart, Vertex)
            AdjustedVertex(Grid%nd+1:) = 1
            do m = 1, Grid%nd
              PrincipalCoords(m) = Grid%coords(m)%values(AdjustedVertex(1),AdjustedVertex(2), &
                AdjustedVertex(3))
            end do
            VertexCoords(:,l) = ovkPeriodicExtend(Grid%cart, Grid%periodic_length, &
              Vertex, PrincipalCoords)
            l = l + 1
          end do
        end do
      end do
    end if

  end subroutine GetCubicCellVertexCoords

  function ovkGridCellExists(Grid, Cell) result(Exists)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    logical :: Exists

    integer, dimension(MAX_ND) :: PaddedCell

    PaddedCell(:Grid%nd) = Cell
    PaddedCell(Grid%nd+1:) = 1

    Exists = Grid%cell_mask%values(PaddedCell(1),PaddedCell(2),PaddedCell(3))

  end function ovkGridCellExists

  function ovkGridCellBounds(Grid, Cell) result(CellBounds)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    type(ovk_bbox) :: CellBounds

    integer, dimension(MAX_ND) :: PaddedCell
    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords

    PaddedCell(:Grid%nd) = Cell
    PaddedCell(Grid%nd+1:) = 1

    if (.not. Grid%cell_mask%values(PaddedCell(1),PaddedCell(2),PaddedCell(3))) then
      CellBounds = ovk_bbox_(Grid%nd)
      return
    end if

    call GetCellVertexCoords(Grid, Cell, VertexCoords)

    CellBounds = ovkBBFromPoints(VertexCoords)

  end function ovkGridCellBounds

  function ovkOverlapsGridCell(Grid, Cell, Coords, OverlapTolerance) result(Overlaps)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    real(rk), intent(in) :: OverlapTolerance
    logical :: Overlaps

    integer :: i
    integer, dimension(MAX_ND) :: PaddedCell
    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords
    real(rk), dimension(Grid%nd) :: Centroid

    PaddedCell(:Grid%nd) = Cell
    PaddedCell(Grid%nd+1:) = 1

    if (.not. Grid%cell_mask%values(PaddedCell(1),PaddedCell(2),PaddedCell(3))) then
      Overlaps = .false.
      return
    end if

    call GetCellVertexCoords(Grid, Cell, VertexCoords)

    if (OverlapTolerance > 0._rk) then
      Centroid = sum(VertexCoords,dim=2)/2._rk**Grid%nd
      do i = 1, 2**Grid%nd
        VertexCoords(:,i) = Centroid + (1._rk+OverlapTolerance) * (VertexCoords(:,i)-Centroid)
      end do
    end if

    select case (Grid%geometry_type)
    case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
      select case (Grid%nd)
      case (2)
        Overlaps = ovkOverlapsRectangle(VertexCoords, Coords)
      case (3)
        Overlaps = ovkOverlapsCuboid(VertexCoords, Coords)
      end select
    case (OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN,OVK_GRID_GEOMETRY_ORIENTED_RECTILINEAR)
      select case (Grid%nd)
      case (2)
        Overlaps = ovkOverlapsOrientedRectangle(VertexCoords, Coords)
      case (3)
        Overlaps = ovkOverlapsOrientedCuboid(VertexCoords, Coords)
      end select
    case default
      select case (Grid%nd)
      case (2)
        Overlaps = ovkOverlapsQuad(VertexCoords, Coords)
      case (3)
        Overlaps = ovkOverlapsHexahedron(VertexCoords, Coords)
      end select
    end select

  end function ovkOverlapsGridCell

  function ovkCoordsInGridCell(Grid, Cell, Coords, Guess) result(CoordsInCell)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    real(rk), dimension(Grid%nd), intent(in), optional :: Guess
    real(rk), dimension(Grid%nd) :: CoordsInCell

    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords
    logical :: Success

    call GetCellVertexCoords(Grid, Cell, VertexCoords)

    Success = .true.

    select case (Grid%geometry_type)
    case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
      select case (Grid%nd)
      case (2)
        CoordsInCell = ovkRectangleIsoInverseLinear(VertexCoords, Coords)
      case (3)
        CoordsInCell = ovkCuboidIsoInverseLinear(VertexCoords, Coords)
      end select
    case (OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN,OVK_GRID_GEOMETRY_ORIENTED_RECTILINEAR)
      select case (Grid%nd)
      case (2)
        CoordsInCell = ovkOrientedRectangleIsoInverseLinear(VertexCoords, Coords)
      case (3)
        CoordsInCell = ovkOrientedCuboidIsoInverseLinear(VertexCoords, Coords)
      end select
    case default
      select case (Grid%nd)
      case (2)
        CoordsInCell = ovkQuadIsoInverseLinear(VertexCoords, Coords, Guess=Guess, Success=Success)
      case (3)
        CoordsInCell = ovkHexahedronIsoInverseLinear(VertexCoords, Coords, Guess=Guess, &
          Success=Success)
      end select
    end select

    if (Grid%logger%verbose) then
      if (.not. Success) then
        write (ERROR_UNIT, '(6a)') "WARNING: Failed to compute isoparametric coordinates of ", &
          trim(CoordsToString(Coords)), " in cell ", trim(TupleToString(Cell)), &
          " of grid ", trim(IntToString(Grid%id))
      end if
    end if

  end function ovkCoordsInGridCell

  function ovkCoordsInCubicGridCell(Grid, Cell, Coords, Guess) result(CoordsInCell)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    real(rk), dimension(Grid%nd), intent(in), optional :: Guess
    real(rk), dimension(Grid%nd) :: CoordsInCell

    real(rk), dimension(Grid%nd,4**Grid%nd) :: VertexCoords
    logical :: Success

    call GetCubicCellVertexCoords(Grid, Cell, VertexCoords)

    Success = .true.

    select case (Grid%geometry_type)
    case (OVK_GRID_GEOMETRY_CARTESIAN)
      select case (Grid%nd)
      case (2)
        CoordsInCell = ovkRectangleIsoInverseCubic(VertexCoords, Coords)
      case (3)
        CoordsInCell = ovkCuboidIsoInverseCubic(VertexCoords, Coords)
      end select
    case (OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN)
      select case (Grid%nd)
      case (2)
        CoordsInCell = ovkOrientedRectangleIsoInverseCubic(VertexCoords, Coords)
      case (3)
        CoordsInCell = ovkOrientedCuboidIsoInverseCubic(VertexCoords, Coords)
      end select
    case default
      select case (Grid%nd)
      case (2)
        CoordsInCell = ovkQuadIsoInverseCubic(VertexCoords, Coords, Guess=Guess, &
          Success=Success)
      case (3)
        CoordsInCell = ovkHexahedronIsoInverseCubic(VertexCoords, Coords, Guess=Guess, &
          Success=Success)
      end select
    end select

    if (Grid%logger%verbose) then
      if (.not. Success) then
        write (ERROR_UNIT, '(6a)') "WARNING: Failed to compute isoparametric coordinates of ", &
          trim(CoordsToString(Coords)), " in cell ", trim(TupleToString(Cell)), &
          " of grid ", trim(IntToString(Grid%id))
      end if
    end if

  end function ovkCoordsInCubicGridCell

  function ovkGridResolution(Grid, Cell, CoordsInCell) result(Resolution)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    real(rk), dimension(Grid%nd), intent(in) :: CoordsInCell
    real(rk) :: Resolution

    integer :: i, j, k
    real(rk), dimension(MAX_ND,0:1) :: InterpBasis
    integer, dimension(MAX_ND) :: CellLower, CellUpper
    real(rk) :: InterpolatedCellVolume
    real(rk) :: CellVolume
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: BasisIndex

    InterpBasis(1,:) = ovkInterpBasisLinear(CoordsInCell(1))
    InterpBasis(2,:) = ovkInterpBasisLinear(CoordsInCell(2))
    if (Grid%nd == 3) then
      InterpBasis(3,:) = ovkInterpBasisLinear(CoordsInCell(3))
    else
      InterpBasis(3,:) = ovkInterpBasisLinear(0._rk)
    end if

    CellLower(:Grid%nd) = Cell
    CellLower(Grid%nd+1:) = 1
    CellUpper(:Grid%nd) = Cell+1
    CellUpper(Grid%nd+1:) = 1

    InterpolatedCellVolume = 0._rk

    AwayFromEdge = ovkCartContains(Grid%cart, CellUpper)

    if (AwayFromEdge) then
      do k = CellLower(3), CellUpper(3)
        do j = CellLower(2), CellUpper(2)
          do i = CellLower(1), CellUpper(1)
            Point = [i,j,k]
            BasisIndex(:Grid%nd) = Point(:Grid%nd) - CellLower(:Grid%nd)
            BasisIndex(Grid%nd+1:) = 0
            CellVolume = 1._rk/Grid%resolution%values(i,j,k)
            InterpolatedCellVolume = InterpolatedCellVolume + CellVolume * &
              InterpBasis(1,BasisIndex(1)) * InterpBasis(2,BasisIndex(2)) * &
              InterpBasis(3,BasisIndex(3))
          end do
        end do
      end do
    else
      do k = CellLower(3), CellUpper(3)
        do j = CellLower(2), CellUpper(2)
          do i = CellLower(1), CellUpper(1)
            Point = [i,j,k]
            BasisIndex(:Grid%nd) = Point(:Grid%nd) - CellLower(:Grid%nd)
            BasisIndex(Grid%nd+1:) = 0
            Point(:Grid%nd) = ovkCartPeriodicAdjust(Grid%cart, Point)
            CellVolume = 1._rk/Grid%resolution%values(Point(1),Point(2),Point(3))
            InterpolatedCellVolume = InterpolatedCellVolume + CellVolume * &
              InterpBasis(1,BasisIndex(1)) * InterpBasis(2,BasisIndex(2)) * &
              InterpBasis(3,BasisIndex(3))
          end do
        end do
      end do
    end if

    Resolution = 1._rk/InterpolatedCellVolume

  end function ovkGridResolution

  subroutine ovkGenerateBBOverlapMask(Grid, Bounds, BBOverlapMask)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_bbox), intent(in) :: Bounds
    type(ovk_field_logical), intent(out) :: BBOverlapMask

    integer :: i, j, k, l
    real(rk), dimension(Grid%nd) :: Coords

    BBOverlapMask = ovk_field_logical_(Grid%cart)

    do k = Grid%cart%is(3), Grid%cart%ie(3)
      do j = Grid%cart%is(2), Grid%cart%ie(2)
        do i = Grid%cart%is(1), Grid%cart%ie(1)
          do l = 1, Grid%nd
            Coords(l) = Grid%coords(l)%values(i,j,k)
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

    integer, dimension(Cart%nd) :: PeriodStart, PeriodEnd
    real(rk), dimension(Cart%nd) :: PositiveAdjustment, NegativeAdjustment

    PeriodStart = Cart%is(:Cart%nd)
    if (Cart%periodic_storage == OVK_NO_OVERLAP_PERIODIC) then
      PeriodEnd = Cart%ie(:Cart%nd)
    else
      PeriodEnd = Cart%ie(:Cart%nd)-1
    end if

    PositiveAdjustment = merge(PeriodicLength, 0._rk, Cart%periodic(:Cart%nd) .and. &
      Point > PeriodEnd)
    NegativeAdjustment = merge(PeriodicLength, 0._rk, Cart%periodic(:Cart%nd) .and. &
      Point < PeriodStart)

    ExtendedCoords = Coords + PositiveAdjustment - NegativeAdjustment

  end function ovkPeriodicExtend

  subroutine ovkExportGridCoords(Grid, DimIndex, ExportCart, Coords)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(in) :: DimIndex
    type(ovk_cart), intent(in) :: ExportCart
    type(ovk_field_real), intent(out) :: Coords

    integer :: i, j, k, d
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: AdjustedPoint
    real(rk), dimension(Grid%nd) :: PrincipalCoords
    real(rk), dimension(Grid%nd) :: ExtendedCoords

    if (OVK_DEBUG) then
      if (.not. ovkCartIsCompatible(ExportCart, Grid%cart)) then
        write (ERROR_UNIT, '(a)') "ERROR: Export cart is incompatible with grid."
        stop 1
      end if
    end if

    Coords = ovk_field_real_(ExportCart)

    do k = ExportCart%is(3), ExportCart%ie(3)
      do j = ExportCart%is(2), ExportCart%ie(2)
        do i = ExportCart%is(1), ExportCart%ie(1)
          Point = [i,j,k]
          if (.not. ovkCartContains(Grid%cart, Point)) then
            AdjustedPoint(:Grid%nd) = ovkCartPeriodicAdjust(Grid%cart, Point)
            AdjustedPoint(Grid%nd+1:) = 1
            do d = 1, Grid%nd
              PrincipalCoords(d) = Grid%coords(d)%values(AdjustedPoint(1),AdjustedPoint(2), &
                AdjustedPoint(3))
            end do
            ExtendedCoords = ovkPeriodicExtend(Grid%cart, Grid%periodic_length, Point, &
              PrincipalCoords)
            Coords%values(i,j,k) = ExtendedCoords(DimIndex)
          else
            Coords%values(i,j,k) = Grid%coords(DimIndex)%values(Point(1),Point(2),Point(3))
          end if
        end do
      end do
    end do

  end subroutine ovkExportGridCoords

  function t_grid_edits_(NumDims) result(Edits)

    integer :: NumDims
    type(t_grid_edits) :: Edits

    Edits%coords = .false.
    Edits%mask = .false.
    Edits%boundary_mask = .false.
    Edits%internal_boundary_mask = .false.

  end function t_grid_edits_

end module ovkGrid
