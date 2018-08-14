! Copyright (c) 2018 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkGrid

  use ovkBoundingBox
  use ovkCart
  use ovkField
  use ovkFieldOps
  use ovkGeometryOps
  use ovkGlobal
  use ovkLogger
  implicit none

  private

  ! Public API
  public :: ovk_grid
  public :: ovkGridExists
  public :: ovkGetGridID
  public :: ovkGetGridDimension
  public :: ovkGetGridSize
  public :: ovkGetGridCart
  public :: ovkGetGridPeriodicity
  public :: ovkGetGridPeriodicStorage
  public :: ovkGetGridPeriodicLength
  public :: ovkGetGridGeometryType
  public :: ovkGetGridCoords
  public :: ovkEditGridCoords
  public :: ovkReleaseGridCoords
  public :: ovkGetGridState
  public :: ovkEditGridState
  public :: ovkReleaseGridState
  public :: ovkResetGridState
  public :: ovkFilterGridState
  public :: ovkGetGridBounds
  public :: ovkGridCellBounds
  public :: ovkOverlapsGridCell
  public :: ovkCoordsInGridCell
  public :: ovkPeriodicExtend
  public :: OVK_GEOMETRY_UNIFORM
  public :: OVK_GEOMETRY_ORIENTED_UNIFORM
  public :: OVK_GEOMETRY_RECTILINEAR
  public :: OVK_GEOMETRY_ORIENTED_RECTILINEAR
  public :: OVK_GEOMETRY_CURVILINEAR
  public :: OVK_STATE_GRID
  public :: OVK_STATE_INTERIOR
  public :: OVK_STATE_BOUNDARY
  public :: OVK_STATE_EXTERIOR
  public :: OVK_STATE_DOMAIN_BOUNDARY
  public :: OVK_STATE_INTERNAL_BOUNDARY
  public :: OVK_STATE_OVERLAPPED
  public :: OVK_STATE_INFERRED_DOMAIN_BOUNDARY
  public :: OVK_STATE_BOUNDARY_HOLE
  public :: OVK_STATE_OCCLUDED
  public :: OVK_STATE_FRINGE
  public :: OVK_STATE_OUTER_FRINGE
  public :: OVK_STATE_INNER_FRINGE
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
  public :: OVK_EXTERIOR_POINT

  ! Internal API
  public :: ovk_grid_
  public :: t_grid_edits
  public :: CreateGrid
  public :: DestroyGrid
  public :: GetGridEdits
  public :: ResetGridEdits

  ! For tests only
  public :: EditingCoords
  public :: EditingState

  type t_grid_edits
    logical :: coords
    logical :: mask
    logical :: boundary
  end type t_grid_edits

  type ovk_grid
    type(t_noconstruct) :: noconstruct
    type(t_existence_flag) :: existence_flag
    type(t_logger) :: logger
    integer :: id
    integer :: nd
    integer, dimension(MAX_DIMS) :: npoints
    type(ovk_cart) :: cart
    type(ovk_cart) :: cell_cart
    real(rk), dimension(MAX_DIMS) :: periodic_length
    integer :: geometry_type
    type(ovk_field_real), dimension(:), pointer :: coords
    integer :: coords_edit_ref_count
    type(ovk_field_int), pointer :: state
    type(ovk_field_int), pointer :: prev_state
    integer :: state_edit_ref_count
    type(t_grid_edits), pointer :: edits
    type(ovk_bbox) :: bounds
    type(ovk_field_logical) :: mask
    type(ovk_field_logical) :: cell_mask
    type(ovk_field_real) :: volumes
    type(ovk_field_real) :: cell_volumes
  end type ovk_grid

  integer, parameter :: OVK_GEOMETRY_UNIFORM = 1
  integer, parameter :: OVK_GEOMETRY_ORIENTED_UNIFORM = 2
  integer, parameter :: OVK_GEOMETRY_RECTILINEAR = 3
  integer, parameter :: OVK_GEOMETRY_ORIENTED_RECTILINEAR = 4
  integer, parameter :: OVK_GEOMETRY_CURVILINEAR = 5

  ! User-modifiable states
  integer, parameter :: NUM_USER_STATES = 6
  integer, parameter :: USER_STATE_OFFSET = 0
  integer, parameter :: USER_STATES = ishft(ishft(1,NUM_USER_STATES)-1,USER_STATE_OFFSET)
  integer, parameter :: FIRST_USER_STATE = ishft(1,USER_STATE_OFFSET)

  integer, parameter :: OVK_STATE_GRID = ishft(FIRST_USER_STATE,0)
  integer, parameter :: OVK_STATE_INTERIOR = ishft(FIRST_USER_STATE,1)
  integer, parameter :: OVK_STATE_BOUNDARY = ishft(FIRST_USER_STATE,2)
  integer, parameter :: OVK_STATE_EXTERIOR = ishft(FIRST_USER_STATE,3)
  integer, parameter :: OVK_STATE_DOMAIN_BOUNDARY = ishft(FIRST_USER_STATE,4)
  integer, parameter :: OVK_STATE_INTERNAL_BOUNDARY = ishft(FIRST_USER_STATE,5)

  ! Assembly states
  integer, parameter :: NUM_ASSEMBLY_STATES = 10
  integer, parameter :: ASSEMBLY_STATE_OFFSET = NUM_USER_STATES
  integer, parameter :: ASSEMBLY_STATES = ishft(ishft(1,NUM_ASSEMBLY_STATES)-1, ASSEMBLY_STATE_OFFSET)
  integer, parameter :: FIRST_ASSEMBLY_STATE = ishft(1,ASSEMBLY_STATE_OFFSET)

  integer, parameter :: OVK_STATE_OVERLAPPED = ishft(FIRST_ASSEMBLY_STATE,0)
  integer, parameter :: OVK_STATE_INFERRED_DOMAIN_BOUNDARY = ishft(FIRST_ASSEMBLY_STATE,1)
  integer, parameter :: OVK_STATE_BOUNDARY_HOLE = ishft(FIRST_ASSEMBLY_STATE,2)
  integer, parameter :: OVK_STATE_OCCLUDED = ishft(FIRST_ASSEMBLY_STATE,3)
  integer, parameter :: OVK_STATE_FRINGE = ishft(FIRST_ASSEMBLY_STATE,4)
  integer, parameter :: OVK_STATE_OUTER_FRINGE = ishft(FIRST_ASSEMBLY_STATE,5)
  integer, parameter :: OVK_STATE_INNER_FRINGE = ishft(FIRST_ASSEMBLY_STATE,6)
  integer, parameter :: OVK_STATE_OVERLAP_MINIMIZED = ishft(FIRST_ASSEMBLY_STATE,7)
  integer, parameter :: OVK_STATE_RECEIVER = ishft(FIRST_ASSEMBLY_STATE,8)
  integer, parameter :: OVK_STATE_ORPHAN = ishft(FIRST_ASSEMBLY_STATE,9)

  ! Debug states
  integer, parameter :: NUM_DEBUG_STATES = 5
  integer, parameter :: DEBUG_STATE_OFFSET = NUM_USER_STATES + NUM_ASSEMBLY_STATES
  integer, parameter :: DEBUG_STATES = ishft(ishft(1,NUM_DEBUG_STATES)-1, DEBUG_STATE_OFFSET)
  integer, parameter :: FIRST_DEBUG_STATE = ishft(1,DEBUG_STATE_OFFSET)

  integer, parameter :: OVK_STATE_DEBUG1 = ishft(FIRST_DEBUG_STATE,0)
  integer, parameter :: OVK_STATE_DEBUG2 = ishft(FIRST_DEBUG_STATE,1)
  integer, parameter :: OVK_STATE_DEBUG3 = ishft(FIRST_DEBUG_STATE,2)
  integer, parameter :: OVK_STATE_DEBUG4 = ishft(FIRST_DEBUG_STATE,3)
  integer, parameter :: OVK_STATE_DEBUG5 = ishft(FIRST_DEBUG_STATE,4)

  ! State presets for users
  integer, parameter :: OVK_INTERIOR_POINT = ior(OVK_STATE_GRID,OVK_STATE_INTERIOR)
  integer, parameter :: OVK_DOMAIN_BOUNDARY_POINT = ior(OVK_STATE_GRID,ior(OVK_STATE_BOUNDARY, &
    OVK_STATE_DOMAIN_BOUNDARY))
  integer, parameter :: OVK_INTERNAL_BOUNDARY_POINT = ior(OVK_STATE_GRID,ior(OVK_STATE_BOUNDARY, &
    OVK_STATE_INTERNAL_BOUNDARY))
  integer, parameter :: OVK_EXTERIOR_POINT = OVK_STATE_EXTERIOR

contains

  function ovk_grid_() result(Grid)

    type(ovk_grid) :: Grid

    Grid%logger = t_logger_()
    Grid%id = 0
    Grid%nd = 2
    Grid%npoints = [0,0,1]
    Grid%cart = ovk_cart_()
    Grid%cell_cart = ovk_cart_()
    Grid%periodic_length = 0._rk
    Grid%geometry_type = OVK_GEOMETRY_CURVILINEAR
    nullify(Grid%coords)
    Grid%coords_edit_ref_count = 0
    nullify(Grid%state)
    Grid%state_edit_ref_count = 0
    nullify(Grid%prev_state)
    nullify(Grid%edits)
    Grid%bounds = ovk_bbox_()
    Grid%mask = ovk_field_logical_()
    Grid%cell_mask = ovk_field_logical_()
    Grid%volumes = ovk_field_real_()
    Grid%cell_volumes = ovk_field_real_()

    call SetExists(Grid%existence_flag, .false.)

  end function ovk_grid_

  subroutine CreateGrid(Grid, Logger, ID, Cart, PeriodicLength, GeometryType)

    type(ovk_grid), intent(out) :: Grid
    type(t_logger), intent(in) :: Logger
    integer, intent(in) :: ID
    type(ovk_cart), intent(in) :: Cart
    real(rk), dimension(Cart%nd), intent(in), optional :: PeriodicLength
    integer, intent(in), optional :: GeometryType

    integer :: d, i, j, k
    real(rk), dimension(MAX_DIMS) :: PeriodicLength_
    integer :: GeometryType_
    integer, dimension(MAX_DIMS) :: Point
    type(ovk_cart) :: CellEdgeDistCart

    Grid%logger = Logger

    if (present(PeriodicLength)) then
      PeriodicLength_(:Cart%nd) = PeriodicLength
      PeriodicLength_(Cart%nd+1:) = 0._rk
    else
      PeriodicLength_ = 0._rk
    end if

    if (present(GeometryType)) then
      GeometryType_ = GeometryType
    else
      GeometryType_ = OVK_GEOMETRY_CURVILINEAR
    end if

    Grid%id = ID
    Grid%nd = Cart%nd
    Grid%npoints(:Cart%nd) = ovkCartSize(Cart)
    Grid%npoints(Cart%nd+1:) = 1
    Grid%cart = Cart
    Grid%cell_cart = ovkCartPointToCell(Grid%cart)
    Grid%periodic_length = PeriodicLength_
    Grid%geometry_type = GeometryType_

    allocate(Grid%coords(Grid%nd))
    do d = 1, Grid%nd
      Grid%coords(d) = ovk_field_real_(Grid%cart)
      do k = Grid%cart%is(3), Grid%cart%ie(3)
        do j = Grid%cart%is(2), Grid%cart%ie(2)
          do i = Grid%cart%is(1), Grid%cart%ie(1)
            Point = [i,j,k]
            Grid%coords(d)%values(i,j,k) = real(Point(d)-1,kind=rk)
          end do
        end do
      end do
    end do

    Grid%coords_edit_ref_count = 0

    allocate(Grid%state)
    Grid%state = ovk_field_int_(Grid%cart, OVK_INTERIOR_POINT)

    nullify(Grid%prev_state)

    Grid%state_edit_ref_count = 0

    allocate(Grid%edits)
    Grid%edits = t_grid_edits_(Cart%nd)

    Grid%bounds = ovk_bbox_(Grid%nd)
    do d = 1, Grid%nd
      Grid%bounds%b(d) = real(Grid%cart%is(d)-1,kind=rk)
      if (Grid%cart%periodic(d) .and. Grid%cart%periodic_storage == OVK_PERIODIC_STORAGE_UNIQUE) &
        then
        Grid%bounds%e(d) = real(Grid%cart%ie(d),kind=rk)
      else
        Grid%bounds%e(d) = real(Grid%cart%ie(d)-1,kind=rk)
      end if
    end do

    Grid%mask = ovk_field_logical_(Grid%cart, .true.)
    Grid%cell_mask = ovk_field_logical_(Grid%cell_cart, .true.)
    Grid%volumes = ovk_field_real_(Grid%cart, 1._rk)
    Grid%cell_volumes = ovk_field_real_(Grid%cell_cart, 1._rk)

    CellEdgeDistCart = Grid%cell_cart
    CellEdgeDistCart%is(:Grid%nd) = CellEdgeDistCart%is(:Grid%nd) - merge(0, 1, &
      Grid%cart%periodic(:Grid%nd))
    CellEdgeDistCart%ie(:Grid%nd) = CellEdgeDistCart%ie(:Grid%nd) + merge(0, 1, &
      Grid%cart%periodic(:Grid%nd))

    call SetExists(Grid%existence_flag, .true.)

  end subroutine CreateGrid

  subroutine DestroyGrid(Grid)

    type(ovk_grid), intent(inout) :: Grid

    if (.not. ovkGridExists(Grid)) return

    call SetExists(Grid%existence_flag, .false.)

    Grid%mask = ovk_field_logical_()
    Grid%cell_mask = ovk_field_logical_()
    Grid%volumes = ovk_field_real_()
    Grid%cell_volumes = ovk_field_real_()

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

  subroutine ovkGetGridCart(Grid, Cart)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_cart), intent(out) :: Cart

    Cart = Grid%cart

  end subroutine ovkGetGridCart

  subroutine ovkGetGridPeriodicity(Grid, Periodic)

    type(ovk_grid), intent(in) :: Grid
    logical, dimension(Grid%nd), intent(out) :: Periodic

    Periodic = Grid%cart%periodic(:Grid%nd)

  end subroutine ovkGetGridPeriodicity

  subroutine ovkGetGridPeriodicStorage(Grid, PeriodicStorage)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(out) :: PeriodicStorage

    PeriodicStorage = Grid%cart%periodic_storage

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
          if (.not. EditingState(Grid)) then
            call UpdateBounds(Grid)
          end if
          call UpdateVolumes(Grid)
          Grid%edits%coords = .true.
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
    logical :: ModifiedBoundary
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

          ModifiedBoundary = .false.
          L2: &
          do k = Grid%cart%is(3), Grid%cart%ie(3)
            do j = Grid%cart%is(2), Grid%cart%ie(2)
              do i = Grid%cart%is(1), Grid%cart%ie(1)
                if (State%values(i,j,k) /= PrevState%values(i,j,k)) then
                  StateValue = iand(State%values(i,j,k), ior(OVK_STATE_DOMAIN_BOUNDARY, &
                    OVK_STATE_INTERNAL_BOUNDARY))
                  PrevStateValue = iand(PrevState%values(i,j,k), ior(OVK_STATE_DOMAIN_BOUNDARY, &
                    OVK_STATE_INTERNAL_BOUNDARY))
                  if (StateValue /= PrevStateValue) then
                    ModifiedBoundary = .true.
                    exit L2
                  end if
                end if
              end do
            end do
          end do L2

          if (ModifiedMask) then
            call UpdateMasks(Grid)
            if (.not. EditingCoords(Grid)) then
              call UpdateBounds(Grid)
            end if
            Grid%edits%mask = .true.
          end if

          if (ModifiedBoundary) then
            Grid%edits%boundary = .true.
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

    call ovkEditGridState(Grid, State)

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
            State%values(i,j,k) = iand(State%values(i,j,k), not(OVK_STATE_EXTERIOR))
            State%values(i,j,k) = ior(State%values(i,j,k), OVK_STATE_GRID)
          end if
          State%values(i,j,k) = iand(State%values(i,j,k),not(ASSEMBLY_STATES))
          State%values(i,j,k) = iand(State%values(i,j,k),not(DEBUG_STATES))
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
    case (OVK_NOT_ALL)
      do k = State%cart%is(3), State%cart%ie(3)
        do j = State%cart%is(2), State%cart%ie(2)
          do i = State%cart%is(1), State%cart%ie(1)
            MatchingMask%values(i,j,k) = iand(State%values(i,j,k), StateBits) /= StateBits
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

  subroutine ovkGetGridBounds(Grid, Bounds)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_bbox), intent(out) :: Bounds

    Bounds = Grid%bounds

  end subroutine ovkGetGridBounds

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
    end if

  end subroutine TryReleaseState

  subroutine UpdateBounds(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: i, j, d, l
    integer :: idir, jdir
    integer, dimension(MAX_DIMS) :: Point, AdjustedPoint
    real(rk), dimension(Grid%nd) :: PeriodicCoords

    Grid%bounds = ovk_bbox_(Grid%nd)
    do d = 1, Grid%nd
      Grid%bounds%b(d) = minval(Grid%coords(d)%values)
      Grid%bounds%e(d) = maxval(Grid%coords(d)%values)
    end do

    ! Add contribution from periodic points
    if (Grid%cart%periodic_storage == OVK_PERIODIC_STORAGE_UNIQUE .and. &
      any(Grid%cart%periodic .and. Grid%periodic_length > 0._rk)) then
      do d = 1, Grid%nd
        if (Grid%cart%periodic(d)) then
          idir = modulo((d+1)-1,MAX_DIMS) + 1
          jdir = modulo((d+2)-1,MAX_DIMS) + 1
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

  subroutine UpdateMasks(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: m, n, o
    integer(lk) :: l
    integer(lk) :: NumCells
    integer, dimension(MAX_DIMS) :: Cell
    integer, dimension(MAX_DIMS) :: VertexLower, VertexUpper
    integer, dimension(MAX_DIMS) :: Vertex
    logical :: AwayFromEdge

    call ovkFilterGridState(Grid, OVK_STATE_GRID, OVK_ALL, Grid%mask)

    NumCells = ovkCartCount(Grid%cell_cart)

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(NONE) &
!$OMP&  PRIVATE(m, n, o, l, Cell, VertexLower, VertexUpper, AwayFromEdge, Vertex) &
!$OMP&  FIRSTPRIVATE(NumCells) &
!$OMP&  SHARED(Grid)
    do l = 1_lk, NumCells
      Cell(:Grid%nd) = ovkCartIndexToTuple(Grid%cell_cart, l)
      Cell(Grid%nd+1:) = 1
      VertexLower = Cell
      VertexUpper(:Grid%nd) = Cell(:Grid%nd) + 1
      VertexUpper(Grid%nd+1:) = 1
      AwayFromEdge = ovkCartContains(Grid%cart, VertexUpper)
      if (AwayFromEdge) then
        Grid%cell_mask%values(Cell(1),Cell(2),Cell(3)) = .true.
        L1: &
        do o = VertexLower(3), VertexUpper(3)
          do n = VertexLower(2), VertexUpper(2)
            do m = VertexLower(1), VertexUpper(1)
              if (.not. Grid%mask%values(m,n,o)) then
                Grid%cell_mask%values(Cell(1),Cell(2),Cell(3)) = .false.
                exit L1
              end if
            end do
          end do
        end do L1
      else
        Grid%cell_mask%values(Cell(1),Cell(2),Cell(3)) = .true.
        L2: &
        do o = VertexLower(3), VertexUpper(3)
          do n = VertexLower(2), VertexUpper(2)
            do m = VertexLower(1), VertexUpper(1)
              Vertex = [m,n,o]
              Vertex(:Grid%nd) = ovkCartPeriodicAdjust(Grid%cart, Vertex)
              if (.not. Grid%mask%values(Vertex(1),Vertex(2),Vertex(3))) then
                Grid%cell_mask%values(Cell(1),Cell(2),Cell(3)) = .false.
                exit L2
              end if
            end do
          end do
        end do L2
      end if
    end do
!$OMP END PARALLEL DO

  end subroutine UpdateMasks

  subroutine UpdateVolumes(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: m, n, o
    integer(lk) :: l
    integer(lk) :: NumPoints
    integer(lk) :: NumCells
    integer, dimension(MAX_DIMS) :: Cell
    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords
    integer, dimension(MAX_DIMS) :: Point
    integer, dimension(MAX_DIMS) :: NeighborCellLower, NeighborCellUpper
    integer, dimension(MAX_DIMS) :: Neighbor
    integer :: NumNeighbors
    logical :: AwayFromEdge

    NumPoints = ovkCartCount(Grid%cart)
    NumCells = ovkCartCount(Grid%cell_cart)

!$OMP PARALLEL DO &
!$OMP&  DEFAULT(NONE) &
!$OMP&  PRIVATE(l, Cell, VertexCoords) &
!$OMP&  FIRSTPRIVATE(NumCells) &
!$OMP&  SHARED(Grid)
    do l = 1_lk, NumCells
      Cell(:Grid%nd) = ovkCartIndexToTuple(Grid%cell_cart, l)
      Cell(Grid%nd+1:) = 1
      call GetCellVertexCoordsLinear(Grid, Cell, VertexCoords)
      select case (Grid%geometry_type)
      case (OVK_GEOMETRY_UNIFORM,OVK_GEOMETRY_RECTILINEAR)
        select case (Grid%nd)
        case (2)
          Grid%cell_volumes%values(Cell(1),Cell(2),Cell(3)) = ovkRectangleSize(VertexCoords)
        case (3)
          Grid%cell_volumes%values(Cell(1),Cell(2),Cell(3)) = ovkCuboidSize(VertexCoords)
        end select
      case (OVK_GEOMETRY_ORIENTED_UNIFORM,OVK_GEOMETRY_ORIENTED_RECTILINEAR)
        select case (Grid%nd)
        case (2)
          Grid%cell_volumes%values(Cell(1),Cell(2),Cell(3)) = ovkOrientedRectangleSize(VertexCoords)
        case (3)
          Grid%cell_volumes%values(Cell(1),Cell(2),Cell(3)) = ovkOrientedCuboidSize(VertexCoords)
        end select
      case default
        select case (Grid%nd)
        case (2)
          Grid%cell_volumes%values(Cell(1),Cell(2),Cell(3)) = ovkQuadSize(VertexCoords)
        case (3)
          Grid%cell_volumes%values(Cell(1),Cell(2),Cell(3)) = ovkHexahedronSize(VertexCoords)
        end select
      end select
    end do
!$OMP END PARALLEL DO

    ! Compute the volume at each point by averaging the volumes of neighboring cells
!$OMP PARALLEL DO &
!$OMP&  DEFAULT(NONE) &
!$OMP&  PRIVATE(m, n, o, l, Point, NeighborCellLower, NeighborCellUpper, NumNeighbors) &
!$OMP&  PRIVATE(AwayFromEdge, Neighbor) &
!$OMP&  FIRSTPRIVATE(NumPoints) &
!$OMP&  SHARED(Grid)
    do l = 1_lk, NumPoints
      Point(:Grid%nd) = ovkCartIndexToTuple(Grid%cart, l)
      Point(Grid%nd+1:) = 1
      NeighborCellLower(:Grid%nd) = Point(:Grid%nd)-1
      NeighborCellLower(Grid%nd+1:) = 1
      NeighborCellUpper(:Grid%nd) = Point(:Grid%nd)
      NeighborCellUpper(Grid%nd+1:) = 1
      Grid%volumes%values(Point(1),Point(2),Point(3)) = 0._rk
      NumNeighbors = 0
      AwayFromEdge = ovkCartContains(Grid%cell_cart, NeighborCellLower) .and. &
        ovkCartContains(Grid%cell_cart, NeighborCellUpper)
      if (AwayFromEdge) then
        do o = NeighborCellLower(3), NeighborCellUpper(3)
          do n = NeighborCellLower(2), NeighborCellUpper(2)
            do m = NeighborCellLower(1), NeighborCellUpper(1)
              Grid%volumes%values(Point(1),Point(2),Point(3)) = Grid%volumes%values(Point(1),&
                Point(2),Point(3)) + Grid%cell_volumes%values(m,n,o)
              NumNeighbors = NumNeighbors + 1
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
                Grid%volumes%values(Point(1),Point(2),Point(3)) = Grid%volumes%values(Point(1), &
                  Point(2),Point(3)) + Grid%cell_volumes%values(Neighbor(1),Neighbor(2),Neighbor(3))
                NumNeighbors = NumNeighbors + 1
              end if
            end do
          end do
        end do
      end if
      Grid%volumes%values(Point(1),Point(2),Point(3)) = Grid%volumes%values(Point(1),Point(2), &
        Point(3))/real(max(NumNeighbors,1),kind=rk)
    end do
!$OMP END PARALLEL DO

  end subroutine UpdateVolumes

  subroutine GetGridEdits(Grid, Edits)

    type(ovk_grid), intent(in) :: Grid
    type(t_grid_edits), pointer, intent(out) :: Edits

    Edits => Grid%edits

  end subroutine GetGridEdits

  subroutine ResetGridEdits(Grid)

    type(ovk_grid), intent(inout) :: Grid

    Grid%edits%coords = .false.
    Grid%edits%mask = .false.
    Grid%edits%boundary = .false.

  end subroutine ResetGridEdits

  subroutine GetCellVertexCoordsLinear(Grid, Cell, VertexCoords)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    real(rk), dimension(Grid%nd,2**Grid%nd), intent(out) :: VertexCoords

    integer :: i, j, k, l, m
    integer, dimension(MAX_DIMS) :: VertexLower, VertexUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_DIMS) :: Vertex
    integer, dimension(MAX_DIMS) :: AdjustedVertex
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

  end subroutine GetCellVertexCoordsLinear

  subroutine GetCellVertexCoordsCubic(Grid, Cell, VertexCoords)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell
    real(rk), dimension(Grid%nd,4**Grid%nd), intent(out) :: VertexCoords

    integer :: i, j, k, l, m
    integer, dimension(MAX_DIMS) :: VertexLower, VertexUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_DIMS) :: Vertex
    integer, dimension(MAX_DIMS) :: AdjustedVertex
    real(rk), dimension(Grid%nd) :: PrincipalCoords

    VertexLower(:Grid%nd) = Cell-1
    VertexLower(Grid%nd+1:) = 1
    VertexUpper(:Grid%nd) = Cell+2
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

  end subroutine GetCellVertexCoordsCubic

  function ovkGridCellBounds(Grid, Cell_) result(CellBounds)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell_
    type(ovk_bbox) :: CellBounds

    integer, dimension(MAX_DIMS) :: Cell
    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords

    Cell(:Grid%nd) = Cell_
    Cell(Grid%nd+1:) = 1

    call GetCellVertexCoordsLinear(Grid, Cell, VertexCoords)

    CellBounds = ovkBBFromPoints(VertexCoords)

  end function ovkGridCellBounds

  function ovkOverlapsGridCell(Grid, Cell_, Coords, OverlapTolerance) result(Overlaps)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell_
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    real(rk), intent(in) :: OverlapTolerance
    logical :: Overlaps

    integer :: i
    integer, dimension(MAX_DIMS) :: Cell
    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords
    real(rk), dimension(Grid%nd) :: Centroid

    Cell(:Grid%nd) = Cell_
    Cell(Grid%nd+1:) = 1

    if (.not. Grid%cell_mask%values(Cell(1),Cell(2),Cell(3))) then
      Overlaps = .false.
      return
    end if

    call GetCellVertexCoordsLinear(Grid, Cell, VertexCoords)

    if (OverlapTolerance > 0._rk) then
      Centroid = sum(VertexCoords,dim=2)/2._rk**Grid%nd
      do i = 1, 2**Grid%nd
        VertexCoords(:,i) = Centroid + (1._rk+OverlapTolerance) * (VertexCoords(:,i)-Centroid)
      end do
    end if

    select case (Grid%geometry_type)
    case (OVK_GEOMETRY_UNIFORM,OVK_GEOMETRY_RECTILINEAR)
      select case (Grid%nd)
      case (2)
        Overlaps = ovkOverlapsRectangle(VertexCoords, Coords)
      case (3)
        Overlaps = ovkOverlapsCuboid(VertexCoords, Coords)
      end select
    case (OVK_GEOMETRY_ORIENTED_UNIFORM,OVK_GEOMETRY_ORIENTED_RECTILINEAR)
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

  function ovkCoordsInGridCell(Grid, Cell_, Coords, Success) result(CoordsInCell)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%nd), intent(in) :: Cell_
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    logical, intent(out), optional :: Success
    real(rk), dimension(Grid%nd) :: CoordsInCell

    logical :: Success_
    integer, dimension(MAX_DIMS) :: Cell

    Success_ = .true.

    Cell(:Grid%nd) = Cell_
    Cell(Grid%nd+1:) = 1

    select case (Grid%geometry_type)
    case (OVK_GEOMETRY_UNIFORM)
      CoordsInCell = CoordsInGridCellLinear(Grid, Cell, Coords)
    case (OVK_GEOMETRY_ORIENTED_UNIFORM)
      CoordsInCell = CoordsInGridCellLinearOriented(Grid, Cell, Coords)
    case default
      CoordsInCell = CoordsInGridCellCubic(Grid, Cell, Coords, Success_)
    end select

    if (present(Success)) then
      Success = Success_
    end if

  end function ovkCoordsInGridCell

  function CoordsInGridCellLinear(Grid, Cell, Coords) result(CoordsInCell)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(:), intent(in) :: Cell
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    real(rk), dimension(Grid%nd) :: CoordsInCell

    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords

    call GetCellVertexCoordsLinear(Grid, Cell, VertexCoords)

    select case (Grid%nd)
    case (2)
      CoordsInCell = ovkRectangleIsoInverseLinear(VertexCoords, Coords)
    case (3)
      CoordsInCell = ovkCuboidIsoInverseLinear(VertexCoords, Coords)
    end select

  end function CoordsInGridCellLinear

  function CoordsInGridCellLinearOriented(Grid, Cell, Coords) result(CoordsInCell)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(:), intent(in) :: Cell
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    real(rk), dimension(Grid%nd) :: CoordsInCell

    real(rk), dimension(Grid%nd,2**Grid%nd) :: VertexCoords

    call GetCellVertexCoordsLinear(Grid, Cell, VertexCoords)

    select case (Grid%nd)
    case (2)
      CoordsInCell = ovkOrientedRectangleIsoInverseLinear(VertexCoords, Coords)
    case (3)
      CoordsInCell = ovkOrientedCuboidIsoInverseLinear(VertexCoords, Coords)
    end select

  end function CoordsInGridCellLinearOriented

  function CoordsInGridCellCubic(Grid, Cell, Coords, Success) result(CoordsInCell)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(:), intent(in) :: Cell
    real(rk), dimension(Grid%nd), intent(in) :: Coords
    logical, intent(out) :: Success
    real(rk), dimension(Grid%nd) :: CoordsInCell

    integer :: d
    integer, dimension(MAX_DIMS) :: StencilShift
    integer, dimension(MAX_DIMS) :: ShiftedCell
    real(rk), dimension(Grid%nd,4**Grid%nd) :: VertexCoords

    StencilShift = 0
    do d = 1, Grid%nd
      if (.not. Grid%cart%periodic(d)) then
        if (Cell(d) == Grid%cell_cart%is(d)) then
          StencilShift(d) = StencilShift(d) + 1
        else if (Cell(d) == Grid%cell_cart%ie(d)) then
          StencilShift(d) = StencilShift(d) - 1
        end if
      end if
    end do

    ShiftedCell = Cell + StencilShift

    call GetCellVertexCoordsCubic(Grid, ShiftedCell, VertexCoords)

    select case (Grid%nd)
    case (2)
      CoordsInCell = ovkQuadIsoInverseCubic(VertexCoords, Coords, Success=Success)
    case (3)
      CoordsInCell = ovkHexahedronIsoInverseCubic(VertexCoords, Coords, Success=Success)
    end select

    CoordsInCell = CoordsInCell + real(StencilShift(:Grid%nd),kind=rk)

  end function CoordsInGridCellCubic

  pure function ovkPeriodicExtend(Cart, PeriodicLength, Point, PrincipalCoords) &
    result(ExtendedCoords)

    type(ovk_cart), intent(in) :: Cart
    real(rk), dimension(Cart%nd), intent(in) :: PeriodicLength
    integer, dimension(Cart%nd), intent(in) :: Point
    real(rk), dimension(Cart%nd), intent(in) :: PrincipalCoords
    real(rk), dimension(Cart%nd) :: ExtendedCoords

    integer, dimension(Cart%nd) :: PeriodStart, PeriodEnd
    real(rk), dimension(Cart%nd) :: PositiveAdjustment, NegativeAdjustment

    PeriodStart = Cart%is(:Cart%nd)
    if (Cart%periodic_storage == OVK_PERIODIC_STORAGE_UNIQUE) then
      PeriodEnd = Cart%ie(:Cart%nd)
    else
      PeriodEnd = Cart%ie(:Cart%nd)-1
    end if

    PositiveAdjustment = merge(PeriodicLength, 0._rk, Cart%periodic(:Cart%nd) .and. &
      Point > PeriodEnd)
    NegativeAdjustment = merge(PeriodicLength, 0._rk, Cart%periodic(:Cart%nd) .and. &
      Point < PeriodStart)

    ExtendedCoords = PrincipalCoords + PositiveAdjustment - NegativeAdjustment

  end function ovkPeriodicExtend

  function t_grid_edits_(NumDims) result(Edits)

    integer :: NumDims
    type(t_grid_edits) :: Edits

    Edits%coords = .false.
    Edits%mask = .false.
    Edits%boundary = .false.

  end function t_grid_edits_

end module ovkGrid
