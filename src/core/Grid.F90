! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkGrid

  use ovkBoundingBox
  use ovkCart
  use ovkField
  use ovkFieldOps
  use ovkGeometry
  use ovkGlobal
  implicit none

  private

  ! API
  public :: ovk_grid
  public :: ovk_grid_
  public :: ovk_grid_properties
  public :: ovk_grid_properties_
  public :: ovk_grid_event_flags
  public :: ovk_grid_event_flags_
  public :: ovk_grid_description
  public :: ovk_grid_description_
  public :: ovkCreateGrid
  public :: ovkDestroyGrid
  public :: ovkGetGridProperties
  public :: ovkEditGridProperties
  public :: ovkReleaseGridProperties
  public :: ovkGetGridDescription
  public :: ovkGetGridCart
  public :: ovkGetGridCoords
  public :: ovkEditGridCoords
  public :: ovkReleaseGridCoords
  public :: ovkGetGridMask
  public :: ovkEditGridMask
  public :: ovkReleaseGridMask
  public :: ovkGetGridBoundaryMask
  public :: ovkEditGridBoundaryMask
  public :: ovkReleaseGridBoundaryMask
  public :: ovkGetGridInternalBoundaryMask
  public :: ovkEditGridInternalBoundaryMask
  public :: ovkReleaseGridInternalBoundaryMask
  public :: ovkGetCellVertexData
  public :: ovkOverlapsCell
  public :: ovkCoordsInCell
  public :: ovkGridResolution
  public :: ovkGenerateBBOverlapMask
  public :: ovkPeriodicExtend
  public :: ovkExportGridCoords
  public :: ovkGetGridPropertyID
  public :: ovkGetGridPropertyDimension
  public :: ovkGetGridPropertySize
  public :: ovkGetGridPropertyPeriodicity
  public :: ovkGetGridPropertyPeriodicStorage
  public :: ovkGetGridPropertyPeriodicLength
  public :: ovkGetGridPropertyGeometryType
  public :: ovkGetGridPropertyVerbose
  public :: ovkSetGridPropertyVerbose
  public :: ovkGetGridPropertyMaxEdgeDistance
  public :: ovkSetGridPropertyMaxEdgeDistance
  public :: ovkResetGridEventFlags
  public :: OVK_GRID_GEOMETRY_CARTESIAN
  public :: OVK_GRID_GEOMETRY_CARTESIAN_ROTATED
  public :: OVK_GRID_GEOMETRY_RECTILINEAR
  public :: OVK_GRID_GEOMETRY_RECTILINEAR_ROTATED
  public :: OVK_GRID_GEOMETRY_CURVILINEAR

  type ovk_grid_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: id
    integer :: nd
    integer, dimension(MAX_ND) :: npoints
    logical, dimension(MAX_ND) :: periodic
    integer :: periodic_storage
    real(rk), dimension(MAX_ND) :: periodic_length
    integer :: geometry_type
    ! Read/write
    logical :: verbose
    integer :: max_edge_dist
  end type ovk_grid_properties

  type t_grid_editor
    integer :: properties_ref_count
    integer, dimension(:), allocatable :: xyz_ref_count
    integer :: grid_mask_ref_count
    integer :: boundary_mask_ref_count
    integer :: internal_boundary_mask_ref_count
  end type t_grid_editor

  type ovk_grid_event_flags
    type(t_noconstruct) :: noconstruct
    logical :: modified_xyz
    logical :: modified_grid_mask
    logical :: modified_boundary_mask
    logical :: modified_internal_boundary_mask
  end type ovk_grid_event_flags

  type ovk_grid_description
    type(t_noconstruct) :: noconstruct
    integer :: nd
    integer, dimension(MAX_ND) :: npoints
    logical, dimension(MAX_ND) :: periodic
    integer :: periodic_storage
    real(rk), dimension(MAX_ND) :: periodic_length
    integer :: geometry_type
  end type ovk_grid_description

  type ovk_grid
    type(t_noconstruct) :: noconstruct
    type(ovk_grid_properties), pointer :: properties
    type(ovk_grid_properties) :: prev_properties
    type(ovk_cart) :: cart
    type(ovk_cart) :: cell_cart
    type(ovk_field_real), dimension(:), pointer :: xyz
    type(ovk_field_logical), pointer :: grid_mask
    type(ovk_field_logical), pointer :: boundary_mask
    type(ovk_field_logical), pointer :: internal_boundary_mask
    type(ovk_bbox) :: bounds
    type(ovk_field_logical) :: cell_grid_mask
    type(ovk_field_real) :: resolution
    type(ovk_field_int) :: edge_dist
    type(ovk_field_int) :: cell_edge_dist
    type(t_grid_editor) :: editor
    type(ovk_grid_event_flags), pointer :: event_flags
  end type ovk_grid

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_grid_
    module procedure ovk_grid_Default
  end interface ovk_grid_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_grid_description_
    module procedure ovk_grid_description_Default
    module procedure ovk_grid_description_Assigned
  end interface ovk_grid_description_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_grid_properties_
    module procedure ovk_grid_properties_Default
    module procedure ovk_grid_properties_Assigned
  end interface ovk_grid_properties_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_grid_editor_
    module procedure t_grid_editor_Default
    module procedure t_grid_editor_Assigned
  end interface t_grid_editor_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_grid_event_flags_
    module procedure ovk_grid_event_flags_Default
    module procedure ovk_grid_event_flags_Assigned
  end interface ovk_grid_event_flags_

  integer, parameter :: OVK_GRID_GEOMETRY_CARTESIAN = 1
  integer, parameter :: OVK_GRID_GEOMETRY_CARTESIAN_ROTATED = 2
  integer, parameter :: OVK_GRID_GEOMETRY_RECTILINEAR = 3
  integer, parameter :: OVK_GRID_GEOMETRY_RECTILINEAR_ROTATED = 4
  integer, parameter :: OVK_GRID_GEOMETRY_CURVILINEAR = 5

contains

  function ovk_grid_Default() result(Grid)

    type(ovk_grid) :: Grid

    nullify(Grid%properties)
    Grid%prev_properties = ovk_grid_properties_()
    Grid%cart = ovk_cart_()
    Grid%cell_cart = ovk_cart_()
    nullify(Grid%xyz)
    nullify(Grid%grid_mask)
    nullify(Grid%boundary_mask)
    nullify(Grid%internal_boundary_mask)
    Grid%bounds = ovk_bbox_()
    Grid%cell_grid_mask = ovk_field_logical_()
    Grid%resolution = ovk_field_real_()
    Grid%edge_dist = ovk_field_int_()
    Grid%cell_edge_dist = ovk_field_int_()
    Grid%editor = t_grid_editor_()
    nullify(Grid%event_flags)

  end function ovk_grid_Default

  subroutine ovkCreateGrid(Grid, ID, NumDims, NumPoints, Periodic, PeriodicStorage, &
    PeriodicLength, GeometryType, Verbose, EventFlags)

    type(ovk_grid), intent(out) :: Grid
    integer, intent(in) :: ID
    integer, intent(in) :: NumDims
    integer, dimension(NumDims), intent(in) :: NumPoints
    logical, dimension(NumDims), intent(in), optional :: Periodic
    integer, intent(in), optional :: PeriodicStorage
    real(rk), dimension(NumDims), intent(in), optional :: PeriodicLength
    integer, intent(in), optional :: GeometryType
    logical, intent(in), optional :: Verbose
    type(ovk_grid_event_flags), target, intent(inout), optional :: EventFlags

    logical, dimension(MAX_ND) :: Periodic_
    integer :: PeriodicStorage_
    real(rk), dimension(MAX_ND) :: PeriodicLength_
    integer :: GeometryType_
    logical :: Verbose_
    integer :: d

    if (present(Periodic)) then
      Periodic_(:NumDims) = Periodic
      Periodic_(NumDims+1:) = .false.
    else
      Periodic_ = .false.
    end if

    if (present(PeriodicStorage)) then
      PeriodicStorage_ = PeriodicStorage
    else
      PeriodicStorage_ = OVK_NO_OVERLAP_PERIODIC
    end if

    if (present(PeriodicLength)) then
      PeriodicLength_(:NumDims) = PeriodicLength
      PeriodicLength_(NumDims+1:) = 0._rk
    else
      PeriodicLength_ = 0._rk
    end if

    if (present(GeometryType)) then
      GeometryType_ = GeometryType
    else
      GeometryType_ = OVK_GRID_GEOMETRY_CURVILINEAR
    end if

    if (present(Verbose)) then
      Verbose_ = Verbose
    else
      Verbose_ = .false.
    end if

    allocate(Grid%properties)
    Grid%properties = ovk_grid_properties_(NumDims)
    Grid%properties%id = ID
    Grid%properties%npoints(:NumDims) = NumPoints
    Grid%properties%periodic = Periodic_
    Grid%properties%periodic_storage = PeriodicStorage_
    Grid%properties%periodic_length = PeriodicLength_
    Grid%properties%geometry_type = GeometryType_
    Grid%properties%verbose = Verbose_
    Grid%properties%max_edge_dist = 1

    Grid%prev_properties = Grid%properties

    Grid%cart = ovk_cart_(NumDims, NumPoints, Periodic_, PeriodicStorage_)
    Grid%cart = ovkCartConvertPeriodicStorage(Grid%cart, OVK_NO_OVERLAP_PERIODIC)

    Grid%cell_cart = ovkCartPointToCell(Grid%cart)

    allocate(Grid%xyz(NumDims))
    do d = 1, NumDims
      Grid%xyz(d) = ovk_field_real_(Grid%cart, 0._rk)
    end do

    allocate(Grid%grid_mask)
    Grid%grid_mask = ovk_field_logical_(Grid%cart, .true.)

    allocate(Grid%boundary_mask)
    Grid%boundary_mask = ovk_field_logical_(Grid%cart, .false.)

    allocate(Grid%internal_boundary_mask)
    Grid%internal_boundary_mask = ovk_field_logical_(Grid%cart, .false.)

    Grid%bounds = ovk_bbox_(NumDims)
    if (ovkCartCount(Grid%cart) > 0) then
      Grid%bounds%b(:NumDims) = 0._rk
      Grid%bounds%e(:NumDims) = 0._rk
    end if

    Grid%cell_grid_mask = ovk_field_logical_(Grid%cell_cart, .true.)
    Grid%resolution = ovk_field_real_(Grid%cart, 0._rk)

    Grid%edge_dist = ovk_field_int_(Grid%cart, 1)
    Grid%cell_edge_dist = ovk_field_int_(Grid%cell_cart, 1)

    Grid%editor = t_grid_editor_(NumDims)

    if (present(EventFlags)) then
      Grid%event_flags => EventFlags
      Grid%event_flags = ovk_grid_event_flags_(NumDims)
    else
      nullify(Grid%event_flags)
    end if

  end subroutine ovkCreateGrid

  subroutine ovkDestroyGrid(Grid)

    type(ovk_grid), intent(inout) :: Grid

    if (associated(Grid%properties)) deallocate(Grid%properties)
    Grid%prev_properties = ovk_grid_properties_()

    if (associated(Grid%xyz)) deallocate(Grid%xyz)
    if (associated(Grid%grid_mask)) deallocate(Grid%grid_mask)
    if (associated(Grid%boundary_mask)) deallocate(Grid%boundary_mask)
    if (associated(Grid%internal_boundary_mask)) deallocate(Grid%internal_boundary_mask)

    Grid%cell_grid_mask = ovk_field_logical_()
    Grid%resolution = ovk_field_real_()
    Grid%edge_dist = ovk_field_int_()
    Grid%cell_edge_dist = ovk_field_int_()

    Grid%editor = t_grid_editor_()

  end subroutine ovkDestroyGrid

  subroutine ovkGetGridProperties(Grid, Properties)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_grid_properties), pointer, intent(out) :: Properties

    Properties => Grid%properties

  end subroutine ovkGetGridProperties

  subroutine ovkEditGridProperties(Grid, Properties)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_grid_properties), pointer, intent(out) :: Properties

    integer :: d
    logical :: EditingCoords
    logical :: EditingGridMask
    logical :: EditingBoundaryMask
    logical :: EditingInternalBoundaryMask
    logical :: Success

    EditingCoords = .false.
    do d = 1, Grid%cart%nd
      EditingCoords = EditingCoords .or. Grid%editor%xyz_ref_count(d) > 0
    end do
    EditingGridMask = Grid%editor%grid_mask_ref_count > 0
    EditingBoundaryMask = Grid%editor%boundary_mask_ref_count > 0
    EditingInternalBoundaryMask = Grid%editor%internal_boundary_mask_ref_count > 0

    Success = &
      .not. EditingCoords .and. &
      .not. EditingGridMask .and. &
      .not. EditingBoundaryMask .and. &
      .not. EditingInternalBoundaryMask

    if (Success) then

      if (Grid%editor%properties_ref_count == 0) then
        Grid%prev_properties = Grid%properties
      end if

      Grid%editor%properties_ref_count = Grid%editor%properties_ref_count + 1

      Properties => Grid%properties

    else

      if (OVK_DEBUG) then
        if (EditingCoords) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid properties while editing coordinates."
        else if (EditingGridMask) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid properties while editing grid mask."
        else if (EditingBoundaryMask) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid properties while editing boundary mask."
        else if (EditingInternalBoundaryMask) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid properties while editing internal boundary mask."
        end if
        stop 1
      end if

      nullify(Properties)

    end if

  end subroutine ovkEditGridProperties

  subroutine ovkReleaseGridProperties(Grid, Properties)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_grid_properties), pointer, intent(inout) :: Properties

    if (associated(Properties, Grid%properties)) then

      nullify(Properties)

      if (Grid%editor%properties_ref_count > 0) then

        Grid%editor%properties_ref_count = Grid%editor%properties_ref_count - 1

        if (Grid%editor%properties_ref_count == 0) then
          if (Grid%properties%max_edge_dist /= Grid%prev_properties%max_edge_dist) then
            call UpdateEdgeDistance(Grid)
          end if
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release grid properties; invalid pointer."
        stop 1
      end if

    end if

  end subroutine ovkReleaseGridProperties

  subroutine ovkGetGridDescription(Grid, Description)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_grid_description), intent(out) :: Description

    Description = ovk_grid_description_(Grid%properties%nd)
    Description%npoints = Grid%properties%npoints
    Description%periodic = Grid%properties%periodic
    Description%periodic_storage = Grid%properties%periodic_storage
    Description%periodic_length = Grid%properties%periodic_length
    Description%geometry_type = Grid%properties%geometry_type

  end subroutine ovkGetGridDescription

  subroutine ovkGetGridCart(Grid, Cart)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_cart), intent(out) :: Cart

    Cart = Grid%cart

  end subroutine ovkGetGridCart

  subroutine ovkGetGridCoords(Grid, DimIndex, Coords)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(in) :: DimIndex
    type(ovk_field_real), pointer, intent(out) :: Coords

    Coords => Grid%xyz(DimIndex)

  end subroutine ovkGetGridCoords

  subroutine ovkEditGridCoords(Grid, DimIndex, Coords)

    type(ovk_grid), intent(inout) :: Grid
    integer, intent(in) :: DimIndex
    type(ovk_field_real), pointer, intent(out) :: Coords

    logical :: EditingProperties
    logical :: Success

    EditingProperties = Grid%editor%properties_ref_count > 0

    Success = .not. EditingProperties

    if (Success) then

      Grid%editor%xyz_ref_count(DimIndex) = Grid%editor%xyz_ref_count(DimIndex) + 1

      Coords => Grid%xyz(DimIndex)

    else

      if (OVK_DEBUG) then
        if (EditingProperties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid coordinates while editing properties."
        end if
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

    DimIndex = 0
    do d = 1, Grid%cart%nd
      if (associated(Coords, Grid%xyz(d))) then
        DimIndex = d
        exit
      end if
    end do

    if (DimIndex /= 0) then

      nullify(Coords)

      if (Grid%editor%xyz_ref_count(DimIndex) > 0) then

        Grid%editor%xyz_ref_count(DimIndex) = Grid%editor%xyz_ref_count(DimIndex) - 1

        if (all(Grid%editor%xyz_ref_count == 0)) then
          if (associated(Grid%event_flags)) then
            Grid%event_flags%modified_xyz = .true.
          end if
          call AttemptUpdateBounds(Grid)
          call AttemptUpdateResolution(Grid)
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release grid coordinates; invalid pointer."
        stop 1
      end if

    end if

  end subroutine ovkReleaseGridCoords

  subroutine ovkGetGridMask(Grid, GridMask)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_field_logical), pointer, intent(out) :: GridMask

    GridMask => Grid%grid_mask

  end subroutine ovkGetGridMask

  subroutine ovkEditGridMask(Grid, GridMask)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_logical), pointer, intent(out) :: GridMask

    logical :: EditingProperties
    logical :: Success

    EditingProperties = Grid%editor%properties_ref_count > 0

    Success = .not. EditingProperties

    if (Success) then

      Grid%editor%grid_mask_ref_count = Grid%editor%grid_mask_ref_count + 1

      GridMask => Grid%grid_mask

    else

      if (OVK_DEBUG) then
        if (EditingProperties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid mask while editing properties."
        end if
        stop 1
      end if

      nullify(GridMask)

    end if

  end subroutine ovkEditGridMask

  subroutine ovkReleaseGridMask(Grid, GridMask)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_logical), pointer, intent(inout) :: GridMask

    if (associated(GridMask, Grid%grid_mask)) then

      nullify(GridMask)

      if (Grid%editor%grid_mask_ref_count > 0) then

        Grid%editor%grid_mask_ref_count = Grid%editor%grid_mask_ref_count - 1

        if (Grid%editor%grid_mask_ref_count == 0) then
          if (associated(Grid%event_flags)) then
            Grid%event_flags%modified_grid_mask = .true.
          end if
          call UpdateCellGridMask(Grid)
          call UpdateEdgeDistance(Grid)
          call AttemptUpdateBounds(Grid)
          call AttemptUpdateResolution(Grid)
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release grid mask; invalid pointer."
        stop 1
      end if

    end if

  end subroutine ovkReleaseGridMask

  subroutine ovkGetGridBoundaryMask(Grid, BoundaryMask)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_field_logical), pointer, intent(out) :: BoundaryMask

    BoundaryMask => Grid%boundary_mask

  end subroutine ovkGetGridBoundaryMask

  subroutine ovkEditGridBoundaryMask(Grid, BoundaryMask)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_logical), pointer, intent(out) :: BoundaryMask

    logical :: EditingProperties
    logical :: Success

    EditingProperties = Grid%editor%properties_ref_count > 0

    Success = .not. EditingProperties

    if (Success) then

      Grid%editor%boundary_mask_ref_count = Grid%editor%boundary_mask_ref_count + 1

      BoundaryMask => Grid%boundary_mask

    else

      if (OVK_DEBUG) then
        if (EditingProperties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid boundary mask while editing properties."
        end if
        stop 1
      end if

      nullify(BoundaryMask)

    end if

  end subroutine ovkEditGridBoundaryMask

  subroutine ovkReleaseGridBoundaryMask(Grid, BoundaryMask)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_logical), pointer, intent(inout) :: BoundaryMask

    if (associated(BoundaryMask, Grid%boundary_mask)) then

      nullify(BoundaryMask)

      if (Grid%editor%boundary_mask_ref_count > 0) then

        Grid%editor%boundary_mask_ref_count = Grid%editor%boundary_mask_ref_count - 1

        if (Grid%editor%boundary_mask_ref_count == 0) then
          if (associated(Grid%event_flags)) then
            Grid%event_flags%modified_boundary_mask = .true.
          end if
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release boundary mask; invalid pointer."
        stop 1
      end if

    end if

  end subroutine ovkReleaseGridBoundaryMask

  subroutine ovkGetGridInternalBoundaryMask(Grid, InternalBoundaryMask)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_field_logical), pointer, intent(out) :: InternalBoundaryMask

    InternalBoundaryMask => Grid%internal_boundary_mask

  end subroutine ovkGetGridInternalBoundaryMask

  subroutine ovkEditGridInternalBoundaryMask(Grid, InternalBoundaryMask)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_logical), pointer, intent(out) :: InternalBoundaryMask

    logical :: EditingProperties
    logical :: Success

    EditingProperties = Grid%editor%properties_ref_count > 0

    Success = .not. EditingProperties

    if (Success) then

      Grid%editor%internal_boundary_mask_ref_count = Grid%editor%internal_boundary_mask_ref_count + 1

      InternalBoundaryMask => Grid%internal_boundary_mask

    else

      if (OVK_DEBUG) then
        if (EditingProperties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid internal boundary mask while editing properties."
        end if
        stop 1
      end if

      nullify(InternalBoundaryMask)

    end if

  end subroutine ovkEditGridInternalBoundaryMask

  subroutine ovkReleaseGridInternalBoundaryMask(Grid, InternalBoundaryMask)

    type(ovk_grid), intent(inout) :: Grid
    type(ovk_field_logical), pointer, intent(inout) :: InternalBoundaryMask

    if (associated(InternalBoundaryMask, Grid%internal_boundary_mask)) then

      nullify(InternalBoundaryMask)

      if (Grid%editor%internal_boundary_mask_ref_count > 0) then

        Grid%editor%internal_boundary_mask_ref_count = &
          Grid%editor%internal_boundary_mask_ref_count - 1

        if (Grid%editor%internal_boundary_mask_ref_count == 0) then
          if (associated(Grid%event_flags)) then
            Grid%event_flags%modified_internal_boundary_mask = .true.
          end if
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release internal boundary mask; invalid pointer."
        stop 1
      end if

    end if

  end subroutine ovkReleaseGridInternalBoundaryMask

  subroutine AttemptUpdateBounds(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: i, j, d, l
    logical :: EditingCoords
    logical :: EditingGridMask
    integer :: idir, jdir
    integer, dimension(MAX_ND) :: Point, AdjustedPoint
    real(rk), dimension(Grid%cart%nd) :: PeriodicCoords

    if (ovkCartCount(Grid%cart) == 0) return

    EditingCoords = .false.
    do d = 1, Grid%cart%nd
      EditingCoords = EditingCoords .or. Grid%editor%xyz_ref_count(d) > 0
    end do
    EditingGridMask = Grid%editor%grid_mask_ref_count > 0

    if (EditingCoords .or. EditingGridMask) return

    Grid%bounds = ovk_bbox_(Grid%cart%nd)
    do d = 1, Grid%cart%nd
      Grid%bounds%b(d) = minval(Grid%xyz(d)%values)
      Grid%bounds%e(d) = maxval(Grid%xyz(d)%values)
    end do

    ! Add contribution from periodic points
    if (any(Grid%cart%periodic .and. Grid%properties%periodic_length > 0._rk)) then
      do d = 1, Grid%cart%nd
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
              do l = 1, Grid%cart%nd
                PeriodicCoords(l) = Grid%xyz(l)%values(AdjustedPoint(1),AdjustedPoint(2), &
                  AdjustedPoint(3))
              end do
              PeriodicCoords = ovkPeriodicExtend(Grid%cart, Grid%properties%periodic_length, &
                Point, PeriodicCoords)
              Grid%bounds = ovkBBExtend(Grid%bounds, PeriodicCoords)
            end do
          end do
        end if
      end do
    end if

  end subroutine AttemptUpdateBounds

  subroutine UpdateCellGridMask(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: i, j, k, m, n, o
    integer, dimension(MAX_ND) :: VertexLower, VertexUpper
    integer, dimension(MAX_ND) :: Vertex
    logical :: AwayFromEdge

    if (ovkCartCount(Grid%cart) == 0) return

    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          VertexLower = [i,j,k]
          VertexUpper(:Grid%cart%nd) = VertexLower(:Grid%Cart%nd) + 1
          VertexUpper(Grid%cart%nd+1:) = 1
          AwayFromEdge = ovkCartContains(Grid%cart, VertexUpper)
          if (AwayFromEdge) then
            Grid%cell_grid_mask%values(i,j,k) = .true.
            L1: &
            do o = VertexLower(3), VertexUpper(3)
              do n = VertexLower(2), VertexUpper(2)
                do m = VertexLower(1), VertexUpper(1)
                  if (.not. Grid%grid_mask%values(m,n,o)) then
                    Grid%cell_grid_mask%values(i,j,k) = .false.
                    exit L1
                  end if
                end do
              end do
            end do L1
          else
            Grid%cell_grid_mask%values(i,j,k) = .true.
            L2: &
            do o = VertexLower(3), VertexUpper(3)
              do n = VertexLower(2), VertexUpper(2)
                do m = VertexLower(1), VertexUpper(1)
                  Vertex = [m,n,o]
                  Vertex(:Grid%cart%nd) = ovkCartPeriodicAdjust(Grid%cart, Vertex)
                  if (.not. Grid%grid_mask%values(Vertex(1),Vertex(2),Vertex(3))) then
                    Grid%cell_grid_mask%values(i,j,k) = .false.
                    exit L2
                  end if
                end do
              end do
            end do L2
          end if
        end do
      end do
    end do

  end subroutine UpdateCellGridMask

  subroutine AttemptUpdateResolution(Grid)

    type(ovk_grid), intent(inout) :: Grid

    integer :: i, j, k, d, m, n, o
    logical :: EditingCoords
    logical :: EditingGridMask
    type(ovk_field_real) :: CellSizes
    integer, dimension(MAX_ND) :: Cell
    real(rk), dimension(Grid%cart%nd,2**Grid%cart%nd) :: VertexCoords
    logical, dimension(2**Grid%cart%nd) :: VertexGridMaskValues
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: NeighborCellLower, NeighborCellUpper
    integer, dimension(MAX_ND) :: Neighbor
    real(rk) :: AvgCellSize
    integer :: NumCells
    logical :: AwayFromEdge

    if (ovkCartCount(Grid%cart) == 0) return

    EditingCoords = .false.
    do d = 1, Grid%cart%nd
      EditingCoords = EditingCoords .or. Grid%editor%xyz_ref_count(d) > 0
    end do
    EditingGridMask = Grid%editor%grid_mask_ref_count > 0

    if (EditingCoords .or. EditingGridMask) return

    CellSizes = ovk_field_real_(Grid%cell_cart)

    do k = Grid%cell_cart%is(3), Grid%cell_cart%ie(3)
      do j = Grid%cell_cart%is(2), Grid%cell_cart%ie(2)
        do i = Grid%cell_cart%is(1), Grid%cell_cart%ie(1)
          Cell = [i,j,k]
          call ovkGetCellVertexData(Grid, Cell, VertexCoords=VertexCoords, &
            VertexGridMaskValues=VertexGridMaskValues)
          if (all(VertexGridMaskValues)) then
            select case (Grid%cart%nd)
            case (2)
              select case (Grid%properties%geometry_type)
              case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
                CellSizes%values(i,j,k) = ovkRectangleSize(VertexCoords)
              case default
                CellSizes%values(i,j,k) = ovkQuadSize(VertexCoords)
              end select
            case (3)
              select case (Grid%properties%geometry_type)
              case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
                CellSizes%values(i,j,k) = ovkCuboidSize(VertexCoords)
              case default
                CellSizes%values(i,j,k) = ovkHexahedronSize(VertexCoords)
              end select
            end select
          else
            CellSizes%values(i,j,k) = 0._rk
          end if
        end do
      end do
    end do

    ! Compute the grid resolution at each point by averaging the sizes of neighboring cells
    Grid%resolution = ovk_field_real_(Grid%cart)

    do k = Grid%cart%is(3), Grid%cart%ie(3)
      do j = Grid%cart%is(2), Grid%cart%ie(2)
        do i = Grid%cart%is(1), Grid%cart%ie(1)
          Point = [i,j,k]
          NeighborCellLower(:Grid%cart%nd) = Point(:Grid%cart%nd)-1
          NeighborCellLower(Grid%cart%nd+1:) = 1
          NeighborCellUpper(:Grid%cart%nd) = Point(:Grid%cart%nd)
          NeighborCellUpper(Grid%cart%nd+1:) = 1
          AvgCellSize = 0._rk
          NumCells = 0
          AwayFromEdge = ovkCartContains(Grid%cell_cart, NeighborCellLower) .and. &
            ovkCartContains(Grid%cell_cart, NeighborCellUpper)
          if (AwayFromEdge) then
            do o = NeighborCellLower(3), NeighborCellUpper(3)
              do n = NeighborCellLower(2), NeighborCellUpper(2)
                do m = NeighborCellLower(1), NeighborCellUpper(1)
                  if (Grid%cell_grid_mask%values(m,n,o)) then
                    AvgCellSize = AvgCellSize + CellSizes%values(m,n,o)
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
                  Neighbor(:Grid%cart%nd) = ovkCartPeriodicAdjust(Grid%cell_cart, Neighbor)
                  if (ovkCartContains(Grid%cell_cart, Neighbor)) then
                    if (Grid%cell_grid_mask%values(Neighbor(1),Neighbor(2),Neighbor(3))) then
                      AvgCellSize = AvgCellSize + CellSizes%values(Neighbor(1),Neighbor(2), &
                        Neighbor(3))
                      NumCells = NumCells + 1
                    end if
                  end if
                end do
              end do
            end do
          end if
          AvgCellSize = AvgCellSize/real(max(NumCells,1), kind=rk)
          Grid%resolution%values(i,j,k) = AvgCellSize
        end do
      end do
    end do

  end subroutine AttemptUpdateResolution

  subroutine UpdateEdgeDistance(Grid)

    type(ovk_grid), intent(inout) :: Grid

    type(ovk_field_logical) :: NotMask

    if (ovkCartCount(Grid%cart) == 0) return

    NotMask = ovk_field_logical_(Grid%cart)
    NotMask%values = .not. Grid%grid_mask%values

    call ovkDistanceField(NotMask, Grid%properties%max_edge_dist, OVK_TRUE, Grid%edge_dist)

    NotMask = ovk_field_logical_(Grid%cell_cart)
    NotMask%values = .not. Grid%cell_grid_mask%values

    call ovkDistanceField(NotMask, Grid%properties%max_edge_dist, OVK_TRUE, Grid%cell_edge_dist)

  end subroutine UpdateEdgeDistance

  subroutine ovkGetCellVertexData(Grid, Cell, VertexCoords, VertexGridMaskValues)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%cart%nd), intent(in) :: Cell
    real(rk), dimension(Grid%cart%nd,2**Grid%cart%nd), intent(out), optional :: VertexCoords
    logical, dimension(2**Grid%cart%nd), intent(out), optional :: VertexGridMaskValues

    integer :: i, j, k, l, m
    integer, dimension(MAX_ND) :: VertexLower, VertexUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Vertex
    integer, dimension(MAX_ND) :: AdjustedVertex
    real(rk), dimension(Grid%cart%nd) :: PrincipalCoords

    VertexLower(:Grid%cart%nd) = Cell
    VertexLower(Grid%cart%nd+1:) = 1
    VertexUpper(:Grid%cart%nd) = Cell+1
    VertexUpper(Grid%cart%nd+1:) = 1

    AwayFromEdge = ovkCartContains(Grid%cart, VertexLower) .and. &
      ovkCartContains(Grid%cart, VertexUpper)

    if (AwayFromEdge) then
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
              VertexCoords(:,l) = ovkPeriodicExtend(Grid%cart, Grid%properties%periodic_length, &
                Vertex, PrincipalCoords)
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

  function ovkOverlapsCell(Grid, Cell, Coords, OverlapTolerance) result(Overlaps)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%cart%nd), intent(in) :: Cell
    real(rk), dimension(Grid%cart%nd), intent(in) :: Coords
    real(rk), intent(in), optional :: OverlapTolerance
    logical :: Overlaps

    real(rk) :: OverlapTolerance_
    integer :: i
    real(rk), dimension(Grid%cart%nd,2**Grid%cart%nd) :: VertexCoords
    logical, dimension(2**Grid%cart%nd) :: VertexGridMaskValues
    real(rk), dimension(Grid%cart%nd) :: Centroid

    if (present(OverlapTolerance)) then
      OverlapTolerance_ = OverlapTolerance
    else
      OverlapTolerance_ = 0._rk
    end if

    call ovkGetCellVertexData(Grid, Cell, VertexCoords=VertexCoords, &
      VertexGridMaskValues=VertexGridMaskValues)

    if (.not. all(VertexGridMaskValues)) then
      Overlaps = .false.
      return
    end if

    if (OverlapTolerance_ > 0._rk) then
      Centroid = sum(VertexCoords,dim=2)/2._rk**Grid%cart%nd
      do i = 1, 2**Grid%cart%nd
        VertexCoords(:,i) = Centroid + (1._rk+OverlapTolerance_) * (VertexCoords(:,i)-Centroid)
      end do
    end if

    select case (Grid%cart%nd)
    case (2)
      select case (Grid%properties%geometry_type)
      case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
        Overlaps = ovkOverlapsRectangle(VertexCoords, Coords)
      case default
        Overlaps = ovkOverlapsQuad(VertexCoords, Coords)
      end select
    case (3)
      select case (Grid%properties%geometry_type)
      case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
        Overlaps = ovkOverlapsCuboid(VertexCoords, Coords)
      case default
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
      select case (Grid%properties%geometry_type)
      case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
        CoordsInCell = ovkRectangleIsoInverseLinear(VertexCoords, Coords)
      case default
        CoordsInCell = ovkQuadIsoInverseLinear(VertexCoords, Coords, Success=Success)
      end select
    case (3)
      select case (Grid%properties%geometry_type)
      case (OVK_GRID_GEOMETRY_CARTESIAN,OVK_GRID_GEOMETRY_RECTILINEAR)
        CoordsInCell = ovkCuboidIsoInverseLinear(VertexCoords, Coords)
      case default
        CoordsInCell = ovkHexahedronIsoInverseLinear(VertexCoords, Coords, Success=Success)
      end select
    end select

    if (OVK_VERBOSE) then
      if (.not. Success) then
        write (ERROR_UNIT, '(6a)') "WARNING: Failed to compute isoparametric coordinates of ", &
          trim(CoordsToString(Coords)), " in cell ", trim(TupleToString(Cell)), &
          " of grid ", trim(IntToString(Grid%properties%id))
      end if
    end if

  end function ovkCoordsInCell

  function ovkGridResolution(Grid, Cell, CoordsInCell) result(Resolution)

    type(ovk_grid), intent(in) :: Grid
    integer, dimension(Grid%cart%nd), intent(in) :: Cell
    real(rk), dimension(Grid%cart%nd), intent(in) :: CoordsInCell
    real(rk) :: Resolution

    integer, dimension(MAX_ND) :: PaddedCell
    integer :: i, j, k
    real(rk), dimension(MAX_ND,0:1) :: InterpBasis
    integer, dimension(MAX_ND) :: CellLower, CellUpper
    logical :: AwayFromEdge
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: BasisIndex

    PaddedCell(:Grid%cart%nd) = Cell
    PaddedCell(Grid%cart%nd+1:) = 1

    InterpBasis(1,:) = ovkInterpBasisLinear(CoordsInCell(1))
    InterpBasis(2,:) = ovkInterpBasisLinear(CoordsInCell(2))
    if (Grid%cart%nd == 3) then
      InterpBasis(3,:) = ovkInterpBasisLinear(CoordsInCell(3))
    else
      InterpBasis(3,:) = ovkInterpBasisLinear(0._rk)
    end if

    CellLower(:Grid%cart%nd) = Cell
    CellLower(Grid%cart%nd+1:) = 1
    CellUpper(:Grid%cart%nd) = Cell+1
    CellUpper(Grid%cart%nd+1:) = 1

    Resolution = 0._rk

    AwayFromEdge = ovkCartContains(Grid%cart, CellUpper)

    if (AwayFromEdge) then
      do k = CellLower(3), CellUpper(3)
        do j = CellLower(2), CellUpper(2)
          do i = CellLower(1), CellUpper(1)
            Point = [i,j,k]
            BasisIndex(:Grid%cart%nd) = Point(:Grid%cart%nd) - CellLower(:Grid%cart%nd)
            BasisIndex(Grid%cart%nd+1:) = 0
            Resolution = Resolution + Grid%resolution%values(i,j,k) * &
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
            BasisIndex(:Grid%cart%nd) = Point(:Grid%cart%nd) - CellLower(:Grid%cart%nd)
            BasisIndex(Grid%cart%nd+1:) = 0
            Point(:Grid%cart%nd) = ovkCartPeriodicAdjust(Grid%cart, Point)
            Resolution = Resolution + Grid%resolution%values(Point(1),Point(2),Point(3)) * &
              InterpBasis(1,BasisIndex(1)) * InterpBasis(2,BasisIndex(2)) * &
              InterpBasis(3,BasisIndex(3))
          end do
        end do
      end do
    end if

  end function ovkGridResolution

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

  subroutine ovkExportGridCoords(Grid, DimIndex, ExportCart, Coords)

    type(ovk_grid), intent(in) :: Grid
    integer, intent(in) :: DimIndex
    type(ovk_cart), intent(in) :: ExportCart
    type(ovk_field_real), intent(out) :: Coords

    integer :: i, j, k, d
    integer, dimension(MAX_ND) :: Point
    integer, dimension(MAX_ND) :: AdjustedPoint
    real(rk), dimension(Grid%cart%nd) :: PrincipalCoords
    real(rk), dimension(Grid%cart%nd) :: ExtendedCoords

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
            AdjustedPoint(:Grid%cart%nd) = ovkCartPeriodicAdjust(Grid%cart, Point)
            AdjustedPoint(Grid%cart%nd+1:) = 1
            do d = 1, Grid%cart%nd
              PrincipalCoords(d) = Grid%xyz(d)%values(AdjustedPoint(1),AdjustedPoint(2), &
                AdjustedPoint(3))
            end do
            ExtendedCoords = ovkPeriodicExtend(Grid%cart, Grid%properties%periodic_length, Point, &
              PrincipalCoords)
            Coords%values(i,j,k) = ExtendedCoords(DimIndex)
          else
            Coords%values(i,j,k) = Grid%xyz(DimIndex)%values(Point(1),Point(2),Point(3))
          end if
        end do
      end do
    end do

  end subroutine ovkExportGridCoords

  function ovk_grid_description_Default() result(Description)

    type(ovk_grid_description) :: Description

    Description = ovk_grid_description_Assigned(2)

  end function ovk_grid_description_Default

  function ovk_grid_description_Assigned(NumDims) result(Description)

    integer, intent(in) :: NumDims
    type(ovk_grid_description) :: Description

    Description%nd = NumDims
    Description%npoints(:NumDims) = 0
    Description%npoints(NumDims+1:) = 1
    Description%periodic = .false.
    Description%periodic_storage = OVK_NO_OVERLAP_PERIODIC
    Description%periodic_length = 0._rk
    Description%geometry_type = OVK_GRID_GEOMETRY_CURVILINEAR

  end function ovk_grid_description_Assigned

  function ovk_grid_properties_Default() result(Properties)

    type(ovk_grid_properties) :: Properties

    Properties = ovk_grid_properties_Assigned(2)

  end function ovk_grid_properties_Default

  function ovk_grid_properties_Assigned(NumDims) result(Properties)

    integer, intent(in) :: NumDims
    type(ovk_grid_properties) :: Properties

    Properties%id = 0
    Properties%nd = NumDims
    Properties%npoints(:NumDims) = 0
    Properties%npoints(NumDims+1:) = 1
    Properties%periodic = .false.
    Properties%periodic_storage = OVK_NO_OVERLAP_PERIODIC
    Properties%periodic_length = 0._rk
    Properties%geometry_type = OVK_GRID_GEOMETRY_CURVILINEAR
    Properties%verbose = .false.
    Properties%max_edge_dist = 1

  end function ovk_grid_properties_Assigned

  subroutine ovkGetGridPropertyID(Properties, ID)

    type(ovk_grid_properties), intent(in) :: Properties
    integer, intent(out) :: ID

    ID = Properties%id

  end subroutine ovkGetGridPropertyID

  subroutine ovkGetGridPropertyDimension(Properties, NumDims)

    type(ovk_grid_properties), intent(in) :: Properties
    integer, intent(out) :: NumDims

    NumDims = Properties%nd

  end subroutine ovkGetGridPropertyDimension

  subroutine ovkGetGridPropertySize(Properties, NumPoints)

    type(ovk_grid_properties), intent(in) :: Properties
    integer, dimension(Properties%nd), intent(out) :: NumPoints

    NumPoints = Properties%npoints(:Properties%nd)

  end subroutine ovkGetGridPropertySize

  subroutine ovkGetGridPropertyPeriodicity(Properties, Periodic)

    type(ovk_grid_properties), intent(in) :: Properties
    logical, dimension(Properties%nd), intent(out) :: Periodic

    Periodic = Properties%periodic(:Properties%nd)

  end subroutine ovkGetGridPropertyPeriodicity

  subroutine ovkGetGridPropertyPeriodicStorage(Properties, PeriodicStorage)

    type(ovk_grid_properties), intent(in) :: Properties
    integer, intent(out) :: PeriodicStorage

    PeriodicStorage = Properties%periodic_storage

  end subroutine ovkGetGridPropertyPeriodicStorage

  subroutine ovkGetGridPropertyPeriodicLength(Properties, PeriodicLength)

    type(ovk_grid_properties), intent(in) :: Properties
    real(rk), dimension(Properties%nd), intent(out) :: PeriodicLength

    PeriodicLength = Properties%periodic_length(:Properties%nd)

  end subroutine ovkGetGridPropertyPeriodicLength

  subroutine ovkGetGridPropertyGeometryType(Properties, GeometryType)

    type(ovk_grid_properties), intent(in) :: Properties
    integer, intent(out) :: GeometryType

    GeometryType = Properties%geometry_type

  end subroutine ovkGetGridPropertyGeometryType

  subroutine ovkGetGridPropertyVerbose(Properties, Verbose)

    type(ovk_grid_properties), intent(in) :: Properties
    logical, intent(out) :: Verbose

    Verbose = Properties%verbose

  end subroutine ovkGetGridPropertyVerbose

  subroutine ovkSetGridPropertyVerbose(Properties, Verbose)

    type(ovk_grid_properties), intent(inout) :: Properties
    logical, intent(in) :: Verbose

    Properties%verbose = Verbose

  end subroutine ovkSetGridPropertyVerbose

  subroutine ovkGetGridPropertyMaxEdgeDistance(Properties, MaxEdgeDist)

    type(ovk_grid_properties), intent(in) :: Properties
    integer, intent(out) :: MaxEdgeDist

    MaxEdgeDist = Properties%max_edge_dist

  end subroutine ovkGetGridPropertyMaxEdgeDistance

  subroutine ovkSetGridPropertyMaxEdgeDistance(Properties, MaxEdgeDist)

    type(ovk_grid_properties), intent(inout) :: Properties
    integer, intent(in) :: MaxEdgeDist

    Properties%max_edge_dist = MaxEdgeDist

  end subroutine ovkSetGridPropertyMaxEdgeDistance

  function t_grid_editor_Default() result(Editor)

    type(t_grid_editor) :: Editor

    Editor%properties_ref_count = 0
    Editor%grid_mask_ref_count = 0
    Editor%boundary_mask_ref_count = 0
    Editor%internal_boundary_mask_ref_count = 0

  end function t_grid_editor_Default

  function t_grid_editor_Assigned(NumDims) result(Editor)

    integer :: NumDims
    type(t_grid_editor) :: Editor

    Editor%properties_ref_count = 0
    allocate(Editor%xyz_ref_count(NumDims))
    Editor%xyz_ref_count = 0
    Editor%grid_mask_ref_count = 0
    Editor%boundary_mask_ref_count = 0
    Editor%internal_boundary_mask_ref_count = 0

  end function t_grid_editor_Assigned

  function ovk_grid_event_flags_Default() result(EventFlags)

    type(ovk_grid_event_flags) :: EventFlags

    EventFlags = ovk_grid_event_flags_Assigned(2)

  end function ovk_grid_event_flags_Default

  function ovk_grid_event_flags_Assigned(NumDims) result(EventFlags)

    integer, intent(in) :: NumDims
    type(ovk_grid_event_flags) :: EventFlags

    call ovkResetGridEventFlags(EventFlags)

  end function ovk_grid_event_flags_Assigned

  subroutine ovkResetGridEventFlags(EventFlags)

    type(ovk_grid_event_flags), intent(inout) :: EventFlags

    EventFlags%modified_xyz = .false.
    EventFlags%modified_grid_mask = .false.
    EventFlags%modified_boundary_mask = .false.
    EventFlags%modified_internal_boundary_mask = .false.

  end subroutine ovkResetGridEventFlags

end module ovkGrid
