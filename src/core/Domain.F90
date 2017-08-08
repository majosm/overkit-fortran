! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkDomain

  use ovkGlobal
  use ovkGrid
  implicit none

  private

  ! API
  public :: ovk_domain
  public :: ovk_domain_
  public :: ovk_domain_properties
  public :: ovk_domain_properties_
  public :: ovk_domain_event_flags
  public :: ovk_domain_event_flags_
  public :: ovkCreateDomain
  public :: ovkDestroyDomain
  public :: ovkGetDomainProperties
  public :: ovkEditDomainProperties
  public :: ovkReleaseDomainProperties
  public :: ovkCreateDomainGrid
  public :: ovkDestroyDomainGrid
  public :: ovkGetDomainGrid
  public :: ovkEditDomainGrid
  public :: ovkReleaseDomainGrid
  public :: ovkGetDomainPropertyDimension
  public :: ovkGetDomainPropertyGridCount
  public :: ovkGetDomainPropertyVerbose
  public :: ovkSetDomainPropertyVerbose
  public :: ovkResetDomainEventFlags

  type ovk_domain_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: nd
    integer :: ngrids
    ! Read/write
    logical :: verbose
  end type ovk_domain_properties

  type t_domain_editor
    integer :: properties_ref_count
    integer, dimension(:), allocatable :: grid_ref_count
  end type t_domain_editor

  type ovk_domain_event_flags
    type(t_noconstruct) :: noconstruct
    logical, dimension(:), allocatable :: modified_cart
    logical, dimension(:), allocatable :: modified_xyz
    logical, dimension(:), allocatable :: modified_grid_mask
    logical, dimension(:), allocatable :: modified_boundary_mask
    logical, dimension(:), allocatable :: modified_internal_boundary_mask
  end type ovk_domain_event_flags

  type ovk_domain
    type(t_noconstruct) :: noconstruct
    type(ovk_domain_properties), pointer :: properties
    type(ovk_domain_properties) :: prev_properties
    type(ovk_grid), dimension(:), pointer :: grids
    type(ovk_grid_event_flags), dimension(:), allocatable :: grid_event_flags
    type(t_domain_editor) :: editor
    type(ovk_domain_event_flags), pointer :: event_flags
  end type ovk_domain

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_domain_
    module procedure ovk_domain_Default
  end interface ovk_domain_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_domain_properties_
    module procedure ovk_domain_properties_Default
    module procedure ovk_domain_properties_Allocated
  end interface ovk_domain_properties_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_domain_editor_
    module procedure t_domain_editor_Default
    module procedure t_domain_editor_Allocated
  end interface t_domain_editor_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_domain_event_flags_
    module procedure ovk_domain_event_flags_Default
    module procedure ovk_domain_event_flags_Allocated
  end interface ovk_domain_event_flags_

contains

  function ovk_domain_Default() result(Domain)

    type(ovk_domain) :: Domain

    nullify(Domain%properties)
    Domain%prev_properties = ovk_domain_properties_()
    nullify(Domain%grids)
    Domain%editor = t_domain_editor_()
    nullify(Domain%event_flags)

  end function ovk_domain_Default

  subroutine ovkCreateDomain(Domain, NumDims, NumGrids, Verbose, EventFlags)

    type(ovk_domain), intent(out) :: Domain
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    logical, intent(in), optional :: Verbose
    type(ovk_domain_event_flags), target, intent(inout), optional :: EventFlags

    logical :: Verbose_
    integer :: m
    integer, dimension(MAX_ND) :: ZeroPoints

    if (present(Verbose)) then
      Verbose_ = Verbose
    else
      Verbose_ = .false.
    end if

    allocate(Domain%properties)
    Domain%properties = ovk_domain_properties_(NumDims, NumGrids)
    Domain%properties%verbose = Verbose_

    Domain%prev_properties = Domain%properties

    ZeroPoints(:NumDims) = 0
    ZeroPoints(NumDims+1:) = 1

    allocate(Domain%grids(NumGrids))
    allocate(Domain%grid_event_flags(NumGrids))
    do m = 1, NumGrids
      call ovkCreateGrid(Domain%grids(m), m, NumDims, ZeroPoints, Verbose=Verbose_, &
        EventFlags=Domain%grid_event_flags(m))
    end do

    Domain%editor = t_domain_editor_(NumGrids)

    if (present(EventFlags)) then
      Domain%event_flags => EventFlags
      Domain%event_flags = ovk_domain_event_flags_(NumDims, NumGrids)
    else
      nullify(Domain%event_flags)
    end if

  end subroutine ovkCreateDomain

  subroutine ovkDestroyDomain(Domain)

    type(ovk_domain), intent(inout) :: Domain

    integer :: m

    if (associated(Domain%properties)) deallocate(Domain%properties)
    Domain%prev_properties = ovk_domain_properties_()

    if (associated(Domain%grids)) then
      do m = 1, size(Domain%grids)
        call ovkDestroyGrid(Domain%grids(m))
      end do
      deallocate(Domain%grids)
      deallocate(Domain%grid_event_flags)
    end if

    Domain%editor = t_domain_editor_()

  end subroutine ovkDestroyDomain

  subroutine ovkGetDomainProperties(Domain, Properties)

    type(ovk_domain), intent(in) :: Domain
    type(ovk_domain_properties), pointer, intent(out) :: Properties

    Properties => Domain%properties

  end subroutine ovkGetDomainProperties

  subroutine ovkEditDomainProperties(Domain, Properties)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_domain_properties), pointer, intent(out) :: Properties

    integer :: m
    logical :: EditingGrids
    logical :: Success

    EditingGrids = .false.
    do m = 1, Domain%properties%ngrids
      EditingGrids = EditingGrids .or. Domain%editor%grid_ref_count(m) > 0
    end do

    Success = .not. EditingGrids

    if (Success) then

      if (Domain%editor%properties_ref_count == 0) then
        Domain%prev_properties = Domain%properties
      end if

      Domain%editor%properties_ref_count = Domain%editor%properties_ref_count + 1

      Properties => Domain%properties

    else

      if (OVK_DEBUG) then
        if (EditingGrids) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit domain properties while editing grids."
        end if
        stop 1
      end if

      nullify(Properties)

    end if

  end subroutine ovkEditDomainProperties

  subroutine ovkReleaseDomainProperties(Domain, Properties)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_domain_properties), pointer, intent(inout) :: Properties

    if (associated(Properties, Domain%properties)) then

      nullify(Properties)

      if (Domain%editor%properties_ref_count > 0) then

        Domain%editor%properties_ref_count = Domain%editor%properties_ref_count - 1

        if (Domain%editor%properties_ref_count == 0) then
          if (Domain%properties%verbose .neqv. Domain%prev_properties%verbose) then
            call UpdateVerbose(Domain)
          end if
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release domain properties; invalid pointer."
        stop 1
      end if

    end if

  end subroutine ovkReleaseDomainProperties

  subroutine UpdateVerbose(Domain)

    type(ovk_domain), intent(inout) :: Domain

    integer :: m
    type(ovk_grid_properties), pointer :: GridProperties

    do m = 1, Domain%properties%ngrids
      call ovkEditGridProperties(Domain%grids(m), GridProperties)
      call ovkSetGridPropertyVerbose(GridProperties, Domain%properties%verbose)
      call ovkReleaseGridProperties(Domain%grids(m), GridProperties)
    end do

  end subroutine UpdateVerbose

  subroutine ovkCreateDomainGrid(Domain, GridID, NumPoints, Periodic, PeriodicStorage, &
    PeriodicLength, GeometryType)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    integer, dimension(Domain%properties%nd), intent(in) :: NumPoints
    logical, dimension(Domain%properties%nd), intent(in), optional :: Periodic
    integer, intent(in), optional :: PeriodicStorage
    real(rk), dimension(Domain%properties%nd), intent(in), optional :: PeriodicLength
    integer, intent(in), optional :: GeometryType

    type(ovk_grid_properties), pointer :: GridProperties
    type(ovk_grid_properties) :: SavedProperties

    call ovkGetGridProperties(Domain%grids(GridID), GridProperties)
    SavedProperties = GridProperties

    call ovkDestroyGrid(Domain%grids(GridID))

    call ovkCreateGrid(Domain%grids(GridID), GridID, Domain%properties%nd, NumPoints, &
      Periodic=Periodic, PeriodicStorage=PeriodicStorage, PeriodicLength=PeriodicLength, &
      GeometryType=GeometryType, Verbose=SavedProperties%verbose, &
      EventFlags=Domain%grid_event_flags(GridID))

    call ovkEditGridProperties(Domain%grids(GridID), GridProperties)
    call ovkSetGridPropertyMaxEdgeDistance(GridProperties, SavedProperties%max_edge_dist)
    call ovkReleaseGridProperties(Domain%grids(GridID), GridProperties)

    if (associated(Domain%event_flags)) then
      Domain%event_flags%modified_cart(GridID) = .true.
      Domain%event_flags%modified_xyz(GridID) = .true.
      Domain%event_flags%modified_grid_mask(GridID) = .true.
      Domain%event_flags%modified_boundary_mask(GridID) = .true.
      Domain%event_flags%modified_internal_boundary_mask(GridID) = .true.
    end if

  end subroutine ovkCreateDomainGrid

  subroutine ovkDestroyDomainGrid(Domain, GridID)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID

    type(ovk_grid_properties), pointer :: GridProperties
    type(ovk_grid_properties) :: SavedProperties
    integer, dimension(MAX_ND) :: ZeroPoints

    call ovkGetGridProperties(Domain%grids(GridID), GridProperties)
    SavedProperties = GridProperties

    call ovkDestroyGrid(Domain%grids(GridID))

    ZeroPoints(:Domain%properties%nd) = 0
    ZeroPoints(Domain%properties%nd+1:) = 1

    call ovkCreateGrid(Domain%grids(GridID), GridID, Domain%properties%nd, ZeroPoints, &
      Verbose=SavedProperties%verbose, EventFlags=Domain%grid_event_flags(GridID))

    call ovkEditGridProperties(Domain%grids(GridID), GridProperties)
    call ovkSetGridPropertyMaxEdgeDistance(GridProperties, SavedProperties%max_edge_dist)
    call ovkReleaseGridProperties(Domain%grids(GridID), GridProperties)

    if (associated(Domain%event_flags)) then
      Domain%event_flags%modified_cart(GridID) = .true.
      Domain%event_flags%modified_xyz(GridID) = .true.
      Domain%event_flags%modified_grid_mask(GridID) = .true.
      Domain%event_flags%modified_boundary_mask(GridID) = .true.
      Domain%event_flags%modified_internal_boundary_mask(GridID) = .true.
    end if

  end subroutine ovkDestroyDomainGrid

  subroutine ovkGetDomainGrid(Domain, GridID, Grid)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: GridID
    type(ovk_grid), pointer, intent(out) :: Grid

    Grid => Domain%grids(GridID)

  end subroutine ovkGetDomainGrid

  subroutine ovkEditDomainGrid(Domain, GridID, Grid)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    type(ovk_grid), pointer, intent(out) :: Grid

    logical :: EditingProperties
    logical :: Success

    EditingProperties = Domain%editor%properties_ref_count > 0

    Success = .not. EditingProperties

    if (Success) then

      Domain%editor%grid_ref_count(GridID) = Domain%editor%grid_ref_count(GridID) + 1

      Grid => Domain%grids(GridID)

    else

      if (OVK_DEBUG) then
        if (EditingProperties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid while editing domain properties."
        end if
        stop 1
      end if

      nullify(Grid)

    end if

  end subroutine ovkEditDomainGrid

  subroutine ovkReleaseDomainGrid(Domain, Grid)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_grid), pointer, intent(inout) :: Grid

    integer :: m
    integer :: GridID

    GridID = 0
    do m = 1, Domain%properties%ngrids
      if (associated(Grid, Domain%grids(m))) then
        GridID = m
        exit
      end if
    end do

    if (GridID /= 0) then

      nullify(Grid)

      if (Domain%editor%grid_ref_count(GridID) > 0) then

        Domain%editor%grid_ref_count(GridID) = Domain%editor%grid_ref_count(GridID) - 1

        if (Domain%editor%grid_ref_count(GridID) == 0) then
          call PropagateGridEventFlags(Domain, GridID)
          call ovkResetGridEventFlags(Domain%grid_event_flags(GridID))
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release grid; invalid pointer."
        stop 1
      end if

    end if

  end subroutine ovkReleaseDomainGrid

  subroutine PropagateGridEventFlags(Domain, GridID)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID

    if (associated(Domain%event_flags)) then
      Domain%event_flags%modified_xyz(GridID) = &
        Domain%grid_event_flags(GridID)%modified_xyz
      Domain%event_flags%modified_grid_mask(GridID) = &
        Domain%grid_event_flags(GridID)%modified_grid_mask
      Domain%event_flags%modified_boundary_mask(GridID) = &
        Domain%grid_event_flags(GridID)%modified_boundary_mask
      Domain%event_flags%modified_internal_boundary_mask(GridID) = &
        Domain%grid_event_flags(GridID)%modified_internal_boundary_mask
    end if

  end subroutine PropagateGridEventFlags

  function ovk_domain_properties_Default() result(Properties)

    type(ovk_domain_properties) :: Properties

    Properties%nd = 2
    Properties%ngrids = 0
    Properties%verbose = .false.

  end function ovk_domain_properties_Default

  function ovk_domain_properties_Allocated(NumDims, NumGrids) result(Properties)

    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    type(ovk_domain_properties) :: Properties

    Properties%nd = NumDims
    Properties%ngrids = NumGrids
    Properties%verbose = .false.

  end function ovk_domain_properties_Allocated

  subroutine ovkGetDomainPropertyDimension(Properties, NumDims)

    type(ovk_domain_properties), intent(in) :: Properties
    integer, intent(out) :: NumDims

    NumDims = Properties%nd

  end subroutine ovkGetDomainPropertyDimension

  subroutine ovkGetDomainPropertyGridCount(Properties, NumGrids)

    type(ovk_domain_properties), intent(in) :: Properties
    integer, intent(out) :: NumGrids

    NumGrids = Properties%ngrids

  end subroutine ovkGetDomainPropertyGridCount

  subroutine ovkGetDomainPropertyVerbose(Properties, Verbose)

    type(ovk_domain_properties), intent(in) :: Properties
    logical, intent(out) :: Verbose

    Verbose = Properties%verbose

  end subroutine ovkGetDomainPropertyVerbose

  subroutine ovkSetDomainPropertyVerbose(Properties, Verbose)

    type(ovk_domain_properties), intent(inout) :: Properties
    logical, intent(in) :: Verbose

    Properties%verbose = Verbose

  end subroutine ovkSetDomainPropertyVerbose

  function t_domain_editor_Default() result(Editor)

    type(t_domain_editor) :: Editor

    Editor%properties_ref_count = 0

  end function t_domain_editor_Default

  function t_domain_editor_Allocated(NumGrids) result(Editor)

    integer, intent(in) :: NumGrids
    type(t_domain_editor) :: Editor

    Editor%properties_ref_count = 0

    allocate(Editor%grid_ref_count(NumGrids))
    Editor%grid_ref_count = 0

  end function t_domain_editor_Allocated

  function ovk_domain_event_flags_Default() result(EventFlags)

    type(ovk_domain_event_flags) :: EventFlags

  end function ovk_domain_event_flags_Default

  function ovk_domain_event_flags_Allocated(NumDims, NumGrids) result(EventFlags)

    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    type(ovk_domain_event_flags) :: EventFlags

    allocate(EventFlags%modified_cart(NumGrids))
    allocate(EventFlags%modified_xyz(NumGrids))
    allocate(EventFlags%modified_grid_mask(NumGrids))
    allocate(EventFlags%modified_boundary_mask(NumGrids))
    allocate(EventFlags%modified_internal_boundary_mask(NumGrids))

    call ovkResetDomainEventFlags(EventFlags)

  end function ovk_domain_event_flags_Allocated

  subroutine ovkResetDomainEventFlags(EventFlags)

    type(ovk_domain_event_flags), intent(inout) :: EventFlags

    EventFlags%modified_cart = .false.
    EventFlags%modified_xyz = .false.
    EventFlags%modified_grid_mask = .false.
    EventFlags%modified_boundary_mask = .false.
    EventFlags%modified_internal_boundary_mask = .false.

  end subroutine ovkResetDomainEventFlags

end module ovkDomain
