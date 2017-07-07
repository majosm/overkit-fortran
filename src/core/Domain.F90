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
  public :: ovkCreateDomain
  public :: ovkDestroyDomain
  public :: ovkUpdateDomain
  public :: ovkGetDomainProperties
  public :: ovkEditDomainProperties
  public :: ovkReleaseDomainProperties
  public :: ovkCreateDomainGrid
  public :: ovkDestroyDomainGrid
  public :: ovkResetDomainGrid
  public :: ovkGetDomainGrid
  public :: ovkEditDomainGrid
  public :: ovkReleaseDomainGrid
  public :: ovkGetDomainPropertyDimension
  public :: ovkGetDomainPropertyGridCount
  public :: ovkGetDomainPropertyVerbose
  public :: ovkSetDomainPropertyVerbose
  public :: ovkGetDomainPropertyMaxEdgeDistance
  public :: ovkSetDomainPropertyMaxEdgeDistance

  type t_grid_properties
    integer :: max_edge_dist
  end type t_grid_properties

  type ovk_domain_properties
    ! Read-only
    integer :: nd
    integer :: ngrids
    ! Read/write
    logical :: verbose
    type(t_grid_properties), dimension(:), allocatable :: grids
  end type ovk_domain_properties

  type ovk_domain
    type(ovk_domain_properties), pointer :: properties
    type(ovk_domain_properties) :: prev_properties
    type(ovk_grid), dimension(:), pointer :: grids
    logical, dimension(:), allocatable :: grid_exists
    logical :: editing_properties
    logical, dimension(:), allocatable :: editing_grid
    logical, dimension(:), allocatable :: changed_grid
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
  interface t_grid_properties_
    module procedure t_grid_properties_Default
  end interface t_grid_properties_

contains

  function ovk_domain_Default() result(Domain)

    type(ovk_domain) :: Domain

    nullify(Domain%properties)
    Domain%prev_properties = ovk_domain_properties_(2)
    nullify(Domain%grids)
    Domain%editing_properties = .false.

  end function ovk_domain_Default

  subroutine ovkCreateDomain(Domain, NumDims, NumGrids, Verbose)

    type(ovk_domain), intent(out) :: Domain
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    logical, intent(in), optional :: Verbose

    logical :: Verbose_
    integer :: m

    if (present(Verbose)) then
      Verbose_ = Verbose
    else
      Verbose_ = .false.
    end if

    allocate(Domain%properties)
    Domain%properties = ovk_domain_properties_(NumDims, NumGrids)
    Domain%properties%verbose = Verbose_

    Domain%prev_properties = Domain%properties

    allocate(Domain%grids(NumGrids))
    do m = 1, NumGrids
      Domain%grids(m) = ovk_grid_()
    end do

    allocate(Domain%grid_exists(NumGrids))
    Domain%grid_exists = .false.

    Domain%editing_properties = .false.

    allocate(Domain%editing_grid(NumGrids))
    Domain%editing_grid = .false.

    allocate(Domain%changed_grid(NumGrids))
    Domain%changed_grid = .false.

  end subroutine ovkCreateDomain

  subroutine ovkDestroyDomain(Domain)

    type(ovk_domain), intent(inout) :: Domain

    integer :: m

    if (associated(Domain%properties)) deallocate(Domain%properties)
    Domain%prev_properties = ovk_domain_properties_(2)

    if (associated(Domain%grids)) then
      do m = 1, size(Domain%grids)
        if (Domain%grid_exists(m)) then
          call ovkDestroyGrid(Domain%grids(m))
        end if
      end do
      deallocate(Domain%grids)
    end if

    if (allocated(Domain%editing_grid)) deallocate(Domain%editing_grid)
    if (allocated(Domain%changed_grid)) deallocate(Domain%changed_grid)

  end subroutine ovkDestroyDomain

  subroutine ovkUpdateDomain(Domain)

    type(ovk_domain), intent(inout) :: Domain

    logical :: CannotUpdate

    CannotUpdate = &
      Domain%editing_properties .or. &
      any(Domain%editing_grid)

    if (CannotUpdate) then

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Cannot update domain; still being edited."
        stop 1
      end if

      return

    else

      Domain%changed_grid = .false.

    end if

  end subroutine ovkUpdateDomain

  subroutine ovkGetDomainProperties(Domain, Properties)

    type(ovk_domain), intent(in) :: Domain
    type(ovk_domain_properties), pointer, intent(out) :: Properties

    Properties => Domain%properties

  end subroutine ovkGetDomainProperties

  subroutine ovkEditDomainProperties(Domain, Properties)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_domain_properties), pointer, intent(out) :: Properties

    logical :: CannotEdit

    CannotEdit = &
      any(Domain%editing_grid)

    if (CannotEdit) then

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Cannot edit domain properties while editing grids."
        stop 1
      end if

      nullify(Properties)

    else

      if (.not. Domain%editing_properties) then
        Domain%prev_properties = Domain%properties
      end if

      Properties => Domain%properties
      Domain%editing_properties = .true.

    end if

  end subroutine ovkEditDomainProperties

  subroutine ovkReleaseDomainProperties(Domain, Properties)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_domain_properties), pointer, intent(inout) :: Properties

    integer :: m
    integer :: NumGrids
    type(ovk_grid_properties), pointer :: GridProperties

    if (.not. associated(Properties, Domain%properties)) return
    if (.not. Domain%editing_properties) then
      nullify(Properties)
      return
    end if

    NumGrids = Domain%properties%ngrids

    ! Propagate property changes to grids
    do m = 1, NumGrids

      if (Domain%grid_exists(m)) then

        call ovkEditGridProperties(Domain%grids(m), GridProperties)

        ! Verbose
        if (Domain%properties%verbose .neqv. Domain%prev_properties%verbose) then
          call ovkSetGridPropertyVerbose(GridProperties, Domain%properties%verbose)
        end if

        ! Maximum edge distance
        if (Domain%properties%grids(m)%max_edge_dist /= &
          Domain%prev_properties%grids(m)%max_edge_dist) then
          call ovkSetGridPropertyMaxEdgeDistance(GridProperties, &
            Domain%properties%grids(m)%max_edge_dist)
        end if

        call ovkReleaseGridProperties(Domain%grids(m), GridProperties)

      end if

    end do

    Domain%prev_properties = Domain%properties

    nullify(Properties)
    Domain%editing_properties = .false.

  end subroutine ovkReleaseDomainProperties

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

    call ovkCreateGrid(Domain%grids(GridID), GridID, Domain%properties%nd, NumPoints, &
      Periodic=Periodic, PeriodicStorage=PeriodicStorage, PeriodicLength=PeriodicLength, &
      GeometryType=GeometryType, Verbose=Domain%properties%verbose)

    call ovkEditGridProperties(Domain%grids(GridID), GridProperties)
    call ovkSetGridPropertyMaxEdgeDistance(GridProperties, &
      Domain%properties%grids(GridID)%max_edge_dist)
    call ovkReleaseGridProperties(Domain%grids(GridID), GridProperties)

    Domain%grid_exists(GridID) = .true.
    Domain%changed_grid(GridID) = .true.

  end subroutine ovkCreateDomainGrid

  subroutine ovkDestroyDomainGrid(Domain, GridID)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID

    call ovkDestroyGrid(Domain%grids(GridID))

    Domain%grid_exists(GridID) = .false.
    Domain%changed_grid(GridID) = .true.

  end subroutine ovkDestroyDomainGrid

  subroutine ovkResetDomainGrid(Domain, GridID)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID

    call ovkResetGrid(Domain%grids(GridID))
    Domain%changed_grid(GridID) = .true.

  end subroutine ovkResetDomainGrid

  subroutine ovkGetDomainGrid(Domain, GridID, Grid)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: GridID
    type(ovk_grid), pointer, intent(out) :: Grid

    if (.not. Domain%grid_exists(GridID)) then

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Grid does not exist."
        stop 1
      end if

      nullify(Grid)

    else

      Grid => Domain%grids(GridID)

    end if

  end subroutine ovkGetDomainGrid

  subroutine ovkEditDomainGrid(Domain, GridID, Grid)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    type(ovk_grid), pointer, intent(out) :: Grid

    logical :: CannotEdit

    CannotEdit = &
      .not. Domain%grid_exists(GridID) .or. &
      Domain%editing_properties

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Domain%editing_properties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid while editing domain properties."
        else
          write (ERROR_UNIT, '(a)') "ERROR: Grid does not exist."
        end if
        stop 1
      end if

      nullify(Grid)

    else

      Grid => Domain%grids(GridID)
      Domain%editing_grid(GridID) = .true.

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

    if (GridID == 0) return
    if (.not. Domain%editing_grid(GridID)) then
      nullify(Grid)
      return
    end if

    call ovkUpdateGrid(Domain%grids(GridID))

    nullify(Grid)
    Domain%editing_grid(GridID) = .false.
    Domain%changed_grid(GridID) = .true.

  end subroutine ovkReleaseDomainGrid

  function ovk_domain_properties_Default(NumDims) result(Properties)

    integer, intent(in) :: NumDims
    type(ovk_domain_properties) :: Properties

    Properties%nd = NumDims
    Properties%ngrids = 0
    Properties%verbose = .false.

  end function ovk_domain_properties_Default

  function ovk_domain_properties_Allocated(NumDims, NumGrids) result(Properties)

    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    type(ovk_domain_properties) :: Properties

    integer :: m

    Properties%nd = NumDims
    Properties%ngrids = NumGrids
    Properties%verbose = .false.
    allocate(Properties%grids(NumGrids))
    do m = 1, NumGrids
      Properties%grids(m) = t_grid_properties_(NumDims)
    end do

  end function ovk_domain_properties_Allocated

  function t_grid_properties_Default(NumDims) result(Properties)

    integer, intent(in) :: NumDims
    type(t_grid_properties) :: Properties

    Properties%max_edge_dist = 1

  end function t_grid_properties_Default

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

  subroutine ovkGetDomainPropertyMaxEdgeDistance(Properties, GridID, MaxEdgeDist)

    type(ovk_domain_properties), intent(in) :: Properties
    integer, intent(in) :: GridID
    integer, intent(out) :: MaxEdgeDist

    MaxEdgeDist = Properties%grids(GridID)%max_edge_dist

  end subroutine ovkGetDomainPropertyMaxEdgeDistance

  subroutine ovkSetDomainPropertyMaxEdgeDistance(Properties, GridID, MaxEdgeDist)

    type(ovk_domain_properties), intent(inout) :: Properties
    integer, intent(in) :: GridID
    integer, intent(in) :: MaxEdgeDist

    Properties%grids(GridID)%max_edge_dist = MaxEdgeDist

  end subroutine ovkSetDomainPropertyMaxEdgeDistance

end module ovkDomain
