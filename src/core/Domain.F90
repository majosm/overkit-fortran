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
  public :: ovkGetDomainGrid
  public :: ovkEditDomainGrid
  public :: ovkReleaseDomainGrid
  public :: ovkGetDomainPropertyDimension
  public :: ovkGetDomainPropertyGridCount
  public :: ovkGetDomainPropertyVerbose
  public :: ovkSetDomainPropertyVerbose

  type ovk_domain_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: nd
    integer :: ngrids
    ! Read/write
    logical :: verbose
  end type ovk_domain_properties

  type ovk_domain
    type(t_noconstruct) :: noconstruct
    type(ovk_domain_properties), pointer :: properties
    type(ovk_domain_properties) :: prev_properties
    type(ovk_grid), dimension(:), pointer :: grids
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
    do m = 1, NumGrids
      call ovkCreateGrid(Domain%grids(m), m, NumDims, ZeroPoints, Verbose=Verbose_)
    end do

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
        call ovkDestroyGrid(Domain%grids(m))
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

      call ovkEditGridProperties(Domain%grids(m), GridProperties)

      ! Verbose
      if (Domain%properties%verbose .neqv. Domain%prev_properties%verbose) then
        call ovkSetGridPropertyVerbose(GridProperties, Domain%properties%verbose)
      end if

      call ovkReleaseGridProperties(Domain%grids(m), GridProperties)

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
    type(ovk_grid_properties) :: SavedProperties

    call ovkGetGridProperties(Domain%grids(GridID), GridProperties)
    SavedProperties = GridProperties

    call ovkDestroyGrid(Domain%grids(GridID))

    call ovkCreateGrid(Domain%grids(GridID), GridID, Domain%properties%nd, NumPoints, &
      Periodic=Periodic, PeriodicStorage=PeriodicStorage, PeriodicLength=PeriodicLength, &
      GeometryType=GeometryType, Verbose=SavedProperties%verbose)

    call ovkEditGridProperties(Domain%grids(GridID), GridProperties)
    call ovkSetGridPropertyMaxEdgeDistance(GridProperties, SavedProperties%max_edge_dist)
    call ovkReleaseGridProperties(Domain%grids(GridID), GridProperties)

    Domain%changed_grid(GridID) = .true.

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
      Verbose=SavedProperties%verbose)

    call ovkEditGridProperties(Domain%grids(GridID), GridProperties)
    call ovkSetGridPropertyMaxEdgeDistance(GridProperties, SavedProperties%max_edge_dist)
    call ovkReleaseGridProperties(Domain%grids(GridID), GridProperties)

    Domain%changed_grid(GridID) = .true.

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

    logical :: CannotEdit

    CannotEdit = Domain%editing_properties

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Domain%editing_properties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid while editing domain properties."
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

end module ovkDomain
