! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkDomain

  use ovkAssemblyOptions
  use ovkCart
  use ovkConnectivity
  use ovkField
  use ovkGlobal
  use ovkGrid
  use ovkLogger
  use ovkOverlap
  implicit none

  private

  ! API
  public :: ovk_domain
  public :: ovk_domain_properties
  public :: ovkCreateDomain
  public :: ovkDestroyDomain
  public :: ovkGetDomainProperties
  public :: ovkEditDomainProperties
  public :: ovkReleaseDomainProperties
  public :: ovkCreateGrid
  public :: ovkDestroyGrid
  public :: ovkGridExists
  public :: ovkGetGrid
  public :: ovkEditGrid
  public :: ovkReleaseGrid
  public :: ovkOverlapExists
  public :: ovkGetOverlap
  public :: ovkEditOverlap
  public :: ovkReleaseOverlap
  public :: ovkConnectivityExists
  public :: ovkGetConnectivity
  public :: ovkEditConnectivity
  public :: ovkReleaseConnectivity
  public :: ovkGetDomainPropertyDimension
  public :: ovkGetDomainPropertyGridCount
  public :: ovkGetDomainPropertyVerbose
  public :: ovkSetDomainPropertyVerbose

  ! Internal
  public :: ovk_domain_
  public :: t_domain_edits
  public :: DomainExists
  public :: GetDomainEdits
  public :: ResetDomainEdits
  public :: PrintDomainSummary

  type ovk_domain_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: nd
    integer :: ngrids
    ! Read/write
    logical :: verbose
  end type ovk_domain_properties

  type t_domain_edits
    logical, dimension(:,:), allocatable :: overlap_dependencies
    logical, dimension(:), allocatable :: boundary_hole_dependencies
    logical, dimension(:,:), allocatable :: connectivity_dependencies
  end type t_domain_edits

  type ovk_domain
    type(t_noconstruct) :: noconstruct
    type(ovk_domain_properties), pointer :: properties
    type(ovk_domain_properties), pointer :: prev_properties
    integer :: properties_edit_ref_count
    type(t_logger), pointer :: logger
    type(t_domain_edits), pointer :: edits
    type(ovk_grid), dimension(:), pointer :: grid
    integer, dimension(:), allocatable :: grid_edit_ref_counts
    type(ovk_overlap), dimension(:,:), pointer :: overlap
    integer, dimension(:,:), allocatable :: overlap_edit_ref_counts
    type(ovk_connectivity), dimension(:,:), pointer :: connectivity
    integer, dimension(:,:), allocatable :: connectivity_edit_ref_counts
    type(ovk_assembly_options) :: cached_assembly_options
  end type ovk_domain

contains

  function ovk_domain_() result(Domain)

    type(ovk_domain) :: Domain

    nullify(Domain%properties)
    nullify(Domain%prev_properties)
    Domain%properties_edit_ref_count = 0
    nullify(Domain%logger)
    nullify(Domain%edits)
    nullify(Domain%grid)
    nullify(Domain%overlap)
    nullify(Domain%connectivity)
    Domain%cached_assembly_options = ovk_assembly_options_()

  end function ovk_domain_

  subroutine ovkCreateDomain(Domain, NumDims, NumGrids, Verbose)

    type(ovk_domain), intent(out) :: Domain
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    logical, intent(in), optional :: Verbose

    logical :: Verbose_
    integer :: m, n

    if (present(Verbose)) then
      Verbose_ = Verbose
    else
      Verbose_ = .false.
    end if

    allocate(Domain%properties)
    Domain%properties = ovk_domain_properties_(NumDims, NumGrids)
    Domain%properties%verbose = Verbose_

    nullify(Domain%prev_properties)

    Domain%properties_edit_ref_count = 0

    allocate(Domain%logger)
    Domain%logger = t_logger_(Verbose_)

    allocate(Domain%edits)
    Domain%edits = t_domain_edits_(NumDims, NumGrids)

    allocate(Domain%grid(NumGrids))
    do m = 1, NumGrids
      Domain%grid(m) = ovk_grid_()
    end do

    allocate(Domain%grid_edit_ref_counts(NumGrids))
    Domain%grid_edit_ref_counts = 0

    allocate(Domain%overlap(NumGrids,NumGrids))
    do n = 1, NumGrids
      do m = 1, NumGrids
        Domain%overlap(m,n) = ovk_overlap_()
      end do
    end do

    allocate(Domain%overlap_edit_ref_counts(NumGrids,NumGrids))
    Domain%overlap_edit_ref_counts = 0

    allocate(Domain%connectivity(NumGrids,NumGrids))
    do n = 1, NumGrids
      do m = 1, NumGrids
        Domain%connectivity(m,n) = ovk_connectivity_()
      end do
    end do

    allocate(Domain%connectivity_edit_ref_counts(NumGrids,NumGrids))
    Domain%connectivity_edit_ref_counts = 0

    Domain%cached_assembly_options = ovk_assembly_options_(NumDims, NumGrids)

  end subroutine ovkCreateDomain

  subroutine ovkDestroyDomain(Domain)

    type(ovk_domain), intent(inout) :: Domain

    integer :: m, n

    if (.not. DomainExists(Domain)) return

    do m = 1, Domain%properties%ngrids
      if (GridExists(Domain%grid(m))) then
        call DestroyGrid(Domain%grid(m))
      end if
    end do
    deallocate(Domain%grid)

    deallocate(Domain%grid_edit_ref_counts)

    do n = 1, Domain%properties%ngrids
      do m = 1, Domain%properties%ngrids
        if (OverlapExists(Domain%overlap(m,n))) then
          call DestroyOverlap(Domain%overlap(m,n))
        end if
      end do
    end do
    deallocate(Domain%overlap)

    deallocate(Domain%overlap_edit_ref_counts)

    do n = 1, Domain%properties%ngrids
      do m = 1, Domain%properties%ngrids
        if (ConnectivityExists(Domain%connectivity(m,n))) then
          call DestroyConnectivity(Domain%connectivity(m,n))
        end if
      end do
    end do
    deallocate(Domain%connectivity)

    deallocate(Domain%connectivity_edit_ref_counts)

    deallocate(Domain%edits)

    deallocate(Domain%logger)

    deallocate(Domain%properties)
    if (associated(Domain%prev_properties)) deallocate(Domain%prev_properties)

    Domain%cached_assembly_options = ovk_assembly_options_()

  end subroutine ovkDestroyDomain

  function DomainExists(Domain) result(Exists)

    type(ovk_domain), intent(in) :: Domain
    logical :: Exists

    Exists = associated(Domain%properties)

  end function DomainExists

  subroutine ovkGetDomainProperties(Domain, Properties)

    type(ovk_domain), intent(in) :: Domain
    type(ovk_domain_properties), pointer, intent(out) :: Properties

    Properties => Domain%properties

  end subroutine ovkGetDomainProperties

  subroutine ovkEditDomainProperties(Domain, Properties)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_domain_properties), pointer, intent(out) :: Properties

    logical :: Success, StartEdit

    call TryEditProperties(Domain, Success, StartEdit)

    if (Success) then
      if (StartEdit) then
        allocate(Domain%prev_properties)
        Domain%prev_properties = Domain%properties
      end if
      Properties => Domain%properties
    else
      if (OVK_DEBUG) then
        if (EditingGrid(Domain)) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit properties while editing grids."
        else if (EditingOverlap(Domain)) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit properties while editing overlap."
        else if (EditingConnectivity(Domain)) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit properties while editing connectivity."
        end if
        stop 1
      end if
      nullify(Properties)
    end if

  end subroutine ovkEditDomainProperties

  subroutine ovkReleaseDomainProperties(Domain, Properties)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_domain_properties), pointer, intent(inout) :: Properties

    integer :: NumGrids
    logical :: Success, EndEdit
    type(ovk_domain_properties), pointer :: PrevProperties

    if (associated(Properties, Domain%properties)) then

      call TryReleaseProperties(Domain, Success, EndEdit)

      if (Success) then

        if (EndEdit) then

          NumGrids = Domain%properties%ngrids

          PrevProperties => Domain%prev_properties
          nullify(Domain%prev_properties)

          if (Properties%verbose .neqv. PrevProperties%verbose) then
            Domain%logger%verbose = Domain%properties%verbose
          end if

        end if

      else

        if (OVK_DEBUG) then
          write (ERROR_UNIT, '(2a)') "ERROR: Unable to release properties; not ", &
            "currently being edited."
          stop 1
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release properties; invalid pointer."
        stop 1
      end if

    end if

    nullify(Properties)

  end subroutine ovkReleaseDomainProperties

  subroutine ovkCreateGrid(Domain, GridID, NumPoints, Periodic, PeriodicStorage, PeriodicLength, &
    GeometryType)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    integer, dimension(Domain%properties%nd), intent(in) :: NumPoints
    logical, dimension(Domain%properties%nd), intent(in), optional :: Periodic
    integer, intent(in), optional :: PeriodicStorage
    real(rk), dimension(Domain%properties%nd), intent(in), optional :: PeriodicLength
    integer, intent(in), optional :: GeometryType

    logical, dimension(Domain%properties%nd) :: Periodic_
    integer :: PeriodicStorage_
    real(rk), dimension(Domain%properties%nd) :: PeriodicLength_
    integer :: GeometryType_
    logical :: Success
    type(ovk_cart) :: Cart
    type(t_domain_edits), pointer :: Edits

    if (present(Periodic)) then
      Periodic_ = Periodic
    else
      Periodic_ = .false.
    end if

    if (present(PeriodicStorage)) then
      PeriodicStorage_ = PeriodicStorage
    else
      PeriodicStorage_ = OVK_NO_OVERLAP_PERIODIC
    end if

    if (present(PeriodicLength)) then
      PeriodicLength_ = PeriodicLength
    else
      PeriodicLength_ = 0._rk
    end if

    if (present(GeometryType)) then
      GeometryType_ = GeometryType
    else
      GeometryType_ = OVK_GRID_GEOMETRY_CURVILINEAR
    end if

    if (ValidID(Domain, GridID)) then

      Success = &
        .not. EditingProperties(Domain) .and. &
        .not. EditingGrid(Domain) .and. &
        .not. EditingOverlap(Domain) .and. &
        .not. EditingConnectivity(Domain) .and. &
        .not. ovkGridExists(Domain, GridID)

      if (Success) then

        Cart = ovk_cart_(Domain%properties%nd, NumPoints, Periodic_, PeriodicStorage_)

        call CreateGrid(Domain%grid(GridID), GridID, Domain%logger, Cart, PeriodicLength_, &
          GeometryType_)

        Edits => Domain%edits

        Edits%overlap_dependencies(GridID,:) = .true.
        Edits%overlap_dependencies(:,GridID) = .true.
        Edits%boundary_hole_dependencies(:) = .true.
        Edits%connectivity_dependencies(:,:) = .true.

      else 

        if (OVK_DEBUG) then
          if (EditingProperties(Domain) .or. EditingGrid(Domain) .or. EditingOverlap(Domain) .or. &
            EditingConnectivity(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot create grid while editing."
          else if (ovkGridExists(Domain, GridID)) then
            write (ERROR_UNIT, '(a)') "ERROR: Grid already exists."
          end if
          stop 1
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to create grid; invalid ID."
        stop 1
      end if

    end if

  end subroutine ovkCreateGrid

  subroutine ovkDestroyGrid(Domain, GridID)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID

    logical :: Success
    type(t_domain_edits), pointer :: Edits

    if (ValidID(Domain, GridID)) then

      Success = &
        .not. EditingProperties(Domain) .and. &
        .not. EditingGrid(Domain) .and. &
        .not. EditingOverlap(Domain) .and. &
        .not. EditingConnectivity(Domain) .and. &
        ovkGridExists(Domain, GridID)

      if (Success) then

        call DestroyGrid(Domain%grid(GridID))

        Edits => Domain%edits

        Edits%overlap_dependencies(GridID,:) = .true.
        Edits%overlap_dependencies(:,GridID) = .true.
        Edits%boundary_hole_dependencies(:) = .true.
        Edits%connectivity_dependencies(:,:) = .true.

      else

        if (OVK_DEBUG) then
          if (EditingProperties(Domain) .or. EditingGrid(Domain) .or. EditingOverlap(Domain) .or. &
            EditingConnectivity(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot destroy grid while editing."
          else if (.not. ovkGridExists(Domain, GridID)) then
            write (ERROR_UNIT, '(a)') "ERROR: Grid does not exist."
          end if
          stop 1
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to destroy grid; invalid ID."
        stop 1
      end if

    end if

  end subroutine ovkDestroyGrid

  function ovkGridExists(Domain, GridID) result(Exists)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: GridID
    logical :: Exists

    if (ValidID(Domain, GridID)) then
      Exists = GridExists(Domain%grid(GridID))
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      Exists = .false.
    end if

  end function ovkGridExists

  subroutine ovkGetGrid(Domain, GridID, Grid)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: GridID
    type(ovk_grid), pointer, intent(out) :: Grid

    if (ovkGridExists(Domain, GridID)) then
      Grid => Domain%grid(GridID)
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to retrieve grid; does not exist."
        stop 1
      end if
      nullify(Grid)
    end if

  end subroutine ovkGetGrid

  subroutine ovkEditGrid(Domain, GridID, Grid)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    type(ovk_grid), pointer, intent(out) :: Grid

    logical :: Success, StartEdit

    if (ovkGridExists(Domain, GridID)) then

      call TryEditGrid(Domain, GridID, Success, StartEdit)

      if (Success) then
        Grid => Domain%grid(GridID)
      else
        if (OVK_DEBUG) then
          if (EditingProperties(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid while editing properties."
          else if (EditingOverlap(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid while editing overlap."
          else if (EditingConnectivity(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit grid while editing connectivity."
          end if
          stop 1
        end if
        nullify(Grid)
      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to edit grid; does not exist."
        stop 1
      end if

      nullify(Grid)

    end if

  end subroutine ovkEditGrid

  subroutine ovkReleaseGrid(Domain, Grid)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_grid), pointer, intent(inout) :: Grid

    integer :: m
    integer :: GridID
    logical :: Success, EndEdit
    type(t_domain_edits), pointer :: Edits
    type(t_grid_edits), pointer :: GridEdits

    GridID = 0
    do m = 1, Domain%properties%ngrids
      if (associated(Grid, Domain%grid(m))) then
        GridID = m
        exit
      end if
    end do

    if (GridID /= 0) then

      call TryReleaseGrid(Domain, GridID, Success, EndEdit)

      if (Success) then

        if (EndEdit) then

          Edits => Domain%edits

          call GetGridEdits(Domain%grid(GridID), GridEdits)

          if (GridEdits%coords .or. GridEdits%mask) then
            Edits%overlap_dependencies(GridID,:) = .true.
            Edits%overlap_dependencies(:,GridID) = .true.
            Edits%boundary_hole_dependencies(:) = .true.
            Edits%connectivity_dependencies(:,:) = .true.
          end if

          if (GridEdits%boundary_mask .or. GridEdits%internal_boundary_mask) then
            Edits%boundary_hole_dependencies(:) = .true.
            Edits%connectivity_dependencies(:,:) = .true.
          end if

        end if

      else

        if (OVK_DEBUG) then
          write (ERROR_UNIT, '(a)') "ERROR: Unable to release grid; not currently being edited."
          stop 1
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release grid; invalid pointer."
        stop 1
      end if

    end if

    nullify(Grid)

  end subroutine ovkReleaseGrid

  function ovkOverlapExists(Domain, OverlappingGridID, OverlappedGridID) result(Exists)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical :: Exists

    if (ValidID(Domain, OverlappingGridID) .and. ValidID(Domain, OverlappedGridID)) then
      Exists = OverlapExists(Domain%overlap(OverlappingGridID,OverlappedGridID))
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID(s)."
        stop 1
      end if
      Exists = .false.
    end if

  end function ovkOverlapExists

  subroutine ovkGetOverlap(Domain, OverlappingGridID, OverlappedGridID, Overlap)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    type(ovk_overlap), pointer, intent(out) :: Overlap

    if (ovkOverlapExists(Domain, OverlappingGridID, OverlappedGridID)) then
      Overlap => Domain%overlap(OverlappingGridID,OverlappedGridID)
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to retrieve overlap; does not exist."
        stop 1
      end if
      nullify(Overlap)
    end if

  end subroutine ovkGetOverlap

  subroutine ovkEditOverlap(Domain, OverlappingGridID, OverlappedGridID, Overlap)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    type(ovk_overlap), pointer, intent(out) :: Overlap

    logical :: Success, StartEdit

    if (ovkOverlapExists(Domain, OverlappingGridID, OverlappedGridID)) then

      call TryEditOverlap(Domain, OverlappingGridID, OverlappedGridID, Success, StartEdit)

      if (Success) then
        Overlap => Domain%overlap(OverlappingGridID,OverlappedGridID)
      else
        if (OVK_DEBUG) then
          if (EditingProperties(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit overlap while editing properties."
          else if (EditingGrid(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit overlap while editing grid."
          else if (EditingConnectivity(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit overlap while editing connectivity."
          end if
          stop 1
        end if
        nullify(Overlap)
      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to edit overlap; does not exist."
        stop 1
      end if

      nullify(Overlap)

    end if

  end subroutine ovkEditOverlap

  subroutine ovkReleaseOverlap(Domain, Overlap)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_overlap), pointer, intent(inout) :: Overlap

    integer :: m, n
    integer :: OverlappingGridID, OverlappedGridID
    logical :: Success, EndEdit

    OverlappingGridID = 0
    OverlappedGridID = 0
    do n = 1, Domain%properties%ngrids
      do m = 1, Domain%properties%ngrids
        if (associated(Overlap, Domain%overlap(m,n))) then
          OverlappingGridID = m
          OverlappedGridID = n
          exit
        end if
      end do
    end do

    if (OverlappingGridID /= 0 .and. OverlappedGridID /= 0) then

      call TryReleaseOverlap(Domain, OverlappingGridID, OverlappedGridID, Success, EndEdit)

      if (Success) then

        ! Nothing to do here at the moment

      else

        if (OVK_DEBUG) then
          write (ERROR_UNIT, '(a)') "ERROR: Unable to release overlap; not currently being edited."
          stop 1
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release overlap; invalid pointer."
        stop 1
      end if

    end if

    nullify(Overlap)

  end subroutine ovkReleaseOverlap

  function ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID) result(Exists)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: DonorGridID, ReceiverGridID
    logical :: Exists

    if (ValidID(Domain, DonorGridID) .and. ValidID(Domain, ReceiverGridID)) then
      Exists = ConnectivityExists(Domain%connectivity(DonorGridID,ReceiverGridID))
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID(s)."
        stop 1
      end if
      Exists = .false.
    end if

  end function ovkConnectivityExists

  subroutine ovkGetConnectivity(Domain, DonorGridID, ReceiverGridID, Connectivity)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: DonorGridID, ReceiverGridID
    type(ovk_connectivity), pointer, intent(out) :: Connectivity

    if (ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID)) then
      Connectivity => Domain%connectivity(DonorGridID,ReceiverGridID)
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to retrieve connectivity; does not exist."
        stop 1
      end if
      nullify(Connectivity)
    end if

  end subroutine ovkGetConnectivity

  subroutine ovkEditConnectivity(Domain, DonorGridID, ReceiverGridID, Connectivity)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: DonorGridID, ReceiverGridID
    type(ovk_connectivity), pointer, intent(out) :: Connectivity

    logical :: Success, StartEdit

    if (ovkConnectivityExists(Domain, DonorGridID, ReceiverGridID)) then

      call TryEditConnectivity(Domain, DonorGridID, ReceiverGridID, Success, StartEdit)

      if (Success) then
        Connectivity => Domain%connectivity(DonorGridID,ReceiverGridID)
      else
        if (OVK_DEBUG) then
          if (EditingProperties(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity while editing properties."
          else if (EditingOverlap(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity while editing overlap."
          else if (EditingGrid(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity while editing grid."
          end if
          stop 1
        end if
        nullify(Connectivity)
      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to edit connectivity; does not exist."
        stop 1
      end if

      nullify(Connectivity)

    end if

  end subroutine ovkEditConnectivity

  subroutine ovkReleaseConnectivity(Domain, Connectivity)

    type(ovk_domain), intent(inout) :: Domain
    type(ovk_connectivity), pointer, intent(inout) :: Connectivity

    integer :: m, n
    integer :: DonorGridID, ReceiverGridID
    logical :: Success, EndEdit

    DonorGridID = 0
    ReceiverGridID = 0
    do n = 1, Domain%properties%ngrids
      do m = 1, Domain%properties%ngrids
        if (associated(Connectivity, Domain%connectivity(m,n))) then
          DonorGridID = m
          ReceiverGridID = n
          exit
        end if
      end do
    end do

    if (DonorGridID /= 0 .and. ReceiverGridID /= 0) then

      call TryReleaseConnectivity(Domain, DonorGridID, ReceiverGridID, Success, EndEdit)

      if (Success) then

        ! Nothing to do here at the moment

      else

        if (OVK_DEBUG) then
          write (ERROR_UNIT, '(a)') "ERROR: Unable to release connectivity; not currently being edited."
          stop 1
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release connectivity; invalid pointer."
        stop 1
      end if

    end if

    nullify(Connectivity)

  end subroutine ovkReleaseConnectivity

  function EditingProperties(Domain) result(Editing)

    type(ovk_domain), intent(in) :: Domain
    logical :: Editing

    Editing = Domain%properties_edit_ref_count > 0

  end function EditingProperties

  function EditingGrid(Domain, GridID) result(Editing)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in), optional :: GridID
    logical :: Editing

    integer :: m

    if (present(GridID)) then
      Editing = Domain%grid_edit_ref_counts(GridID) > 0
    else
      Editing = .false.
      do m = 1, Domain%properties%ngrids
        if (Domain%grid_edit_ref_counts(m) > 0) then
          Editing = .true.
          exit
        end if
      end do
    end if

  end function EditingGrid

  function EditingOverlap(Domain, OverlappingGridID, OverlappedGridID) result(Editing)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in), optional :: OverlappingGridID, OverlappedGridID
    logical :: Editing

    integer :: m, n

    if (present(OverlappingGridID) .and. present(OverlappedGridID)) then
      Editing = Domain%overlap_edit_ref_counts(OverlappingGridID,OverlappedGridID) > 0
    else if (present(OverlappingGridID)) then
      Editing = .false.
      do n = 1, Domain%properties%ngrids
        if (Domain%overlap_edit_ref_counts(OverlappingGridID,n) > 0) then
          Editing = .true.
          exit
        end if
      end do
    else if (present(OverlappedGridID)) then
      Editing = .false.
      do m = 1, Domain%properties%ngrids
        if (Domain%overlap_edit_ref_counts(m,OverlappedGridID) > 0) then
          Editing = .true.
          exit
        end if
      end do
    else
      Editing = .false.
      do n = 1, Domain%properties%ngrids
        do m = 1, Domain%properties%ngrids
          if (Domain%overlap_edit_ref_counts(m,n) > 0) then
            Editing = .true.
            exit
          end if
        end do
      end do
    end if

  end function EditingOverlap

  function EditingConnectivity(Domain, DonorGridID, ReceiverGridID) result(Editing)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in), optional :: DonorGridID, ReceiverGridID
    logical :: Editing

    integer :: m, n

    if (present(DonorGridID) .and. present(ReceiverGridID)) then
      Editing = Domain%connectivity_edit_ref_counts(DonorGridID,ReceiverGridID) > 0
    else if (present(DonorGridID)) then
      Editing = .false.
      do n = 1, Domain%properties%ngrids
        if (Domain%connectivity_edit_ref_counts(DonorGridID,n) > 0) then
          Editing = .true.
          exit
        end if
      end do
    else if (present(ReceiverGridID)) then
      Editing = .false.
      do m = 1, Domain%properties%ngrids
        if (Domain%connectivity_edit_ref_counts(m,ReceiverGridID) > 0) then
          Editing = .true.
          exit
        end if
      end do
    else
      Editing = .false.
      do n = 1, Domain%properties%ngrids
        do m = 1, Domain%properties%ngrids
          if (Domain%connectivity_edit_ref_counts(m,n) > 0) then
            Editing = .true.
            exit
          end if
        end do
      end do
    end if

  end function EditingConnectivity

  subroutine TryEditProperties(Domain, Success, StartEdit)

    type(ovk_domain), intent(inout) :: Domain
    logical, intent(out) :: Success
    logical, intent(out) :: StartEdit

    Success = &
      .not. EditingGrid(Domain) .and. &
      .not. EditingOverlap(Domain) .and. &
      .not. EditingConnectivity(Domain)

    if (Success) then
      StartEdit = Domain%properties_edit_ref_count == 0
      Domain%properties_edit_ref_count = Domain%properties_edit_ref_count + 1
    else
      StartEdit = .false.
    end if

  end subroutine TryEditProperties

  subroutine TryReleaseProperties(Domain, Success, EndEdit)

    type(ovk_domain), intent(inout) :: Domain
    logical, intent(out) :: Success
    logical, intent(out) :: EndEdit

    Success = EditingProperties(Domain)

    if (Success) then
      Domain%properties_edit_ref_count = Domain%properties_edit_ref_count - 1
      EndEdit = Domain%properties_edit_ref_count == 0
    else
      EndEdit = .false.
    end if

  end subroutine TryReleaseProperties

  subroutine TryEditGrid(Domain, GridID, Success, StartEdit)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    logical, intent(out) :: Success
    logical, intent(out) :: StartEdit

    Success = &
      .not. EditingProperties(Domain) .and. &
      .not. EditingOverlap(Domain) .and. &
      .not. EditingConnectivity(Domain)

    if (Success) then
      StartEdit = Domain%grid_edit_ref_counts(GridID) == 0
      Domain%grid_edit_ref_counts(GridID) = Domain%grid_edit_ref_counts(GridID) + 1
    else
      StartEdit = .false.
    end if

  end subroutine TryEditGrid

  subroutine TryReleaseGrid(Domain, GridID, Success, EndEdit)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    logical, intent(out) :: Success
    logical, intent(out) :: EndEdit

    Success = EditingGrid(Domain, GridID=GridID)

    if (Success) then
      Domain%grid_edit_ref_counts(GridID) = Domain%grid_edit_ref_counts(GridID) - 1
      EndEdit = Domain%grid_edit_ref_counts(GridID) == 0
    else
      EndEdit = .false.
    end if

  end subroutine TryReleaseGrid

  subroutine TryEditOverlap(Domain, OverlappingGridID, OverlappedGridID, Success, StartEdit)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(out) :: Success
    logical, intent(out) :: StartEdit

    Success = &
      .not. EditingProperties(Domain) .and. &
      .not. EditingGrid(Domain) .and. &
      .not. EditingConnectivity(Domain)

    if (Success) then
      StartEdit = Domain%overlap_edit_ref_counts(OverlappingGridID,OverlappedGridID) == 0
      Domain%overlap_edit_ref_counts(OverlappingGridID,OverlappedGridID) = &
        Domain%overlap_edit_ref_counts(OverlappingGridID,OverlappedGridID) + 1
    else
      StartEdit = .false.
    end if

  end subroutine TryEditOverlap

  subroutine TryReleaseOverlap(Domain, OverlappingGridID, OverlappedGridID, Success, EndEdit)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(out) :: Success
    logical, intent(out) :: EndEdit

    Success = EditingOverlap(Domain, OverlappingGridID=OverlappingGridID, &
      OverlappedGridID=OverlappedGridID)

    if (Success) then
      Domain%overlap_edit_ref_counts(OverlappingGridID,OverlappedGridID) = &
        Domain%overlap_edit_ref_counts(OverlappingGridID,OverlappedGridID) - 1
      EndEdit = Domain%overlap_edit_ref_counts(OverlappingGridID,OverlappedGridID) == 0
    else
      EndEdit = .false.
    end if

  end subroutine TryReleaseOverlap

  subroutine TryEditConnectivity(Domain, DonorGridID, ReceiverGridID, Success, StartEdit)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: DonorGridID, ReceiverGridID
    logical, intent(out) :: Success
    logical, intent(out) :: StartEdit

    Success = &
      .not. EditingProperties(Domain) .and. &
      .not. EditingGrid(Domain) .and. &
      .not. EditingOverlap(Domain)

    if (Success) then
      StartEdit = Domain%connectivity_edit_ref_counts(DonorGridID,ReceiverGridID) == 0
      Domain%connectivity_edit_ref_counts(DonorGridID,ReceiverGridID) = &
        Domain%connectivity_edit_ref_counts(DonorGridID,ReceiverGridID) + 1
    else
      StartEdit = .false.
    end if

  end subroutine TryEditConnectivity

  subroutine TryReleaseConnectivity(Domain, DonorGridID, ReceiverGridID, Success, EndEdit)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: DonorGridID, ReceiverGridID
    logical, intent(out) :: Success
    logical, intent(out) :: EndEdit

    Success = EditingConnectivity(Domain, DonorGridID=DonorGridID, ReceiverGridID=ReceiverGridID)

    if (Success) then
      Domain%connectivity_edit_ref_counts(DonorGridID,ReceiverGridID) = &
        Domain%connectivity_edit_ref_counts(DonorGridID,ReceiverGridID) - 1
      EndEdit = Domain%connectivity_edit_ref_counts(DonorGridID,ReceiverGridID) == 0
    else
      EndEdit = .false.
    end if

  end subroutine TryReleaseConnectivity

  subroutine GetDomainEdits(Domain, Edits)

    type(ovk_domain), intent(in) :: Domain
    type(t_domain_edits), pointer, intent(out) :: Edits

    Edits => Domain%edits

  end subroutine GetDomainEdits

  subroutine ResetDomainEdits(Domain)

    type(ovk_domain), intent(inout) :: Domain

    integer :: m

    Domain%edits%overlap_dependencies = .false.
    Domain%edits%boundary_hole_dependencies = .false.
    Domain%edits%connectivity_dependencies = .false.

    do m = 1, Domain%properties%ngrids
      if (GridExists(Domain%grid(m))) then
        call ResetGridEdits(Domain%grid(m))
      end if
    end do

!     do n = 1, Domain%properties%ngrids
!       do m = 1, Domain%properties%ngrids
!         if (Domain%properties%overlappable(m,n)) then
!           call ResetOverlapEdits(Domain%overlap(m,n))
!         end if
!       end do
!     end do

!     do n = 1, Domain%properties%ngrids
!       do m = 1, Domain%properties%ngrids
!         if (Domain%properties%connection_type(m,n) /= OVK_CONNECTION_NONE) then
!           call ResetConnectivityEdits(Domain%connectivity(m,n))
!         end if
!       end do
!     end do

  end subroutine ResetDomainEdits

  subroutine PrintDomainSummary(Domain)

    type(ovk_domain), intent(in) :: Domain

    integer :: n
    integer(lk) :: TotalPoints
    type(ovk_grid), pointer :: Grid
    character(len=STRING_LENGTH) :: TotalPointsString
    character(len=STRING_LENGTH) :: iSString, iEString, jSString, jEString, kSString, kEString

    if (Domain%logger%verbose) then
      write (*, '(a)') "Domain info:"
      write (*, '(3a)') "* Dimension: ", trim(IntToString(Domain%properties%nd)), "D"
      write (*, '(2a)') "* Number of grids: ", trim(IntToString(Domain%properties%ngrids))
      TotalPoints = 0_lk
      do n = 1, Domain%properties%ngrids
        Grid => Domain%grid(n)
        TotalPoints = TotalPoints + ovkCartCount(Grid%cart)
      end do
      TotalPointsString = LargeIntToString(TotalPoints)
      write (*, '(2a)') "* Total number of grid points: ", trim(TotalPointsString)
      do n = 1, Domain%properties%ngrids
        Grid => Domain%grid(n)
        if (.not. ovkCartIsEmpty(Grid%cart)) then
          TotalPointsString = LargeIntToString(ovkCartCount(Grid%cart))
          iSString = IntToString(Grid%cart%is(1))
          iEString = IntToString(Grid%cart%ie(1))
          jSString = IntToString(Grid%cart%is(2))
          jEString = IntToString(Grid%cart%ie(2))
          kSString = IntToString(Grid%cart%is(3))
          kEString = IntToString(Grid%cart%ie(3))
          write (*, '(3a)', advance="no") "* Grid ", trim(IntToString(Grid%properties%id)), ": "
          write (*, '(2a)', advance="no") trim(TotalPointsString), " points "
          select case (Domain%properties%nd)
          case (2)
            write (*, '(9a)', advance="no") "(i=", trim(iSString), ":", trim(iEString), &
              ", j=", trim(jSString), ":", trim(jEString), ")"
          case (3)
            write (*, '(13a)', advance="no") "(i=", trim(iSString), ":", trim(iEString), &
              ", j=", trim(jSString), ":", trim(jEString), ", k=", trim(kSString), ":", &
              trim(kEString), ")"
          end select
          write (*, '(a)') ""
        end if
      end do
    end if

  end subroutine PrintDomainSummary

  function ovk_domain_properties_(NumDims, NumGrids) result(Properties)

    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    type(ovk_domain_properties) :: Properties

    Properties%nd = NumDims
    Properties%ngrids = NumGrids
    Properties%verbose = .false.

  end function ovk_domain_properties_

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

  function t_domain_edits_(NumDims, NumGrids) result(Edits)

    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    type(t_domain_edits) :: Edits

    allocate(Edits%overlap_dependencies(NumGrids,NumGrids))
    Edits%overlap_dependencies = .false.

    allocate(Edits%boundary_hole_dependencies(NumGrids))
    Edits%boundary_hole_dependencies = .false.

    allocate(Edits%connectivity_dependencies(NumGrids,NumGrids))
    Edits%connectivity_dependencies = .false.

  end function t_domain_edits_

  function ValidID(Domain, GridID)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: GridID
    logical :: ValidID

    ValidID = GridID >= 1 .and. GridID <= Domain%properties%ngrids

  end function ValidID

end module ovkDomain
