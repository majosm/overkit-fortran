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
  public :: ovkCreateDomain
  public :: ovkDestroyDomain
  public :: ovkDomainExists
  public :: ovkGetDomainDimension
  public :: ovkGetDomainGridCount
  public :: ovkGetDomainVerbose
  public :: ovkSetDomainVerbose
  public :: ovkCreateGrid
  public :: ovkDestroyGrid
  public :: ovkHasGrid
  public :: ovkGetGrid
  public :: ovkEditGrid
  public :: ovkReleaseGrid
  public :: ovkHasOverlap
  public :: ovkGetOverlap
  public :: ovkEditOverlap
  public :: ovkReleaseOverlap
  public :: ovkHasConnectivity
  public :: ovkGetConnectivity
  public :: ovkEditConnectivity
  public :: ovkReleaseConnectivity

  ! Internal
  public :: ovk_domain_
  public :: t_domain_edits
  public :: GetDomainEdits
  public :: ResetDomainEdits
  public :: PrintDomainSummary

  type t_domain_edits
    logical, dimension(:,:), allocatable :: overlap_dependencies
    logical, dimension(:), allocatable :: boundary_hole_dependencies
    logical, dimension(:,:), allocatable :: connectivity_dependencies
  end type t_domain_edits

  type ovk_domain
    type(t_noconstruct) :: noconstruct
    type(t_existence_flag) :: existence_flag
    type(t_logger), pointer :: logger
    integer :: nd
    integer :: ngrids
    logical :: verbose
    type(ovk_grid), dimension(:), pointer :: grid
    integer, dimension(:), allocatable :: grid_edit_ref_counts
    type(ovk_overlap), dimension(:,:), pointer :: overlap
    integer, dimension(:,:), allocatable :: overlap_edit_ref_counts
    type(ovk_connectivity), dimension(:,:), pointer :: connectivity
    integer, dimension(:,:), allocatable :: connectivity_edit_ref_counts
    type(t_domain_edits), pointer :: edits
    type(ovk_assembly_options) :: cached_assembly_options
  end type ovk_domain

contains

  function ovk_domain_() result(Domain)

    type(ovk_domain) :: Domain

    nullify(Domain%logger)
    Domain%nd = 2
    Domain%ngrids = 0
    Domain%verbose = .false.
    nullify(Domain%grid)
    nullify(Domain%overlap)
    nullify(Domain%connectivity)
    nullify(Domain%edits)
    Domain%cached_assembly_options = ovk_assembly_options_()

    call SetExists(Domain%existence_flag, .false.)

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

    allocate(Domain%logger)
    Domain%logger = t_logger_(Verbose_)

    Domain%nd = NumDims
    Domain%ngrids = NumGrids
    Domain%verbose = Verbose_

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

    allocate(Domain%edits)
    Domain%edits = t_domain_edits_(NumDims, NumGrids)

    Domain%cached_assembly_options = ovk_assembly_options_(NumDims, NumGrids)

    call SetExists(Domain%existence_flag, .true.)

  end subroutine ovkCreateDomain

  subroutine ovkDestroyDomain(Domain)

    type(ovk_domain), intent(inout) :: Domain

    integer :: m, n

    if (.not. ovkDomainExists(Domain)) return

    call SetExists(Domain%existence_flag, .false.)

    Domain%cached_assembly_options = ovk_assembly_options_()

    deallocate(Domain%edits)

    do m = 1, Domain%ngrids
      if (ovkGridExists(Domain%grid(m))) then
        call DestroyGrid(Domain%grid(m))
      end if
    end do
    deallocate(Domain%grid)

    deallocate(Domain%grid_edit_ref_counts)

    do n = 1, Domain%ngrids
      do m = 1, Domain%ngrids
        if (ovkOverlapExists(Domain%overlap(m,n))) then
          call DestroyOverlap(Domain%overlap(m,n))
        end if
      end do
    end do
    deallocate(Domain%overlap)

    deallocate(Domain%overlap_edit_ref_counts)

    do n = 1, Domain%ngrids
      do m = 1, Domain%ngrids
        if (ovkConnectivityExists(Domain%connectivity(m,n))) then
          call DestroyConnectivity(Domain%connectivity(m,n))
        end if
      end do
    end do
    deallocate(Domain%connectivity)

    deallocate(Domain%connectivity_edit_ref_counts)

    deallocate(Domain%logger)

  end subroutine ovkDestroyDomain

  function ovkDomainExists(Domain) result(Exists)

    type(ovk_domain), intent(in) :: Domain
    logical :: Exists

    Exists = CheckExists(Domain%existence_flag)

  end function ovkDomainExists

  subroutine ovkGetDomainDimension(Domain, NumDims)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(out) :: NumDims

    NumDims = Domain%nd

  end subroutine ovkGetDomainDimension

  subroutine ovkGetDomainGridCount(Domain, NumGrids)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(out) :: NumGrids

    NumGrids = Domain%ngrids

  end subroutine ovkGetDomainGridCount

  subroutine ovkGetDomainVerbose(Domain, Verbose)

    type(ovk_domain), intent(in) :: Domain
    logical, intent(out) :: Verbose

    Verbose = Domain%verbose

  end subroutine ovkGetDomainVerbose

  subroutine ovkSetDomainVerbose(Domain, Verbose)

    type(ovk_domain), intent(inout) :: Domain
    logical, intent(in) :: Verbose

    Domain%verbose = Verbose

  end subroutine ovkSetDomainVerbose

  subroutine ovkCreateGrid(Domain, GridID, NumPoints, Periodic, PeriodicStorage, PeriodicLength, &
    GeometryType)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    integer, dimension(Domain%nd), intent(in) :: NumPoints
    logical, dimension(Domain%nd), intent(in), optional :: Periodic
    integer, intent(in), optional :: PeriodicStorage
    real(rk), dimension(Domain%nd), intent(in), optional :: PeriodicLength
    integer, intent(in), optional :: GeometryType

    logical, dimension(Domain%nd) :: Periodic_
    integer :: PeriodicStorage_
    real(rk), dimension(Domain%nd) :: PeriodicLength_
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
        .not. EditingGrid(Domain) .and. &
        .not. EditingOverlap(Domain) .and. &
        .not. EditingConnectivity(Domain) .and. &
        .not. ovkHasGrid(Domain, GridID)

      if (Success) then

        Cart = ovk_cart_(Domain%nd, NumPoints, Periodic_, PeriodicStorage_)

        call CreateGrid(Domain%grid(GridID), GridID, Domain%logger, Cart, PeriodicLength_, &
          GeometryType_)

        Edits => Domain%edits

        Edits%overlap_dependencies(GridID,:) = .true.
        Edits%overlap_dependencies(:,GridID) = .true.
        Edits%boundary_hole_dependencies(:) = .true.
        Edits%connectivity_dependencies(:,:) = .true.

      else 

        if (OVK_DEBUG) then
          if (EditingGrid(Domain) .or. EditingOverlap(Domain) .or. EditingConnectivity(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot create grid while editing."
          else if (ovkHasGrid(Domain, GridID)) then
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
        .not. EditingGrid(Domain) .and. &
        .not. EditingOverlap(Domain) .and. &
        .not. EditingConnectivity(Domain) .and. &
        ovkHasGrid(Domain, GridID)

      if (Success) then

        call DestroyGrid(Domain%grid(GridID))

        Edits => Domain%edits

        Edits%overlap_dependencies(GridID,:) = .true.
        Edits%overlap_dependencies(:,GridID) = .true.
        Edits%boundary_hole_dependencies(:) = .true.
        Edits%connectivity_dependencies(:,:) = .true.

      else

        if (OVK_DEBUG) then
          if (EditingGrid(Domain) .or. EditingOverlap(Domain) .or. EditingConnectivity(Domain)) then
            write (ERROR_UNIT, '(a)') "ERROR: Cannot destroy grid while editing."
          else if (.not. ovkHasGrid(Domain, GridID)) then
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

  function ovkHasGrid(Domain, GridID) result(HasGrid)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: GridID
    logical :: HasGrid

    if (ValidID(Domain, GridID)) then
      HasGrid = ovkGridExists(Domain%grid(GridID))
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      HasGrid = .false.
    end if

  end function ovkHasGrid

  subroutine ovkGetGrid(Domain, GridID, Grid)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: GridID
    type(ovk_grid), pointer, intent(out) :: Grid

    if (ovkHasGrid(Domain, GridID)) then
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

    if (ovkHasGrid(Domain, GridID)) then

      call TryEditGrid(Domain, GridID, Success, StartEdit)

      if (Success) then
        Grid => Domain%grid(GridID)
      else
        if (OVK_DEBUG) then
          if (EditingOverlap(Domain)) then
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
    do m = 1, Domain%ngrids
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

  function ovkHasOverlap(Domain, OverlappingGridID, OverlappedGridID) result(HasOverlap)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical :: HasOverlap

    if (ValidID(Domain, OverlappingGridID) .and. ValidID(Domain, OverlappedGridID)) then
      HasOverlap = ovkOverlapExists(Domain%overlap(OverlappingGridID,OverlappedGridID))
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID(s)."
        stop 1
      end if
      HasOverlap = .false.
    end if

  end function ovkHasOverlap

  subroutine ovkGetOverlap(Domain, OverlappingGridID, OverlappedGridID, Overlap)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    type(ovk_overlap), pointer, intent(out) :: Overlap

    if (ovkHasOverlap(Domain, OverlappingGridID, OverlappedGridID)) then
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

    if (ovkHasOverlap(Domain, OverlappingGridID, OverlappedGridID)) then

      call TryEditOverlap(Domain, OverlappingGridID, OverlappedGridID, Success, StartEdit)

      if (Success) then
        Overlap => Domain%overlap(OverlappingGridID,OverlappedGridID)
      else
        if (OVK_DEBUG) then
          if (EditingGrid(Domain)) then
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
    do n = 1, Domain%ngrids
      do m = 1, Domain%ngrids
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

  function ovkHasConnectivity(Domain, DonorGridID, ReceiverGridID) result(HasConnectivity)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: DonorGridID, ReceiverGridID
    logical :: HasConnectivity

    if (ValidID(Domain, DonorGridID) .and. ValidID(Domain, ReceiverGridID)) then
      HasConnectivity = ovkConnectivityExists(Domain%connectivity(DonorGridID,ReceiverGridID))
    else
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID(s)."
        stop 1
      end if
      HasConnectivity = .false.
    end if

  end function ovkHasConnectivity

  subroutine ovkGetConnectivity(Domain, DonorGridID, ReceiverGridID, Connectivity)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in) :: DonorGridID, ReceiverGridID
    type(ovk_connectivity), pointer, intent(out) :: Connectivity

    if (ovkHasConnectivity(Domain, DonorGridID, ReceiverGridID)) then
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

    if (ovkHasConnectivity(Domain, DonorGridID, ReceiverGridID)) then

      call TryEditConnectivity(Domain, DonorGridID, ReceiverGridID, Success, StartEdit)

      if (Success) then
        Connectivity => Domain%connectivity(DonorGridID,ReceiverGridID)
      else
        if (OVK_DEBUG) then
          if (EditingOverlap(Domain)) then
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
    do n = 1, Domain%ngrids
      do m = 1, Domain%ngrids
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

  function EditingGrid(Domain, GridID) result(Editing)

    type(ovk_domain), intent(in) :: Domain
    integer, intent(in), optional :: GridID
    logical :: Editing

    integer :: m

    if (present(GridID)) then
      Editing = Domain%grid_edit_ref_counts(GridID) > 0
    else
      Editing = .false.
      do m = 1, Domain%ngrids
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
      do n = 1, Domain%ngrids
        if (Domain%overlap_edit_ref_counts(OverlappingGridID,n) > 0) then
          Editing = .true.
          exit
        end if
      end do
    else if (present(OverlappedGridID)) then
      Editing = .false.
      do m = 1, Domain%ngrids
        if (Domain%overlap_edit_ref_counts(m,OverlappedGridID) > 0) then
          Editing = .true.
          exit
        end if
      end do
    else
      Editing = .false.
      do n = 1, Domain%ngrids
        do m = 1, Domain%ngrids
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
      do n = 1, Domain%ngrids
        if (Domain%connectivity_edit_ref_counts(DonorGridID,n) > 0) then
          Editing = .true.
          exit
        end if
      end do
    else if (present(ReceiverGridID)) then
      Editing = .false.
      do m = 1, Domain%ngrids
        if (Domain%connectivity_edit_ref_counts(m,ReceiverGridID) > 0) then
          Editing = .true.
          exit
        end if
      end do
    else
      Editing = .false.
      do n = 1, Domain%ngrids
        do m = 1, Domain%ngrids
          if (Domain%connectivity_edit_ref_counts(m,n) > 0) then
            Editing = .true.
            exit
          end if
        end do
      end do
    end if

  end function EditingConnectivity

  subroutine TryEditGrid(Domain, GridID, Success, StartEdit)

    type(ovk_domain), intent(inout) :: Domain
    integer, intent(in) :: GridID
    logical, intent(out) :: Success
    logical, intent(out) :: StartEdit

    Success = &
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

    do m = 1, Domain%ngrids
      if (ovkGridExists(Domain%grid(m))) then
        call ResetGridEdits(Domain%grid(m))
      end if
    end do

!     do n = 1, Domain%ngrids
!       do m = 1, Domain%ngrids
!         if (ovkOverlapExists(Domain%overlap(m,n))) then
!           call ResetOverlapEdits(Domain%overlap(m,n))
!         end if
!       end do
!     end do

!     do n = 1, Domain%ngrids
!       do m = 1, Domain%ngrids
!         if (ovkConnectivityExists(Domain%connectivity(m,n))) then
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
      write (*, '(3a)') "* Dimension: ", trim(IntToString(Domain%nd)), "D"
      write (*, '(2a)') "* Number of grids: ", trim(IntToString(Domain%ngrids))
      TotalPoints = 0_lk
      do n = 1, Domain%ngrids
        Grid => Domain%grid(n)
        TotalPoints = TotalPoints + ovkCartCount(Grid%cart)
      end do
      TotalPointsString = LargeIntToString(TotalPoints)
      write (*, '(2a)') "* Total number of grid points: ", trim(TotalPointsString)
      do n = 1, Domain%ngrids
        Grid => Domain%grid(n)
        if (.not. ovkCartIsEmpty(Grid%cart)) then
          TotalPointsString = LargeIntToString(ovkCartCount(Grid%cart))
          iSString = IntToString(Grid%cart%is(1))
          iEString = IntToString(Grid%cart%ie(1))
          jSString = IntToString(Grid%cart%is(2))
          jEString = IntToString(Grid%cart%ie(2))
          kSString = IntToString(Grid%cart%is(3))
          kEString = IntToString(Grid%cart%ie(3))
          write (*, '(3a)', advance="no") "* Grid ", trim(IntToString(Grid%id)), ": "
          write (*, '(2a)', advance="no") trim(TotalPointsString), " points "
          select case (Domain%nd)
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

    ValidID = GridID >= 1 .and. GridID <= Domain%ngrids

  end function ValidID

end module ovkDomain
