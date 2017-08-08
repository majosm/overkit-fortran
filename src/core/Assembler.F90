! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkAssembler

  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDomain
  use ovkField
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkInterp
  implicit none

  private

  ! API
  public :: ovk_assembler
  public :: ovk_assembler_
  public :: ovk_assembler_properties
  public :: ovk_assembler_properties_
  public :: ovkCreateAssembler
  public :: ovkDestroyAssembler
  public :: ovkGetAssemblerProperties
  public :: ovkEditAssemblerProperties
  public :: ovkReleaseAssemblerProperties
  public :: ovkGetAssemblerDomain
  public :: ovkEditAssemblerDomain
  public :: ovkReleaseAssemblerDomain
!   public :: ovkGetAssemblerOverlap
!   public :: ovkEditAssemblerOverlap
!   public :: ovkReleaseAssemblerOverlap
  public :: ovkGetAssemblerConnectivity
  public :: ovkEditAssemblerConnectivity
  public :: ovkReleaseAssemblerConnectivity
  public :: ovkGetAssemblerDebugField
  public :: ovkGetAssemblerPropertyDimension
  public :: ovkGetAssemblerPropertyGridCount
  public :: ovkGetAssemblerPropertyVerbose
  public :: ovkSetAssemblerPropertyVerbose
  public :: ovkGetAssemblerPropertyManualPadding
  public :: ovkSetAssemblerPropertyManualPadding
  public :: ovkGetAssemblerPropertyInferBoundaries
  public :: ovkSetAssemblerPropertyInferBoundaries
  public :: ovkGetAssemblerPropertyOverlappable
  public :: ovkSetAssemblerPropertyOverlappable
  public :: ovkGetAssemblerPropertyOverlapTolerance
  public :: ovkSetAssemblerPropertyOverlapTolerance
  public :: ovkGetAssemblerPropertyBoundaryHoleCutting
  public :: ovkSetAssemblerPropertyBoundaryHoleCutting
  public :: ovkGetAssemblerPropertyOverlapHoleCutting
  public :: ovkSetAssemblerPropertyOverlapHoleCutting
  public :: ovkGetAssemblerPropertyConnectionType
  public :: ovkSetAssemblerPropertyConnectionType
  public :: ovkGetAssemblerPropertyDisjointConnection
  public :: ovkSetAssemblerPropertyDisjointConnection
  public :: ovkGetAssemblerPropertyInterpScheme
  public :: ovkSetAssemblerPropertyInterpScheme
  public :: ovkGetAssemblerPropertyFringeSize
  public :: ovkSetAssemblerPropertyFringeSize
  public :: ovkGetAssemblerPropertyFringePadding
  public :: ovkSetAssemblerPropertyFringePadding

  type ovk_assembler_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: nd
    integer :: ngrids
    ! Read/write
    logical :: verbose
    logical :: manual_padding
    logical, dimension(:), allocatable :: infer_boundaries
    logical, dimension(:,:), allocatable :: overlappable
    real(rk), dimension(:,:), allocatable :: overlap_tolerance
    logical, dimension(:,:), allocatable :: boundary_hole_cutting
    logical, dimension(:,:), allocatable :: overlap_hole_cutting
    integer, dimension(:,:), allocatable :: connection_type
    logical, dimension(:,:), allocatable :: disjoint_connection
    integer, dimension(:,:), allocatable :: interp_scheme
    integer, dimension(:,:), allocatable :: fringe_size
    integer, dimension(:,:), allocatable :: fringe_padding
  end type ovk_assembler_properties

  type ovk_assembler
    type(t_noconstruct) :: noconstruct
    type(ovk_assembler_properties), pointer :: properties
    type(ovk_assembler_properties) :: prev_properties
    type(ovk_domain), pointer :: domain
!     type(ovk_overlap), pointer :: overlap
    type(ovk_connectivity), pointer :: connectivity
    type(ovk_field_int), dimension(:), pointer :: debug_fields
    logical :: editing_properties
    logical :: editing_domain
!     logical :: editing_overlap
    logical :: editing_connectivity
  end type ovk_assembler

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_assembler_
    module procedure ovk_assembler_Default
  end interface ovk_assembler_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_assembler_properties_
    module procedure ovk_assembler_properties_Default
    module procedure ovk_assembler_properties_Allocated
  end interface ovk_assembler_properties_

contains

  function ovk_assembler_Default() result(Assembler)

    type(ovk_assembler) :: Assembler

    nullify(Assembler%properties)
    Assembler%prev_properties = ovk_assembler_properties_()
    nullify(Assembler%domain)
!     nullify(Assembler%overlap)
    nullify(Assembler%connectivity)
    nullify(Assembler%debug_fields)
    Assembler%editing_properties = .false.
    Assembler%editing_domain = .false.
!     Assembler%editing_overlap = .false.
    Assembler%editing_connectivity = .false.

  end function ovk_assembler_Default

  subroutine ovkCreateAssembler(Assembler, NumDims, NumGrids, Verbose)

    type(ovk_assembler), intent(out) :: Assembler
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

    allocate(Assembler%properties)
    Assembler%properties = ovk_assembler_properties_(NumDims, NumGrids)
    Assembler%properties%verbose = Verbose_

    Assembler%prev_properties = Assembler%properties

    allocate(Assembler%domain)
    call ovkCreateDomain(Assembler%domain, NumDims, NumGrids, Verbose=Verbose_)

!     allocate(Assembler%overlap)
!     call ovkCreateOverlap(Assembler%overlap, NumDims, NumGrids, Verbose=Verbose_)

    allocate(Assembler%connectivity)
    call ovkCreateConnectivity(Assembler%connectivity, NumDims, NumGrids, Verbose=Verbose_)

    if (OVK_DEBUG) then
      allocate(Assembler%debug_fields(NumGrids))
      do m = 1, NumGrids
        Assembler%debug_fields(m) = ovk_field_int_()
      end do
    else
      nullify(Assembler%debug_fields)
    end if

    Assembler%editing_properties = .false.
    Assembler%editing_domain = .false.
!     Assembler%editing_overlap = .false.
    Assembler%editing_connectivity = .false.

  end subroutine ovkCreateAssembler

  subroutine ovkDestroyAssembler(Assembler)

    type(ovk_assembler), intent(inout) :: Assembler

    if (associated(Assembler%properties)) deallocate(Assembler%properties)
    Assembler%prev_properties = ovk_assembler_properties_()

    if (associated(Assembler%domain)) then
      call ovkDestroyDomain(Assembler%domain)
      deallocate(Assembler%domain)
    end if

!     if (associated(Assembler%overlap)) then
!       call ovkDestroyOverlap(Assembler%overlap)
!       deallocate(Assembler%overlap)
!     end if

    if (associated(Assembler%connectivity)) then
      call ovkDestroyConnectivity(Assembler%connectivity)
      deallocate(Assembler%connectivity)
    end if

    if (OVK_DEBUG) then
      if (associated(Assembler%debug_fields)) then
        deallocate(Assembler%debug_fields)
      end if
    end if

  end subroutine ovkDestroyAssembler

  subroutine ovkGetAssemblerProperties(Assembler, Properties)

    type(ovk_assembler), intent(in) :: Assembler
    type(ovk_assembler_properties), pointer, intent(out) :: Properties

    Properties => Assembler%properties

  end subroutine ovkGetAssemblerProperties

  subroutine ovkEditAssemblerProperties(Assembler, Properties)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_assembler_properties), pointer, intent(out) :: Properties

    logical :: CannotEdit

    CannotEdit = &
      Assembler%editing_domain .or. &
      Assembler%editing_connectivity

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Assembler%editing_domain) then
          write (ERROR_UNIT, '(2a)') "ERROR: Cannot edit assembler properties while editing ", &
            "domain."
        else
          write (ERROR_UNIT, '(2a)') "ERROR: Cannot edit assembler properties while editing ", &
            "connectivity."
        end if
        stop 1
      end if

      nullify(Properties)

    else

      if (.not. Assembler%editing_properties) then
        Assembler%prev_properties = Assembler%properties
      end if

      Properties => Assembler%properties
      Assembler%editing_properties = .true.

    end if

  end subroutine ovkEditAssemblerProperties

  subroutine ovkReleaseAssemblerProperties(Assembler, Properties)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_assembler_properties), pointer, intent(inout) :: Properties

    integer :: m, n
    integer :: NumGrids
    type(ovk_domain_properties), pointer :: DomainProperties
    type(ovk_connectivity_properties), pointer :: ConnectivityProperties
    type(ovk_grid), pointer :: Grid
    type(ovk_grid_properties), pointer :: GridProperties
    integer :: InterpScheme
    integer :: FringeSize
    integer :: MaxEdgeDist, PrevMaxEdgeDist

    if (.not. associated(Properties, Assembler%properties)) return
    if (.not. Assembler%editing_properties) then
      nullify(Properties)
      return
    end if

    NumGrids = Assembler%properties%ngrids

    if (OVK_DEBUG) then

      do m = 1, NumGrids
        if (ovkCartCount(Assembler%domain%grids(m)%cart) > 0) then
          write (ERROR_UNIT, '(2a)') "ERROR: Dynamically updating assembler properties is not ", &
            "yet supported; must set properties before grids are created."
          stop 1
        end if
      end do

      do n = 1, NumGrids
        InterpScheme = OVK_INTERP_LINEAR
        do m = 1, NumGrids
          if (Assembler%properties%connection_type(m,n) /= OVK_CONNECTION_NONE) then
            InterpScheme = Assembler%properties%interp_scheme(m,n)
            exit
          end if
        end do
        do m = 1, NumGrids
          if (Assembler%properties%connection_type(m,n) /= OVK_CONNECTION_NONE) then
            if (Assembler%properties%interp_scheme(m,n) /= InterpScheme) then
              write (ERROR_UNIT, '(2a)') "ERROR: Pairwise interpolation scheme specification is ", &
                "not currently supported; must be set uniformly for each receiver grid."
              stop 1
            end if
          end if
        end do
      end do

      do n = 1, NumGrids
        FringeSize = 0
        do m = 1, NumGrids
          if (Assembler%properties%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
            FringeSize = Assembler%properties%fringe_size(m,n)
            exit
          end if
        end do
        do m = 1, NumGrids
          if (Assembler%properties%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
            if (Assembler%properties%fringe_size(m,n) /= FringeSize) then
              write (ERROR_UNIT, '(2a)') "ERROR: Pairwise fringe size specification is ", &
                "not currently supported; must be set uniformly for each receiver grid."
              stop 1
            end if
          end if
        end do
      end do

    end if

    ! Propagate verbose flag to sub-objects
    call ovkEditDomainProperties(Assembler%domain, DomainProperties)
    call ovkSetDomainPropertyVerbose(DomainProperties, Assembler%properties%verbose)
    call ovkReleaseDomainProperties(Assembler%domain, DomainProperties)
!     call ovkEditOverlapProperties(Assembler%overlap, OverlapProperties)
!     call ovkSetOverlapPropertyVerbose(OverlapProperties, Assembler%properties%verbose)
!     call ovkReleaseOverlapProperties(Assembler%overlap, OverlapProperties)
    call ovkEditConnectivityProperties(Assembler%connectivity, ConnectivityProperties)
    call ovkSetConnectivityPropertyVerbose(ConnectivityProperties, Assembler%properties%verbose)
    call ovkReleaseConnectivityProperties(Assembler%connectivity, ConnectivityProperties)

    do n = 1, NumGrids
      ! No self-intersection support (yet)
      Assembler%properties%overlappable(n,n) = .false.
      ! No overlap implies no hole cutting, no communication
      do m = 1, NumGrids
        if (.not. Assembler%properties%overlappable(m,n)) then
          Assembler%properties%boundary_hole_cutting(m,n) = .false.
          Assembler%properties%overlap_hole_cutting(m,n) = .false.
          Assembler%properties%connection_type(m,n) = OVK_CONNECTION_NONE
        end if
      end do
    end do

    do n = 1, NumGrids
      MaxEdgeDist = GetMaxEdgeDistance(Assembler%properties, n)
      PrevMaxEdgeDist = GetMaxEdgeDistance(Assembler%prev_properties, n)
      if (MaxEdgeDist /= PrevMaxEdgeDist) then
        call ovkEditDomainGrid(Assembler%domain, n, Grid)
        call ovkEditGridProperties(Grid, GridProperties)
        call ovkSetGridPropertyMaxEdgeDistance(GridProperties, MaxEdgeDist)
        call ovkReleaseGridProperties(Grid, GridProperties)
        call ovkReleaseDomainGrid(Assembler%domain, Grid)
      end if
    end do

    Assembler%prev_properties = Assembler%properties

    nullify(Properties)
    Assembler%editing_properties = .false.

  end subroutine ovkReleaseAssemblerProperties

  subroutine ovkGetAssemblerDomain(Assembler, Domain)

    type(ovk_assembler), intent(in) :: Assembler
    type(ovk_domain), pointer, intent(out) :: Domain

    Domain => Assembler%domain

  end subroutine ovkGetAssemblerDomain

  subroutine ovkEditAssemblerDomain(Assembler, Domain)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_domain), pointer, intent(out) :: Domain

    logical :: CannotEdit

    CannotEdit = &
      Assembler%editing_properties .or. &
      Assembler%editing_connectivity

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Assembler%editing_properties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit domain while editing assembler properties."
        else
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit domain while editing connectivity."
        end if
        stop 1
      end if

      nullify(Domain)

    else

      Domain => Assembler%domain
      Assembler%editing_domain = .true.

    end if

  end subroutine ovkEditAssemblerDomain

  subroutine ovkReleaseAssemblerDomain(Assembler, Domain)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_domain), pointer, intent(inout) :: Domain

    integer :: m, n
    integer :: NumGrids
    logical, dimension(:), allocatable :: ChangedGrid
    type(ovk_bbox) :: Bounds

    if (.not. associated(Domain, Assembler%domain)) return
    if (.not. Assembler%editing_domain) then
      nullify(Domain)
      return
    end if

    NumGrids = Assembler%properties%ngrids

    allocate(ChangedGrid(NumGrids))
    ChangedGrid = Domain%changed_grid

    call ovkUpdateDomain(Assembler%domain)

!     do n = 1, NumGrids
!       if (ChangedGrid(n)) then
!         call ovkResetConnectivityDonors(Connectivity, GridID)
!         call ovkResetConnectivityReceivers(Connectivity, GridID)
!       end if
!     end do

    do n = 1, NumGrids
      do m = 1, NumGrids
        if (ChangedGrid(m) .or. ChangedGrid(n)) then
          if (Assembler%properties%overlappable(m,n)) then
            Bounds = ovkBBIntersect(Domain%grids(m)%bounds, Domain%grids(n)%bounds)
            if (ovkBBIsEmpty(Bounds)) then
              Assembler%properties%overlappable(m,n) = .false.
              Assembler%properties%boundary_hole_cutting(m,n) = .false.
              Assembler%properties%overlap_hole_cutting(m,n) = .false.
              Assembler%properties%connection_type(m,n) = OVK_CONNECTION_NONE
            end if
          end if
        end if
      end do
    end do

    if (OVK_DEBUG) then
      do m = 1, NumGrids
        if (ChangedGrid(m)) then
          Assembler%debug_fields(m) = ovk_field_int_(Domain%grids(m)%cart, 0)
        end if
      end do
    end if

    nullify(Domain)
    Assembler%editing_domain = .false.

  end subroutine ovkReleaseAssemblerDomain

  subroutine ovkGetAssemblerConnectivity(Assembler, Connectivity)

    type(ovk_assembler), intent(in) :: Assembler
    type(ovk_connectivity), pointer, intent(out) :: Connectivity

    Connectivity => Assembler%connectivity

  end subroutine ovkGetAssemblerConnectivity

  subroutine ovkEditAssemblerConnectivity(Assembler, Connectivity)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_connectivity), pointer, intent(out) :: Connectivity

    logical :: CannotEdit

    CannotEdit = &
      Assembler%editing_properties .or. &
      Assembler%editing_domain

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Assembler%editing_properties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity while editing assembler properties."
        else
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity while editing domain."
        end if
        stop 1
      end if

      nullify(Connectivity)

    else

      Connectivity => Assembler%connectivity
      Assembler%editing_connectivity = .true.

    end if

  end subroutine ovkEditAssemblerConnectivity

  subroutine ovkReleaseAssemblerConnectivity(Assembler, Connectivity)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_connectivity), pointer, intent(inout) :: Connectivity

    if (.not. associated(Connectivity, Assembler%connectivity)) return
    if (.not. Assembler%editing_connectivity) then
      nullify(Connectivity)
      return
    end if

    call ovkUpdateConnectivity(Assembler%connectivity)

    nullify(Connectivity)
    Assembler%editing_connectivity = .false.

  end subroutine ovkReleaseAssemblerConnectivity

  subroutine ovkGetAssemblerDebugField(Assembler, GridID, DebugField)

    type(ovk_assembler), intent(in) :: Assembler
    integer, intent(in) :: GridID
    type(ovk_field_int), pointer, intent(out) :: DebugField

    if (OVK_DEBUG) then
      DebugField => Assembler%debug_fields(GridID)
    else
      nullify(DebugField)
    end if

  end subroutine ovkGetAssemblerDebugField

  function ovk_assembler_properties_Default() result(Properties)

    type(ovk_assembler_properties) :: Properties

    Properties%nd = 2
    Properties%ngrids = 0
    Properties%verbose = .false.
    Properties%manual_padding = .true.

  end function ovk_assembler_properties_Default

  function ovk_assembler_properties_Allocated(NumDims, NumGrids) result(Properties)

    type(ovk_assembler_properties) :: Properties
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids

    Properties%nd = NumDims
    Properties%ngrids = NumGrids
    Properties%verbose = .false.
    Properties%manual_padding = .true.

    allocate(Properties%infer_boundaries(NumGrids))
    Properties%infer_boundaries = .false.

    allocate(Properties%overlappable(NumGrids,NumGrids))
    Properties%overlappable = .false.

    allocate(Properties%overlap_tolerance(NumGrids,NumGrids))
    Properties%overlap_tolerance = 0._rk

    allocate(Properties%boundary_hole_cutting(NumGrids,NumGrids))
    Properties%boundary_hole_cutting = .false.

    allocate(Properties%overlap_hole_cutting(NumGrids,NumGrids))
    Properties%overlap_hole_cutting = .false.

    allocate(Properties%connection_type(NumGrids,NumGrids))
    Properties%connection_type = OVK_CONNECTION_NONE

    allocate(Properties%disjoint_connection(NumGrids,NumGrids))
    Properties%disjoint_connection = .false.

    allocate(Properties%interp_scheme(NumGrids,NumGrids))
    Properties%interp_scheme = OVK_INTERP_LINEAR

    allocate(Properties%fringe_size(NumGrids,NumGrids))
    Properties%fringe_size = 0

    allocate(Properties%fringe_padding(NumGrids,NumGrids))
    Properties%fringe_padding = 0

  end function ovk_assembler_properties_Allocated

  subroutine ovkGetAssemblerPropertyDimension(Properties, NumDims)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(out) :: NumDims

    NumDims = Properties%nd

  end subroutine ovkGetAssemblerPropertyDimension

  subroutine ovkGetAssemblerPropertyGridCount(Properties, NumGrids)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(out) :: NumGrids

    NumGrids = Properties%ngrids

  end subroutine ovkGetAssemblerPropertyGridCount

  subroutine ovkGetAssemblerPropertyVerbose(Properties, Verbose)

    type(ovk_assembler_properties), intent(in) :: Properties
    logical, intent(out) :: Verbose

    Verbose = Properties%verbose

  end subroutine ovkGetAssemblerPropertyVerbose

  subroutine ovkSetAssemblerPropertyVerbose(Properties, Verbose)

    type(ovk_assembler_properties), intent(inout) :: Properties
    logical, intent(in) :: Verbose

    Properties%verbose = Verbose

  end subroutine ovkSetAssemblerPropertyVerbose

  subroutine ovkGetAssemblerPropertyManualPadding(Properties, ManualPadding)

    type(ovk_assembler_properties), intent(in) :: Properties
    logical, intent(out) :: ManualPadding

    ManualPadding = Properties%manual_padding

  end subroutine ovkGetAssemblerPropertyManualPadding

  subroutine ovkSetAssemblerPropertyManualPadding(Properties, ManualPadding)

    type(ovk_assembler_properties), intent(inout) :: Properties
    logical, intent(in) :: ManualPadding

    Properties%manual_padding = ManualPadding

  end subroutine ovkSetAssemblerPropertyManualPadding

  subroutine ovkGetAssemblerPropertyInferBoundaries(Properties, GridID, InferBoundaries)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: GridID
    logical, intent(out) :: InferBoundaries

    InferBoundaries = Properties%infer_boundaries(GridID)

  end subroutine ovkGetAssemblerPropertyInferBoundaries

  subroutine ovkSetAssemblerPropertyInferBoundaries(Properties, GridID, InferBoundaries)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: GridID
    logical, intent(in) :: InferBoundaries

    integer :: m
    integer :: ms, me

    call GridIDRange(Properties%ngrids, GridID, ms, me)

    do m = ms, me
      Properties%infer_boundaries(m) = InferBoundaries
    end do

  end subroutine ovkSetAssemblerPropertyInferBoundaries

  subroutine ovkGetAssemblerPropertyOverlappable(Properties, OverlappingGridID, OverlappedGridID, &
    Overlappable)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(out) :: Overlappable

    Overlappable = Properties%overlappable(OverlappingGridID,OverlappedGridID)

  end subroutine ovkGetAssemblerPropertyOverlappable

  subroutine ovkSetAssemblerPropertyOverlappable(Properties, OverlappingGridID, OverlappedGridID, &
    Overlappable)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(in) :: Overlappable

    integer :: m, n
    integer :: ms, me, ns, ne

    call GridIDRange(Properties%ngrids, OverlappingGridID, ms, me)
    call GridIDRange(Properties%ngrids, OverlappedGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%overlappable(m,n) = Overlappable
      end do
    end do

  end subroutine ovkSetAssemblerPropertyOverlappable

  subroutine ovkGetAssemblerPropertyOverlapTolerance(Properties, OverlappingGridID, &
    OverlappedGridID, OverlapTolerance)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    real(rk), intent(out) :: OverlapTolerance

    OverlapTolerance = Properties%overlap_tolerance(OverlappingGridID,OverlappedGridID)

  end subroutine ovkGetAssemblerPropertyOverlapTolerance

  subroutine ovkSetAssemblerPropertyOverlapTolerance(Properties, OverlappingGridID, &
    OverlappedGridID, OverlapTolerance)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    real(rk), intent(in) :: OverlapTolerance

    integer :: m, n
    integer :: ms, me, ns, ne

    if (OVK_DEBUG) then
      if (OverlapTolerance < 0._rk) then
        write (ERROR_UNIT, '(a)') "ERROR: Overlap tolerance must be nonnegative."
        stop 1
      end if
    end if

    call GridIDRange(Properties%ngrids, OverlappingGridID, ms, me)
    call GridIDRange(Properties%ngrids, OverlappedGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%overlap_tolerance(m,n) = OverlapTolerance
      end do
    end do

  end subroutine ovkSetAssemblerPropertyOverlapTolerance

  subroutine ovkGetAssemblerPropertyBoundaryHoleCutting(Properties, CuttingGridID, CutGridID, &
    BoundaryHoleCutting)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(out) :: BoundaryHoleCutting

    BoundaryHoleCutting = Properties%boundary_hole_cutting(CuttingGridID,CutGridID)

  end subroutine ovkGetAssemblerPropertyBoundaryHoleCutting

  subroutine ovkSetAssemblerPropertyBoundaryHoleCutting(Properties, CuttingGridID, CutGridID, &
    BoundaryHoleCutting)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(in) :: BoundaryHoleCutting

    integer :: m, n
    integer :: ms, me, ns, ne

    call GridIDRange(Properties%ngrids, CuttingGridID, ms, me)
    call GridIDRange(Properties%ngrids, CutGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%boundary_hole_cutting(m,n) = BoundaryHoleCutting
      end do
    end do

  end subroutine ovkSetAssemblerPropertyBoundaryHoleCutting

  subroutine ovkGetAssemblerPropertyOverlapHoleCutting(Properties, CuttingGridID, CutGridID, &
    OverlapHoleCutting)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(out) :: OverlapHoleCutting

    OverlapHoleCutting = Properties%overlap_hole_cutting(CuttingGridID,CutGridID)

  end subroutine ovkGetAssemblerPropertyOverlapHoleCutting

  subroutine ovkSetAssemblerPropertyOverlapHoleCutting(Properties, CuttingGridID, CutGridID, &
    OverlapHoleCutting)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(in) :: OverlapHoleCutting

    integer :: m, n
    integer :: ms, me, ns, ne

    call GridIDRange(Properties%ngrids, CuttingGridID, ms, me)
    call GridIDRange(Properties%ngrids, CutGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%overlap_hole_cutting(m,n) = OverlapHoleCutting
      end do
    end do

  end subroutine ovkSetAssemblerPropertyOverlapHoleCutting

  subroutine ovkGetAssemblerPropertyConnectionType(Properties, DonorGridID, ReceiverGridID, &
    ConnectionType)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: ConnectionType

    if (OVK_DEBUG) then
      if (ConnectionType == OVK_CONNECTION_FULL_GRID) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_CONNECTION_FULL_GRID is not yet supported."
        stop 1
      end if
    end if

    ConnectionType = Properties%connection_type(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerPropertyConnectionType

  subroutine ovkSetAssemblerPropertyConnectionType(Properties, DonorGridID, ReceiverGridID, &
    ConnectionType)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(in) :: ConnectionType

    logical :: ValidValue
    integer :: m, n
    integer :: ms, me, ns, ne

    if (OVK_DEBUG) then
      ValidValue = &
        ConnectionType == OVK_CONNECTION_NONE .or. &
        ConnectionType == OVK_CONNECTION_FRINGE .or. &
        ConnectionType == OVK_CONNECTION_FULL_GRID
      if (.not. ValidValue) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid connection type."
        stop 1
      end if
    end if

    call GridIDRange(Properties%ngrids, DonorGridID, ms, me)
    call GridIDRange(Properties%ngrids, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%connection_type(m,n) = ConnectionType
      end do
    end do

  end subroutine ovkSetAssemblerPropertyConnectionType

  subroutine ovkGetAssemblerPropertyDisjointConnection(Properties, DonorGridID, ReceiverGridID, &
    DisjointConnection)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    logical, intent(out) :: DisjointConnection

    DisjointConnection = Properties%disjoint_connection(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerPropertyDisjointConnection

  subroutine ovkSetAssemblerPropertyDisjointConnection(Properties, DonorGridID, ReceiverGridID, &
    DisjointConnection)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    logical, intent(in) :: DisjointConnection

    integer :: m, n
    integer :: ms, me, ns, ne

    call GridIDRange(Properties%ngrids, DonorGridID, ms, me)
    call GridIDRange(Properties%ngrids, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%disjoint_connection(m,n) = DisjointConnection
      end do
    end do

  end subroutine ovkSetAssemblerPropertyDisjointConnection

  subroutine ovkGetAssemblerPropertyInterpScheme(Properties, DonorGridID, ReceiverGridID, &
    InterpScheme)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: InterpScheme

    InterpScheme = Properties%interp_scheme(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerPropertyInterpScheme

  subroutine ovkSetAssemblerPropertyInterpScheme(Properties, DonorGridID, ReceiverGridID, &
    InterpScheme)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(in) :: InterpScheme

    logical :: ValidValue
    integer :: m, n
    integer :: ms, me, ns, ne

    if (OVK_DEBUG) then
      ValidValue = &
        InterpScheme == OVK_INTERP_LINEAR .or. &
        InterpScheme == OVK_INTERP_CUBIC
      if (.not. ValidValue) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid interpolation scheme."
        stop 1
      end if
    end if

    call GridIDRange(Properties%ngrids, DonorGridID, ms, me)
    call GridIDRange(Properties%ngrids, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%interp_scheme(m,n) = InterpScheme
      end do
    end do

  end subroutine ovkSetAssemblerPropertyInterpScheme

  subroutine ovkGetAssemblerPropertyFringeSize(Properties, DonorGridID, ReceiverGridID, FringeSize)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: FringeSize

    FringeSize = Properties%fringe_size(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerPropertyFringeSize

  subroutine ovkSetAssemblerPropertyFringeSize(Properties, DonorGridID, ReceiverGridID, FringeSize)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(in) :: FringeSize

    integer :: m, n
    integer :: ms, me, ns, ne

    if (OVK_DEBUG) then
      if (FringeSize < 0) then
        write (ERROR_UNIT, '(a)') "ERROR: Fringe size must be nonnegative."
        stop 1
      end if
    end if

    call GridIDRange(Properties%ngrids, DonorGridID, ms, me)
    call GridIDRange(Properties%ngrids, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%fringe_size(m,n) = FringeSize
      end do
    end do

  end subroutine ovkSetAssemblerPropertyFringeSize

  subroutine ovkGetAssemblerPropertyFringePadding(Properties, DonorGridID, ReceiverGridID, &
    FringePadding)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: FringePadding

    FringePadding = Properties%fringe_padding(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerPropertyFringePadding

  subroutine ovkSetAssemblerPropertyFringePadding(Properties, DonorGridID, ReceiverGridID, &
    FringePadding)

    type(ovk_assembler_properties), intent(inout) :: Properties
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(in) :: FringePadding

    integer :: m, n
    integer :: ms, me, ns, ne

    if (OVK_DEBUG) then
      if (FringePadding < 0) then
        write (ERROR_UNIT, '(a)') "ERROR: Fringe padding must be nonnegative."
        stop 1
      end if
    end if

    call GridIDRange(Properties%ngrids, DonorGridID, ms, me)
    call GridIDRange(Properties%ngrids, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Properties%fringe_padding(m,n) = FringePadding
      end do
    end do

  end subroutine ovkSetAssemblerPropertyFringePadding

  subroutine GridIDRange(NumGrids, GridID, StartID, EndID)

    integer, intent(in) :: NumGrids
    integer, intent(in) :: GridID
    integer, intent(out) :: StartID, EndID

    if (GridID == OVK_ALL_GRIDS) then
      StartID = 1
      EndID = NumGrids
    else
      StartID = GridID
      EndID = GridID
    end if

  end subroutine GridIDRange

  function GetMaxEdgeDistance(Properties, GridID) result(MaxEdgeDist)

    type(ovk_assembler_properties), intent(in) :: Properties
    integer, intent(in) :: GridID
    integer :: MaxEdgeDist

    integer :: m, n

    n = GridID

    MaxEdgeDist = 0

    do m = 1, Properties%ngrids
      if (Properties%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
        if (Properties%disjoint_connection(m,n)) then
          MaxEdgeDist = max(MaxEdgeDist, Properties%fringe_size(m,n) + Properties%fringe_padding(m,n))
        else
          MaxEdgeDist = max(MaxEdgeDist, max(Properties%fringe_size(m,n), Properties%fringe_padding(m,n)))
        end if
      end if
    end do

    ! A little extra
    MaxEdgeDist = MaxEdgeDist + 2

  end function GetMaxEdgeDistance

end module ovkAssembler
