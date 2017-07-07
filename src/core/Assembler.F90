! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkAssembler

  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDonors
  use ovkDomain
  use ovkField
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkInterp
  use ovkMask
  implicit none

  private

  ! API
  public :: ovk_assembler
  public :: ovk_assembler_
  public :: ovk_assembler_properties
  public :: ovk_assembler_properties_
  public :: ovk_assembler_graph
  public :: ovk_assembler_graph_
  public :: ovkCreateAssembler
  public :: ovkDestroyAssembler
  public :: ovkGetAssemblerProperties
  public :: ovkEditAssemblerProperties
  public :: ovkReleaseAssemblerProperties
  public :: ovkGetAssemblerGraph
  public :: ovkEditAssemblerGraph
  public :: ovkReleaseAssemblerGraph
  public :: ovkResetAssemblerGraph
  public :: ovkGetAssemblerDomain
  public :: ovkEditAssemblerDomain
  public :: ovkReleaseAssemblerDomain
!   public :: ovkGetAssemblerOverlap
!   public :: ovkEditAssemblerOverlap
!   public :: ovkReleaseAssemblerOverlap
  public :: ovkGetAssemblerConnectivity
  public :: ovkEditAssemblerConnectivity
  public :: ovkReleaseAssemblerConnectivity
  public :: ovkGetAssemblerPropertyDimension
  public :: ovkGetAssemblerPropertyGridCount
  public :: ovkGetAssemblerPropertyVerbose
  public :: ovkSetAssemblerPropertyVerbose
  public :: ovkGetAssemblerGraphOverlap
  public :: ovkSetAssemblerGraphOverlap
  public :: ovkGetAssemblerGraphOverlapTolerance
  public :: ovkSetAssemblerGraphOverlapTolerance
  public :: ovkGetAssemblerGraphBoundaryHoleCutting
  public :: ovkSetAssemblerGraphBoundaryHoleCutting
  public :: ovkGetAssemblerGraphOverlapHoleCutting
  public :: ovkSetAssemblerGraphOverlapHoleCutting
  public :: ovkGetAssemblerGraphConnectionType
  public :: ovkSetAssemblerGraphConnectionType
  public :: ovkGetAssemblerGraphDisjointConnection
  public :: ovkSetAssemblerGraphDisjointConnection
  public :: ovkGetAssemblerGraphInterpScheme
  public :: ovkSetAssemblerGraphInterpScheme
  public :: ovkGetAssemblerGraphFringeSize
  public :: ovkSetAssemblerGraphFringeSize
  public :: ovkGetAssemblerGraphFringePadding
  public :: ovkSetAssemblerGraphFringePadding

  type ovk_assembler_properties
    ! Read-only
    integer :: nd
    integer :: ngrids
    ! Read/write
    logical :: verbose
  end type ovk_assembler_properties

  type ovk_assembler_graph
    integer :: ngrids
    logical, dimension(:,:), allocatable :: overlap
    real(rk), dimension(:,:), allocatable :: overlap_tolerance
    logical, dimension(:,:), allocatable :: boundary_hole_cutting
    logical, dimension(:,:), allocatable :: overlap_hole_cutting
    integer, dimension(:,:), allocatable :: connection_type
    logical, dimension(:,:), allocatable :: disjoint_connection
    integer, dimension(:,:), allocatable :: interp_scheme
    integer, dimension(:,:), allocatable :: fringe_size
    integer, dimension(:,:), allocatable :: fringe_padding
  end type ovk_assembler_graph

  type ovk_assembler
    type(ovk_assembler_properties), pointer :: properties
    type(ovk_assembler_graph), pointer :: graph
    type(ovk_assembler_graph) :: prev_graph
    type(ovk_domain), pointer :: domain
!     type(ovk_overlap), pointer :: overlap
    type(ovk_connectivity), pointer :: connectivity
    type(ovk_donors), dimension(:), pointer :: donors
    logical :: editing_properties
    logical :: editing_graph
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
  end interface ovk_assembler_properties_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_assembler_graph_
    module procedure ovk_assembler_graph_Default
    module procedure ovk_assembler_graph_Allocated
  end interface ovk_assembler_graph_

contains

  function ovk_assembler_Default() result(Assembler)

    type(ovk_assembler) :: Assembler

    nullify(Assembler%properties)
    nullify(Assembler%graph)
    Assembler%prev_graph = ovk_assembler_graph_()
    nullify(Assembler%domain)
!     nullify(Assembler%overlap)
    nullify(Assembler%donors)
    nullify(Assembler%connectivity)
    Assembler%editing_properties = .false.
    Assembler%editing_graph = .false.
    Assembler%editing_domain = .false.
!     Assembler%editing_overlap = .false.
    Assembler%editing_connectivity = .false.

  end function ovk_assembler_Default

  subroutine ovkCreateAssembler(Assembler, NumDims, NumGrids)

    type(ovk_assembler), intent(out) :: Assembler
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids

    integer :: m

    allocate(Assembler%properties)
    Assembler%properties = ovk_assembler_properties_()

    Assembler%properties%nd = NumDims
    Assembler%properties%ngrids = NumGrids

    allocate(Assembler%graph)
    Assembler%graph = ovk_assembler_graph_(NumGrids)

    Assembler%prev_graph = Assembler%graph

    allocate(Assembler%domain)
    call ovkCreateDomain(Assembler%domain, NumDims, NumGrids, Verbose=.false.)

!     allocate(Assembler%overlap)
!     call ovkCreateOverlap(Assembler%overlap, NumDims, NumGrids, Verbose=.false.)

    allocate(Assembler%donors(NumGrids))
    do m = 1, NumGrids
      Assembler%donors(m) = ovk_donors_()
    end do

    allocate(Assembler%connectivity)
    call ovkCreateConnectivity(Assembler%connectivity, NumDims, NumGrids, Verbose=.false.)

    Assembler%editing_properties = .false.
    Assembler%editing_graph = .false.
    Assembler%editing_domain = .false.
!     Assembler%editing_overlap = .false.
    Assembler%editing_connectivity = .false.

  end subroutine ovkCreateAssembler

  subroutine ovkDestroyAssembler(Assembler)

    type(ovk_assembler), intent(inout) :: Assembler

    integer :: m

    if (associated(Assembler%properties)) deallocate(Assembler%properties)

    if (associated(Assembler%graph)) deallocate(Assembler%graph)
    Assembler%prev_graph = ovk_assembler_graph_()

    if (associated(Assembler%domain)) then
      call ovkDestroyDomain(Assembler%domain)
      deallocate(Assembler%domain)
    end if

!     if (associated(Assembler%overlap)) then
!       call ovkDestroyOverlap(Assembler%overlap)
!       deallocate(Assembler%overlap)
!     end if

    if (associated(Assembler%donors)) then
      do m = 1, size(Assembler%donors)
        call ovkDestroyDonors(Assembler%donors(m))
      end do
      deallocate(Assembler%donors)
    end if

    if (associated(Assembler%connectivity)) then
      call ovkDestroyConnectivity(Assembler%connectivity)
      deallocate(Assembler%connectivity)
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
      Assembler%editing_graph .or. &
      Assembler%editing_domain .or. &
      Assembler%editing_connectivity

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Assembler%editing_graph) then
          write (ERROR_UNIT, '(2a)') "ERROR: Cannot edit assembler properties while editing ", &
            "assembler graph."
        else if (Assembler%editing_domain) then
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

      Properties => Assembler%properties
      Assembler%editing_properties = .true.

    end if

  end subroutine ovkEditAssemblerProperties

  subroutine ovkReleaseAssemblerProperties(Assembler, Properties)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_assembler_properties), pointer, intent(inout) :: Properties

    integer :: NumGrids
    type(ovk_domain_properties), pointer :: DomainProperties
    type(ovk_connectivity_properties), pointer :: ConnectivityProperties

    if (.not. associated(Properties, Assembler%properties)) return
    if (.not. Assembler%editing_properties) then
      nullify(Properties)
      return
    end if

    NumGrids = Assembler%properties%ngrids

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

    nullify(Properties)
    Assembler%editing_properties = .false.

  end subroutine ovkReleaseAssemblerProperties

  subroutine ovkGetAssemblerGraph(Assembler, Graph)

    type(ovk_assembler), intent(in) :: Assembler
    type(ovk_assembler_graph), pointer, intent(out) :: Graph

    Graph => Assembler%graph

  end subroutine ovkGetAssemblerGraph

  subroutine ovkEditAssemblerGraph(Assembler, Graph)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_assembler_graph), pointer, intent(out) :: Graph

    logical :: CannotEdit

    CannotEdit = &
      Assembler%editing_properties .or. &
      Assembler%editing_domain .or. &
      Assembler%editing_connectivity

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Assembler%editing_properties) then
          write (ERROR_UNIT, '(2a)') "ERROR: Cannot edit assembler graph while editing assembler ", &
            "properties."
        else if (Assembler%editing_domain) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit assembler graph while editing domain."
        else
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit assembler graph while editing connectivity."
        end if
        stop 1
      end if

      nullify(Graph)

    else

      if (.not. Assembler%editing_graph) then
        Assembler%prev_graph = Assembler%graph
      end if

      Graph => Assembler%graph
      Assembler%editing_graph = .true.

    end if

  end subroutine ovkEditAssemblerGraph

  subroutine ovkReleaseAssemblerGraph(Assembler, Graph)

    type(ovk_assembler), intent(inout) :: Assembler
    type(ovk_assembler_graph), pointer, intent(inout) :: Graph

    integer :: m, n
    integer :: NumGrids
    integer :: InterpScheme
    integer :: FringeSize
    integer :: FringePadding
    type(ovk_domain_properties), pointer :: DomainProperties
    integer :: MaxEdgeDist, PrevMaxEdgeDist

    if (.not. associated(Graph, Assembler%graph)) return
    if (.not. Assembler%editing_graph) then
      nullify(Graph)
      return
    end if

    NumGrids = Assembler%properties%ngrids

    if (OVK_DEBUG) then

      do m = 1, NumGrids
        if (ovkCartCount(Assembler%domain%grids(m)%cart) > 0) then
          write (ERROR_UNIT, '(2a)') "ERROR: Dynamically updating assembler graph is not yet ", &
            "supported; must be set before grids are created."
          stop 1
        end if
      end do

      if (any(Assembler%graph%connection_type == OVK_CONNECTION_FULL_GRID)) then
        write (ERROR_UNIT, '(a)') "ERROR: OVK_CONNECTION_FULL_GRID is not yet supported."
        stop 1
      end if

      do n = 1, NumGrids
        InterpScheme = OVK_INTERP_LINEAR
        do m = 1, NumGrids
          if (Assembler%graph%connection_type(m,n) /= OVK_CONNECTION_NONE) then
            InterpScheme = Assembler%graph%interp_scheme(m,n)
            exit
          end if
        end do
        do m = 1, NumGrids
          if (Assembler%graph%connection_type(m,n) /= OVK_CONNECTION_NONE) then
            if (Assembler%graph%interp_scheme(m,n) /= InterpScheme) then
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
          if (Assembler%graph%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
            FringeSize = Assembler%graph%fringe_size(m,n)
            exit
          end if
        end do
        do m = 1, NumGrids
          if (Assembler%graph%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
            if (Assembler%graph%fringe_size(m,n) /= FringeSize) then
              write (ERROR_UNIT, '(2a)') "ERROR: Pairwise fringe size specification is ", &
                "not currently supported; must be set uniformly for each receiver grid."
              stop 1
            end if
          end if
        end do
      end do

      do n = 1, NumGrids
        FringePadding = 0
        do m = 1, NumGrids
          if (Assembler%graph%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
            FringePadding = Assembler%graph%fringe_padding(m,n)
            exit
          end if
        end do
        do m = 1, NumGrids
          if (Assembler%graph%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
            if (Assembler%graph%fringe_padding(m,n) /= FringePadding) then
              write (ERROR_UNIT, '(2a)') "ERROR: Pairwise fringe padding specification is ", &
                "not currently supported; must be set uniformly for each receiver grid."
              stop 1
            end if
          end if
        end do
      end do

    end if

    do n = 1, NumGrids
      ! No self-intersection support (yet)
      Assembler%graph%overlap(n,n) = .false.
      ! No overlap implies no hole cutting, no communication
      do m = 1, NumGrids
        if (.not. Assembler%graph%overlap(m,n)) then
          Assembler%graph%boundary_hole_cutting(m,n) = .false.
          Assembler%graph%overlap_hole_cutting(m,n) = .false.
          Assembler%graph%connection_type(m,n) = OVK_CONNECTION_NONE
        end if
      end do
    end do

    do n = 1, NumGrids
      MaxEdgeDist = GetMaxEdgeDistance(Assembler%graph, n)
      PrevMaxEdgeDist = GetMaxEdgeDistance(Assembler%prev_graph, n)
      if (MaxEdgeDist /= PrevMaxEdgeDist) then
        call ovkEditDomainProperties(Assembler%domain, DomainProperties)
        call ovkSetDomainPropertyMaxEdgeDistance(DomainProperties, n, MaxEdgeDist)
        call ovkReleaseDomainProperties(Assembler%domain, DomainProperties)
      end if
    end do

    Assembler%prev_graph = Assembler%graph

    nullify(Graph)
    Assembler%editing_graph = .false.

  end subroutine ovkReleaseAssemblerGraph

  subroutine ovkResetAssemblerGraph(Assembler)

    type(ovk_assembler), intent(inout) :: Assembler

    type(ovk_assembler_graph), pointer :: Graph

    call ovkEditAssemblerGraph(Assembler, Graph)
    Graph = ovk_assembler_graph_(Assembler%properties%ngrids)
    call ovkReleaseAssemblerGraph(Assembler, Graph)

  end subroutine ovkResetAssemblerGraph

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
      Assembler%editing_graph .or. &
      Assembler%editing_connectivity

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Assembler%editing_properties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit domain while editing assembler properties."
        else if (Assembler%editing_graph) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit domain while editing assembler graph."
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
          if (Assembler%graph%overlap(m,n)) then
            Bounds = ovkBBIntersect(Domain%grids(m)%bounds, Domain%grids(n)%bounds)
            if (ovkBBIsEmpty(Bounds)) then
              Assembler%graph%overlap(m,n) = .false.
              Assembler%graph%boundary_hole_cutting(m,n) = .false.
              Assembler%graph%overlap_hole_cutting(m,n) = .false.
              Assembler%graph%connection_type(m,n) = OVK_CONNECTION_NONE
            end if
          end if
        end if
      end do
    end do

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

    CannotEdit = Assembler%editing_properties .or. Assembler%editing_graph .or. &
      Assembler%editing_domain

    if (CannotEdit) then

      if (OVK_DEBUG) then
        if (Assembler%editing_properties) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity while editing assembler properties."
        else if (Assembler%editing_graph) then
          write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity while editing assembler graph."
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

  function ovk_assembler_properties_Default() result(Properties)

    type(ovk_assembler_properties) :: Properties

    Properties%nd = 2
    Properties%ngrids = 0
    Properties%verbose = .false.

  end function ovk_assembler_properties_Default

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

  function ovk_assembler_graph_Default() result(Graph)

    type(ovk_assembler_graph) :: Graph

    Graph%ngrids = 0

  end function ovk_assembler_graph_Default

  function ovk_assembler_graph_Allocated(NumGrids) result(Graph)

    type(ovk_assembler_graph) :: Graph
    integer, intent(in) :: NumGrids

    Graph%ngrids = NumGrids

    allocate(Graph%overlap(NumGrids,NumGrids))
    Graph%overlap = .false.

    allocate(Graph%overlap_tolerance(NumGrids,NumGrids))
    Graph%overlap_tolerance = 0._rk

    allocate(Graph%boundary_hole_cutting(NumGrids,NumGrids))
    Graph%boundary_hole_cutting = .false.

    allocate(Graph%overlap_hole_cutting(NumGrids,NumGrids))
    Graph%overlap_hole_cutting = .false.

    allocate(Graph%connection_type(NumGrids,NumGrids))
    Graph%connection_type = OVK_CONNECTION_NONE

    allocate(Graph%disjoint_connection(NumGrids,NumGrids))
    Graph%disjoint_connection = .false.

    allocate(Graph%interp_scheme(NumGrids,NumGrids))
    Graph%interp_scheme = OVK_INTERP_LINEAR

    allocate(Graph%fringe_size(NumGrids,NumGrids))
    Graph%fringe_size = 0

    allocate(Graph%fringe_padding(NumGrids,NumGrids))
    Graph%fringe_padding = 0

  end function ovk_assembler_graph_Allocated

  subroutine ovkGetAssemblerGraphOverlap(Graph, OverlappingGridID, OverlappedGridID, Overlap)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(out) :: Overlap

    Overlap = Graph%overlap(OverlappingGridID,OverlappedGridID)

  end subroutine ovkGetAssemblerGraphOverlap

  subroutine ovkSetAssemblerGraphOverlap(Graph, OverlappingGridID, OverlappedGridID, Overlap)

    type(ovk_assembler_graph), intent(inout) :: Graph
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(in) :: Overlap

    integer :: m, n
    integer :: ms, me, ns, ne

    call AssemblerGraphIndexRange(Graph, OverlappingGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, OverlappedGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%overlap(m,n) = Overlap
      end do
    end do

  end subroutine ovkSetAssemblerGraphOverlap

  subroutine ovkGetAssemblerGraphOverlapTolerance(Graph, OverlappingGridID, OverlappedGridID, &
    OverlapTolerance)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    real(rk), intent(out) :: OverlapTolerance

    OverlapTolerance = Graph%overlap_tolerance(OverlappingGridID,OverlappedGridID)

  end subroutine ovkGetAssemblerGraphOverlapTolerance

  subroutine ovkSetAssemblerGraphOverlapTolerance(Graph, OverlappingGridID, OverlappedGridID, &
    OverlapTolerance)

    type(ovk_assembler_graph), intent(inout) :: Graph
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

    call AssemblerGraphIndexRange(Graph, OverlappingGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, OverlappedGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%overlap_tolerance(m,n) = OverlapTolerance
      end do
    end do

  end subroutine ovkSetAssemblerGraphOverlapTolerance

  subroutine ovkGetAssemblerGraphBoundaryHoleCutting(Graph, CuttingGridID, CutGridID, &
    BoundaryHoleCutting)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(out) :: BoundaryHoleCutting

    BoundaryHoleCutting = Graph%boundary_hole_cutting(CuttingGridID,CutGridID)

  end subroutine ovkGetAssemblerGraphBoundaryHoleCutting

  subroutine ovkSetAssemblerGraphBoundaryHoleCutting(Graph, CuttingGridID, CutGridID, &
    BoundaryHoleCutting)

    type(ovk_assembler_graph), intent(inout) :: Graph
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(in) :: BoundaryHoleCutting

    integer :: m, n
    integer :: ms, me, ns, ne

    call AssemblerGraphIndexRange(Graph, CuttingGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, CutGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%boundary_hole_cutting(m,n) = BoundaryHoleCutting
      end do
    end do

  end subroutine ovkSetAssemblerGraphBoundaryHoleCutting

  subroutine ovkGetAssemblerGraphOverlapHoleCutting(Graph, CuttingGridID, CutGridID, &
    OverlapHoleCutting)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(out) :: OverlapHoleCutting

    OverlapHoleCutting = Graph%overlap_hole_cutting(CuttingGridID,CutGridID)

  end subroutine ovkGetAssemblerGraphOverlapHoleCutting

  subroutine ovkSetAssemblerGraphOverlapHoleCutting(Graph, CuttingGridID, CutGridID, &
    OverlapHoleCutting)

    type(ovk_assembler_graph), intent(inout) :: Graph
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(in) :: OverlapHoleCutting

    integer :: m, n
    integer :: ms, me, ns, ne

    call AssemblerGraphIndexRange(Graph, CuttingGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, CutGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%overlap_hole_cutting(m,n) = OverlapHoleCutting
      end do
    end do

  end subroutine ovkSetAssemblerGraphOverlapHoleCutting

  subroutine ovkGetAssemblerGraphConnectionType(Graph, DonorGridID, ReceiverGridID, ConnectionType)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: ConnectionType

    ConnectionType = Graph%connection_type(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerGraphConnectionType

  subroutine ovkSetAssemblerGraphConnectionType(Graph, DonorGridID, ReceiverGridID, ConnectionType)

    type(ovk_assembler_graph), intent(inout) :: Graph
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

    call AssemblerGraphIndexRange(Graph, DonorGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%connection_type(m,n) = ConnectionType
      end do
    end do

  end subroutine ovkSetAssemblerGraphConnectionType

  subroutine ovkGetAssemblerGraphDisjointConnection(Graph, DonorGridID, ReceiverGridID, &
    DisjointConnection)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: DonorGridID, ReceiverGridID
    logical, intent(out) :: DisjointConnection

    DisjointConnection = Graph%disjoint_connection(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerGraphDisjointConnection

  subroutine ovkSetAssemblerGraphDisjointConnection(Graph, DonorGridID, ReceiverGridID, &
    DisjointConnection)

    type(ovk_assembler_graph), intent(inout) :: Graph
    integer, intent(in) :: DonorGridID, ReceiverGridID
    logical, intent(in) :: DisjointConnection

    integer :: m, n
    integer :: ms, me, ns, ne

    call AssemblerGraphIndexRange(Graph, DonorGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%disjoint_connection(m,n) = DisjointConnection
      end do
    end do

  end subroutine ovkSetAssemblerGraphDisjointConnection

  subroutine ovkGetAssemblerGraphInterpScheme(Graph, DonorGridID, ReceiverGridID, InterpScheme)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: InterpScheme

    InterpScheme = Graph%interp_scheme(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerGraphInterpScheme

  subroutine ovkSetAssemblerGraphInterpScheme(Graph, DonorGridID, ReceiverGridID, InterpScheme)

    type(ovk_assembler_graph), intent(inout) :: Graph
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

    call AssemblerGraphIndexRange(Graph, DonorGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%interp_scheme(m,n) = InterpScheme
      end do
    end do

  end subroutine ovkSetAssemblerGraphInterpScheme

  subroutine ovkGetAssemblerGraphFringeSize(Graph, DonorGridID, ReceiverGridID, FringeSize)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: FringeSize

    FringeSize = Graph%fringe_size(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerGraphFringeSize

  subroutine ovkSetAssemblerGraphFringeSize(Graph, DonorGridID, ReceiverGridID, FringeSize)

    type(ovk_assembler_graph), intent(inout) :: Graph
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

    call AssemblerGraphIndexRange(Graph, DonorGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%fringe_size(m,n) = FringeSize
      end do
    end do

  end subroutine ovkSetAssemblerGraphFringeSize

  subroutine ovkGetAssemblerGraphFringePadding(Graph, DonorGridID, ReceiverGridID, FringePadding)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: FringePadding

    FringePadding = Graph%fringe_padding(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblerGraphFringePadding

  subroutine ovkSetAssemblerGraphFringePadding(Graph, DonorGridID, ReceiverGridID, FringePadding)

    type(ovk_assembler_graph), intent(inout) :: Graph
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

    call AssemblerGraphIndexRange(Graph, DonorGridID, ms, me)
    call AssemblerGraphIndexRange(Graph, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Graph%fringe_padding(m,n) = FringePadding
      end do
    end do

  end subroutine ovkSetAssemblerGraphFringePadding

  subroutine AssemblerGraphIndexRange(Graph, i, is, ie)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: i
    integer, intent(out) :: is, ie

    if (i == OVK_ALL_GRIDS) then
      is = 1
      ie = Graph%ngrids
    else
      is = i
      ie = i
    end if

  end subroutine AssemblerGraphIndexRange

  function GetMaxEdgeDistance(Graph, GridID) result(MaxEdgeDist)

    type(ovk_assembler_graph), intent(in) :: Graph
    integer, intent(in) :: GridID
    integer :: MaxEdgeDist

    integer :: m, n

    n = GridID

    MaxEdgeDist = 0

    do m = 1, Graph%ngrids
      if (Graph%connection_type(m,n) == OVK_CONNECTION_FRINGE) then
        if (Graph%disjoint_connection(m,n)) then
          MaxEdgeDist = max(MaxEdgeDist, Graph%fringe_size(m,n) + Graph%fringe_padding(m,n))
        else
          MaxEdgeDist = max(MaxEdgeDist, max(Graph%fringe_size(m,n), Graph%fringe_padding(m,n)))
        end if
      end if
    end do

    ! A little extra
    MaxEdgeDist = MaxEdgeDist + 2

  end function GetMaxEdgeDistance

end module ovkAssembler
