! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

program Inlet

  use Overkit
  use ovsCommandLine
  use ovsGlobal
  implicit none

  integer :: i, j, m, n
  character(len=256), dimension(:), allocatable :: RawArguments
  character(len=256) :: Usage
  character(len=256) :: Description
  type(t_cmd_opt), dimension(1) :: Options
  integer :: NumPoints
  integer, dimension(2) :: NumPointsBox, NumPointsInlet
  type(ovk_bbox) :: BoundsBox, BoundsInlet
  real(rk) :: hX, hY
  integer :: iLeftBoundary, iRightBoundary, jBoundary
  type(ovk_domain) :: Domain
  type(ovk_domain_properties), pointer :: Properties
  type(ovk_grid), pointer :: Grid
  type(ovk_field_real), pointer :: X, Y
  type(ovk_field_int), pointer :: State
  real(rk) :: U, V
  integer, dimension(2,2) :: NumPointsAll
  type(ovk_plot3d_grid_file) :: GridFile
  type(ovk_cart) :: Cart
  integer :: ConnectionType
  type(ovk_connectivity), pointer :: Connectivity
  integer, dimension(:,:), pointer :: ReceiverPoints
  integer, dimension(MAX_ND) :: Point
  type(ovk_field_logical) :: Mask
  type(ovk_field_int) :: IBlank

  allocate(RawArguments(command_argument_count()))
  do i = 1, size(RawArguments)
    call get_command_argument(i, RawArguments(i))
  end do

  Usage = "Inlet [<options> ...]"
  Description = "Generates an overset mesh representing a box with an inlet."
  Options(1) = t_cmd_opt_("size", "N", CMD_OPT_INTEGER, "Size of box grid in each " // &
    "direction (inlet grid is proportional) [ Default: 101 ]")

  call ParseArguments(RawArguments, Usage=Usage, Description=Description, Options=Options)

  call GetOptionValue(Options(1), NumPoints, 101)

  ! Initialize the domain
  call ovkCreateDomain(Domain, NumDims=2, NumGrids=2, Verbose=.true.)

  call ovkEditDomainProperties(Domain, Properties)

  ! Indicate which grids can intersect
  call ovkSetDomainPropertyOverlappable(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
  call ovkSetDomainPropertyOverlapTolerance(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 1.e-10_rk)

  ! Automatically define boundaries in non-overlapping regions
  call ovkSetDomainPropertyInferBoundaries(Properties, OVK_ALL_GRIDS, .true.)

  ! Retain some extra overlap between grids
  call ovkSetDomainPropertyEdgePadding(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 2)
  call ovkSetDomainPropertyEdgeSmoothing(Properties, OVK_ALL_GRIDS, 2)

  ! Indicate which grids can communicate and how
  call ovkSetDomainPropertyConnectionType(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_FRINGE)
  call ovkSetDomainPropertyInterpScheme(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_INTERP_CUBIC)
  call ovkSetDomainPropertyFringeSize(Properties, OVK_ALL_GRIDS, 2)
  call ovkSetDomainPropertyOverlapMinimization(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)

  call ovkReleaseDomainProperties(Domain, Properties)

  !==========
  ! Box grid
  !==========

  NumPointsBox = [NumPoints,NumPoints/2]

  BoundsBox = ovk_bbox_(2, [-1._rk, -0.5_rk], [1._rk, 0.5_rk])

  ! Initialize grid data structure for box grid
  ! Set geometry type as a hint for potential performance improvements
  call ovkCreateGrid(Domain, 1, NumPoints=NumPointsBox, GeometryType=OVK_GRID_GEOMETRY_CARTESIAN)
  call ovkEditGrid(Domain, 1, Grid)

  ! Generate coordinates for box grid
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBox(2)
    do i = 1, NumPointsBox(1)
      U = real(i-1,kind=rk)/real(NumPointsBox(1)-1,kind=rk)
      V = real(j-1,kind=rk)/real(NumPointsBox(2)-1,kind=rk)
      X%values(i,j,1) = (1._rk-U)*BoundsBox%b(1) + U*BoundsBox%e(1)
      Y%values(i,j,1) = (1._rk-V)*BoundsBox%b(2) + V*BoundsBox%e(2)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  ! Lower wall boundary (other boundaries can be inferred)
  call ovkEditGridState(Grid, State)
  State%values(:,1,1) = OVK_DOMAIN_BOUNDARY_POINT
  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  !============
  ! Inlet grid
  !============

  NumPointsInlet = [NumPoints,NumPoints]

  iLeftBoundary = NumPointsInlet(1)/4+1
  iRightBoundary = NumPointsInlet(1)-(NumPointsInlet(1)/4+1)
  jBoundary = NumPointsInlet(2)/2+1

  hX = 0.5_rk/real(iRightBoundary-iLeftBoundary+1,kind=rk)
  hY = 0.5_rk/real(jBoundary-1,kind=rk)

  BoundsInlet = ovk_bbox_(2)
  BoundsInlet%b(1) = -0.25_rk - hX*real(iLeftBoundary-1,kind=rk)
  BoundsInlet%e(1) = 0.25_rk + hX*real(NumPointsInlet(1)-iRightBoundary,kind=rk)
  BoundsInlet%b(2) = BoundsBox%b(2) - 0.5_rk
  BoundsInlet%e(2) = BoundsBox%b(2) + hY*real(NumPointsInlet(2)-jBoundary,kind=rk)

  ! Initialize grid data structure for inlet grid
  ! Set geometry type as a hint for potential performance improvements
  call ovkCreateGrid(Domain, 2, NumPoints=NumPointsInlet, GeometryType=OVK_GRID_GEOMETRY_CARTESIAN)
  call ovkEditGrid(Domain, 2, Grid)

  ! Generate coordinates for inlet grid
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsInlet(2)
    do i = 1, NumPointsInlet(1)
      U = real(i-1, kind=rk)/real(NumPointsInlet(1)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPointsInlet(2)-1, kind=rk)
      X%values(i,j,1) = (1._rk-U)*BoundsInlet%b(1) + U*BoundsInlet%e(1)
      Y%values(i,j,1) = (1._rk-V)*BoundsInlet%b(2) + V*BoundsInlet%e(2)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  call ovkEditGridState(Grid, State)

  ! Collar hole regions
  State%values(:iLeftBoundary-1,:jBoundary-1,1) = OVK_HOLE_POINT
  State%values(iRightBoundary+1:,:jBoundary-1,1) = OVK_HOLE_POINT

  ! Inlet boundaries
  State%values(:iLeftBoundary,jBoundary,1) = OVK_DOMAIN_BOUNDARY_POINT
  State%values(iLeftBoundary,:jBoundary,1) = OVK_DOMAIN_BOUNDARY_POINT
  State%values(iRightBoundary,:jBoundary,1) = OVK_DOMAIN_BOUNDARY_POINT
  State%values(iRightBoundary:,jBoundary,1) = OVK_DOMAIN_BOUNDARY_POINT
  State%values(iLeftBoundary:iRightBoundary,1,1) = OVK_DOMAIN_BOUNDARY_POINT

  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  !==================
  ! Overset assembly
  !==================

  call ovkAssemble(Domain)

  !========
  ! Output
  !========

  NumPointsAll(:,1) = NumPointsBox
  NumPointsAll(:,2) = NumPointsInlet

  ! Write a PLOT3D grid file
  call ovkCreateP3D(GridFile, "inlet.xyz", NumDims=2, NumGrids=2, NumPointsAll=NumPointsAll, &
    WithIBlank=.true., Verbose=.true.)

  call ovkGetDomainProperties(Domain, Properties)

  do n = 1, 2

    call ovkGetGrid(Domain, n, Grid)
    call ovkGetGridCart(Grid, Cart)
    call ovkGetGridCoords(Grid, 1, X)
    call ovkGetGridCoords(Grid, 2, Y)

    ! Use IBlank data to visualize status of grid points
    ! IBlank == 1 => Normal
    IBlank = ovk_field_int_(Cart, 1)

    ! IBlank == 0 => Hole
    call ovkFilterGridState(Grid, OVK_STATE_HOLE, OVK_ALL, Mask)
    IBlank%values = merge(0, IBlank%values, Mask%values)

    ! IBlank == -N => Receives from grid N
    do m = 1, 2
      call ovkGetDomainPropertyConnectionType(Properties, m, n, ConnectionType)
      if (ConnectionType /= OVK_CONNECTION_NONE) then
        call ovkGetConnectivity(Domain, m, n, Connectivity)
        call ovkGetConnectivityReceiverPoints(Connectivity, ReceiverPoints)
        do i = 1, size(ReceiverPoints,2)
          Point = ReceiverPoints(:,i)
          IBlank%values(Point(1),Point(2),Point(3)) = -m
        end do
      end if
    end do

    ! IBlank == 7 => Orphan
    call ovkFilterGridState(Grid, OVK_STATE_ORPHAN, OVK_ALL, Mask)
    IBlank%values = merge(7, IBlank%values, Mask%values)

    call ovkWriteP3D(GridFile, n, X, Y, IBlank)

  end do

  ! Finalize the PLOT3D file
  call ovkCloseP3D(GridFile)

  call ovkDestroyDomain(Domain)

end program Inlet
