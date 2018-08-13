! Copyright (c) 2018 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

program Inlet

  use Overkit
  use ovsCommandLine
  use ovsGlobal
  implicit none

  integer :: i, j, k, m, n
  character(len=256), dimension(:), allocatable :: RawArguments
  character(len=256) :: Usage
  character(len=256) :: Description
  type(t_cmd_opt), dimension(1) :: Options
  integer :: NumPoints
  type(ovk_domain) :: Domain
  type(ovk_grid), pointer :: Grid
  type(ovk_field_real), pointer :: X, Y, Z
  type(ovk_field_int), pointer :: State
  integer, dimension(3) :: NumPointsBox, NumPointsInletOuter, NumPointsInletInner
  type(ovk_bbox) :: BoundsBox
  integer :: iBoundary, kBoundary
  real(rk) :: U, V, W
  real(rk) :: RMin, RMax
  real(rk) :: Radius
  real(rk) :: Theta
  real(rk) :: ZShift
  type(ovk_assembly_options) :: AssemblyOptions
  integer, dimension(3,3) :: NumPointsAll
  type(ovk_plot3d_grid_file) :: GridFile
  type(ovk_cart) :: Cart
  type(ovk_connectivity), pointer :: Connectivity
  integer, dimension(:,:), pointer :: ReceiverPoints
  integer, dimension(MAX_DIMS) :: Point
  type(ovk_field_logical) :: Mask
  type(ovk_field_int) :: IBlank

  allocate(RawArguments(command_argument_count()))
  do i = 1, size(RawArguments)
    call get_command_argument(i, RawArguments(i))
  end do

  Usage = "Inlet [<options> ...]"
  Description = "Generates an overset mesh representing a box with a cylindrical inlet."
  Options(1) = t_cmd_opt_("size", "N", CMD_OPT_INTEGER, "Size of box grid in each " // &
    "direction (inlet grids are proportional) [ Default: 81 ]")

  call ParseCommandLineArguments(RawArguments, Usage=Usage, Description=Description, &
    Options=Options)

  call GetCommandLineOptionValue(Options(1), NumPoints, 81)

  ! Initialize the domain
  call ovkCreateDomain(Domain, NumDims=3, NumGrids=3, StatusLogFile=OUTPUT_UNIT, &
      ErrorLogFile=ERROR_UNIT)

  !==========
  ! Box grid
  !==========

  NumPointsBox = [NumPoints,NumPoints,NumPoints]

  BoundsBox = ovk_bbox_(3, [-1._rk, -1._rk, -1._rk], [1._rk, 1._rk, 1._rk])

  ! Initialize grid data structure for box grid
  ! Set geometry type as a hint for potential performance improvements
  call ovkCreateGrid(Domain, 1, NumPoints=NumPointsBox, GeometryType=OVK_GEOMETRY_UNIFORM)
  call ovkEditGrid(Domain, 1, Grid)

  ! Generate coordinates for box grid
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  call ovkEditGridCoords(Grid, 3, Z)
  do k = 1, NumPointsBox(3)
    do j = 1, NumPointsBox(2)
      do i = 1, NumPointsBox(1)
        U = real(i-1,kind=rk)/real(NumPointsBox(1)-1,kind=rk)
        V = real(j-1,kind=rk)/real(NumPointsBox(2)-1,kind=rk)
        W = real(k-1,kind=rk)/real(NumPointsBox(3)-1,kind=rk)
        X%values(i,j,k) = (1._rk-U)*BoundsBox%b(1) + U*BoundsBox%e(1)
        Y%values(i,j,k) = (1._rk-V)*BoundsBox%b(2) + V*BoundsBox%e(2)
        Z%values(i,j,k) = (1._rk-W)*BoundsBox%b(3) + W*BoundsBox%e(3)
      end do
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)
  call ovkReleaseGridCoords(Grid, Z)

  ! Lower wall boundary (other boundaries can be inferred)
  call ovkEditGridState(Grid, State)
  State%values(:,:,1) = OVK_DOMAIN_BOUNDARY_POINT
  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  !==================
  ! Inlet outer grid
  !==================

  NumPointsInletOuter = [NumPoints/3,2*NumPoints,NumPoints]

  iBoundary = (2*NumPointsInletOuter(1))/3+1
  kBoundary = NumPointsInletOuter(3)/2+1

  ! Initialize grid data structure for inlet grid
  call ovkCreateGrid(Domain, 2, NumPoints=NumPointsInletOuter, Periodic=[.false.,.true.,.false.], &
    PeriodicStorage=OVK_PERIODIC_STORAGE_DUPLICATED)

  call ovkEditGrid(Domain, 2, Grid)

  ! Generate coordinates for inlet grid
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  call ovkEditGridCoords(Grid, 3, Z)

  ! Initially place in the center of the box
  do k = 1, NumPointsInletOuter(3)
    do j = 1, NumPointsInletOuter(2)
      do i = 1, NumPointsInletOuter(1)
        U = real(i-1,kind=rk)/real(NumPointsInletOuter(1)-1,kind=rk)
        V = real(j-1,kind=rk)/real(NumPointsInletOuter(2)-1,kind=rk)
        W = real(k-1,kind=rk)/real(NumPointsInletOuter(3)-1,kind=rk)
        RMin = 0.2_rk
        RMax = 0.8_rk
        Radius = (RMin+0.5_rk) * ((RMax+0.5_rk)/(RMin+0.5_rk))**U - 0.5_rk
        Theta = 2._rk * Pi * V
        X%values(i,j,k) = Radius * cos(Theta)
        Y%values(i,j,k) = Radius * sin(Theta)
        Z%values(i,j,k) = (1._rk-W)*BoundsBox%b(3) + W*BoundsBox%e(3)
      end do
    end do
  end do

  ! Then shift down so kBoundary point lines up with box wall
  ZShift = -(Z%values(1,1,kBoundary)+1._rk)
  do k = 1, NumPointsInletOuter(3)
    do j = 1, NumPointsInletOuter(2)
      do i = 1, NumPointsInletOuter(1)
        Z%values(i,j,k) = Z%values(i,j,k) + ZShift
      end do
    end do
  end do

  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)
  call ovkReleaseGridCoords(Grid, Z)

  call ovkEditGridState(Grid, State)

  ! Collar hole regions
  State%values(iBoundary+1:,:,:kBoundary-1) = OVK_EXTERIOR_POINT

  ! Inlet boundaries
  State%values(:iBoundary,:,1) = OVK_DOMAIN_BOUNDARY_POINT
  State%values(iBoundary,:,:kBoundary) = OVK_DOMAIN_BOUNDARY_POINT
  State%values(iBoundary:,:,kBoundary) = OVK_DOMAIN_BOUNDARY_POINT

  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  !==================
  ! Inlet inner grid
  !==================

  NumPointsInletInner = [NumPoints/2,NumPoints/2,NumPoints]

  kBoundary = NumPointsInletInner(3)/2+1

  ! Initialize grid data structure for inlet grid
  ! Set geometry type as a hint for potential performance improvements
  call ovkCreateGrid(Domain, 3, NumPoints=NumPointsInletInner, GeometryType=OVK_GEOMETRY_UNIFORM)
  call ovkEditGrid(Domain, 3, Grid)

  ! Generate coordinates for inlet grid
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  call ovkEditGridCoords(Grid, 3, Z)

  ! Initially place in the center of the box
  do k = 1, NumPointsInletInner(3)
    do j = 1, NumPointsInletInner(2)
      do i = 1, NumPointsInletInner(1)
        U = real(i-1,kind=rk)/real(NumPointsInletInner(1)-1,kind=rk)
        V = real(j-1,kind=rk)/real(NumPointsInletInner(2)-1,kind=rk)
        W = real(k-1,kind=rk)/real(NumPointsInletInner(3)-1,kind=rk)
        X%values(i,j,k) = (1._rk-U)*0.3_rk*BoundsBox%b(1) + U*0.3_rk*BoundsBox%e(1)
        Y%values(i,j,k) = (1._rk-V)*0.3_rk*BoundsBox%b(2) + V*0.3_rk*BoundsBox%e(2)
        Z%values(i,j,k) = (1._rk-W)*BoundsBox%b(3) + W*BoundsBox%e(3)
      end do
    end do
  end do

  ! Then shift down so kBoundary point lines up with box wall
  ZShift = -(Z%values(1,1,kBoundary)+1._rk)
  do k = 1, NumPointsInletInner(3)
    do j = 1, NumPointsInletInner(2)
      do i = 1, NumPointsInletInner(1)
        Z%values(i,j,k) = Z%values(i,j,k) + ZShift
      end do
    end do
  end do

  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)
  call ovkReleaseGridCoords(Grid, Z)

  call ovkEditGridState(Grid, State)

  ! Inlet boundaries
  State%values(:,:,1) = OVK_DOMAIN_BOUNDARY_POINT

  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  !==================
  ! Overset assembly
  !==================

  AssemblyOptions = ovk_assembly_options_(3, 3)

  ! Indicate which grids can intersect
  call ovkSetAssemblyOptionOverlappable(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
  call ovkSetAssemblyOptionOverlapTolerance(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 1.e-10_rk)

  ! Automatically define boundaries in non-overlapping regions
  call ovkSetAssemblyOptionInferBoundaries(AssemblyOptions, OVK_ALL_GRIDS, .true.)

  ! Indicate which grids can cut each other
  call ovkSetAssemblyOptionCutBoundaryHoles(AssemblyOptions, 2, 3, .true.)

  ! Indicate how to treat overlap between grids
  call ovkSetAssemblyOptionOccludes(AssemblyOptions, 2, 1, OVK_OCCLUDES_ALL)
  call ovkSetAssemblyOptionOccludes(AssemblyOptions, 3, 1, OVK_OCCLUDES_ALL)
  call ovkSetAssemblyOptionOccludes(AssemblyOptions, 3, 2, OVK_OCCLUDES_ALL)

  ! Retain some extra overlap between grids
  call ovkSetAssemblyOptionEdgePadding(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 2)
  call ovkSetAssemblyOptionEdgeSmoothing(AssemblyOptions, OVK_ALL_GRIDS, 4)

  ! Indicate which grids can communicate and how
  call ovkSetAssemblyOptionConnectionType(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_CUBIC)
  call ovkSetAssemblyOptionFringeSize(AssemblyOptions, OVK_ALL_GRIDS, 2)
  call ovkSetAssemblyOptionMinimizeOverlap(AssemblyOptions, OVK_ALL_GRIDS, 1, .true.)
  call ovkSetAssemblyOptionMinimizeOverlap(AssemblyOptions, OVK_ALL_GRIDS, 2, .true.)

  call ovkAssemble(Domain, AssemblyOptions)

  !========
  ! Output
  !========

  NumPointsAll(:,1) = NumPointsBox
  NumPointsAll(:,2) = NumPointsInletOuter
  NumPointsAll(:,3) = NumPointsInletInner

  ! Write a PLOT3D grid file
  call ovkCreateP3D(GridFile, "inlet.xyz", NumDims=3, NumGrids=3, NumPointsAll=NumPointsAll, &
    WithIBlank=.true., StatusLogFile=OUTPUT_UNIT, ErrorLogFile=ERROR_UNIT)

  do n = 1, 3

    call ovkGetGrid(Domain, n, Grid)
    call ovkGetGridCart(Grid, Cart)
    call ovkGetGridCoords(Grid, 1, X)
    call ovkGetGridCoords(Grid, 2, Y)
    call ovkGetGridCoords(Grid, 3, Z)

    ! Use IBlank data to visualize status of grid points
    ! IBlank == 1 => Normal
    IBlank = ovk_field_int_(Cart, 1)

    ! IBlank == 0 => Hole
    call ovkFilterGridState(Grid, OVK_STATE_GRID, OVK_NONE, Mask)
    IBlank%values = merge(0, IBlank%values, Mask%values)

    ! IBlank == -N => Receives from grid N
    do m = 1, 3
      if (ovkHasConnectivity(Domain, m, n)) then
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

    call ovkWriteP3D(GridFile, n, X, Y, Z, IBlank)

  end do

  ! Finalize the PLOT3D file
  call ovkCloseP3D(GridFile)

  call ovkDestroyDomain(Domain)

end program Inlet
