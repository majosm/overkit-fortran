! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

program Blobs

  use Overkit
  use ovsCommandLine
  use ovsGlobal
  implicit none

  integer :: i, j, m, n
  character(len=256), dimension(:), allocatable :: RawArguments
  character(len=256) :: Usage
  character(len=256) :: Description
  type(t_cmd_opt), dimension(2) :: Options
  integer :: NumPoints
  character(len=64) :: VisMode
  integer, dimension(2) :: NumPointsBackground
  integer, dimension(2) :: NumPointsBlob
  real(rk) :: SeparationScale
  type(ovk_domain) :: Domain
  type(ovk_grid), pointer :: Grid
  type(ovk_field_real), pointer :: X, Y
  type(ovk_field_int), pointer :: State
  real(rk) :: U, V
  real(rk) :: RMin, RMax
  real(rk) :: Radius
  real(rk) :: Theta
  type(ovk_assembly_options) :: AssemblyOptions
  integer, dimension(2,4) :: NumPointsAll
  type(ovk_plot3d_grid_file) :: GridFile
  type(ovk_cart) :: Cart
  type(ovk_connectivity), pointer :: Connectivity
  integer, dimension(:,:), pointer :: ReceiverPoints
  integer, dimension(MAX_ND) :: Point
  type(ovk_field_logical) :: Mask
  type(ovk_field_int) :: IBlank

  allocate(RawArguments(command_argument_count()))
  do i = 1, size(RawArguments)
    call get_command_argument(i, RawArguments(i))
  end do

  Usage = "Blobs [<options> ...]"
  Description = "Generates an overset mesh representing a rectangular domain with several " // &
    "blob-like holes in it."
  Options(1) = t_cmd_opt_("size", "N", CMD_OPT_INTEGER, "Size of background grid in each " // &
    "direction (other grids are proportional) [ Default: 81 ]")
  Options(2) = t_cmd_opt_("vis-mode", "V", CMD_OPT_STRING, "How to set IBlank of output file " // &
    "for visualization (fringe or state) [ Default: fringe ]")

  call ParseCommandLineArguments(RawArguments, Usage=Usage, Description=Description, &
    Options=Options)

  call GetCommandLineOptionValue(Options(1), NumPoints, 81)

  call GetCommandLineOptionValue(Options(2), VisMode, "fringe")
  if (VisMode /= "fringe" .and. VisMode /= "state") then
    write (ERROR_UNIT, '(3a)') "ERROR: Invalid visualization mode '", trim(VisMode), "'."
    stop 1
  end if

  ! Initialize the domain
  call ovkCreateDomain(Domain, NumDims=2, NumGrids=4, StatusLogFile=OUTPUT_UNIT, &
    ErrorLogFile=ERROR_UNIT)

  !=================
  ! Background grid
  !=================

  NumPointsBackground = [NumPoints,NumPoints]

  ! Initialize grid data structure for background grid
  ! Set geometry type as a hint for potential performance improvements
  call ovkCreateGrid(Domain, 1, NumPoints=NumPointsBackground, GeometryType=OVK_GEOMETRY_CARTESIAN)
  call ovkEditGrid(Domain, 1, Grid)

  ! Generate coordinates for background grid
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBackground(2)
    do i = 1, NumPointsBackground(1)
      U = real(i-1,kind=rk)/real(NumPointsBackground(1)-1,kind=rk)
      V = real(j-1,kind=rk)/real(NumPointsBackground(2)-1,kind=rk)
      X%values(i,j,1) = 2._rk * (U-0.5_rk)
      Y%values(i,j,1) = 2._rk * (V-0.5_rk)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  call ovkReleaseGrid(Domain, Grid)

  !============
  ! Blob grids
  !============

  NumPointsBlob = [NumPoints,2*NumPoints]
  SeparationScale = 0.8_rk

  ! Initialize grid data structure for blob #1
  ! Periodic in the angular direction, with the last set of points being equal to the first
  call ovkCreateGrid(Domain, 2, NumPoints=NumPointsBlob, Periodic=[.false.,.true.], &
    PeriodicStorage=OVK_OVERLAP_PERIODIC)
  call ovkEditGrid(Domain, 2, Grid)

  ! Generate coordinates for blob #1
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBlob(2)
    do i = 1, NumPointsBlob(1)
      U = real(i-1, kind=rk)/real(NumPointsBlob(1)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPointsBlob(2)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(3._rk*Theta) + 0.1_rk*sin(2._rk*Theta+Pi/4._rk))
      RMax = 0.5_rk + RMin
      Radius = RMin * (RMax/RMin)**U
      X%values(i,j,1) = -0.425_rk*SeparationScale + Radius * cos(Theta)
      Y%values(i,j,1) = -0.025_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  ! Inner radial boundary on blob #1
  call ovkEditGridState(Grid, State)
  State%values(1,:,1) = OVK_DOMAIN_BOUNDARY_POINT
  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  ! Initialize grid data structures for blob #2
  ! Periodic in the angular direction, with the last set of points being equal to the first
  call ovkCreateGrid(Domain, 3, NumPoints=NumPointsBlob, Periodic=[.false.,.true.], &
    PeriodicStorage=OVK_OVERLAP_PERIODIC)
  call ovkEditGrid(Domain, 3, Grid)

  ! Generate coordinates for blob #2
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBlob(2)
    do i = 1, NumPointsBlob(1)
      U = real(i-1, kind=rk)/real(NumPointsBlob(1)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPointsBlob(2)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(4._rk*Theta+Pi/4._rk) + 0.1_rk*sin(2._rk*Theta))
      RMax = 0.5_rk + RMin
      Radius = RMin * (RMax/RMin)**U
      X%values(i,j,1) = 0.075_rk*SeparationScale + Radius * cos(Theta)
      Y%values(i,j,1) = 0.425_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  ! Inner radial boundary on blob #2
  call ovkEditGridState(Grid, State)
  State%values(1,:,1) = OVK_DOMAIN_BOUNDARY_POINT
  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  ! Initialize grid data structures for blob #3
  ! Periodic in the angular direction, with the last set of points being equal to the first
  call ovkCreateGrid(Domain, 4, NumPoints=NumPointsBlob, Periodic=[.false.,.true.], &
    PeriodicStorage=OVK_OVERLAP_PERIODIC)
  call ovkEditGrid(Domain, 4, Grid)

  ! Generate coordinates for blob #3
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBlob(2)
    do i = 1, NumPointsBlob(1)
      U = real(i-1, kind=rk)/real(NumPointsBlob(1)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPointsBlob(2)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(5._rk*Theta+Pi/4._rk) + 0.1_rk*sin(3._rk*Theta))
      RMax = 0.5_rk + RMin
      Radius = RMin * (RMax/RMin)**U
      X%values(i,j,1) = 0.375_rk*SeparationScale + Radius * cos(Theta)
      Y%values(i,j,1) = -0.375_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  ! Inner radial boundary on blob #3
  call ovkEditGridState(Grid, State)
  State%values(1,:,1) = OVK_DOMAIN_BOUNDARY_POINT
  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  !==================
  ! Overset assembly
  !==================

  AssemblyOptions = ovk_assembly_options_(2, 4)

  ! Indicate which grids can overlap
  call ovkSetAssemblyOptionOverlappable(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)

  ! Automatically define boundaries in non-overlapping regions
  call ovkSetAssemblyOptionInferBoundaries(AssemblyOptions, OVK_ALL_GRIDS, .true.)
  ! Indicate which grids can cut each other
  call ovkSetAssemblyOptionCutBoundaryHoles(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)

  ! Indicate how to treat overlap between grids
  call ovkSetAssemblyOptionOccludes(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_OCCLUDES_COARSE)
  call ovkSetAssemblyOptionOccludes(AssemblyOptions, OVK_ALL_GRIDS, 1, OVK_OCCLUDES_ALL)

  ! Retain some extra overlap between grids
  call ovkSetAssemblyOptionEdgePadding(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 2)
  call ovkSetAssemblyOptionEdgeSmoothing(AssemblyOptions, OVK_ALL_GRIDS, 2)

  ! Indicate which grids can communicate and how
  call ovkSetAssemblyOptionConnectionType(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_LINEAR)
  ! Set the number of interpolation fringe layers
  call ovkSetAssemblyOptionFringeSize(AssemblyOptions, OVK_ALL_GRIDS, 2)
  ! Cut out non-fringe interpolating points
  call ovkSetAssemblyOptionMinimizeOverlap(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)

  call ovkAssemble(Domain, AssemblyOptions)

  !========
  ! Output
  !========

  NumPointsAll(:,1) = NumPointsBackground
  NumPointsAll(:,2) = NumPointsBlob
  NumPointsAll(:,3) = NumPointsBlob
  NumPointsAll(:,4) = NumPointsBlob

  ! Write a PLOT3D grid file
  call ovkCreateP3D(GridFile, "blobs.xyz", NumDims=2, NumGrids=4, NumPointsAll=NumPointsAll, &
    WithIBlank=.true., StatusLogFile=OUTPUT_UNIT, ErrorLogFile=ERROR_UNIT)

  do n = 1, 4

    call ovkGetGrid(Domain, n, Grid)
    call ovkGetGridCart(Grid, Cart)
    call ovkGetGridCoords(Grid, 1, X)
    call ovkGetGridCoords(Grid, 2, Y)

    ! Use IBlank data to visualize status of grid points

    select case (VisMode)
    case ("fringe")

      ! IBlank == 1 => Normal
      IBlank = ovk_field_int_(Cart, 1)

      ! IBlank == 0 => Hole
      call ovkFilterGridState(Grid, OVK_STATE_GRID, OVK_NONE, Mask)
      IBlank%values = merge(0, IBlank%values, Mask%values)

      ! IBlank == -N => Receives from grid N
      do m = 1, 4
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

    case ("state")

      ! IBlank == 1 => Normal
      IBlank = ovk_field_int_(Cart, 1)

      ! IBlank == 2 => Domain boundary
      call ovkFilterGridState(Grid, OVK_STATE_DOMAIN_BOUNDARY, OVK_ALL, Mask)
      IBlank%values = merge(2, IBlank%values, Mask%values)

      ! IBlank == 3 => Boundary hole
      call ovkFilterGridState(Grid, OVK_STATE_BOUNDARY_HOLE, OVK_ALL, Mask)
      IBlank%values = merge(3, IBlank%values, Mask%values)

      ! IBlank == 4 => Overlap minimized
      call ovkFilterGridState(Grid, OVK_STATE_OVERLAP_MINIMIZED, OVK_ALL, Mask)
      IBlank%values = merge(4, IBlank%values, Mask%values)

      ! IBlank == 5 => Receiver
      call ovkFilterGridState(Grid, OVK_STATE_RECEIVER, OVK_ALL, Mask)
      IBlank%values = merge(5, IBlank%values, Mask%values)

      ! IBlank == 6 => Orphan
      call ovkFilterGridState(Grid, OVK_STATE_ORPHAN, OVK_ALL, Mask)
      IBlank%values = merge(6, IBlank%values, Mask%values)

    end select

    call ovkWriteP3D(GridFile, n, X, Y, IBlank)

  end do

  ! Finalize the PLOT3D file
  call ovkCloseP3D(GridFile)

  call ovkDestroyDomain(Domain)

end program Blobs
