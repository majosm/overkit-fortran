! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

program Circle

  use Overkit
  use ovsCommandLine
  use ovsGlobal
  implicit none

  integer :: i
  character(len=256), dimension(:), allocatable :: RawArguments
  character(len=256) :: Usage
  character(len=256) :: Description
  type(t_cmd_opt), dimension(2) :: Options
  integer :: NumPoints
  character(len=64) :: Variant
  integer, dimension(2) :: NumPointsBackground, NumPointsCylinder
  integer, dimension(2) :: NumPointsBlock
  real(rk), dimension(:,:,:), allocatable :: CoordsBackground, CoordsCylinder
  real(rk), dimension(:,:,:), allocatable :: CoordsEastBlock, CoordsNorthBlock, CoordsWestBlock, &
    CoordsSouthBlock
  integer, dimension(:,:), allocatable :: IBlankBackground

  allocate(RawArguments(command_argument_count()))
  do i = 1, size(RawArguments)
    call get_command_argument(i, RawArguments(i))
  end do

  Usage = "Circle [<options> ...]"
  Description = "Generates an overset mesh representing a rectangular domain with a circular cutout."
  Options(1) = t_cmd_opt_("size", "N", CMD_OPT_INTEGER, "Characteristic size of grids " // &
    "[ Default: 81 ]")
  Options(2) = t_cmd_opt_("variant", "R", CMD_OPT_STRING, "Variant of mesh to generate (" // &
    "cylinder, block, or remap) [ Default: cylinder ]")

  call ParseCommandLineArguments(RawArguments, Usage=Usage, Description=Description, &
    Options=Options)

  call GetCommandLineOptionValue(Options(1), NumPoints, 81)

  call GetCommandLineOptionValue(Options(2), Variant, "cylinder")
  if (Variant /= "cylinder" .and. Variant /= "block" .and. Variant /= "remap") then
    write (ERROR_UNIT, '(3a)') "ERROR: Invalid mesh variant '", trim(Variant), "'."
    stop 1
  end if

  if (Variant == "cylinder") then
    call GenerateCylinder()
  else if (Variant == "block") then
    call GenerateBlock()
  else
    call GenerateCylinder()
    call GenerateBlock()
    call GenerateRemap()
  end if

contains

  subroutine GenerateCylinder()

    integer :: i, j, m, n
    type(ovk_domain) :: Domain
    type(ovk_grid), pointer :: Grid
    type(ovk_field_real), pointer :: X, Y
    type(ovk_field_int), pointer :: State
    real(rk) :: U, V
    real(rk) :: RMin, RMax
    real(rk) :: Radius
    real(rk) :: Theta
    type(ovk_assembly_options) :: AssemblyOptions
    integer, dimension(2,2) :: NumPointsAll
    type(ovk_plot3d_grid_file) :: GridFile
    type(ovk_cart) :: Cart
    type(ovk_connectivity), pointer :: Connectivity
    integer, dimension(:,:), pointer :: ReceiverPoints
    integer, dimension(MAX_ND) :: Point
    type(ovk_field_logical) :: Mask
    type(ovk_field_int) :: IBlank

    !=====================
    ! Data initialization
    !=====================

    NumPointsBackground = [NumPoints,NumPoints]
    NumPointsCylinder = [NumPoints/3,2*NumPoints]

    allocate(CoordsBackground(NumPointsBackground(1),NumPointsBackground(2),2))
    allocate(CoordsCylinder(NumPointsCylinder(1),NumPointsCylinder(2),2))

    ! Generate coordinates for background
    do j = 1, NumPointsBackground(2)
      do i = 1, NumPointsBackground(1)
        U = real(i-1,kind=rk)/real(NumPointsBackground(1)-1,kind=rk)
        V = real(j-1,kind=rk)/real(NumPointsBackground(2)-1,kind=rk)
        CoordsBackground(i,j,1) = 2._rk * (U-0.5_rk)
        CoordsBackground(i,j,2) = 2._rk * (V-0.5_rk)
      end do
    end do

    ! Generate coordinates for cylinder
    do j = 1, NumPointsCylinder(2)
      do i = 1, NumPointsCylinder(1)
        U = real(i-1, kind=rk)/real(NumPointsCylinder(1)-1, kind=rk)
        V = real(j-1, kind=rk)/real(NumPointsCylinder(2)-1, kind=rk)
        RMin = 0.3_rk
        RMax = 0.7_rk
        Radius = RMin * (RMax/RMin)**U
        Theta = 2._rk * Pi * V
        CoordsCylinder(i,j,1) = Radius * cos(Theta)
        CoordsCylinder(i,j,2) = Radius * sin(Theta)
      end do
    end do

    allocate(IBlankBackground(NumPointsBackground(1),NumPointsBackground(2)))

    ! Initialize IBlank to 1 for now
    IBlankBackground = 1

    !========================
    ! Overkit initialization
    !========================

    ! Initialize the domain
    call ovkCreateDomain(Domain, NumDims=2, NumGrids=2, StatusLogFile=OUTPUT_UNIT, &
      ErrorLogFile=ERROR_UNIT)

    !===============
    ! Grid creation
    !===============

    ! Initialize grid data structure for background grid
    ! Set geometry type as a hint for potential performance improvements
    call ovkCreateGrid(Domain, 1, NumPoints=NumPointsBackground, &
      GeometryType=OVK_GRID_GEOMETRY_CARTESIAN)
    call ovkEditGrid(Domain, 1, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBackground(1),1:NumPointsBackground(2),1) = CoordsBackground(:,:,1)
    Y%values(1:NumPointsBackground(1),1:NumPointsBackground(2),1) = CoordsBackground(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for cylinder grid
    ! Periodic in the angular direction, with the last set of points being equal to the first
    call ovkCreateGrid(Domain, 2, NumPoints=NumPointsCylinder, Periodic=[.false.,.true.], &
      PeriodicStorage=OVK_OVERLAP_PERIODIC)
    call ovkEditGrid(Domain, 2, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsCylinder(1),1:NumPointsCylinder(2),1) = CoordsCylinder(:,:,1)
    Y%values(1:NumPointsCylinder(1),1:NumPointsCylinder(2),1) = CoordsCylinder(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    ! Inner radial boundary
    call ovkEditGridState(Grid, State)
    State%values(1,:,1) = OVK_DOMAIN_BOUNDARY_POINT
    call ovkReleaseGridState(Grid, State)

    call ovkReleaseGrid(Domain, Grid)

    !==================
    ! Overset assembly
    !==================

    AssemblyOptions = ovk_assembly_options_(2, 2)

    ! Indicate which grids can intersect
    call ovkSetAssemblyOptionOverlappable(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)

    ! Automatically define boundaries in non-overlapping regions
    call ovkSetAssemblyOptionInferBoundaries(AssemblyOptions, OVK_ALL_GRIDS, .true.)

    ! Indicate which grids can cut each other
    call ovkSetAssemblyOptionCutBoundaryHoles(AssemblyOptions, 2, 1, .true.)

    ! Indicate how to treat overlap between grids
    call ovkSetAssemblyOptionOccludes(AssemblyOptions, 2, 1, OVK_OCCLUDES_ALL)

    ! Retain some extra overlap between grids
    call ovkSetAssemblyOptionEdgePadding(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 2)
    call ovkSetAssemblyOptionEdgeSmoothing(AssemblyOptions, OVK_ALL_GRIDS, 2)

    ! Indicate which grids can communicate and how
    call ovkSetAssemblyOptionConnectionType(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_CUBIC)
    call ovkSetAssemblyOptionFringeSize(AssemblyOptions, OVK_ALL_GRIDS, 2)
    call ovkSetAssemblyOptionMinimizeOverlap(AssemblyOptions, 2, 1, .true.)

    call ovkAssemble(Domain, AssemblyOptions)

    ! Update background IBlank with holes (used later if generating remap)
    call ovkGetGrid(Domain, 1, Grid)
    call ovkFilterGridState(Grid, OVK_STATE_GRID, OVK_ALL, Mask)
    IBlankBackground = merge(IBlankBackground, 0, Mask%values(:,:,1))

    !========
    ! Output
    !========

    NumPointsAll(:,1) = NumPointsBackground
    NumPointsAll(:,2) = NumPointsCylinder

    ! Write a PLOT3D grid file
    call ovkCreateP3D(GridFile, "circle_cylinder.xyz", NumDims=2, NumGrids=2, &
      NumPointsAll=NumPointsAll, WithIBlank=.true., StatusLogFile=OUTPUT_UNIT, &
      ErrorLogFile=ERROR_UNIT)

    do n = 1, 2

      call ovkGetGrid(Domain, n, Grid)
      call ovkGetGridCart(Grid, Cart)
      call ovkGetGridCoords(Grid, 1, X)
      call ovkGetGridCoords(Grid, 2, Y)

      ! Use IBlank data to visualize status of grid points
      ! IBlank == 1 => Normal
      IBlank = ovk_field_int_(Cart, 1)

      ! IBlank == 0 => Hole
      call ovkFilterGridState(Grid, OVK_STATE_GRID, OVK_NONE, Mask)
      IBlank%values = merge(0, IBlank%values, Mask%values)

      ! IBlank == -N => Receives from grid N
      do m = 1, 2
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

      call ovkWriteP3D(GridFile, n, X, Y, IBlank)

    end do

    ! Finalize the PLOT3D file
    call ovkCloseP3D(GridFile)

    call ovkDestroyDomain(Domain)

  end subroutine GenerateCylinder

  subroutine GenerateBlock()

    integer :: i, j, m, n
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

    !=====================
    ! Data initialization
    !=====================

    NumPointsBlock = [NumPoints/2,NumPoints/2]

    allocate(CoordsEastBlock(NumPointsBlock(1),NumPointsBlock(2),2))
    allocate(CoordsNorthBlock(NumPointsBlock(1),NumPointsBlock(2),2))
    allocate(CoordsWestBlock(NumPointsBlock(1),NumPointsBlock(2),2))
    allocate(CoordsSouthBlock(NumPointsBlock(1),NumPointsBlock(2),2))

    ! Generate coordinates for east block
    do j = 1, NumPointsBlock(2)
      do i = 1, NumPointsBlock(1)
        U = real(i-1,kind=rk)/real(NumPointsBlock(1)-1,kind=rk)
        V = real(j-1,kind=rk)/real(NumPointsBlock(2)-1,kind=rk)
        Theta = 0.5_rk * Pi * (V-0.5_rk)
        RMin = 0.3_rk
        RMax = 1._rk/cos(Theta)
        Radius = (RMin+0.5_rk) * ((RMax+0.5_rk)/(RMin+0.5_rk))**U-0.5_rk
        CoordsEastBlock(i,j,1) = Radius * cos(Theta)
        CoordsEastBlock(i,j,2) = Radius * sin(Theta)
      end do
    end do

    ! Generate coordinates for north block
    do j = 1, NumPointsBlock(2)
      do i = 1, NumPointsBlock(1)
        U = real(i-1,kind=rk)/real(NumPointsBlock(1)-1,kind=rk)
        V = real(j-1,kind=rk)/real(NumPointsBlock(2)-1,kind=rk)
        Theta = 0.5_rk * Pi * (V+0.5_rk)
        RMin = 0.3_rk
        RMax = 1._rk/sin(Theta)
        Radius = (RMin+0.5_rk) * ((RMax+0.5_rk)/(RMin+0.5_rk))**U-0.5_rk
        CoordsNorthBlock(i,j,1) = Radius * cos(Theta)
        CoordsNorthBlock(i,j,2) = Radius * sin(Theta)
      end do
    end do

    ! Generate coordinates for west block
    do j = 1, NumPointsBlock(2)
      do i = 1, NumPointsBlock(1)
        U = real(i-1,kind=rk)/real(NumPointsBlock(1)-1,kind=rk)
        V = real(j-1,kind=rk)/real(NumPointsBlock(2)-1,kind=rk)
        Theta = 0.5_rk * Pi * (V+1.5_rk)
        RMin = 0.3_rk
        RMax = -1._rk/cos(Theta)
        Radius = (RMin+0.5_rk) * ((RMax+0.5_rk)/(RMin+0.5_rk))**U-0.5_rk
        CoordsWestBlock(i,j,1) = Radius * cos(Theta)
        CoordsWestBlock(i,j,2) = Radius * sin(Theta)
      end do
    end do

    ! Generate coordinates for south block
    do j = 1, NumPointsBlock(2)
      do i = 1, NumPointsBlock(1)
        U = real(i-1,kind=rk)/real(NumPointsBlock(1)-1,kind=rk)
        V = real(j-1,kind=rk)/real(NumPointsBlock(2)-1,kind=rk)
        Theta = 0.5_rk * Pi * (V+2.5_rk)
        RMin = 0.3_rk
        RMax = -1._rk/sin(Theta)
        Radius = (RMin+0.5_rk) * ((RMax+0.5_rk)/(RMin+0.5_rk))**U-0.5_rk
        CoordsSouthBlock(i,j,1) = Radius * cos(Theta)
        CoordsSouthBlock(i,j,2) = Radius * sin(Theta)
      end do
    end do

    !========================
    ! Overkit initialization
    !========================

    ! Initialize the domain
    call ovkCreateDomain(Domain, NumDims=2, NumGrids=4, StatusLogFile=OUTPUT_UNIT, &
      ErrorLogFile=ERROR_UNIT)

    !===============
    ! Grid creation
    !===============

    ! Initialize grid data structure for east block grid
    call ovkCreateGrid(Domain, 1, NumPoints=NumPointsBlock)
    call ovkEditGrid(Domain, 1, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsEastBlock(:,:,1)
    Y%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsEastBlock(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkEditGridState(Grid, State)
    ! Inner radial boundary
    State%values(1,:,1) = OVK_DOMAIN_BOUNDARY_POINT
    ! Outer box boundary
    State%values(NumPointsBlock(1),:,1) = OVK_DOMAIN_BOUNDARY_POINT
    call ovkReleaseGridState(Grid, State)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for north block grid
    call ovkCreateGrid(Domain, 2, NumPoints=NumPointsBlock)
    call ovkEditGrid(Domain, 2, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsNorthBlock(:,:,1)
    Y%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsNorthBlock(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkEditGridState(Grid, State)
    ! Inner radial boundary
    State%values(1,:,1) = OVK_DOMAIN_BOUNDARY_POINT
    ! Outer box boundary
    State%values(NumPointsBlock(1),:,1) = OVK_DOMAIN_BOUNDARY_POINT
    call ovkReleaseGridState(Grid, State)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for west block grid
    call ovkCreateGrid(Domain, 3, NumPoints=NumPointsBlock)
    call ovkEditGrid(Domain, 3, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsWestBlock(:,:,1)
    Y%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsWestBlock(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkEditGridState(Grid, State)
    ! Inner radial boundary
    State%values(1,:,1) = OVK_DOMAIN_BOUNDARY_POINT
    ! Outer box boundary
    State%values(NumPointsBlock(1),:,1) = OVK_DOMAIN_BOUNDARY_POINT
    call ovkReleaseGridState(Grid, State)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for south block grid
    call ovkCreateGrid(Domain, 4, NumPoints=NumPointsBlock)
    call ovkEditGrid(Domain, 4, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsSouthBlock(:,:,1)
    Y%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsSouthBlock(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkEditGridState(Grid, State)
    ! Inner radial boundary
    State%values(1,:,1) = OVK_DOMAIN_BOUNDARY_POINT
    ! Outer box boundary
    State%values(NumPointsBlock(1),:,1) = OVK_DOMAIN_BOUNDARY_POINT
    call ovkReleaseGridState(Grid, State)

    call ovkReleaseGrid(Domain, Grid)

    !==================
    ! Overset assembly
    !==================

    AssemblyOptions = ovk_assembly_options_(2, 4)

    ! Indicate which grids can intersect
    call ovkSetAssemblyOptionOverlappable(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
    call ovkSetAssemblyOptionOverlapTolerance(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 1e-6_rk)

    ! Indicate which grids can communicate and how
    call ovkSetAssemblyOptionConnectionType(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_CUBIC)
    call ovkSetAssemblyOptionFringeSize(AssemblyOptions, OVK_ALL_GRIDS, 1)

    call ovkAssemble(Domain, AssemblyOptions)

    !========
    ! Output
    !========

    NumPointsAll(:,1) = NumPointsBlock
    NumPointsAll(:,2) = NumPointsBlock
    NumPointsAll(:,3) = NumPointsBlock
    NumPointsAll(:,4) = NumPointsBlock

    ! Write a PLOT3D grid file
    call ovkCreateP3D(GridFile, "circle_block.xyz", NumDims=2, NumGrids=4, &
      NumPointsAll=NumPointsAll, WithIBlank=.true., StatusLogFile=OUTPUT_UNIT, &
      ErrorLogFile=ERROR_UNIT)

    do n = 1, 4

      call ovkGetGrid(Domain, n, Grid)
      call ovkGetGridCart(Grid, Cart)
      call ovkGetGridCoords(Grid, 1, X)
      call ovkGetGridCoords(Grid, 2, Y)

      ! Use IBlank data to visualize status of grid points
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

      call ovkWriteP3D(GridFile, n, X, Y, IBlank)

    end do

    ! Finalize the PLOT3D file
    call ovkCloseP3D(GridFile)

    call ovkDestroyDomain(Domain)

  end subroutine GenerateBlock

  subroutine GenerateRemap()

    integer :: i, j, m, n
    type(ovk_domain) :: Domain
    type(ovk_grid), pointer :: Grid
    type(ovk_field_real), pointer :: X, Y
    type(ovk_field_int), pointer :: State
    type(ovk_assembly_options) :: AssemblyOptions
    integer, dimension(2,6) :: NumPointsAll
    type(ovk_plot3d_grid_file) :: GridFile
    type(ovk_cart) :: Cart
    type(ovk_connectivity), pointer :: Connectivity
    integer, dimension(:,:), pointer :: ReceiverPoints
    integer, dimension(MAX_ND) :: Point
    type(ovk_field_logical) :: Mask
    type(ovk_field_int) :: IBlank

    !========================
    ! Overkit initialization
    !========================

    ! Initialize the domain
    ! Combined domain with both cylinder (1-2) and block (3-6) grids
    call ovkCreateDomain(Domain, NumDims=2, NumGrids=6, StatusLogFile=OUTPUT_UNIT, &
      ErrorLogFile=ERROR_UNIT)

    !===============
    ! Grid creation
    !===============

    ! Initialize grid data structure for background grid
    ! Set geometry type as a hint for potential performance improvements
    call ovkCreateGrid(Domain, 1, NumPoints=NumPointsBackground, &
      GeometryType=OVK_GRID_GEOMETRY_CARTESIAN)
    call ovkEditGrid(Domain, 1, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBackground(1),1:NumPointsBackground(2),1) = CoordsBackground(:,:,1)
    Y%values(1:NumPointsBackground(1),1:NumPointsBackground(2),1) = CoordsBackground(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkEditGridState(Grid, State)
    ! Update with hole-cutting information from IBlank
    do j = 1, NumPointsBackground(2)
      do i = 1, NumPointsBackground(1)
        if (IBlankBackground(i,j) == 0) then
          State%values(i,j,1) = OVK_EXTERIOR_POINT
        end if
      end do
    end do
    call ovkReleaseGridState(Grid, State)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for cylinder grid
    ! Periodic in the angular direction, with the last set of points being equal to the first
    call ovkCreateGrid(Domain, 2, NumPoints=NumPointsCylinder, Periodic=[.false.,.true.], &
      PeriodicStorage=OVK_OVERLAP_PERIODIC)
    call ovkEditGrid(Domain, 2, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsCylinder(1),1:NumPointsCylinder(2),1) = CoordsCylinder(:,:,1)
    Y%values(1:NumPointsCylinder(1),1:NumPointsCylinder(2),1) = CoordsCylinder(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for east block grid
    call ovkCreateGrid(Domain, 3, NumPoints=NumPointsBlock)
    call ovkEditGrid(Domain, 3, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsEastBlock(:,:,1)
    Y%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsEastBlock(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for north block grid
    call ovkCreateGrid(Domain, 4, NumPoints=NumPointsBlock)
    call ovkEditGrid(Domain, 4, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsNorthBlock(:,:,1)
    Y%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsNorthBlock(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for west block grid
    call ovkCreateGrid(Domain, 5, NumPoints=NumPointsBlock)
    call ovkEditGrid(Domain, 5, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsWestBlock(:,:,1)
    Y%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsWestBlock(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkReleaseGrid(Domain, Grid)

    ! Initialize grid data structure for south block grid
    call ovkCreateGrid(Domain, 6, NumPoints=NumPointsBlock)
    call ovkEditGrid(Domain, 6, Grid)

    ! Fill in coordinates
    call ovkEditGridCoords(Grid, 1, X)
    call ovkEditGridCoords(Grid, 2, Y)
    X%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsSouthBlock(:,:,1)
    Y%values(1:NumPointsBlock(1),1:NumPointsBlock(2),1) = CoordsSouthBlock(:,:,2)
    call ovkReleaseGridCoords(Grid, X)
    call ovkReleaseGridCoords(Grid, Y)

    call ovkReleaseGrid(Domain, Grid)

    !==================
    ! Overset assembly
    !==================

    AssemblyOptions = ovk_assembly_options_(2, 6)

    ! Indicate which grids can intersect
    ! Only need to compute overlap of block grids by cylinder grids
    do n = 3, 6
      do m = 1, 2
        call ovkSetAssemblyOptionOverlappable(AssemblyOptions, m, n, .true.)
      end do
    end do
    call ovkSetAssemblyOptionOverlapTolerance(AssemblyOptions, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 1e-6_rk)

    ! Indicate how to treat overlap between grids
    ! Cylinder grids should occlude block grids
    do n = 3, 6
      do m = 1, 2
        call ovkSetAssemblyOptionOccludes(AssemblyOptions, m, n, OVK_OCCLUDES_ALL)
      end do
    end do

    ! Indicate which grids can communicate and how
    ! Full grid interpolation from cylinder grids to block grids
    do n = 3, 6
      do m = 1, 2
        call ovkSetAssemblyOptionConnectionType(AssemblyOptions, m, n, OVK_CONNECTION_CUBIC)
      end do
    end do

    call ovkAssemble(Domain, AssemblyOptions)

    !========
    ! Output
    !========

    NumPointsAll(:,1) = NumPointsBackground
    NumPointsAll(:,2) = NumPointsCylinder
    NumPointsAll(:,3) = NumPointsBlock
    NumPointsAll(:,4) = NumPointsBlock
    NumPointsAll(:,5) = NumPointsBlock
    NumPointsAll(:,6) = NumPointsBlock

    ! Write a PLOT3D grid file
    call ovkCreateP3D(GridFile, "circle_remap.xyz", NumDims=2, NumGrids=6, &
      NumPointsAll=NumPointsAll, WithIBlank=.true., StatusLogFile=OUTPUT_UNIT, &
      ErrorLogFile=ERROR_UNIT)

    do n = 1, 6

      call ovkGetGrid(Domain, n, Grid)
      call ovkGetGridCart(Grid, Cart)
      call ovkGetGridCoords(Grid, 1, X)
      call ovkGetGridCoords(Grid, 2, Y)

      ! Use IBlank data to visualize status of grid points
      ! IBlank == 1 => Normal
      IBlank = ovk_field_int_(Cart, 1)

      ! IBlank == 0 => Hole
      call ovkFilterGridState(Grid, OVK_STATE_GRID, OVK_NONE, Mask)
      IBlank%values = merge(0, IBlank%values, Mask%values)

      ! IBlank == -N => Receives from grid N
      do m = 1, 6
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

      call ovkWriteP3D(GridFile, n, X, Y, IBlank)

    end do

    ! Finalize the PLOT3D file
    call ovkCloseP3D(GridFile)

    call ovkDestroyDomain(Domain)

  end subroutine GenerateRemap

end program Circle
