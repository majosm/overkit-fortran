! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

program Bump

  use Overkit
  use ovsCommandLine
  use ovsGlobal
  implicit none

  integer :: i, j, k, d, m
  character(len=256), dimension(:), allocatable :: RawArguments
  character(len=256) :: Usage
  character(len=256) :: Description
  type(t_cmd_opt), dimension(3) :: Options
  integer :: N
  integer :: NumDims
  character(len=32) :: InterpSchemeString
  integer :: InterpScheme
  integer, dimension(MAX_ND) :: NumPointsBackground
  integer, dimension(MAX_ND) :: NumPointsBump
  real(rk), dimension(MAX_ND) :: Length
  type(ovk_domain) :: Domain
  type(ovk_domain_properties), pointer :: Properties
  type(ovk_grid), pointer :: Grid
  real(rk), dimension(:,:,:,:), allocatable :: XYZ
  type(ovk_field_real), pointer :: Coords
  type(ovk_field_int), pointer :: State
  integer, dimension(MAX_ND) :: Point
  real(rk) :: U
  real(rk) :: R
  real(rk) :: BumpHeight
  real(rk) :: MinHeight, MaxHeight
  real(rk) :: Shift
  integer, dimension(MAX_ND,2) :: NumPointsAll
  type(ovk_plot3d_grid_file) :: GridFile
  type(ovk_cart) :: Cart
  integer :: ConnectionType
  type(ovk_connectivity), pointer :: Connectivity
  integer, dimension(:,:), pointer :: ReceiverPoints
  type(ovk_field_logical) :: Mask
  type(ovk_field_real), pointer :: X, Y, Z
  type(ovk_field_int) :: IBlank

  allocate(RawArguments(command_argument_count()))
  do i = 1, size(RawArguments)
    call get_command_argument(i, RawArguments(i))
  end do

  Usage = "Bump [<options> ...]"
  Description = "Generates an overset mesh representing a domain with a wall that has a bump in it."
  Options(1) = t_cmd_opt_("dimension", "d", CMD_OPT_INTEGER, "Grid dimension (2 or 3) [ Default: 2 ]")
  Options(2) = t_cmd_opt_("size", "N", CMD_OPT_INTEGER, "Size of background grid in " // &
    "wall-parallel direction(s) (other directions and grids are proportional) [ Default: 81 ]")
  Options(3) = t_cmd_opt_("interp", "i", CMD_OPT_STRING, "Interpolation scheme (linear or " // &
    "cubic) [ Default: cubic ]")

  call ParseArguments(RawArguments, Usage=Usage, Description=Description, Options=Options)

  call GetOptionValue(Options(1), NumDims, 2)
  if (NumDims /= 2 .and. NumDims /= 3) then
    write (ERROR_UNIT, '(a,i0,a)') "ERROR: Invalid grid dimension ", NumDims, "."
    stop 1
  end if

  call GetOptionValue(Options(2), N, 81)

  call GetOptionValue(Options(3), InterpSchemeString, "cubic")
  select case (InterpSchemeString)
  case ("linear")
    InterpScheme = OVK_INTERP_LINEAR
  case ("cubic")
    InterpScheme = OVK_INTERP_CUBIC
  case default
    write (ERROR_UNIT, '(3a)') "ERROR: Unrecognized interpolation scheme '", &
      trim(InterpSchemeString), "'."
    stop 1
  end select

  ! Initialize the problem
  call ovkCreateDomain(Domain, NumDims=NumDims, NumGrids=2, Verbose=.true.)

  call ovkEditDomainProperties(Domain, Properties)

  ! Indicate which grids can intersect
  call ovkSetDomainPropertyOverlappable(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
  call ovkSetDomainPropertyOverlapTolerance(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 0.1_rk)

  ! Automatically define boundaries in non-overlapping regions
  call ovkSetDomainPropertyInferBoundaries(Properties, OVK_ALL_GRIDS, .true.)

  ! Indicate which grids can cut each other
  call ovkSetDomainPropertyBoundaryHoleCutting(Properties, 2, 1, .true.)

  ! Indicate which grids can communicate and how
  call ovkSetDomainPropertyConnectionType(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_FRINGE)
  call ovkSetDomainPropertyInterpScheme(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, InterpScheme)
  call ovkSetDomainPropertyFringeSize(Properties, OVK_ALL_GRIDS, 2)
  call ovkSetDomainPropertyEdgePadding(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 6)
  call ovkSetDomainPropertyOverlapMinimization(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)

  call ovkReleaseDomainProperties(Domain, Properties)

  !=================
  ! Background grid
  !=================

  NumPointsBackground(:NumDims-1) = N
  NumPointsBackground(NumDims) = N/2
  NumPointsBackground(NumDims+1:) = 1

  Length(:NumDims-1) = 2._rk
  Length(NumDims) = 1._rk

  ! Initialize grid data structure for background grid
  ! Set geometry type as a hint for potential performance improvements
  call ovkCreateGrid(Domain, 1, NumPoints=NumPointsBackground, &
    GeometryType=OVK_GRID_GEOMETRY_CARTESIAN)
  call ovkEditGrid(Domain, 1, Grid)

  allocate(XYZ(NumPointsBackground(1),NumPointsBackground(2),NumPointsBackground(3),NumDims))

  ! Generate coordinates for background grid
  do k = 1, NumPointsBackground(3)
    do j = 1, NumPointsBackground(2)
      do i = 1, NumPointsBackground(1)
        Point = [i,j,k]
        do d = 1, NumDims-1
          U = real(Point(d)-1,kind=rk)/real(NumPointsBackground(d)-1,kind=rk)
          XYZ(i,j,k,d) = Length(d) * (U-0.5_rk)
        end do
          U = real(Point(NumDims)-1,kind=rk)/real(NumPointsBackground(NumDims)-1,kind=rk)
          XYZ(i,j,k,NumDims) = Length(NumDims) * U
      end do
    end do
  end do

  ! Feed them to Overkit
  do d = 1, NumDims
    call ovkEditGridCoords(Grid, d, Coords)
    Coords%values = XYZ(:,:,:,d)
    call ovkReleaseGridCoords(Grid, Coords)
  end do

  deallocate(XYZ)

  ! Lower edge boundary on background grid
  call ovkEditGridState(Grid, State)
  select case (NumDims)
  case (2)
    State%values(:,1,1) = OVK_DOMAIN_BOUNDARY_POINT
  case (3)
    State%values(:,:,1) = OVK_DOMAIN_BOUNDARY_POINT
  end select
  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  !===========
  ! Bump grid
  !===========

  NumPointsBump(:NumDims-1) = N
  NumPointsBump(NumDims) = (3*N)/4
  NumPointsBump(NumDims+1:) = 1

  ! Initialize grid data structure for bump grid
  call ovkCreateGrid(Domain, 2, NumPoints=NumPointsBump)
  call ovkEditGrid(Domain, 2, Grid)

  allocate(XYZ(NumPointsBump(1),NumPointsBump(2),NumPointsBump(3),NumDims))

  ! Generate coordinates for bump grid
  do k = 1, NumPointsBump(3)
    do j = 1, NumPointsBump(2)
      do i = 1, NumPointsBump(1)
        Point = [i,j,k]
        do d = 1, NumDims-1
          U = real(Point(d)-1,kind=rk)/real(NumPointsBump(d)-1,kind=rk)
          XYZ(i,j,k,d) = 0.65_rk*Length(d) * (U-0.5_rk)
        end do
        R = sqrt(sum(XYZ(i,j,k,:NumDims-1)**2))
        if (R <= 0.499_rk) then
          BumpHeight = 0.25_rk*exp(1._rk-1._rk/(1._rk-(2._rk*R)**2))
        else
          BumpHeight = 0._rk
        end if
        U = real(Point(NumDims)-1,kind=rk)/real(NumPointsBump(NumDims)-1,kind=rk)
        MinHeight = BumpHeight
        MaxHeight = 0.5_rk * Length(NumDims)
        Shift = BumpHeight-0.1_rk
        XYZ(i,j,k,NumDims) = (MinHeight-Shift)*((MaxHeight-Shift)/(MinHeight-Shift))**U + Shift
      end do
    end do
  end do

  ! Feed them to Overkit
  do d = 1, NumDims
    call ovkEditGridCoords(Grid, d, Coords)
    Coords%values = XYZ(:,:,:,d)
    call ovkReleaseGridCoords(Grid, Coords)
  end do

  deallocate(XYZ)

  ! Lower edge boundary on bump grid
  call ovkEditGridState(Grid, State)
  select case (NumDims)
  case (2)
    State%values(:,1,1) = OVK_DOMAIN_BOUNDARY_POINT
  case (3)
    State%values(:,:,1) = OVK_DOMAIN_BOUNDARY_POINT
  end select
  call ovkReleaseGridState(Grid, State)

  call ovkReleaseGrid(Domain, Grid)

  !==================
  ! Overset assembly
  !==================

  call ovkAssemble(Domain)

  !========
  ! Output
  !========

  NumPointsAll(:,1) = NumPointsBackground
  NumPointsAll(:,2) = NumPointsBump

  ! Write a PLOT3D grid file with IBlank to visualize the result
  call ovkCreateP3D(GridFile, "bump.xyz", NumDims=NumDims, NumGrids=2, NumPointsAll=NumPointsAll, &
    WithIBlank=.true.)

  call ovkGetDomainProperties(Domain, Properties)

  do n = 1, 2

    call ovkGetGrid(Domain, n, Grid)
    call ovkGetGridCart(Grid, Cart)
    call ovkGetGridCoords(Grid, 1, X)
    call ovkGetGridCoords(Grid, 2, Y)
    if (NumDims == 3) then
      call ovkGetGridCoords(Grid, 3, Z)
    end if

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

    if (NumDims == 2) then
      call ovkWriteP3D(GridFile, n, X, Y, IBlank)
    else
      call ovkWriteP3D(GridFile, n, X, Y, Z, IBlank)
    end if

  end do

  ! Finalize the PLOT3D file
  call ovkCloseP3D(GridFile)

  call ovkDestroyDomain(Domain)

end program Bump
