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
  type(ovk_assembler) :: Assembler
  type(ovk_assembler_properties), pointer :: AssemblerProperties
  type(ovk_assembler_graph), pointer :: Graph
  type(ovk_domain), pointer :: Domain
  type(ovk_grid), pointer :: Grid
  real(rk), dimension(:,:,:,:), allocatable :: XYZ
  type(ovk_field_real), pointer :: Coords
  type(ovk_field_logical), pointer :: BoundaryMask
  integer, dimension(MAX_ND) :: Point
  real(rk) :: U
  real(rk) :: R
  real(rk) :: BumpHeight
  real(rk) :: MinHeight, MaxHeight
  real(rk) :: Shift
  integer, dimension(MAX_ND,2) :: NumPointsAll
  type(ovk_plot3d_grid_file) :: GridFile
  type(ovk_connectivity), pointer :: Connectivity
  type(ovk_interp), pointer :: InterpData
  type(ovk_cart) :: Cart
  type(ovk_cart) :: CartOverlap
  type(ovk_field_logical), pointer :: GridMask
  type(ovk_field_logical), pointer :: ReceiverMask
  type(ovk_field_logical), pointer :: OrphanMask
  type(ovk_field_int), pointer :: DonorGridIDs
  type(ovk_field_real) :: XOverlap, YOverlap, ZOverlap
  type(ovk_field_int) :: IBlank
  type(ovk_field_int) :: IBlankOverlap

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
    write (ERROR_UNIT, '(3a)') "ERROR: Invalid grid dimension ", NumDims, "."
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

  NumPointsBackground(:NumDims-1) = N
  NumPointsBackground(NumDims) = N/2
  NumPointsBackground(NumDims+1:) = 1

  NumPointsBump(:NumDims-1) = N
  NumPointsBump(NumDims) = (3*N)/4
  NumPointsBump(NumDims+1:) = 1

  Length(:NumDims-1) = 2._rk
  Length(NumDims) = 1._rk

  ! Initialize the problem
  call ovkCreateAssembler(Assembler, NumDims=NumDims, NumGrids=2)

  ! Enable verbose command line output
  call ovkEditAssemblerProperties(Assembler, AssemblerProperties)
  call ovkSetAssemblerPropertyVerbose(AssemblerProperties, .true.)
  call ovkReleaseAssemblerProperties(Assembler, AssemblerProperties)

  ! Indicate which grids can intersect, cut, communicate, etc.
  call ovkEditAssemblerGraph(Assembler, Graph)
  call ovkSetAssemblerGraphOverlap(Graph, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
  call ovkSetAssemblerGraphOverlapTolerance(Graph, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 0.1_rk)
  call ovkSetAssemblerGraphBoundaryHoleCutting(Graph, 2, 1, .true.)
  call ovkSetAssemblerGraphOverlapHoleCutting(Graph, 2, 1, .true.)
  call ovkSetAssemblerGraphConnectionType(Graph, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_FRINGE)
  call ovkSetAssemblerGraphDisjointConnection(Graph, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
  call ovkSetAssemblerGraphInterpScheme(Graph, OVK_ALL_GRIDS, OVK_ALL_GRIDS, InterpScheme)
  call ovkSetAssemblerGraphFringeSize(Graph, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 2)
  call ovkSetAssemblerGraphFringePadding(Graph, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 6)
  call ovkReleaseAssemblerGraph(Assembler, Graph)

  ! Set up the domain
  call ovkEditAssemblerDomain(Assembler, Domain)

  !=================
  ! Background grid
  !=================

  ! Initialize grid data structure for background grid
  ! Set geometry type as a hint for potential performance improvements
  call ovkCreateDomainGrid(Domain, 1, NumPoints=NumPointsBackground, &
    GeometryType=OVK_GRID_GEOMETRY_CARTESIAN)
  call ovkEditDomainGrid(Domain, 1, Grid)

  ! Generate coordinates for background grid
  allocate(XYZ(NumPointsBackground(1),NumPointsBackground(2),NumPointsBackground(3),NumDims))
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

  ! Outer edge boundaries on background grid
  call ovkEditGridBoundaryMask(Grid, BoundaryMask)
  select case (NumDims)
  case (2)
    BoundaryMask%values(1,:,1) = .true.
    BoundaryMask%values(NumPointsBackground(1),:,1) = .true.
    BoundaryMask%values(:,1,1) = .true.
    BoundaryMask%values(:,NumPointsBackground(2),1) = .true.
  case (3)
    BoundaryMask%values(1,:,:) = .true.
    BoundaryMask%values(NumPointsBackground(1),:,:) = .true.
    BoundaryMask%values(:,1,:) = .true.
    BoundaryMask%values(:,NumPointsBackground(2),:) = .true.
    BoundaryMask%values(:,:,1) = .true.
    BoundaryMask%values(:,:,NumPointsBackground(3)) = .true.
  end select
  call ovkReleaseGridBoundaryMask(Grid, BoundaryMask)

  deallocate(XYZ)

  call ovkReleaseDomainGrid(Domain, Grid)

  !===========
  ! Bump grid
  !===========

  ! Initialize grid data structure for bump grid
  call ovkCreateDomainGrid(Domain, 2, NumPoints=NumPointsBump)
  call ovkEditDomainGrid(Domain, 2, Grid)

  ! Generate coordinates for bump grid
  allocate(XYZ(NumPointsBump(1),NumPointsBump(2),NumPointsBump(3),NumDims))
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

  ! Lower edge boundary on bump grid
  call ovkEditGridBoundaryMask(Grid, BoundaryMask)
  select case (NumDims)
  case (2)
    BoundaryMask%values(:,1,1) = .true.
  case (3)
    BoundaryMask%values(:,:,1) = .true.
  end select
  call ovkReleaseGridBoundaryMask(Grid, BoundaryMask)

  deallocate(XYZ)

  call ovkReleaseDomainGrid(Domain, Grid)

  call ovkReleaseAssemblerDomain(Assembler, Domain)

  !==================
  ! Overset assembly
  !==================

  ! Perform overset assembly
  call ovkAssemble(Assembler)
!   call ovkFindOverlap(Assembler)
!   call ovkCutHoles(Assembler)
!   call ovkConnectGrids(Assembler)

  !========
  ! Output
  !========

  NumPointsAll(:,1) = NumPointsBackground
  NumPointsAll(:,2) = NumPointsBump

  call ovkGetAssemblerDomain(Assembler, Domain)
  call ovkGetAssemblerConnectivity(Assembler, Connectivity)

  ! Write a PLOT3D grid file with IBlank to visualize the result
  call ovkCreateP3D(GridFile, "grid.xyz", NumDims=NumDims, NumGrids=2, NumPointsAll=NumPointsAll, &
    WithIBlank=.true.)

  do m = 1, 2

    call ovkGetDomainGrid(Domain, m, Grid)
    call ovkGetConnectivityInterpData(Connectivity, m, InterpData)

    ! At the moment Overkit converts everything to no-overlap periodic internally, so
    ! we will need to convert back
    call ovkGetGridCart(Grid, Cart)
    CartOverlap = ovkCartConvertPeriodicStorage(Cart, OVK_OVERLAP_PERIODIC)

    call ovkExportGridCoords(Grid, 1, CartOverlap, XOverlap)
    call ovkExportGridCoords(Grid, 2, CartOverlap, YOverlap)
    if (NumDims == 3) then
      call ovkExportGridCoords(Grid, 3, CartOverlap, ZOverlap)
    end if

    call ovkGetGridMask(Grid, GridMask)
    call ovkGetInterpDataReceiverMask(InterpData, ReceiverMask)
    call ovkGetInterpDataOrphanMask(InterpData, OrphanMask)
    call ovkGetInterpDataDonorGridIDs(InterpData, DonorGridIDs)

    ! IBlank values are set as follows:
    !   1 => normal point
    !   0 => hole
    !  -N => receives from grid N
    !   7 => orphan
    IBlank = ovk_field_int_(Cart)
    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          if (GridMask%values(i,j,k)) then
            if (ReceiverMask%values(i,j,k)) then
              IBlank%values(i,j,k) = -DonorGridIDs%values(i,j,k)
            else if (OrphanMask%values(i,j,k)) then
              IBlank%values(i,j,k) = 7
            else
              IBlank%values(i,j,k) = 1
            end if
          else
            IBlank%values(i,j,k) = 0
          end if
        end do
      end do
    end do
    call ovkExportField(IBlank, CartOverlap, IBlankOverlap)

    if (NumDims == 2) then
      call ovkWriteP3D(GridFile, m, XOverlap, YOverlap, IBlankOverlap)
    else
      call ovkWriteP3D(GridFile, m, XOverlap, YOverlap, ZOverlap, IBlankOverlap)
    end if

  end do

  ! Finalize the PLOT3D file
  call ovkCloseP3D(GridFile)

  call ovkDestroyAssembler(Assembler)

end program Bump
