! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

program Blobs

  use Overkit
  use ovsCommandLine
  use ovsGlobal
  implicit none

  integer :: i, j, m
  character(len=256), dimension(:), allocatable :: RawArguments
  character(len=256) :: Usage
  character(len=256) :: Description
  type(t_cmd_opt), dimension(1) :: Options
  integer :: N
  integer, dimension(2) :: NumPointsBackground
  integer, dimension(2) :: NumPointsBlob
  real(rk), dimension(2) :: Length
  real(rk) :: SeparationScale
  type(ovk_assembler) :: Assembler
  type(ovk_assembler_properties), pointer :: Properties
  type(ovk_domain), pointer :: Domain
  type(ovk_grid), pointer :: Grid
  type(ovk_field_real), pointer :: X, Y
  type(ovk_field_logical), pointer :: BoundaryMask
  real(rk) :: U, V
  real(rk) :: RMin, RMax
  real(rk) :: Radius
  real(rk) :: Theta
  integer, dimension(2,4) :: NumPointsAll
  type(ovk_plot3d_grid_file) :: GridFile
  type(ovk_connectivity), pointer :: Connectivity
  type(ovk_interp), pointer :: InterpData
  type(ovk_cart) :: Cart
  type(ovk_cart) :: CartOverlap
  type(ovk_field_logical), pointer :: GridMask
  type(ovk_field_logical), pointer :: ReceiverMask
  type(ovk_field_logical), pointer :: OrphanMask
  type(ovk_field_int), pointer :: DonorGridIDs
  type(ovk_field_real) :: XOverlap, YOverlap
  type(ovk_field_int) :: IBlank
  type(ovk_field_int) :: IBlankOverlap

  allocate(RawArguments(command_argument_count()))
  do i = 1, size(RawArguments)
    call get_command_argument(i, RawArguments(i))
  end do

  Usage = "Blobs [<options> ...]"
  Description = "Generates an overset mesh representing a rectangular domain with several " // &
    "blob-like holes in it."
  Options(1) = t_cmd_opt_("size", "N", CMD_OPT_INTEGER, "Size of background grid in each " // &
    "direction (other grids are proportional) [ Default: 81 ]")

  call ParseArguments(RawArguments, Usage=Usage, Description=Description, Options=Options)

  call GetOptionValue(Options(1), N, 81)

  NumPointsBackground = [N,N]
  NumPointsBlob = [N,2*N]

  Length = 2._rk
  SeparationScale = 0.8_rk

  ! Initialize the problem
  call ovkCreateAssembler(Assembler, NumDims=2, NumGrids=4, Verbose=.true.)

  call ovkEditAssemblerProperties(Assembler, Properties)

  ! Automatically define boundaries in non-overlapping regions
  call ovkSetAssemblerPropertyInferBoundaries(Properties, OVK_ALL_GRIDS, .true.)

  ! Indicate which grids can intersect
  call ovkSetAssemblerPropertyOverlap(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
  call ovkSetAssemblerPropertyOverlapTolerance(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 0._rk)

  ! Indicate which grids can cut each other
  call ovkSetAssemblerPropertyBoundaryHoleCutting(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
  call ovkSetAssemblerPropertyOverlapHoleCutting(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)

  ! Indicate which grids can communicate and how
  call ovkSetAssemblerPropertyConnectionType(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_CONNECTION_FRINGE)
  call ovkSetAssemblerPropertyDisjointConnection(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, .true.)
  call ovkSetAssemblerPropertyInterpScheme(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, OVK_INTERP_LINEAR)
  call ovkSetAssemblerPropertyFringeSize(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 2)
  call ovkSetAssemblerPropertyFringePadding(Properties, OVK_ALL_GRIDS, OVK_ALL_GRIDS, 8)

  call ovkReleaseAssemblerProperties(Assembler, Properties)

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
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBackground(2)
    do i = 1, NumPointsBackground(1)
      U = real(i-1,kind=rk)/real(NumPointsBackground(1)-1,kind=rk)
      V = real(j-1,kind=rk)/real(NumPointsBackground(2)-1,kind=rk)
      X%values(i,j,1) = Length(1) * (U-0.5_rk)
      Y%values(i,j,1) = Length(2) * (V-0.5_rk)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  call ovkReleaseDomainGrid(Domain, Grid)

  !============
  ! Blob grids
  !============

  ! Initialize grid data structure for blob #1
  ! Periodic in the angular direction, with the last set of points being equal to the first
  call ovkCreateDomainGrid(Domain, 2, NumPoints=NumPointsBlob, Periodic=[.false.,.true.], &
    PeriodicStorage=OVK_OVERLAP_PERIODIC)
  call ovkEditDomainGrid(Domain, 2, Grid)

  ! Generate coordinates for blob #1
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBlob(2)-1
    do i = 1, NumPointsBlob(1)
      U = real(i-1, kind=rk)/real(NumPointsBlob(1)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPointsBlob(2)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(3._rk*Theta) + 0.1_rk*sin(2._rk*Theta+Pi/4._rk))
      RMax = 0.3_rk + 2._rk * RMin
      Radius = RMin * (RMax/RMin)**U
      X%values(i,j,1) = -0.5_rk*SeparationScale + Radius * cos(Theta)
      Y%values(i,j,1) = -0.2_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  ! Inner radial boundary on blob #1
  call ovkEditGridBoundaryMask(Grid, BoundaryMask)
  BoundaryMask%values(1,1:NumPointsBlob(2)-1,1) = .true.
  call ovkReleaseGridBoundaryMask(Grid, BoundaryMask)

  call ovkReleaseDomainGrid(Domain, Grid)

  ! Initialize grid data structures for blob #2
  ! Periodic in the angular direction, with the last set of points being equal to the first
  call ovkCreateDomainGrid(Domain, 3, NumPoints=NumPointsBlob, Periodic=[.false.,.true.], &
    PeriodicStorage=OVK_OVERLAP_PERIODIC)
  call ovkEditDomainGrid(Domain, 3, Grid)

  ! Generate coordinates for blob #2
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBlob(2)-1
    do i = 1, NumPointsBlob(1)
      U = real(i-1, kind=rk)/real(NumPointsBlob(1)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPointsBlob(2)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(4._rk*Theta+Pi/4._rk) + 0.1_rk*sin(2._rk*Theta))
      RMax = 0.4_rk + 2_rk * RMin
      Radius = RMin * (RMax/RMin)**U
      X%values(i,j,1) = -0.1_rk*SeparationScale + Radius * cos(Theta)
      Y%values(i,j,1) = 0.35_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  ! Inner radial boundary on blob #2
  call ovkEditGridBoundaryMask(Grid, BoundaryMask)
  BoundaryMask%values(1,1:NumPointsBlob(2)-1,1) = .true.
  call ovkReleaseGridBoundaryMask(Grid, BoundaryMask)

  call ovkReleaseDomainGrid(Domain, Grid)

  ! Initialize grid data structures for blob #3
  ! Periodic in the angular direction, with the last set of points being equal to the first
  call ovkCreateDomainGrid(Domain, 4, NumPoints=NumPointsBlob, Periodic=[.false.,.true.], &
    PeriodicStorage=OVK_OVERLAP_PERIODIC)
  call ovkEditDomainGrid(Domain, 4, Grid)

  ! Generate coordinates for blob #3
  call ovkEditGridCoords(Grid, 1, X)
  call ovkEditGridCoords(Grid, 2, Y)
  do j = 1, NumPointsBlob(2)-1
    do i = 1, NumPointsBlob(1)
      U = real(i-1, kind=rk)/real(NumPointsBlob(1)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPointsBlob(2)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(5._rk*Theta+Pi/4._rk) + 0.1_rk*sin(3._rk*Theta))
      RMax = 0.4_rk + 2_rk * RMin
      Radius = RMin * (RMax/RMin)**U
      X%values(i,j,1) = 0.3_rk*SeparationScale + Radius * cos(Theta)
      Y%values(i,j,1) = -0.3_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do
  call ovkReleaseGridCoords(Grid, X)
  call ovkReleaseGridCoords(Grid, Y)

  ! Inner radial boundary on blob #3
  call ovkEditGridBoundaryMask(Grid, BoundaryMask)
  BoundaryMask%values(1,1:NumPointsBlob(2)-1,1) = .true.
  call ovkReleaseGridBoundaryMask(Grid, BoundaryMask)

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
  NumPointsAll(:,2) = NumPointsBlob
  NumPointsAll(:,3) = NumPointsBlob
  NumPointsAll(:,4) = NumPointsBlob

  call ovkGetAssemblerDomain(Assembler, Domain)
  call ovkGetAssemblerConnectivity(Assembler, Connectivity)

  ! Write a PLOT3D grid file with IBlank to visualize the result
  call ovkCreateP3D(GridFile, "grid.xyz", NumDims=2, NumGrids=4, NumPointsAll=NumPointsAll, &
    WithIBlank=.true.)

  do m = 1, 4

    call ovkGetDomainGrid(Domain, m, Grid)
    call ovkGetConnectivityInterpData(Connectivity, m, InterpData)

    ! At the moment Overkit converts everything to no-overlap periodic internally, so
    ! we will need to convert back
    call ovkGetGridCart(Grid, Cart)
    CartOverlap = ovkCartConvertPeriodicStorage(Cart, OVK_OVERLAP_PERIODIC)

    call ovkExportGridCoords(Grid, 1, CartOverlap, XOverlap)
    call ovkExportGridCoords(Grid, 2, CartOverlap, YOverlap)

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
    do j = Cart%is(2), Cart%ie(2)
      do i = Cart%is(1), Cart%ie(1)
        if (GridMask%values(i,j,1)) then
          if (ReceiverMask%values(i,j,1)) then
            IBlank%values(i,j,1) = -DonorGridIDs%values(i,j,1)
          else if (OrphanMask%values(i,j,1)) then
            IBlank%values(i,j,1) = 7
          else
            IBlank%values(i,j,1) = 1
          end if
        else
          IBlank%values(i,j,1) = 0
        end if
      end do
    end do
    call ovkExportField(IBlank, CartOverlap, IBlankOverlap)

    call ovkWriteP3D(GridFile, m, XOverlap, YOverlap, IBlankOverlap)

  end do

  ! Finalize the PLOT3D file
  call ovkCloseP3D(GridFile)

  call ovkDestroyAssembler(Assembler)

end program Blobs
