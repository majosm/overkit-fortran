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
  integer, dimension(MAX_ND,4) :: NumPoints
  type(ovk_cart), dimension(4) :: Carts
  type(ovk_field_real), dimension(2,4) :: XYZ
  real(rk) :: U, V
  real(rk), dimension(2) :: Length
  real(rk) :: RMin, RMax
  real(rk) :: Radius
  real(rk) :: Theta
  real(rk) :: SeparationScale
  integer, dimension(4) :: GridType
  type(ovk_field_logical), dimension(4) :: BoundaryMasks
  type(ovk_grid), dimension(4) :: Grids
  integer, dimension(4) :: FringeSize
  integer, dimension(4) :: InterpScheme
  type(ovk_interp), dimension(4) :: InterpData
  type(ovk_field_logical), dimension(4) :: HoleMasks
  type(ovk_field_int) :: IBlank
  type(ovk_p3d_grid_file) :: GridFile
  type(ovk_pegasus) :: PegasusData

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

  NumPoints = 1
  NumPoints(:2,1) = N
  NumPoints(1,2) = N
  NumPoints(2,2) = 2*N
  NumPoints(1,3) = N
  NumPoints(2,3) = 2*N
  NumPoints(1,4) = N
  NumPoints(2,4) = 2*N

  Length = 2._rk

  ! Specify the grids' structural properties (dimension, size, periodicity, etc.)
  Carts(1) = ovk_cart_(2, NumPoints(:,1))
  Carts(2) = ovk_cart_(2, NumPoints(:,2), [.false.,.true.], OVK_OVERLAP_PERIODIC)
  Carts(3) = ovk_cart_(2, NumPoints(:,3), [.false.,.true.], OVK_OVERLAP_PERIODIC)
  Carts(4) = ovk_cart_(2, NumPoints(:,4), [.false.,.true.], OVK_OVERLAP_PERIODIC)

  ! Create the coordinate arrays
  XYZ(:,1) = ovk_field_real_(Carts(1))
  XYZ(:,2) = ovk_field_real_(Carts(2))
  XYZ(:,3) = ovk_field_real_(Carts(3))
  XYZ(:,4) = ovk_field_real_(Carts(4))

  ! Grid 1 is the Cartesian background grid
  do j = 1, NumPoints(2,1)
    do i = 1, NumPoints(1,1)
      U = real(i-1,kind=rk)/real(NumPoints(1,1)-1,kind=rk)
      V = real(j-1,kind=rk)/real(NumPoints(2,1)-1,kind=rk)
      XYZ(1,1)%values(i,j,1) = Length(1) * (U-0.5_rk)
      XYZ(2,1)%values(i,j,1) = Length(2) * (V-0.5_rk)
    end do
  end do

  SeparationScale = 0.8_rk

  ! Grid 2 wraps around a blob
  do j = 1, NumPoints(2,2)-1
    do i = 1, NumPoints(1,2)
      U = real(i-1, kind=rk)/real(NumPoints(1,2)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPoints(2,2)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(3._rk*Theta) + 0.1_rk*sin(2._rk*Theta+Pi/4._rk))
      RMax = 0.3_rk + 2._rk * RMin
      Radius = RMin * (RMax/RMin)**U
      XYZ(1,2)%values(i,j,1) = -0.5_rk*SeparationScale + Radius * cos(Theta)
      XYZ(2,2)%values(i,j,1) = -0.2_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do

  ! Grid 3 wraps around another blob
  do j = 1, NumPoints(2,3)-1
    do i = 1, NumPoints(1,3)
      U = real(i-1, kind=rk)/real(NumPoints(1,3)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPoints(2,3)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(4._rk*Theta+Pi/4._rk) + 0.1_rk*sin(2._rk*Theta))
      RMax = 0.4_rk + 2_rk * RMin
      Radius = RMin * (RMax/RMin)**U
      XYZ(1,3)%values(i,j,1) = -0.1_rk*SeparationScale + Radius * cos(Theta)
      XYZ(2,3)%values(i,j,1) = 0.35_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do

  ! Grid 4 wraps around yet another blob
  do j = 1, NumPoints(2,4)-1
    do i = 1, NumPoints(1,4)
      U = real(i-1, kind=rk)/real(NumPoints(1,4)-1, kind=rk)
      V = real(j-1, kind=rk)/real(NumPoints(2,4)-1, kind=rk)
      Theta = 2._rk * Pi * V
      RMin = 0.1_rk * (1._rk + 0.2_rk*sin(5._rk*Theta+Pi/4._rk) + 0.1_rk*sin(3._rk*Theta))
      RMax = 0.4_rk + 2_rk * RMin
      Radius = RMin * (RMax/RMin)**U
      XYZ(1,4)%values(i,j,1) = 0.3_rk*SeparationScale + Radius * cos(Theta)
      XYZ(2,4)%values(i,j,1) = -0.3_rk*SeparationScale + Radius * sin(Theta)
    end do
  end do

  ! Information about grid type helps optimize performance
  GridType(1) = OVK_GRID_TYPE_CARTESIAN
  GridType(2) = OVK_GRID_TYPE_CURVILINEAR
  GridType(3) = OVK_GRID_TYPE_CURVILINEAR
  GridType(4) = OVK_GRID_TYPE_CURVILINEAR

  ! Outer edge boundaries on background grid
  BoundaryMasks(1) = ovk_field_logical_(Carts(1), .false.)
  BoundaryMasks(1)%values(:,1,1) = .true.
  BoundaryMasks(1)%values(:,NumPoints(2,1),1) = .true.
  BoundaryMasks(1)%values(1,:,1) = .true.
  BoundaryMasks(1)%values(NumPoints(1,1),:,1) = .true.

  ! Inner radial boundary on blob #1
  BoundaryMasks(2) = ovk_field_logical_(Carts(2), .false.)
  BoundaryMasks(2)%values(1,:,1) = .true.

  ! Inner radial boundary on blob #2
  BoundaryMasks(3) = ovk_field_logical_(Carts(3), .false.)
  BoundaryMasks(3)%values(1,:,1) = .true.

  ! Inner radial boundary on blob #3
  BoundaryMasks(4) = ovk_field_logical_(Carts(4), .false.)
  BoundaryMasks(4)%values(1,:,1) = .true.

  ! Assemble the grid data structure
  do m = 1, 4
    call ovkMakeGrid(Grids(m), Carts(m), XYZ(:,m), GridType=GridType(m), &
      BoundaryMask=BoundaryMasks(m))
  end do

  ! Set parameters for overset assembly
  FringeSize = 2
  InterpScheme = OVK_INTERP_LINEAR

  ! Perform overset assembly; results are stored in InterpData
  call ovkAssembleOverset(Grids, InterpData, FringeSize=FringeSize, InterpScheme=InterpScheme, &
    HoleMasks=HoleMasks)

  ! Write PLOT3D grid file
  ! IBlank values are set as follows:
  !   1 => normal point
  !   0 => hole
  !  -N => receives from grid N
  call ovkP3DCreate(GridFile, "grid.xyz", NumGrids=4, Carts=Carts, WithIBlank=.true.)
  do m = 1, 4
    IBlank = ovk_field_int_(Carts(m), 1)
    call ovkDonorGridIDToIBlank(InterpData(m), IBlank, Multiplier=-1)
    call ovkMaskToIBlank(HoleMasks(m), IBlank, TrueValue=0)
    call ovkP3DWrite(GridFile, m, Grids(m)%xyz, IBlank)
  end do
  call ovkP3DClose(GridFile)

  ! Write the interpolation data
  call ovkMakePegasusData(Grids, InterpData, Carts, PegasusData)
  call ovkWritePegasusData(PegasusData, "XINTOUT.HO.2D", "XINTOUT.X.2D")

end program Blobs
