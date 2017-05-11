! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

program Bump

  use Overkit
  use ovsCommandLine
  use ovsGlobal
  implicit none

  integer :: i, j, k, l, m
  character(len=256), dimension(:), allocatable :: RawArguments
  character(len=256) :: Usage
  character(len=256) :: Description
  type(t_cmd_opt), dimension(3) :: Options
  integer :: N
  integer :: NumDims
  character(len=32) :: InterpSchemeString
  integer, dimension(MAX_ND,2) :: NumPoints
  real(rk), dimension(MAX_ND) :: Length
  type(ovk_cart), dimension(2) :: Carts
  type(ovk_field_real), dimension(:,:), allocatable :: XYZ
  integer, dimension(MAX_ND) :: Point
  real(rk) :: U
  real(rk), dimension(MAX_ND) :: Coords
  real(rk) :: R
  real(rk) :: BumpHeight
  real(rk) :: MinHeight, MaxHeight
  real(rk) :: Offset
  integer, dimension(2) :: GridType
  type(ovk_field_logical), dimension(2) :: BoundaryMasks
  type(ovk_grid), dimension(2) :: Grids
  integer, dimension(2) :: FringeSize
  integer, dimension(2,2) :: FringePadding
  integer, dimension(2) :: InterpScheme
  real(rk), dimension(2,2) :: OverlapTolerance
  type(ovk_interp), dimension(2) :: InterpData
  type(ovk_field_logical), dimension(2) :: HoleMasks
  type(ovk_field_int) :: IBlank
  type(ovk_p3d_grid_file) :: GridFile
  type(ovk_pegasus) :: PegasusData

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

  NumPoints = 1
  NumPoints(:NumDims-1,1) = N
  NumPoints(NumDims,1) = N/2
  NumPoints(:NumDims-1,2) = N
  NumPoints(NumDims,2) = (3*N)/4

  Length(:NumDims-1) = 2._rk
  Length(NumDims) = 1._rk

  ! Specify the grids' structural properties (dimension, size, periodicity, etc.)
  Carts(1) = ovk_cart_(NumDims, NumPoints(:NumDims,1))
  Carts(2) = ovk_cart_(NumDims, NumPoints(:NumDims,2))

  allocate(XYZ(NumDims,2))

  ! Create the coordinate arrays
  do m = 1, 2
    XYZ(:,m) = ovk_field_real_(Carts(m))
  end do

  ! Grid 1 is the Cartesian background grid
  do k = 1, NumPoints(3,1)
    do j = 1, NumPoints(2,1)
      do i = 1, NumPoints(1,1)
        Point = [i,j,k]
        do l = 1, NumDims-1
          U = real(Point(l)-1,kind=rk)/real(NumPoints(l,1)-1,kind=rk)
          XYZ(l,1)%values(i,j,k) = Length(l) * (U-0.5_rk)
        end do
        U = real(Point(NumDims)-1,kind=rk)/real(NumPoints(NumDims,1)-1,kind=rk)
        XYZ(NumDims,1)%values(i,j,k) = Length(NumDims) * U
      end do
    end do
  end do

  ! Grid 2 represents the bump
  do k = 1, NumPoints(3,2)
    do j = 1, NumPoints(2,2)
      do i = 1, NumPoints(1,2)
        Point = [i,j,k]
        do l = 1, NumDims-1
          U = real(Point(l)-1,kind=rk)/real(NumPoints(l,2)-1,kind=rk)
          Coords(l) = 0.65_rk*Length(l) * (U-0.5_rk)
        end do
        R = sqrt(sum(Coords(:NumDims-1)**2))
        if (R <= 0.499_rk) then
          BumpHeight = 0.25_rk*exp(1._rk-1._rk/(1._rk-(2._rk*R)**2))
        else
          BumpHeight = 0._rk
        end if
        U = real(Point(NumDims)-1,kind=rk)/real(NumPoints(NumDims,2)-1,kind=rk)
        MinHeight = BumpHeight
        MaxHeight = 0.5_rk * Length(NumDims)
        Offset = BumpHeight-0.1_rk
        Coords(NumDims) = (MinHeight-Offset)*((MaxHeight-Offset)/(MinHeight-Offset))**U + Offset
        do l = 1, NumDims
          XYZ(l,2)%values(i,j,k) = Coords(l)
        end do
      end do
    end do
  end do

  ! Information about grid type helps optimize performance
  GridType(1) = OVK_GRID_TYPE_CARTESIAN
  GridType(2) = OVK_GRID_TYPE_CURVILINEAR

  ! Lower edge boundary on background grid (other boundaries don't need to be specified as there
  ! is no overlap there)
  BoundaryMasks(1) = ovk_field_logical_(Carts(1), .false.)
  do k = 1, NumPoints(3,1)
    do j = 1, NumPoints(2,1)
      do i = 1, NumPoints(1,1)
        Point = [i,j,k]
        do l = 1, NumDims
          Coords(l) = XYZ(l,1)%values(i,j,k)
        end do
        R = sqrt(sum(Coords(:NumDims-1)**2))
        if (R >= 0.5_rk .and. Point(NumDims) == 1) then
          BoundaryMasks(1)%values(i,j,k) = .true.
        end if
      end do
    end do
  end do

  ! Lower edge boundary on bump grid
  BoundaryMasks(2) = ovk_field_logical_(Carts(2), .false.)
  select case (NumDims)
  case (2)
    BoundaryMasks(2)%values(:,1,1) = .true.
  case (3)
    BoundaryMasks(2)%values(:,:,1) = .true.
  end select

  ! Assemble the grid data structure
  do m = 1, 2
    call ovkMakeGrid(Grids(m), Carts(m), XYZ(:,m), GridType=GridType(m), &
      BoundaryMask=BoundaryMasks(m))
  end do

  ! Set parameters for overset assembly
  FringeSize = 2
  FringePadding = 2
  OverlapTolerance = 0.1_rk

  ! Perform overset assembly; results are stored in InterpData
  call ovkAssembleOverset(Grids, InterpData, FringeSize=FringeSize, FringePadding=FringePadding, &
    InterpScheme=InterpScheme, OverlapTolerance=OverlapTolerance, HoleMasks=HoleMasks)

  ! Write PLOT3D grid file
  ! IBlank values are set as follows:
  !   1 => normal point
  !   0 => hole
  !  -N => receives from grid N
  call ovkP3DCreate(GridFile, "grid.xyz", NumGrids=2, Carts=Carts, WithIBlank=.true.)
  do m = 1, 2
    IBlank = ovk_field_int_(Carts(m), 1)
    call ovkDonorGridIDToIBlank(InterpData(m), IBlank, Multiplier=-1)
    call ovkMaskToIBlank(HoleMasks(m), IBlank, TrueValue=0)
    call ovkP3DWrite(GridFile, m, Grids(m)%xyz, IBlank)
  end do
  call ovkP3DClose(GridFile)

  ! Write the interpolation data
  call ovkMakePegasusData(Grids, InterpData, Carts, PegasusData)
  call ovkWritePegasusData(PegasusData, "XINTOUT.HO.2D", "XINTOUT.X.2D")

end program Bump
