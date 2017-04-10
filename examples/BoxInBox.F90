! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

program BoxInBox

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
  integer :: nDims
  character(len=32) :: InterpSchemeString
  integer, dimension(MAX_ND,2) :: nPoints
  real(rk), dimension(MAX_ND,2) :: Length
  type(ovk_cart), dimension(2) :: Carts
  type(ovk_field_real), dimension(:,:), allocatable :: XYZ
  integer, dimension(MAX_ND) :: Point
  real(rk) :: U
  integer, dimension(2) :: GridType
  type(ovk_grid), dimension(2) :: Grids
  integer, dimension(2) :: FringeSize
  integer, dimension(2,2) :: FringePadding
  integer, dimension(2) :: InterpScheme
  type(ovk_interp), dimension(2) :: InterpData
  type(ovk_field_logical), dimension(2) :: HoleMasks
  type(ovk_field_int) :: IBlank
  type(ovk_p3d_grid_file) :: GridFile
  type(ovk_pegasus) :: PegasusData

  allocate(RawArguments(command_argument_count()))
  do i = 1, size(RawArguments)
    call get_command_argument(i, RawArguments(i))
  end do

  Usage = "BoxInBox [<options> ...]"
  Description = "Generates an overset mesh consisting of two uniform cartesian grids."
  Options(1) = t_cmd_opt_("dimension", "d", CMD_OPT_INTEGER, "Grid dimension (2 or 3) [ Default: 2 ]")
  Options(2) = t_cmd_opt_("size", "N", CMD_OPT_INTEGER, "Grid size in each direction [ Default: 41 ]")
  Options(3) = t_cmd_opt_("interp", "i", CMD_OPT_STRING, "Interpolation scheme (linear or cubic) [ Default: cubic ]")

  call ParseArguments(RawArguments, Usage=Usage, Description=Description, Options=Options)

  call GetOptionValue(Options(1), nDims, 2)
  if (nDims /= 2 .and. nDims /= 3) then
    write (ERROR_UNIT, '(3a)') "ERROR: Invalid grid dimension ", nDims, "."
    stop 1
  end if

  call GetOptionValue(Options(2), N, 41)

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

  nPoints = 1
  nPoints(:nDims,1) = N
  nPoints(:nDims,2) = N

  Length = 0._rk
  Length(:nDims,1) = 2._rk
  Length(:nDims,2) = 1.25_rk

  ! Specify the grids' structural properties (dimension, size, periodicity, etc.)
  Carts(1) = ovk_cart_(nDims, nPoints(:nDims,1))
  Carts(2) = ovk_cart_(nDims, nPoints(:nDims,2))

  allocate(XYZ(nDims,2))

  ! Create the coordinate arrays
  do m = 1, 2
    XYZ(:,m) = ovk_field_real_(Carts(m))
    do k = 1, nPoints(3,m)
      do j = 1, nPoints(2,m)
        do i = 1, nPoints(1,m)
          Point = [i,j,k]
          do l = 1, nDims
            U = real(Point(l)-1,kind=rk)/real(nPoints(l,m)-1,kind=rk)
            XYZ(l,m)%values(i,j,k) = Length(l,m) * (U-0.5_rk)
          end do
        end do
      end do
    end do
  end do

  ! Information about grid type helps optimize performance
  GridType = OVK_GRID_TYPE_CARTESIAN

  ! Assemble the grid data structure
  do m = 1, 2
    call ovkMakeGrid(Grids(m), Carts(m), XYZ(:,m), GridType=GridType(m))
  end do

  ! Set parameters for overset assembly
  FringeSize = 2
  FringePadding = 2

  ! Perform overset assembly; results are stored in InterpData
  call ovkAssembleOverset(Grids, InterpData, FringeSize=FringeSize, FringePadding=FringePadding, &
    InterpScheme=InterpScheme, HoleMasks=HoleMasks)

  ! Write PLOT3D grid file
  ! IBlank values are set as follows:
  !   1 => normal point
  !   0 => hole
  !  -N => receives from grid N
  call ovkP3DCreate(GridFile, "grid.xyz", nGrids=2, Carts=Carts, WithIBlank=.true.)
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

end program BoxInBox
