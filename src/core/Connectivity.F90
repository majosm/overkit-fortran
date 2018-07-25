! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkConnectivity

  use ovkBoundingBox
  use ovkCart
  use ovkOverlap
  use ovkField
  use ovkFieldOps
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkLogger
  implicit none

  private

  ! API
  public :: ovk_connectivity
  public :: ovkConnectivityExists
  public :: ovkGetConnectivityDonorGridID
  public :: ovkGetConnectivityReceiverGridID
  public :: ovkGetConnectivityDimension
  public :: ovkGetConnectivityMaxDonorSize
  public :: ovkGetConnectivityCount
  public :: ovkGetConnectivityDonorExtents
  public :: ovkGetConnectivityDonorCoords
  public :: ovkGetConnectivityDonorInterpCoefs
  public :: ovkGetConnectivityReceiverPoints

  ! Internal
  public :: ovk_connectivity_
  public :: t_donor_grid_info
  public :: t_donor_grid_info_
  public :: CreateConnectivity
  public :: DestroyConnectivity
  public :: ResizeConnectivity
  public :: CreateDonorGridInfo
  public :: FillConnectivity

  type ovk_connectivity
    type(t_noconstruct) :: noconstruct
    type(t_existence_flag) :: existence_flag
    type(t_logger), pointer :: logger
    integer :: donor_grid_id
    integer :: receiver_grid_id
    integer :: nd
    integer :: max_donor_size
    integer(lk) :: nconnections
    integer, dimension(:,:,:), pointer :: donor_extents
    real(rk), dimension(:,:), pointer :: donor_coords
    real(rk), dimension(:,:,:), pointer :: donor_interp_coefs
    integer, dimension(:,:), pointer :: receiver_points
  end type ovk_connectivity

  type t_donor_grid_info
    integer :: max_stencil_ncells
    integer :: max_stencil_cell_offset
    type(ovk_field_int) :: cell_edge_dists
    type(ovk_field_int), dimension(:), allocatable :: cell_edge_dist_gradients
  end type t_donor_grid_info

contains

  pure function ovk_connectivity_() result(Connectivity)

    type(ovk_connectivity) :: Connectivity

    nullify(Connectivity%logger)
    Connectivity%donor_grid_id = 0
    Connectivity%receiver_grid_id = 0
    Connectivity%nd = 2
    Connectivity%max_donor_size = 1
    Connectivity%nconnections = 0
    nullify(Connectivity%donor_extents)
    nullify(Connectivity%donor_coords)
    nullify(Connectivity%donor_interp_coefs)
    nullify(Connectivity%receiver_points)

    call SetExists(Connectivity%existence_flag, .false.)

  end function ovk_connectivity_

  subroutine CreateConnectivity(Connectivity, DonorGridID, ReceiverGridID, Logger, NumDims)

    type(ovk_connectivity), intent(out) :: Connectivity
    integer, intent(in) :: DonorGridID, ReceiverGridID
    type(t_logger), pointer, intent(in) :: Logger
    integer, intent(in) :: NumDims

    Connectivity%logger => Logger

    Connectivity%donor_grid_id = DonorGridID
    Connectivity%receiver_grid_id = ReceiverGridID
    Connectivity%nd = NumDims
    Connectivity%max_donor_size = 1
    Connectivity%nconnections = 0

    allocate(Connectivity%donor_extents(MAX_ND,2,0))
    allocate(Connectivity%donor_coords(NumDims,0))
    allocate(Connectivity%donor_interp_coefs(1,NumDims,0))
    allocate(Connectivity%receiver_points(MAX_ND,0))

    call SetExists(Connectivity%existence_flag, .true.)

  end subroutine CreateConnectivity

  subroutine DestroyConnectivity(Connectivity)

    type(ovk_connectivity), intent(inout) :: Connectivity

    if (.not. ovkConnectivityExists(Connectivity)) return

    call SetExists(Connectivity%existence_flag, .false.)

    deallocate(Connectivity%donor_extents)
    deallocate(Connectivity%donor_coords)
    deallocate(Connectivity%donor_interp_coefs)
    deallocate(Connectivity%receiver_points)

  end subroutine DestroyConnectivity

  function ovkConnectivityExists(Connectivity) result(Exists)

    type(ovk_connectivity), intent(in) :: Connectivity
    logical :: Exists

    Exists = CheckExists(Connectivity%existence_flag)

  end function ovkConnectivityExists

  subroutine ovkGetConnectivityDonorGridID(Connectivity, DonorGridID)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer, intent(out) :: DonorGridID

    DonorGridID = Connectivity%donor_grid_id

  end subroutine ovkGetConnectivityDonorGridID

  subroutine ovkGetConnectivityReceiverGridID(Connectivity, ReceiverGridID)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer, intent(out) :: ReceiverGridID

    ReceiverGridID = Connectivity%receiver_grid_id

  end subroutine ovkGetConnectivityReceiverGridID

  subroutine ovkGetConnectivityDimension(Connectivity, NumDims)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer, intent(out) :: NumDims

    NumDims = Connectivity%nd

  end subroutine ovkGetConnectivityDimension

  subroutine ovkGetConnectivityMaxDonorSize(Connectivity, MaxDonorSize)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer, intent(out) :: MaxDonorSize

    MaxDonorSize = Connectivity%max_donor_size

  end subroutine ovkGetConnectivityMaxDonorSize

  subroutine ovkGetConnectivityCount(Connectivity, NumConnections)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer(lk), intent(out) :: NumConnections

    NumConnections = Connectivity%nconnections

  end subroutine ovkGetConnectivityCount

  subroutine ovkGetConnectivityDonorExtents(Connectivity, DonorExtents)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer, dimension(:,:,:), pointer, intent(out) :: DonorExtents

    DonorExtents => Connectivity%donor_extents

  end subroutine ovkGetConnectivityDonorExtents

  subroutine ovkGetConnectivityDonorCoords(Connectivity, DonorCoords)

    type(ovk_connectivity), intent(in) :: Connectivity
    real(rk), dimension(:,:), pointer, intent(out) :: DonorCoords

    DonorCoords => Connectivity%donor_coords

  end subroutine ovkGetConnectivityDonorCoords

  subroutine ovkGetConnectivityDonorInterpCoefs(Connectivity, DonorInterpCoefs)

    type(ovk_connectivity), intent(in) :: Connectivity
    real(rk), dimension(:,:,:), pointer, intent(out) :: DonorInterpCoefs

    DonorInterpCoefs => Connectivity%donor_interp_coefs

  end subroutine ovkGetConnectivityDonorInterpCoefs

  subroutine ovkGetConnectivityReceiverPoints(Connectivity, ReceiverPoints)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer, dimension(:,:), pointer, intent(out) :: ReceiverPoints

    ReceiverPoints => Connectivity%receiver_points

  end subroutine ovkGetConnectivityReceiverPoints

  subroutine ResizeConnectivity(Connectivity, NumConnections, MaxDonorSize)

    type(ovk_connectivity), intent(inout) :: Connectivity
    integer(lk), intent(in) :: NumConnections
    integer, intent(in) :: MaxDonorSize

    deallocate(Connectivity%donor_extents)
    deallocate(Connectivity%donor_coords)
    deallocate(Connectivity%donor_interp_coefs)
    deallocate(Connectivity%receiver_points)

    allocate(Connectivity%donor_extents(MAX_ND,2,NumConnections))
    allocate(Connectivity%donor_coords(Connectivity%nd,NumConnections))
    allocate(Connectivity%donor_interp_coefs(MaxDonorSize,Connectivity%nd,NumConnections))
    allocate(Connectivity%receiver_points(MAX_ND,NumConnections))

    Connectivity%max_donor_size = MaxDonorSize
    Connectivity%nconnections = NumConnections

  end subroutine ResizeConnectivity

  pure function t_donor_grid_info_() result(DonorGridInfo)

    type(t_donor_grid_info) :: DonorGridInfo

    DonorGridInfo%max_stencil_ncells = 0
    DonorGridInfo%max_stencil_cell_offset = 0

  end function t_donor_grid_info_

  subroutine CreateDonorGridInfo(DonorGrid, NumGrids, ConnectionTypes, DonorGridInfo)

    type(ovk_grid), intent(in) :: DonorGrid
    integer, intent(in) :: NumGrids
    integer, dimension(:), intent(in) :: ConnectionTypes
    type(t_donor_grid_info), intent(out) :: DonorGridInfo

    integer :: d, i, j, k, n
    integer :: NumDims
    integer, dimension(MAX_ND) :: StencilStart, StencilEnd
    integer, dimension(MAX_ND) :: StencilCellStart, StencilCellEnd
    integer :: NumCells
    type(ovk_cart) :: CellCart
    type(ovk_cart) :: ExtendedCellCart
    type(ovk_field_logical) :: NotMask
    integer :: CellEdgeDistance
    integer, dimension(MAX_ND) :: PrevOffset, NextOffset
    integer, dimension(MAX_ND) :: PrevCell, NextCell
    integer :: PrevDistance, NextDistance

    NumDims = DonorGrid%nd

    DonorGridInfo%max_stencil_ncells = 0
    DonorGridInfo%max_stencil_cell_offset = 0
    do n = 1, NumGrids
      if (ConnectionTypes(n) /= OVK_CONNECTION_NONE) then
        StencilStart = 0
        StencilEnd = 0
        call GetDonorStencilExtents(NumDims, ConnectionTypes(n), StencilStart, StencilEnd)
        StencilCellStart = StencilStart
        StencilCellEnd = StencilEnd
        StencilCellEnd(:NumDims) = StencilCellEnd(:NumDims)-1
        NumCells = 1
        do d = 1, NumDims
          NumCells = NumCells * (StencilCellEnd(d)-StencilCellStart(d)+1)
        end do
        DonorGridInfo%max_stencil_ncells = max(DonorGridInfo%max_stencil_ncells, NumCells)
        DonorGridInfo%max_stencil_cell_offset = max(DonorGridInfo%max_stencil_cell_offset, &
          max(-minval(StencilCellStart), maxval(StencilCellEnd)))
      end if
    end do

    if (DonorGridInfo%max_stencil_cell_offset > 0) then

      CellCart = DonorGrid%cell_cart

      ExtendedCellCart = CellCart
      ExtendedCellCart%is(:NumDims) = ExtendedCellCart%is(:NumDims) - merge(0, 1, &
        DonorGrid%cart%periodic(:NumDims))
      ExtendedCellCart%ie(:NumDims) = ExtendedCellCart%ie(:NumDims) + merge(0, 1, &
        DonorGrid%cart%periodic(:NumDims))

      NotMask = ovk_field_logical_(ExtendedCellCart, .true.)
      NotMask%values(CellCart%is(1):CellCart%ie(1),CellCart%is(2):CellCart%ie(2), &
        CellCart%is(3):CellCart%ie(3)) = .not. DonorGrid%cell_mask%values

      call ovkDistanceField(NotMask, OVK_TRUE, DonorGridInfo%cell_edge_dists)

      allocate(DonorGridInfo%cell_edge_dist_gradients(NumDims))

      do d = 1, NumDims
        DonorGridInfo%cell_edge_dist_gradients(d) = ovk_field_int_(CellCart, 0)
        PrevOffset = 0
        PrevOffset(d) = -1
        NextOffset = 0
        NextOffset(d) = 1
        do k = CellCart%is(3), CellCart%ie(3)
          do j = CellCart%is(2), CellCart%ie(2)
            do i = CellCart%is(1), CellCart%ie(1)
              CellEdgeDistance = DonorGridInfo%cell_edge_dists%values(i,j,k)
              if (CellEdgeDistance > 0 .and. CellEdgeDistance-DonorGridInfo%max_stencil_cell_offset &
                < 1) then
                PrevCell = [i+PrevOffset(1),j+PrevOffset(2),k+PrevOffset(3)]
                PrevCell(:NumDims) = ovkCartPeriodicAdjust(DonorGridInfo%cell_edge_dists%cart, PrevCell)
                NextCell = [i+NextOffset(1),j+NextOffset(2),k+NextOffset(3)]
                NextCell(:NumDims) = ovkCartPeriodicAdjust(DonorGridInfo%cell_edge_dists%cart, NextCell)
                PrevDistance = DonorGridInfo%cell_edge_dists%values(PrevCell(1),PrevCell(2),PrevCell(3))
                NextDistance = DonorGridInfo%cell_edge_dists%values(NextCell(1),NextCell(2),NextCell(3))
                DonorGridInfo%cell_edge_dist_gradients(d)%values(i,j,k) = NextDistance-PrevDistance
              end if
            end do
          end do
        end do
      end do

    end if

  end subroutine CreateDonorGridInfo

  subroutine FillConnectivity(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ConnectionType, &
    ReceiverMask, Connectivity)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_grid_info), intent(in) :: DonorGridInfo
    integer, intent(in) :: ConnectionType
    type(ovk_field_logical), intent(in) :: ReceiverMask
    type(ovk_connectivity), intent(inout) :: Connectivity

    integer :: i, j, k
    integer(lk) :: l, p
    integer :: NumDims
    integer(lk) :: NumConnections
    integer, dimension(MAX_ND) :: StencilLower, StencilUpper
    integer, dimension(MAX_ND) :: StencilSize
    integer :: MaxDonorSize
    integer, dimension(MAX_ND,2) :: DonorExtents
    real(rk), dimension(DonorGrid%nd) :: DonorCoords
    real(rk), dimension(:,:), allocatable :: DonorInterpCoefs
    integer, dimension(MAX_ND) :: ReceiverPoint

    NumDims = DonorGrid%nd

    NumConnections = ovkCountMask(ReceiverMask)

    if (NumConnections > 0_lk) then

      StencilLower = 0
      StencilUpper = 0
      call GetDonorStencilExtents(NumDims, ConnectionType, StencilLower, StencilUpper)
      StencilSize = StencilUpper - StencilLower + 1
      MaxDonorSize = maxval(StencilSize)

      call ResizeConnectivity(Connectivity, NumConnections, MaxDonorSize)

      allocate(DonorInterpCoefs(MaxDonorSize,NumDims))

      l = 1_lk
      p = 1_lk
      DonorExtents = 1
      do k = ReceiverGrid%cart%is(3), ReceiverGrid%cart%ie(3)
        do j = ReceiverGrid%cart%is(2), ReceiverGrid%cart%ie(2)
          do i = ReceiverGrid%cart%is(1), ReceiverGrid%cart%ie(1)
            if (ReceiverMask%values(i,j,k)) then
              ReceiverPoint = [i,j,k]
              call FindDonor(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ReceiverPoint, l, &
                ConnectionType, DonorExtents, DonorCoords, DonorInterpCoefs)
              Connectivity%donor_extents(:,:,p) = DonorExtents
              Connectivity%donor_coords(:,p) = DonorCoords
              Connectivity%donor_interp_coefs(:,:,p) = DonorInterpCoefs
              Connectivity%receiver_points(:,p) = ReceiverPoint
              p = p + 1_lk
            end if
            if (Overlap%mask%values(i,j,k)) then
              l = l + 1_lk
            end if
          end do
        end do
      end do

      deallocate(DonorInterpCoefs)

      if (Connectivity%logger%verbose) then
        write (*, '(7a)') "* ", trim(LargeIntToString(NumConnections)), &
          " donor/receiver pairs between grid ", trim(IntToString(DonorGrid%id)), &
          " and grid ", trim(IntToString(ReceiverGrid%id)), "."
      end if

    end if

  end subroutine FillConnectivity

  subroutine GetDonorStencilExtents(NumDims, ConnectionType, StencilLower, StencilUpper)

    integer, intent(in) :: NumDims
    integer, intent(in) :: ConnectionType
    integer, dimension(NumDims), intent(out) :: StencilLower, StencilUpper

    select case (ConnectionType)
    case (OVK_CONNECTION_NEAREST)
      StencilLower = 0
      StencilUpper = 0
    case (OVK_CONNECTION_LINEAR)
      StencilLower = 0
      StencilUpper = 1
    case (OVK_CONNECTION_CUBIC)
      StencilLower = -1
      StencilUpper = 2
    case default
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid connection type."
        stop 1
      end if
    end select

  end subroutine GetDonorStencilExtents

  subroutine FindDonor(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ReceiverPoint, &
    OverlapIndex, ConnectionType, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_grid_info), intent(in) :: DonorGridInfo
    integer, dimension(:), intent(in) :: ReceiverPoint
    integer(lk), intent(in) :: OverlapIndex
    integer, intent(in) :: ConnectionType
    integer, dimension(:,:), intent(out) :: DonorExtents
    real(rk), dimension(:), intent(out) :: DonorCoords
    real(rk), dimension(:,:), intent(out) :: DonorInterpCoefs

    select case (ConnectionType)
    case (OVK_CONNECTION_NEAREST)
      call FindDonorNearest(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ReceiverPoint, &
        OverlapIndex, DonorExtents, DonorCoords, DonorInterpCoefs)
    case (OVK_CONNECTION_LINEAR)
      call FindDonorLinear(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ReceiverPoint, &
        OverlapIndex, DonorExtents, DonorCoords, DonorInterpCoefs)
    case (OVK_CONNECTION_CUBIC)
      call FindDonorCubic(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ReceiverPoint, &
        OverlapIndex, DonorExtents, DonorCoords, DonorInterpCoefs)
    case default
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid connection type."
        stop 1
      end if
    end select

  end subroutine FindDonor

  subroutine FindDonorNearest(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ReceiverPoint, &
    OverlapIndex, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_grid_info), intent(in) :: DonorGridInfo
    integer, dimension(:), intent(in) :: ReceiverPoint
    integer(lk), intent(in) :: OverlapIndex
    integer, dimension(:,:), intent(out) :: DonorExtents
    real(rk), dimension(:), intent(out) :: DonorCoords
    real(rk), dimension(:,:), intent(out) :: DonorInterpCoefs

    integer :: d
    integer :: NumDims
    integer, dimension(MAX_ND) :: NearestPoint

    NumDims = DonorGrid%nd

    NearestPoint = 1
    do d = 1, NumDims
      NearestPoint(d) = Overlap%cells(d,OverlapIndex) + int(Overlap%coords(d,OverlapIndex) + 0.5_rk)
    end do
    if (.not. ovkCartContains(DonorGrid%cart, NearestPoint)) then
      NearestPoint = ovkCartPeriodicAdjust(DonorGrid%cart, NearestPoint)
    end if

    DonorExtents(:NumDims,1) = NearestPoint(:NumDims)
    DonorExtents(:NumDims,2) = NearestPoint(:NumDims)

    DonorCoords = 0._rk

    DonorInterpCoefs(:,:NumDims) = 0._rk
    do d = 1, NumDims
      DonorInterpCoefs(1,d) = 1._rk
    end do

  end subroutine FindDonorNearest

  subroutine FindDonorLinear(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ReceiverPoint, &
    OverlapIndex, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_grid_info), intent(in) :: DonorGridInfo
    integer, dimension(:), intent(in) :: ReceiverPoint
    integer(lk), intent(in) :: OverlapIndex
    integer, dimension(:,:), intent(out) :: DonorExtents
    real(rk), dimension(:), intent(out) :: DonorCoords
    real(rk), dimension(:,:), intent(out) :: DonorInterpCoefs

    integer :: d
    integer :: NumDims

    NumDims = DonorGrid%nd

    DonorExtents(:NumDims,1) = Overlap%cells(:NumDims,OverlapIndex)
    DonorExtents(:NumDims,2) = Overlap%cells(:NumDims,OverlapIndex) + 1

    DonorCoords = Overlap%coords(:,OverlapIndex)

    DonorInterpCoefs(:,:NumDims) = 0._rk
    do d = 1, NumDims
      DonorInterpCoefs(:2,d) = ovkInterpBasisLinear(DonorCoords(d))
    end do

  end subroutine FindDonorLinear

  subroutine FindDonorCubic(DonorGrid, ReceiverGrid, Overlap, DonorGridInfo, ReceiverPoint, &
    OverlapIndex, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_grid_info), intent(in) :: DonorGridInfo
    integer, dimension(:), intent(in) :: ReceiverPoint
    integer(lk), intent(in) :: OverlapIndex
    integer, dimension(:,:), intent(out) :: DonorExtents
    real(rk), dimension(:), intent(out) :: DonorCoords
    real(rk), dimension(:,:), intent(out) :: DonorInterpCoefs

    integer :: d
    integer :: NumDims
    real(rk), dimension(DonorGrid%nd) :: ReceiverCoords
    integer, dimension(MAX_ND) :: OverlapCell
    real(rk), dimension(DonorGrid%nd) :: OverlapCoords
    integer, dimension(MAX_ND) :: StencilLower, StencilUpper
    integer, dimension(MAX_ND) :: StencilSize
    integer, dimension(MAX_ND) :: StencilCellLower, StencilCellUpper
    logical, dimension(3**DonorGrid%nd) :: StencilCellMask
    integer, dimension(MAX_ND) :: StencilShift
    logical :: Success
    integer :: NumWarnings
    character(len=STRING_LENGTH) :: DonorCellString
    character(len=STRING_LENGTH) :: DonorGridIDString
    character(len=STRING_LENGTH) :: ReceiverPointString
    character(len=STRING_LENGTH) :: ReceiverGridIDString

    NumDims = DonorGrid%nd

    do d = 1, NumDims
      ReceiverCoords(d) = ReceiverGrid%coords(d)%values(ReceiverPoint(1),ReceiverPoint(2), &
        ReceiverPoint(3))
    end do

    OverlapCell = Overlap%cells(:,OverlapIndex)
    OverlapCoords = Overlap%coords(:,OverlapIndex)

    Success = .false.
    NumWarnings = 0

    StencilLower(:NumDims) = -1
    StencilLower(NumDims+1:) = 0
    StencilUpper(:NumDims) = 2
    StencilUpper(NumDims+1:) = 0
    StencilSize = StencilUpper - StencilLower + 1

    if (DonorGridInfo%cell_edge_dists%values(OverlapCell(1),OverlapCell(2),OverlapCell(3)) > 1) then
      StencilShift = 0
      Success = .true.
    else
      StencilCellLower(:NumDims) = -1
      StencilCellLower(NumDims+1:) = 0
      StencilCellUpper(:NumDims) = 1
      StencilCellUpper(NumDims+1:) = 0
      StencilCellMask = .true.
      call ShiftStencil(DonorGrid, DonorGridInfo, StencilCellLower, StencilCellUpper, &
        StencilCellMask, OverlapCell, StencilShift, Success)
    end if

    if (Success) then
      DonorExtents(:NumDims,1) = OverlapCell(:NumDims) + StencilShift(:NumDims) + &
        StencilLower(:NumDims)
      DonorExtents(:NumDims,1) = ovkCartPeriodicAdjust(DonorGrid%cart, DonorExtents(:,1))
      DonorExtents(:NumDims,2) = DonorExtents(:NumDims,1) + StencilSize(:NumDims)
      DonorCoords = OverlapCoords - real(StencilShift(:NumDims),kind=rk)
      DonorCoords = ovkCoordsInCubicGridCell(DonorGrid, DonorExtents(:,1), ReceiverCoords, &
        Guess=DonorCoords)
      DonorInterpCoefs(:,:NumDims) = 0._rk
      do d = 1, NumDims
        DonorInterpCoefs(:,d) = ovkInterpBasisCubic(DonorCoords(d))
      end do
    else
      if (DonorGrid%logger%verbose) then
        DonorCellString = TupleToString(DonorExtents(:NumDims,1))
        DonorGridIDString = IntToString(DonorGrid%id)
        ReceiverPointString = TupleToString(ReceiverPoint(:NumDims))
        ReceiverGridIDString = IntToString(ReceiverGrid%id)
        write (ERROR_UNIT, '(10a)') "WARNING: Could not use cubic ", &
          "interpolation for donor cell ", trim(DonorCellString), " of grid ", &
          trim(DonorGridIDString), " corresponding to receiver point ", &
          trim(ReceiverPointString), " of grid ", trim(ReceiverGridIDString), "; using linear instead."
        if (NumWarnings == 100) then
          write (ERROR_UNIT, '(a)') "WARNING: Further warnings suppressed."
        end if
        NumWarnings = NumWarnings + 1
      end if
      DonorExtents(:NumDims,1) = OverlapCell(:NumDims)
      DonorExtents(:NumDims,2) = OverlapCell(:NumDims) + 1
      DonorCoords = OverlapCoords
      DonorInterpCoefs(:,:NumDims) = 0._rk
      do d = 1, NumDims
        DonorInterpCoefs(:2,d) = ovkInterpBasisLinear(DonorCoords(d))
      end do
    end if

  end subroutine FindDonorCubic

  subroutine ShiftStencil(DonorGrid, DonorGridInfo, StencilCellLower_, StencilCellUpper_, &
    StencilCellMask, OverlapCell_, StencilShift, Success)

    type(ovk_grid), intent(in) :: DonorGrid
    type(t_donor_grid_info), intent(in) :: DonorGridInfo
    integer, dimension(DonorGrid%nd), intent(in) :: StencilCellLower_, StencilCellUpper_
    logical, dimension(:), intent(in) :: StencilCellMask
    integer, dimension(DonorGrid%nd), intent(in) :: OverlapCell_
    integer, dimension(DonorGrid%nd), intent(out) :: StencilShift
    logical, intent(out) :: Success

    integer :: d, i, j, k, l
    integer :: NumDims
    integer, dimension(MAX_ND) :: StencilCellLower, StencilCellUpper
    integer, dimension(MAX_ND) :: OverlapCell
    real(rk), dimension(DonorGrid%nd) :: Gradient
    real(rk), dimension(DonorGrid%nd) :: GradientNormalized
    real(rk) :: DirectionQuality, DistanceQuality
    logical, dimension(DonorGridInfo%max_stencil_ncells) :: ShiftAllowed
    real(rk), dimension(DonorGridInfo%max_stencil_ncells) :: ShiftQuality
    integer, dimension(MAX_ND) :: Shift
    real(rk), dimension(DonorGrid%nd) :: Displacement
    real(rk), dimension(DonorGrid%nd) :: DisplacementNormalized
    integer, dimension(MAX_ND) :: BestShift
    integer :: BestShiftIndex
    real(rk) :: BestShiftQuality
    integer, dimension(MAX_ND) :: ShiftedStencilCellStart
    integer, dimension(MAX_ND) :: ShiftedStencilCellEnd
    logical :: AwayFromBoundary
    integer, dimension(MAX_ND) :: Cell
    logical :: CellExists

    NumDims = DonorGrid%nd

    StencilCellLower = 0
    StencilCellLower(:NumDims) = StencilCellLower_
    StencilCellUpper = 0
    StencilCellUpper(:NumDims) = StencilCellUpper_
    OverlapCell = 1
    OverlapCell(:NumDims) = OverlapCell_

    do d = 1, NumDims
      Gradient(d) = real(DonorGridInfo%cell_edge_dist_gradients(d)%values(OverlapCell(1), &
        OverlapCell(2),OverlapCell(3)), kind=rk)
    end do
    GradientNormalized = Gradient/max(sqrt(sum(Gradient**2)),1._rk)

    l = 1
    do k = StencilCellLower(3), StencilCellUpper(3)
      do j = StencilCellLower(2), StencilCellUpper(2)
        do i = StencilCellLower(1), StencilCellUpper(1)
          ShiftAllowed(l) = StencilCellMask(l)
          if (ShiftAllowed(l)) then
            Shift = [-i,-j,-k]
            Displacement = real(Shift(:NumDims), kind=rk)
            DisplacementNormalized = Displacement/max(sqrt(sum(Displacement**2)),1._rk)
            DirectionQuality = 0.5_rk*(1._rk + dot_product(GradientNormalized,DisplacementNormalized))
            ! No particularly great justification for this coefficient other than to make
            !   qual(dot-1,dist) = qual(dot,dist+1) + some small amount (e.g. 1.e-12)
            ! so that 0 shift is preferred over distance 1 shift in gradient direction
            DistanceQuality = (0.5_rk+1.e-12_rk)*(DonorGridInfo%max_stencil_cell_offset-maxval(abs(Shift)))
            ShiftQuality(l) = DirectionQuality + DistanceQuality
          end if
          l = l + 1
        end do
      end do
    end do

    Success = .false.
    StencilShift = 0
    do
      BestShift = 0
      BestShiftIndex = 0
      BestShiftQuality = -huge(0._rk)
      l = 1
      do k = StencilCellLower(3), StencilCellUpper(3)
        do j = StencilCellLower(2), StencilCellUpper(2)
          do i = StencilCellLower(1), StencilCellUpper(1)
            if (ShiftAllowed(l) .and. ShiftQuality(l) > BestShiftQuality) then
              BestShift = [-i,-j,-k]
              BestShiftIndex = l
              BestShiftQuality = ShiftQuality(l)
            end if
            l = l + 1
          end do
        end do
      end do
      if (BestShiftIndex > 0) then
        ShiftedStencilCellStart = OverlapCell + StencilCellLower + BestShift
        ShiftedStencilCellEnd = OverlapCell + StencilCellUpper + BestShift
        AwayFromBoundary = ovkCartContains(DonorGrid%cell_cart, ShiftedStencilCellStart) .and. &
          ovkCartContains(DonorGrid%cell_cart, ShiftedStencilCellEnd)
        if (AwayFromBoundary) then
          l = 1
      L1: do k = ShiftedStencilCellStart(3), ShiftedStencilCellEnd(3)
            do j = ShiftedStencilCellStart(2), ShiftedStencilCellEnd(2)
              do i = ShiftedStencilCellStart(1), ShiftedStencilCellEnd(1)
                if (StencilCellMask(l)) then
                  Cell = [i,j,k]
                  CellExists = DonorGrid%cell_mask%values(Cell(1),Cell(2),Cell(3))
                  if (.not. CellExists) then
                    ShiftAllowed(BestShiftIndex) = .false.
                    exit L1
                  end if
                end if
                l = l + 1
              end do
            end do
          end do L1
        else
          l = 1
      L2: do k = ShiftedStencilCellStart(3), ShiftedStencilCellEnd(3)
            do j = ShiftedStencilCellStart(2), ShiftedStencilCellEnd(2)
              do i = ShiftedStencilCellStart(1), ShiftedStencilCellEnd(1)
                if (StencilCellMask(l)) then
                  Cell = [i,j,k]
                  Cell(:NumDims) = ovkCartPeriodicAdjust(DonorGrid%cell_cart, Cell)
                  if (ovkCartContains(DonorGrid%cell_cart, Cell)) then
                    CellExists = DonorGrid%cell_mask%values(Cell(1),Cell(2),Cell(3))
                  else
                    CellExists = .false.
                  end if
                  if (.not. CellExists) then
                    ShiftAllowed(BestShiftIndex) = .false.
                    exit L2
                  end if
                end if
                l = l + 1
              end do
            end do
          end do L2
        end if
        if (ShiftAllowed(BestShiftIndex)) then
          Success = .true.
          StencilShift = BestShift(:NumDims)
          exit
        end if
      else
        exit
      end if
    end do

  end subroutine ShiftStencil

end module ovkConnectivity
