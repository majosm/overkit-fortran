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
  public :: ovk_connectivity_properties
  public :: ovkConnectivityExists
  public :: ovkGetConnectivityProperties
  public :: ovkEditConnectivityProperties
  public :: ovkReleaseConnectivityProperties
  public :: ovkGetConnectivityReceiverPoints
  public :: ovkGetConnectivityDonorExtents
  public :: ovkGetConnectivityDonorCoords
  public :: ovkGetConnectivityDonorInterpCoefs
  public :: ovkGetConnectivityPropertyDonorGridID
  public :: ovkGetConnectivityPropertyReceiverGridID
  public :: ovkGetConnectivityPropertyDimension
  public :: ovkGetConnectivityPropertyMaxDonorSize
  public :: ovkGetConnectivityPropertyConnectionCount
  public :: ovkDonorSize
  public :: ovkFindDonor

  ! Internal
  public :: ovk_connectivity_
  public :: CreateConnectivity
  public :: DestroyConnectivity
  public :: ResizeConnectivity
  public :: SetConnectivityPropertyMaxDonorSize

  type ovk_connectivity_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: donor_grid_id
    integer :: receiver_grid_id
    integer :: nd
    integer :: max_donor_size
    integer(lk) :: nconnections
  end type ovk_connectivity_properties

  type ovk_connectivity
    type(t_noconstruct) :: noconstruct
    type(t_existence_flag) :: existence_flag
    type(ovk_connectivity_properties), pointer :: properties
    type(ovk_connectivity_properties), pointer :: prev_properties
    integer :: properties_edit_ref_count
    type(t_logger), pointer :: logger
    integer, dimension(:,:), pointer :: receiver_points
    integer, dimension(:,:,:), pointer :: donor_extents
    real(rk), dimension(:,:), pointer :: donor_coords
    real(rk), dimension(:,:,:), pointer :: donor_interp_coefs
  end type ovk_connectivity

contains

  pure function ovk_connectivity_() result(Connectivity)

    type(ovk_connectivity) :: Connectivity

    nullify(Connectivity%properties)
    nullify(Connectivity%prev_properties)
    Connectivity%properties_edit_ref_count = 0
    nullify(Connectivity%logger)
    nullify(Connectivity%receiver_points)
    nullify(Connectivity%donor_extents)
    nullify(Connectivity%donor_coords)
    nullify(Connectivity%donor_interp_coefs)

    call SetExists(Connectivity%existence_flag, .false.)

  end function ovk_connectivity_

  subroutine CreateConnectivity(Connectivity, DonorGridID, ReceiverGridID, Logger, NumDims, &
    MaxDonorSize)

    type(ovk_connectivity), intent(out) :: Connectivity
    integer, intent(in) :: DonorGridID, ReceiverGridID
    type(t_logger), pointer, intent(in) :: Logger
    integer, intent(in) :: NumDims
    integer, intent(in) :: MaxDonorSize

    allocate(Connectivity%properties)
    Connectivity%properties = ovk_connectivity_properties_(NumDims)
    Connectivity%properties%donor_grid_id = DonorGridID
    Connectivity%properties%receiver_grid_id = ReceiverGridID
    Connectivity%properties%max_donor_size = MaxDonorSize

    nullify(Connectivity%prev_properties)

    Connectivity%properties_edit_ref_count = 0

    Connectivity%logger => Logger

    allocate(Connectivity%receiver_points(MAX_ND,0))
    allocate(Connectivity%donor_extents(MAX_ND,2,0))
    allocate(Connectivity%donor_coords(NumDims,0))
    allocate(Connectivity%donor_interp_coefs(MaxDonorSize,NumDims,0))

    call SetExists(Connectivity%existence_flag, .true.)

  end subroutine CreateConnectivity

  subroutine DestroyConnectivity(Connectivity)

    type(ovk_connectivity), intent(inout) :: Connectivity

    if (.not. ovkConnectivityExists(Connectivity)) return

    call SetExists(Connectivity%existence_flag, .false.)

    deallocate(Connectivity%receiver_points)
    deallocate(Connectivity%donor_extents)
    deallocate(Connectivity%donor_coords)
    deallocate(Connectivity%donor_interp_coefs)

    deallocate(Connectivity%properties)
    if (associated(Connectivity%prev_properties)) deallocate(Connectivity%prev_properties)

  end subroutine DestroyConnectivity

  function ovkConnectivityExists(Connectivity) result(Exists)

    type(ovk_connectivity), intent(in) :: Connectivity
    logical :: Exists

    Exists = CheckExists(Connectivity%existence_flag)

  end function ovkConnectivityExists

  subroutine ovkGetConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(in) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(out) :: Properties

    Properties => Connectivity%properties

  end subroutine ovkGetConnectivityProperties

  subroutine ovkEditConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(inout) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(out) :: Properties

    logical :: Success, StartEdit

    call TryEditProperties(Connectivity, Success, StartEdit)

    if (Success) then
      if (StartEdit) then
        allocate(Connectivity%prev_properties)
        Connectivity%prev_properties = Connectivity%properties
      end if
      Properties => Connectivity%properties
    else
      ! Can't fail right now
    end if

  end subroutine ovkEditConnectivityProperties

  subroutine ovkReleaseConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(inout) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(inout) :: Properties

    logical :: Success, EndEdit
    type(ovk_connectivity_properties), pointer :: PrevProperties

    if (associated(Properties, Connectivity%properties)) then

      call TryReleaseProperties(Connectivity, Success, EndEdit)

      if (Success) then
        PrevProperties => Connectivity%prev_properties
        nullify(Connectivity%prev_properties)
        if (Properties%max_donor_size /= PrevProperties%max_donor_size) then
          call ResizeConnectivity(Connectivity, 0_lk)
        end if
        deallocate(PrevProperties)
      else
        if (OVK_DEBUG) then
          write (ERROR_UNIT, '(2a)') "ERROR: Unable to release properties; not ", &
            "currently being edited."
          stop 1
        end if
      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release properties; invalid pointer."
        stop 1
      end if

    end if

    nullify(Properties)

  end subroutine ovkReleaseConnectivityProperties

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

  function EditingProperties(Connectivity) result(Editing)

    type(ovk_connectivity), intent(in) :: Connectivity
    logical :: Editing

    Editing = Connectivity%properties_edit_ref_count > 0

  end function EditingProperties

  subroutine TryEditProperties(Connectivity, Success, StartEdit)

    type(ovk_connectivity), intent(inout) :: Connectivity
    logical, intent(out) :: Success
    logical, intent(out) :: StartEdit

    Success = .true.

    if (Success) then
      StartEdit = Connectivity%properties_edit_ref_count == 0
      Connectivity%properties_edit_ref_count = Connectivity%properties_edit_ref_count + 1
    else
      StartEdit = .false.
    end if

  end subroutine TryEditProperties

  subroutine TryReleaseProperties(Connectivity, Success, EndEdit)

    type(ovk_connectivity), intent(inout) :: Connectivity
    logical, intent(out) :: Success
    logical, intent(out) :: EndEdit

    Success = EditingProperties(Connectivity)

    if (Success) then
      Connectivity%properties_edit_ref_count = Connectivity%properties_edit_ref_count - 1
      EndEdit = Connectivity%properties_edit_ref_count == 0
    else
      EndEdit = .false.
    end if

  end subroutine TryReleaseProperties

  subroutine ResizeConnectivity(Connectivity, NumConnections)

    type(ovk_connectivity), intent(inout) :: Connectivity
    integer(lk), intent(in) :: NumConnections

    deallocate(Connectivity%receiver_points)
    deallocate(Connectivity%donor_extents)
    deallocate(Connectivity%donor_coords)
    deallocate(Connectivity%donor_interp_coefs)

    allocate(Connectivity%receiver_points(MAX_ND,NumConnections))
    allocate(Connectivity%donor_extents(MAX_ND,2,NumConnections))
    allocate(Connectivity%donor_coords(Connectivity%properties%nd,NumConnections))
    allocate(Connectivity%donor_interp_coefs(Connectivity%properties%max_donor_size, &
      Connectivity%properties%nd,NumConnections))

    Connectivity%properties%nconnections = NumConnections

  end subroutine ResizeConnectivity

  function ovk_connectivity_properties_(NumDims) result(Properties)

    integer, intent(in) :: NumDims
    type(ovk_connectivity_properties) :: Properties

    Properties%donor_grid_id = 0
    Properties%receiver_grid_id = 0
    Properties%nd = NumDims
    Properties%max_donor_size = 1
    Properties%nconnections = 0

  end function ovk_connectivity_properties_

  subroutine ovkGetConnectivityPropertyDonorGridID(Properties, DonorGridID)

    type(ovk_connectivity_properties), intent(in) :: Properties
    integer, intent(out) :: DonorGridID

    DonorGridID = Properties%donor_grid_id

  end subroutine ovkGetConnectivityPropertyDonorGridID

  subroutine ovkGetConnectivityPropertyReceiverGridID(Properties, ReceiverGridID)

    type(ovk_connectivity_properties), intent(in) :: Properties
    integer, intent(out) :: ReceiverGridID

    ReceiverGridID = Properties%receiver_grid_id

  end subroutine ovkGetConnectivityPropertyReceiverGridID

  subroutine ovkGetConnectivityPropertyDimension(Properties, NumDims)

    type(ovk_connectivity_properties), intent(in) :: Properties
    integer, intent(out) :: NumDims

    NumDims = Properties%nd

  end subroutine ovkGetConnectivityPropertyDimension

  subroutine ovkGetConnectivityPropertyMaxDonorSize(Properties, MaxDonorSize)

    type(ovk_connectivity_properties), intent(in) :: Properties
    integer, intent(out) :: MaxDonorSize

    MaxDonorSize = Properties%max_donor_size

  end subroutine ovkGetConnectivityPropertyMaxDonorSize

  subroutine SetConnectivityPropertyMaxDonorSize(Properties, MaxDonorSize)

    type(ovk_connectivity_properties), intent(inout) :: Properties
    integer, intent(in) :: MaxDonorSize

    Properties%max_donor_size = MaxDonorSize

  end subroutine SetConnectivityPropertyMaxDonorSize

  subroutine ovkGetConnectivityPropertyConnectionCount(Properties, ConnectionCount)

    type(ovk_connectivity_properties), intent(in) :: Properties
    integer(lk), intent(out) :: ConnectionCount

    ConnectionCount = Properties%nconnections

  end subroutine ovkGetConnectivityPropertyConnectionCount

  function ovkDonorSize(NumDims, ConnectionType) result(DonorSize)

    integer, intent(in) :: NumDims
    integer, intent(in) :: ConnectionType
    integer, dimension(NumDims) :: DonorSize

    select case (ConnectionType)
    case (OVK_CONNECTION_NEAREST)
      DonorSize = 1
    case (OVK_CONNECTION_LINEAR)
      DonorSize = 2
    case (OVK_CONNECTION_CUBIC)
      DonorSize = 4
    case default
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid connection type."
        stop 1
      end if
    end select

  end function ovkDonorSize

  subroutine ovkFindDonor(DonorGrid, ReceiverGrid, Overlap, ReceiverPoint, OverlapIndex, &
    ConnectionType, DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, dimension(:), intent(in) :: ReceiverPoint
    integer(lk), intent(in) :: OverlapIndex
    integer, intent(in) :: ConnectionType
    integer, dimension(:,:), intent(out) :: DonorExtents
    real(rk), dimension(:), intent(out) :: DonorCoords
    real(rk), dimension(:,:), intent(out) :: DonorInterpCoefs

    select case (ConnectionType)
    case (OVK_CONNECTION_NEAREST)
      call FindDonorNearest(DonorGrid, ReceiverGrid, Overlap, ReceiverPoint, OverlapIndex, &
        DonorExtents, DonorCoords, DonorInterpCoefs)
    case (OVK_CONNECTION_LINEAR)
      call FindDonorLinear(DonorGrid, ReceiverGrid, Overlap, ReceiverPoint, OverlapIndex, &
        DonorExtents, DonorCoords, DonorInterpCoefs)
    case (OVK_CONNECTION_CUBIC)
      call FindDonorCubic(DonorGrid, ReceiverGrid, Overlap, ReceiverPoint, OverlapIndex, &
        DonorExtents, DonorCoords, DonorInterpCoefs)
    case default
      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid connection type."
        stop 1
      end if
    end select

  end subroutine ovkFindDonor

  subroutine FindDonorNearest(DonorGrid, ReceiverGrid, Overlap, ReceiverPoint, OverlapIndex, &
    DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
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

    DonorExtents(:,1) = NearestPoint
    DonorExtents(:,2) = NearestPoint

    DonorCoords = 0._rk

    DonorInterpCoefs = 0._rk
    do d = 1, NumDims
      DonorInterpCoefs(1,d) = 1._rk
    end do

  end subroutine FindDonorNearest

  subroutine FindDonorLinear(DonorGrid, ReceiverGrid, Overlap, ReceiverPoint, OverlapIndex, &
    DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, dimension(:), intent(in) :: ReceiverPoint
    integer(lk), intent(in) :: OverlapIndex
    integer, dimension(:,:), intent(out) :: DonorExtents
    real(rk), dimension(:), intent(out) :: DonorCoords
    real(rk), dimension(:,:), intent(out) :: DonorInterpCoefs

    integer :: d
    integer :: NumDims

    NumDims = DonorGrid%nd

    DonorExtents(:,1) = Overlap%cells(:,OverlapIndex)
    DonorExtents(:NumDims,2) = Overlap%cells(:NumDims,OverlapIndex) + 1
    DonorExtents(NumDims+1:,2) = 1

    DonorCoords = Overlap%coords(:,OverlapIndex)

    DonorInterpCoefs = 0._rk
    do d = 1, NumDims
      DonorInterpCoefs(:2,d) = ovkInterpBasisLinear(DonorCoords(d))
    end do

  end subroutine FindDonorLinear

  subroutine FindDonorCubic(DonorGrid, ReceiverGrid, Overlap, ReceiverPoint, OverlapIndex, &
    DonorExtents, DonorCoords, DonorInterpCoefs)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_overlap), intent(in) :: Overlap
    integer, dimension(:), intent(in) :: ReceiverPoint
    integer(lk), intent(in) :: OverlapIndex
    integer, dimension(:,:), intent(out) :: DonorExtents
    real(rk), dimension(:), intent(out) :: DonorCoords
    real(rk), dimension(:,:), intent(out) :: DonorInterpCoefs

    integer :: d, i, j, k
    integer :: NumDims
    real(rk), dimension(DonorGrid%nd) :: ReceiverCoords
    integer, dimension(MAX_ND) :: OverlapCell
    real(rk), dimension(DonorGrid%nd) :: OverlapCoords
    logical :: Success
    integer :: NumWarnings
    integer, dimension(MAX_ND) :: CellShift
    integer, dimension(MAX_ND) :: PrevCell, NextCell
    integer :: PrevDistance, NextDistance
    real(rk), dimension(DonorGrid%nd) :: Gradient
    integer, dimension(MAX_ND) :: Offset
    integer, dimension(MAX_ND) :: ShiftedCell
    integer, dimension(MAX_ND) :: BestOffset
    real(rk) :: OffsetQuality, BestOffsetQuality
    integer, dimension(MAX_ND) :: NeighborLower, NeighborUpper
    integer, dimension(MAX_ND) :: Neighbor
    logical :: AwayFromBoundary
    real(rk), dimension(DonorGrid%nd) :: Displacement
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

    if (DonorGrid%cell_edge_dist%values(OverlapCell(1),OverlapCell(2),OverlapCell(3)) > 1) then
      CellShift = 0
      Success = .true.
    else
      ! Choose a neighboring cell in the direction that most closely matches the edge
      ! distance gradient
      do d = 1, NumDims
        PrevCell = OverlapCell
        PrevCell(d) = PrevCell(d) - 1
        PrevCell(:NumDims) = ovkCartPeriodicAdjust(DonorGrid%cell_edge_dist%cart, PrevCell)
        NextCell = OverlapCell
        NextCell(d) = NextCell(d) + 1
        NextCell(:NumDims) = ovkCartPeriodicAdjust(DonorGrid%cell_edge_dist%cart, NextCell)
        PrevDistance = DonorGrid%cell_edge_dist%values(PrevCell(1),PrevCell(2),PrevCell(3))
        NextDistance = DonorGrid%cell_edge_dist%values(NextCell(1),NextCell(2),NextCell(3))
        Gradient(d) = real(NextDistance-PrevDistance,kind=rk)/2._rk
      end do
      Offset(:NumDims) = ClosestOffset(Gradient)
      Offset(NumDims+1:) = 0
      ShiftedCell = OverlapCell + Offset
      ShiftedCell(:NumDims) = ovkCartPeriodicAdjust(DonorGrid%cell_cart, ShiftedCell)
      if (DonorGrid%cell_edge_dist%values(ShiftedCell(1),ShiftedCell(2),ShiftedCell(3)) > 1) then
        CellShift = Offset
        Success = .true.
      else
        ! Fall back to iterating over all neighbors
        BestOffset = 0
        BestOffsetQuality = -huge(0._rk)
        NeighborLower(:NumDims) = OverlapCell(:NumDims)-1
        NeighborLower(NumDims+1:) = 1
        NeighborUpper(:NumDims) = OverlapCell(:NumDims)+1
        NeighborUpper(NumDims+1:) = 1
        AwayFromBoundary = ovkCartContains(DonorGrid%cell_cart, NeighborLower) .and. &
          ovkCartContains(DonorGrid%cell_cart, NeighborUpper)
        if (AwayFromBoundary) then
          do k = NeighborLower(3), NeighborUpper(3)
            do j = NeighborLower(2), NeighborUpper(2)
              do i = NeighborLower(1), NeighborUpper(1)
                Neighbor = [i,j,k]
                if (DonorGrid%cell_edge_dist%values(Neighbor(1),Neighbor(2),Neighbor(3)) > 1) then
                  Offset = [i-OverlapCell(1),j-OverlapCell(2),k-OverlapCell(3)]
                  Displacement = real(Offset(:NumDims),kind=rk)
                  OffsetQuality = dot_product(Gradient, Displacement)/sqrt(sum(Displacement**2))
                  if (OffsetQuality > BestOffsetQuality) then
                    BestOffset = Offset
                    BestOffsetQuality = OffsetQuality
                  end if
                end if
              end do
            end do
          end do
        else
          do k = NeighborLower(3), NeighborUpper(3)
            do j = NeighborLower(2), NeighborUpper(2)
              do i = NeighborLower(1), NeighborUpper(1)
                Neighbor = [i,j,k]
                Neighbor(:NumDims) = ovkCartPeriodicAdjust(DonorGrid%cell_cart, Neighbor)
                if (ovkCartContains(DonorGrid%cell_cart, Neighbor)) then
                  if (DonorGrid%cell_edge_dist%values(Neighbor(1),Neighbor(2),Neighbor(3)) > 1) then
                    Offset = [i-OverlapCell(1),j-OverlapCell(2),k-OverlapCell(3)]
                    Displacement = real(Offset(:NumDims),kind=rk)
                    OffsetQuality = dot_product(Gradient, Displacement)/sqrt(sum(Displacement**2))
                    if (OffsetQuality > BestOffsetQuality) then
                      BestOffset = Offset
                      BestOffsetQuality = OffsetQuality
                    end if
                  end if
                end if
              end do
            end do
          end do
        end if
        if (any(BestOffset /= 0)) then
          CellShift = BestOffset
          Success = .true.
        end if
      end if
    end if

    if (Success) then
      DonorExtents(:NumDims,1) = OverlapCell(:NumDims) + CellShift(:NumDims) - 1
      DonorExtents(:NumDims,1) = ovkCartPeriodicAdjust(DonorGrid%cart, DonorExtents(:,1))
      DonorExtents(:NumDims,2) = DonorExtents(:NumDims,1) + 3
      DonorCoords = OverlapCoords - real(CellShift(:NumDims),kind=rk)
      DonorCoords = ovkCoordsInCubicGridCell(DonorGrid, DonorExtents(:,1), ReceiverCoords, &
        Guess=DonorCoords)
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
      call FindDonorLinear(DonorGrid, ReceiverGrid, Overlap, ReceiverPoint, OverlapIndex, &
        DonorExtents, DonorCoords, DonorInterpCoefs)
    end if

  end subroutine FindDonorCubic

  function ClosestOffset(Vector) result(Offset)

    real(rk), dimension(:), intent(in) :: Vector
    integer, dimension(size(Vector)) :: Offset

    real(rk), dimension(size(Vector)) :: AbsVector
    integer, dimension(size(Vector)) :: VectorSign
    integer, dimension(size(Vector)) :: O

    real(rk), parameter :: Alpha = sqrt(2._rk) - 1._rk
    real(rk), parameter :: Beta = sqrt(1.5_rk) - 1._rk

    AbsVector = abs(Vector)
    VectorSign = merge(1, merge(-1, 0, Vector < 0._rk), Vector > 0._rk)

    call VectorComponentOrder(AbsVector, O)

    Offset = 0

    Offset(O(1)) = VectorSign(O(1))

    if (AbsVector(O(2)) >= Alpha*AbsVector(O(1))) then
      Offset(O(2)) = VectorSign(O(2))
    end if

    if (size(Vector) == 3) then
      if (AbsVector(O(3)) >= Beta*(AbsVector(O(1))+AbsVector(O(2)))) then
        Offset(O(3)) = VectorSign(O(3))
      end if
    end if

  end function ClosestOffset

  subroutine VectorComponentOrder(Vector, Order)

    real(rk), dimension(:), intent(in) :: Vector
    integer, dimension(:), intent(out) :: Order

    integer :: N

    N = size(Vector)

    select case (N)
    case (2)
      Order = [1,2]
      if (Vector(Order(1)) < Vector(Order(2))) call Swap(Order(1),Order(2))
    case (3)
      Order = [1,2,3]
      if (Vector(Order(1)) < Vector(Order(2))) call Swap(Order(1),Order(2))
      if (Vector(Order(1)) < Vector(Order(3))) call Swap(Order(1),Order(3))
      if (Vector(Order(2)) < Vector(Order(3))) call Swap(Order(2),Order(3))
    end select

  end subroutine VectorComponentOrder

  subroutine Swap(X, Y)
    integer, intent(inout) :: X, Y
    integer :: Z
    Z = X
    X = Y
    Y = Z
  end subroutine Swap

end module ovkConnectivity
