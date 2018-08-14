! Copyright (c) 2018 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkConnectivity

  use ovkCart
  use ovkDonorStencil
  use ovkField
  use ovkFieldOps
  use ovkGlobal
  use ovkGrid
  use ovkLogger
  use ovkOverlap
  implicit none

  private

  ! Public API
  public :: ovk_connectivity
  public :: ovkConnectivityExists
  public :: ovkGetConnectivityDonorGrid
  public :: ovkGetConnectivityReceiverGrid
  public :: ovkGetConnectivityDimension
  public :: ovkGetConnectivityMaxDonorSize
  public :: ovkGetConnectivityCount
  public :: ovkResetConnectivity
  public :: ovkGetConnectivityDonorExtents
  public :: ovkGetConnectivityDonorCoords
  public :: ovkGetConnectivityDonorInterpCoefs
  public :: ovkGetConnectivityReceiverPoints

  ! Internal API
  public :: ovk_connectivity_
  public :: CreateConnectivity
  public :: DestroyConnectivity
  public :: FillConnectivity

  type ovk_connectivity
    type(t_noconstruct) :: noconstruct
    type(t_existence_flag) :: existence_flag
    type(t_logger) :: logger
    type(ovk_grid), pointer :: donor_grid
    type(ovk_grid), pointer :: receiver_grid
    integer :: nd
    integer :: max_donor_size
    integer(lk) :: nconnections
    integer, dimension(:,:,:), pointer :: donor_extents
    real(rk), dimension(:,:), pointer :: donor_coords
    real(rk), dimension(:,:,:), pointer :: donor_interp_coefs
    integer, dimension(:,:), pointer :: receiver_points
  end type ovk_connectivity

contains

  pure function ovk_connectivity_() result(Connectivity)

    type(ovk_connectivity) :: Connectivity

    Connectivity%logger = t_logger_()
    nullify(Connectivity%donor_grid)
    nullify(Connectivity%receiver_grid)
    Connectivity%nd = 2
    Connectivity%max_donor_size = 1
    Connectivity%nconnections = 0_lk
    nullify(Connectivity%donor_extents)
    nullify(Connectivity%donor_coords)
    nullify(Connectivity%donor_interp_coefs)
    nullify(Connectivity%receiver_points)

    call SetExists(Connectivity%existence_flag, .false.)

  end function ovk_connectivity_

  subroutine CreateConnectivity(Connectivity, Logger, DonorGrid, ReceiverGrid)

    type(ovk_connectivity), intent(out) :: Connectivity
    type(t_logger), intent(in) :: Logger
    type(ovk_grid), pointer, intent(in) :: DonorGrid, ReceiverGrid

    integer :: NumDims

    NumDims = DonorGrid%nd

    Connectivity%logger = Logger
    Connectivity%donor_grid => DonorGrid
    Connectivity%receiver_grid => ReceiverGrid
    Connectivity%nd = NumDims
    Connectivity%max_donor_size = 1
    Connectivity%nconnections = 0_lk

    allocate(Connectivity%donor_extents(MAX_DIMS,2,0))
    allocate(Connectivity%donor_coords(NumDims,0))
    allocate(Connectivity%donor_interp_coefs(1,NumDims,0))
    allocate(Connectivity%receiver_points(MAX_DIMS,0))

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

  subroutine ovkGetConnectivityDonorGrid(Connectivity, DonorGrid)

    type(ovk_connectivity), intent(in) :: Connectivity
    type(ovk_grid), pointer, intent(out) :: DonorGrid

    DonorGrid => Connectivity%donor_grid

  end subroutine ovkGetConnectivityDonorGrid

  subroutine ovkGetConnectivityReceiverGrid(Connectivity, ReceiverGrid)

    type(ovk_connectivity), intent(in) :: Connectivity
    type(ovk_grid), pointer, intent(out) :: ReceiverGrid

    ReceiverGrid => Connectivity%receiver_grid

  end subroutine ovkGetConnectivityReceiverGrid

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

  subroutine ovkResetConnectivity(Connectivity, NumConnections, MaxDonorSize)

    type(ovk_connectivity), intent(inout) :: Connectivity
    integer(lk), intent(in) :: NumConnections
    integer, intent(in) :: MaxDonorSize

    deallocate(Connectivity%donor_extents)
    deallocate(Connectivity%donor_coords)
    deallocate(Connectivity%donor_interp_coefs)
    deallocate(Connectivity%receiver_points)

    Connectivity%max_donor_size = MaxDonorSize
    Connectivity%nconnections = NumConnections

    allocate(Connectivity%donor_extents(MAX_DIMS,2,NumConnections))
    allocate(Connectivity%donor_coords(Connectivity%nd,NumConnections))
    allocate(Connectivity%donor_interp_coefs(MaxDonorSize,Connectivity%nd,NumConnections))
    allocate(Connectivity%receiver_points(MAX_DIMS,NumConnections))

    Connectivity%donor_extents(Connectivity%nd+1:,:,:) = 1
    Connectivity%donor_interp_coefs = 0._rk
    Connectivity%receiver_points(Connectivity%nd+1:,:) = 1

  end subroutine ovkResetConnectivity

  subroutine FillConnectivity(Connectivity, Overlap, DonorStencil, ReceiverMask)

    type(ovk_connectivity), intent(inout) :: Connectivity
    type(ovk_overlap), intent(in) :: Overlap
    type(t_donor_stencil), intent(in) :: DonorStencil
    type(ovk_field_logical), intent(in) :: ReceiverMask

    integer :: i, j, k
    integer(lk) :: l, p
    integer(lk) :: NumConnections
    integer, dimension(MAX_DIMS) :: DonorSize
    integer(lk), dimension(:), allocatable :: OverlapIndices

    NumConnections = ovkCountMask(ReceiverMask)

    if (NumConnections > 0_lk) then

      DonorSize = 1
      call GetDonorStencilSize(DonorStencil, DonorSize)

      call ovkResetConnectivity(Connectivity, NumConnections, maxval(DonorSize))

      allocate(OverlapIndices(NumConnections))

      l = 1_lk
      p = 1_lk
      do k = Connectivity%receiver_grid%cart%is(3), Connectivity%receiver_grid%cart%ie(3)
        do j = Connectivity%receiver_grid%cart%is(2), Connectivity%receiver_grid%cart%ie(2)
          do i = Connectivity%receiver_grid%cart%is(1), Connectivity%receiver_grid%cart%ie(1)
            if (ReceiverMask%values(i,j,k)) then
              Connectivity%receiver_points(:,l) = [i,j,k]
              OverlapIndices(l) = p
              l = l + 1_lk
            end if
            if (Overlap%mask%values(i,j,k)) then
              p = p + 1_lk
            end if
          end do
        end do
      end do

      call FindDonors(Connectivity%donor_grid, Connectivity%receiver_grid, Overlap, DonorStencil, &
        NumConnections, Connectivity%receiver_points, OverlapIndices, Connectivity%donor_extents, &
        Connectivity%donor_coords, Connectivity%donor_interp_coefs)

    end if

  end subroutine FillConnectivity

end module ovkConnectivity
