! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkInterp

  use ovkCart
  use ovkDonors
  use ovkField
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  implicit none

  private

  ! API
  public :: ovk_interp
  public :: ovk_interp_
  public :: ovk_interp_properties
  public :: ovk_interp_properties_
  public :: ovkCreateInterpData
  public :: ovkDestroyInterpData
  public :: ovkUpdateInterpData
  public :: ovkFillInterpData
  public :: ovkGetInterpDataReceiverMask
  public :: ovkGetInterpDataOrphanMask
  public :: ovkGetInterpDataDonorGridIDs
  public :: ovkGetInterpDataDonorCells
  public :: ovkGetInterpDataDonorCellCoords
  public :: ovkGetInterpDataSchemes
  public :: ovkGetInterpDataCoefs
  public :: ovkGetInterpDataPropertyVerbose
  public :: ovkSetInterpDataPropertyVerbose

  type ovk_interp_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: id
    integer :: nd
    integer, dimension(MAX_ND) :: npoints
    logical, dimension(MAX_ND) :: periodic
    integer :: periodic_storage
    integer :: stencil_size
    ! Read/write
    logical :: verbose
  end type ovk_interp_properties

  type ovk_interp
    type(t_noconstruct) :: noconstruct
    type(ovk_interp_properties), pointer :: properties
    type(ovk_cart) :: cart
    type(ovk_field_logical), pointer :: receiver_mask
    type(ovk_field_logical), pointer :: orphan_mask
    type(ovk_field_int), pointer :: donor_grid_ids
    type(ovk_field_int), dimension(:), pointer :: donor_cells
    type(ovk_field_real), dimension(:), pointer :: donor_cell_coords
    type(ovk_field_int), pointer :: schemes
    type(ovk_field_real), dimension(:,:), pointer :: coefs
    logical :: editing_properties
  end type ovk_interp

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_interp_
    module procedure ovk_interp_Default
  end interface ovk_interp_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_interp_properties_
    module procedure ovk_interp_properties_Default
  end interface ovk_interp_properties_

contains

  pure function ovk_interp_Default() result(InterpData)

    type(ovk_interp) :: InterpData

    nullify(InterpData%properties)
    InterpData%cart = ovk_cart_(2)
    nullify(InterpData%receiver_mask)
    nullify(InterpData%orphan_mask)
    nullify(InterpData%donor_grid_ids)
    nullify(InterpData%donor_cells)
    nullify(InterpData%donor_cell_coords)
    nullify(InterpData%schemes)
    nullify(InterpData%coefs)
    InterpData%editing_properties = .false.

  end function ovk_interp_Default

  subroutine ovkCreateInterpData(InterpData, Grid, StencilSize, Verbose)

    type(ovk_interp), intent(out) :: InterpData
    type(ovk_grid), intent(in) :: Grid
    integer, intent(in) :: StencilSize
    logical, intent(in), optional :: Verbose

    logical :: Verbose_
    integer :: d, l

    if (present(Verbose)) then
      Verbose_ = Verbose
    else
      Verbose_ = .false.
    end if

    allocate(InterpData%properties)
    InterpData%properties = ovk_interp_properties_()
    InterpData%properties%id = Grid%properties%id
    InterpData%properties%nd = Grid%properties%nd
    InterpData%properties%npoints = Grid%properties%npoints
    InterpData%properties%periodic = Grid%properties%periodic
    InterpData%properties%periodic_storage = Grid%properties%periodic_storage
    InterpData%properties%stencil_size = StencilSize
    InterpData%properties%verbose = Verbose_

    InterpData%cart = Grid%cart

    allocate(InterpData%receiver_mask)
    InterpData%receiver_mask = ovk_field_logical_(InterpData%cart, .false.)

    allocate(InterpData%orphan_mask)
    InterpData%orphan_mask = ovk_field_logical_(InterpData%cart, .false.)

    allocate(InterpData%donor_grid_ids)
    InterpData%donor_grid_ids = ovk_field_int_(InterpData%cart, 0)

    allocate(InterpData%donor_cells(InterpData%cart%nd))
    do d = 1, InterpData%cart%nd
      InterpData%donor_cells(d) = ovk_field_int_(InterpData%cart, 0)
    end do

    allocate(InterpData%donor_cell_coords(InterpData%cart%nd))
    do d = 1, InterpData%cart%nd
      InterpData%donor_cell_coords(d) = ovk_field_real_(InterpData%cart, 0._rk)
    end do

    allocate(InterpData%schemes)
    InterpData%schemes = ovk_field_int_(InterpData%cart, OVK_INTERP_LINEAR)

    allocate(InterpData%coefs(StencilSize,InterpData%cart%nd))
    do d = 1, InterpData%cart%nd
      do l = 1, StencilSize
        InterpData%coefs(l,d) = ovk_field_real_(InterpData%cart, 0._rk)
      end do
    end do

    InterpData%editing_properties = .false.

  end subroutine ovkCreateInterpData

  subroutine ovkDestroyInterpData(InterpData)

    type(ovk_interp), intent(inout) :: InterpData

    if (associated(InterpData%properties)) deallocate(InterpData%properties)

    if (associated(InterpData%receiver_mask)) deallocate(InterpData%receiver_mask)
    if (associated(InterpData%orphan_mask)) deallocate(InterpData%orphan_mask)
    if (associated(InterpData%donor_grid_ids)) deallocate(InterpData%donor_grid_ids)
    if (associated(InterpData%donor_cells)) deallocate(InterpData%donor_cells)
    if (associated(InterpData%donor_cell_coords)) deallocate(InterpData%donor_cell_coords)
    if (associated(InterpData%schemes)) deallocate(InterpData%schemes)
    if (associated(InterpData%coefs)) deallocate(InterpData%coefs)

  end subroutine ovkDestroyInterpData

  subroutine ovkUpdateInterpData(InterpData)

    type(ovk_interp), intent(inout) :: InterpData

    logical :: CannotUpdate

    CannotUpdate = InterpData%editing_properties

    if (CannotUpdate) then

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Cannot update interpolation data; still being edited."
        stop 1
      end if

      return

    else

      ! Nothing here at the moment

    end if

  end subroutine ovkUpdateInterpData

  subroutine ovkFillInterpData(InterpData, Donors, OrphanMask, InterpScheme)

    type(ovk_interp), intent(inout) :: InterpData
    type(ovk_donors), intent(in) :: Donors
    type(ovk_field_logical), intent(in) :: OrphanMask
    integer, intent(in), optional :: InterpScheme

    integer :: InterpScheme_
    integer :: i, j, k, d, l
    real(rk), dimension(4) :: Basis

    if (present(InterpScheme)) then
      InterpScheme_ = InterpScheme
    else
      InterpScheme_ = OVK_INTERP_LINEAR
    end if

    InterpData%receiver_mask%values = Donors%valid_mask%values
    InterpData%orphan_mask%values = OrphanMask%values

    select case (InterpScheme_)
    case (OVK_INTERP_LINEAR)
      do k = Donors%cart%is(3), Donors%cart%ie(3)
        do j = Donors%cart%is(2), Donors%cart%ie(2)
          do i = Donors%cart%is(1), Donors%cart%ie(1)
            if (InterpData%receiver_mask%values(i,j,k)) then
              InterpData%donor_grid_ids%values(i,j,k) = Donors%grid_ids%values(i,j,k)
              do d = 1, InterpData%cart%nd
                InterpData%donor_cells(d)%values(i,j,k) = Donors%cells(d)%values(i,j,k)
                InterpData%donor_cell_coords(d)%values(i,j,k) = Donors%cell_coords(d)%values(i,j,k)
              end do
              InterpData%schemes%values(i,j,k) = OVK_INTERP_LINEAR
              do d = 1, InterpData%cart%nd
                Basis(:2) = ovkInterpBasisLinear(Donors%cell_coords(d)%values(i,j,k))
                do l = 1, 2
                  InterpData%coefs(l,d)%values(i,j,k) = Basis(l)
                end do
              end do
            end if
          end do
        end do
      end do
    case (OVK_INTERP_CUBIC)
      do k = Donors%cart%is(3), Donors%cart%ie(3)
        do j = Donors%cart%is(2), Donors%cart%ie(2)
          do i = Donors%cart%is(1), Donors%cart%ie(1)
            if (InterpData%receiver_mask%values(i,j,k)) then
              InterpData%donor_grid_ids%values(i,j,k) = Donors%grid_ids%values(i,j,k)
              do d = 1, Donors%cart%nd
                InterpData%donor_cells(d)%values(i,j,k) = Donors%cells(d)%values(i,j,k)
                InterpData%donor_cell_coords(d)%values(i,j,k) = Donors%cell_coords(d)%values(i,j,k)
              end do
              if (Donors%cell_extents%values(i,j,k) == 4) then
                InterpData%schemes%values(i,j,k) = OVK_INTERP_CUBIC
                do d = 1, InterpData%cart%nd
                  Basis = ovkInterpBasisCubic(Donors%cell_coords(d)%values(i,j,k))
                  do l = 1, 4
                    InterpData%coefs(l,d)%values(i,j,k) = Basis(l)
                  end do
                end do
              else if (Donors%cell_extents%values(i,j,k) == 2) then
                InterpData%schemes%values(i,j,k) = OVK_INTERP_LINEAR
                do d = 1, InterpData%cart%nd
                  Basis(:2) = ovkInterpBasisLinear(Donors%cell_coords(d)%values(i,j,k))
                  do l = 1, 2
                    InterpData%coefs(l,d)%values(i,j,k) = Basis(l)
                  end do
                  do l = 3, 4
                    InterpData%coefs(l,d)%values(i,j,k) = 0._rk
                  end do
                end do
              end if
            end if
          end do
        end do
      end do
    case default
      write (ERROR_UNIT, '(a)') "ERROR: Unrecognized interpolation scheme."
      stop 1
    end select

  end subroutine ovkFillInterpData

  subroutine ovkGetInterpDataReceiverMask(InterpData, ReceiverMask)

    type(ovk_interp), intent(in) :: InterpData
    type(ovk_field_logical), pointer, intent(out) :: ReceiverMask

    ReceiverMask => InterpData%receiver_mask

  end subroutine ovkGetInterpDataReceiverMask

  subroutine ovkGetInterpDataOrphanMask(InterpData, OrphanMask)

    type(ovk_interp), intent(in) :: InterpData
    type(ovk_field_logical), pointer, intent(out) :: OrphanMask

    OrphanMask => InterpData%orphan_mask

  end subroutine ovkGetInterpDataOrphanMask

  subroutine ovkGetInterpDataDonorGridIDs(InterpData, DonorGridIDs)

    type(ovk_interp), intent(in) :: InterpData
    type(ovk_field_int), pointer, intent(out) :: DonorGridIDs

    DonorGridIDs => InterpData%donor_grid_ids

  end subroutine ovkGetInterpDataDonorGridIDs

  subroutine ovkGetInterpDataDonorCells(InterpData, DonorCells)

    type(ovk_interp), intent(in) :: InterpData
    type(ovk_field_int), dimension(:), pointer, intent(out) :: DonorCells

    DonorCells => InterpData%donor_cells

  end subroutine ovkGetInterpDataDonorCells

  subroutine ovkGetInterpDataDonorCellCoords(InterpData, DonorCellCoords)

    type(ovk_interp), intent(in) :: InterpData
    type(ovk_field_real), dimension(:), pointer, intent(out) :: DonorCellCoords

    DonorCellCoords => InterpData%donor_cell_coords

  end subroutine ovkGetInterpDataDonorCellCoords

  subroutine ovkGetInterpDataSchemes(InterpData, InterpSchemes)

    type(ovk_interp), intent(in) :: InterpData
    type(ovk_field_int), pointer, intent(out) :: InterpSchemes

    InterpSchemes => InterpData%schemes

  end subroutine ovkGetInterpDataSchemes

  subroutine ovkGetInterpDataCoefs(InterpData, InterpCoefs)

    type(ovk_interp), intent(in) :: InterpData
    type(ovk_field_real), dimension(:,:), pointer, intent(out) :: InterpCoefs

    InterpCoefs => InterpData%coefs

  end subroutine ovkGetInterpDataCoefs

  function ovk_interp_properties_Default() result(Properties)

    type(ovk_interp_properties) :: Properties

    Properties%id = 0
    Properties%nd = 2
    Properties%npoints = [0,0,1]
    Properties%periodic = .false.
    Properties%periodic_storage = OVK_NO_OVERLAP_PERIODIC
    Properties%stencil_size = 2
    Properties%verbose = .false.

  end function ovk_interp_properties_Default

  subroutine ovkGetInterpDataPropertyVerbose(Properties, Verbose)

    type(ovk_interp_properties), intent(in) :: Properties
    logical, intent(out) :: Verbose

    Verbose = Properties%verbose

  end subroutine ovkGetInterpDataPropertyVerbose

  subroutine ovkSetInterpDataPropertyVerbose(Properties, Verbose)

    type(ovk_interp_properties), intent(inout) :: Properties
    logical, intent(in) :: Verbose

    Properties%verbose = Verbose

  end subroutine ovkSetInterpDataPropertyVerbose

end module ovkInterp
