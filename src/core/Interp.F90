! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkInterp

  use ovkCart
  use ovkDonors
  use ovkGeometry
  use ovkGlobal
  use ovkField
  implicit none

  private

  ! API
  public :: ovk_interp
  public :: ovk_interp_
  public :: ovkMakeInterpData
  public :: ovkDestroyInterpData
  public :: ovkDonorGridIDToIBlank
  public :: ovkGenerateInterpData

  type ovk_interp
    type(ovk_cart) :: cart
    integer :: ncoefs
    type(ovk_field_logical) :: valid_mask
    type(ovk_field_int) :: donor_grid_ids
    type(ovk_field_int), dimension(:), allocatable :: donor_cells
    type(ovk_field_real), dimension(:), allocatable :: donor_cell_coords
    type(ovk_field_int) :: schemes
    type(ovk_field_real), dimension(:,:), allocatable :: coefs
  end type ovk_interp

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_interp_
    module procedure ovk_interp_Default
  end interface ovk_interp_

contains

  pure function ovk_interp_Default() result(InterpData)

    type(ovk_interp) :: InterpData

    InterpData%cart = ovk_cart_(2)
    InterpData%ncoefs = 0
    InterpData%valid_mask = ovk_field_logical_(2)
    InterpData%donor_grid_ids = ovk_field_int_(2)
    InterpData%schemes = ovk_field_int_(2)

  end function ovk_interp_Default

  subroutine ovkMakeInterpData(InterpData, Cart, InterpScheme)

    type(ovk_interp), intent(out) :: InterpData
    type(ovk_cart), intent(in) :: Cart
    integer, intent(in), optional :: InterpScheme

    integer :: i, j

    InterpData%cart = ovkCartConvertPeriodicStorage(Cart, OVK_NO_OVERLAP_PERIODIC)

    if (present(InterpScheme)) then
      select case (InterpScheme)
      case (OVK_INTERP_LINEAR)
        InterpData%ncoefs = 2
      case (OVK_INTERP_CUBIC)
        InterpData%ncoefs = 4
      end select
    else
      InterpData%ncoefs = 2
    end if

    InterpData%valid_mask = ovk_field_logical_(InterpData%cart)
    InterpData%donor_grid_ids = ovk_field_int_(InterpData%cart)

    allocate(InterpData%donor_cells(InterpData%cart%nd))
    do i = 1, InterpData%cart%nd
      InterpData%donor_cells(i) = ovk_field_int_(InterpData%cart)
    end do

    allocate(InterpData%donor_cell_coords(InterpData%cart%nd))
    do i = 1, InterpData%cart%nd
      InterpData%donor_cell_coords(i) = ovk_field_real_(InterpData%cart)
    end do

    InterpData%schemes = ovk_field_int_(InterpData%cart)

    allocate(InterpData%coefs(InterpData%ncoefs,InterpData%cart%nd))
    do j = 1, InterpData%cart%nd
      do i = 1, InterpData%ncoefs
        InterpData%coefs(i,j) = ovk_field_real_(InterpData%cart)
      end do
    end do

  end subroutine ovkMakeInterpData

  subroutine ovkDestroyInterpData(InterpData)

    type(ovk_interp), intent(inout) :: InterpData

    InterpData%valid_mask = ovk_field_logical_(2)
    InterpData%donor_grid_ids = ovk_field_int_(2)

    if (allocated(InterpData%donor_cells)) deallocate(InterpData%donor_cells)
    if (allocated(InterpData%donor_cell_coords)) deallocate(InterpData%donor_cell_coords)

    InterpData%schemes = ovk_field_int_(2)

    if (allocated(InterpData%coefs)) deallocate(InterpData%coefs)

  end subroutine ovkDestroyInterpData

  subroutine ovkDonorGridIDToIBlank(InterpData, IBlank, Multiplier, Offset)

    type(ovk_interp), intent(in) :: InterpData
    type(ovk_field_int), intent(inout) :: IBlank
    integer, intent(in), optional :: Multiplier
    integer, intent(in), optional :: Offset

    integer :: i, j, k
    integer :: Multiplier_
    integer :: Offset_

    if (present(Multiplier)) then
      Multiplier_ = Multiplier
    else
      Multiplier_ = 1
    end if

    if (present(Offset)) then
      Offset_ = Offset
    else
      Offset_ = 0
    end if

    do k = InterpData%cart%is(3), InterpData%cart%ie(3)
      do j = InterpData%cart%is(2), InterpData%cart%ie(2)
        do i = InterpData%cart%is(1), InterpData%cart%ie(1)
          if (InterpData%valid_mask%values(i,j,k)) then
            IBlank%values(i,j,k) = Multiplier_ * InterpData%donor_grid_ids%values(i,j,k) + Offset_
          end if
        end do
      end do
    end do

  end subroutine ovkDonorGridIDToIBlank

  subroutine ovkGenerateInterpData(Donors, InterpData, InterpScheme)

    type(ovk_donors), intent(in) :: Donors
    type(ovk_interp), intent(out) :: InterpData
    integer, intent(in), optional :: InterpScheme

    integer :: InterpScheme_
    integer :: i, j, k, l, m
    real(rk), dimension(4) :: Basis

    if (present(InterpScheme)) then
      InterpScheme_ = InterpScheme
    else
      InterpScheme_ = OVK_INTERP_LINEAR
    end if

    call ovkMakeInterpData(InterpData, Donors%cart, InterpScheme=InterpScheme)
    InterpData%valid_mask%values = .false.

    select case (InterpScheme_)
    case (OVK_INTERP_LINEAR)
      do k = Donors%cart%is(3), Donors%cart%ie(3)
        do j = Donors%cart%is(2), Donors%cart%ie(2)
          do i = Donors%cart%is(1), Donors%cart%ie(1)
            if (Donors%valid_mask%values(i,j,k)) then
              InterpData%valid_mask%values(i,j,k) = .true.
              InterpData%donor_grid_ids%values(i,j,k) = Donors%grid_ids%values(i,j,k)
              do l = 1, InterpData%cart%nd
                InterpData%donor_cells(l)%values(i,j,k) = Donors%cells(l)%values(i,j,k)
                InterpData%donor_cell_coords(l)%values(i,j,k) = Donors%cell_coords(l)%values(i,j,k)
              end do
              InterpData%schemes%values(i,j,k) = OVK_INTERP_LINEAR
              do l = 1, InterpData%cart%nd
                Basis(:2) = ovkInterpBasisLinear(Donors%cell_coords(l)%values(i,j,k))
                do m = 1, 2
                  InterpData%coefs(m,l)%values(i,j,k) = Basis(m)
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
            if (Donors%valid_mask%values(i,j,k)) then
              InterpData%valid_mask%values(i,j,k) = .true.
              InterpData%donor_grid_ids%values(i,j,k) = Donors%grid_ids%values(i,j,k)
              do l = 1, Donors%cart%nd
                InterpData%donor_cells(l)%values(i,j,k) = Donors%cells(l)%values(i,j,k)
                InterpData%donor_cell_coords(l)%values(i,j,k) = Donors%cell_coords(l)%values(i,j,k)
              end do
              if (Donors%cell_extents%values(i,j,k) == 4) then
                InterpData%schemes%values(i,j,k) = OVK_INTERP_CUBIC
                do l = 1, InterpData%cart%nd
                  Basis = ovkInterpBasisCubic(Donors%cell_coords(l)%values(i,j,k))
                  do m = 1, 4
                    InterpData%coefs(m,l)%values(i,j,k) = Basis(m)
                  end do
                end do
              else if (Donors%cell_extents%values(i,j,k) == 2) then
                InterpData%schemes%values(i,j,k) = OVK_INTERP_LINEAR
                do l = 1, InterpData%cart%nd
                  Basis(:2) = ovkInterpBasisLinear(Donors%cell_coords(l)%values(i,j,k))
                  do m = 1, 2
                    InterpData%coefs(m,l)%values(i,j,k) = Basis(m)
                  end do
                  do m = 3, 4
                    InterpData%coefs(m,l)%values(i,j,k) = 0._rk
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

  end subroutine ovkGenerateInterpData

end module ovkInterp
