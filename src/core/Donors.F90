! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkDonors

  use ovkBoundingBox
  use ovkCart
  use ovkDonorAccel
  use ovkField
  use ovkGlobal
  use ovkGrid
  use ovkHashGrid
  use ovkMask
  implicit none

  private

  ! API
  public :: ovk_donors
  public :: ovk_donors_
  public :: ovkMakeDonors
  public :: ovkDestroyDonors
  public :: ovkFindDonors
  public :: ovkChooseDonors
  public :: ovkMergeDonors
  public :: ovkPrintDonors
  public :: ovkGenerateReceiverMask
  public :: ovkGenerateDonorMask
  public :: ovkGenerateOverlapMask
  public :: ovkGenerateCoarseToFineMask
  public :: ovkGenerateCyclicReceiverMask
  public :: ovkGenerateNearCrossoverMask
  public :: ovkGenerateOrphanMask

  type ovk_donors
    type(ovk_cart) :: cart
    type(ovk_field_logical) :: valid_mask
    type(ovk_field_int) :: grid_ids
    type(ovk_field_int), dimension(:), allocatable :: cells
    type(ovk_field_int) :: cell_extents
    type(ovk_field_real), dimension(:), allocatable :: cell_coords
    type(ovk_field_real) :: cell_diff_params
  end type ovk_donors

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_donors_
    module procedure ovk_donors_Default
  end interface ovk_donors_

contains

  pure function ovk_donors_Default(nDims) result(Donors)

    integer, intent(in) :: nDims
    type(ovk_donors) :: Donors

    Donors%cart = ovk_cart_(nDims)
    Donors%valid_mask = ovk_field_logical_(nDims)
    Donors%grid_ids = ovk_field_int_(nDims)
    Donors%cell_extents = ovk_field_int_(nDims)
    Donors%cell_diff_params = ovk_field_real_(nDims)

  end function ovk_donors_Default

  subroutine ovkMakeDonors(Donors, Cart)

    type(ovk_donors), intent(out) :: Donors
    type(ovk_cart), intent(in) :: Cart

    integer :: i

    Donors%cart = ovkCartConvertPeriodicStorage(Cart, OVK_NO_OVERLAP_PERIODIC)

    Donors%valid_mask = ovk_field_logical_(Donors%cart)
    Donors%grid_ids = ovk_field_int_(Donors%cart)

    allocate(Donors%cells(Donors%cart%nd))
    do i = 1, Donors%cart%nd
      Donors%cells(i) = ovk_field_int_(Donors%cart)
    end do

    Donors%cell_extents = ovk_field_int_(Donors%cart)

    allocate(Donors%cell_coords(Donors%cart%nd))
    do i = 1, Donors%cart%nd
      Donors%cell_coords(i) = ovk_field_real_(Donors%cart)
    end do

    Donors%cell_diff_params = ovk_field_real_(Donors%cart)

  end subroutine ovkMakeDonors

  subroutine ovkDestroyDonors(Donors)

    type(ovk_donors), intent(inout) :: Donors

    Donors%valid_mask = ovk_field_logical_(Donors%cart%nd)
    Donors%grid_ids = ovk_field_int_(Donors%cart%nd)

    deallocate(Donors%cells)

    Donors%cell_extents = ovk_field_int_(Donors%cart%nd)

    deallocate(Donors%cell_coords)

    Donors%cell_diff_params = ovk_field_real_(Donors%cart%nd)

  end subroutine ovkDestroyDonors

  subroutine ovkFindDonors(DonorGrid, ReceiverGrid, DonorAccel, Donors, ReceiverSubset)

    type(ovk_grid), intent(in) :: DonorGrid, ReceiverGrid
    type(ovk_donor_accel), intent(in) :: DonorAccel
    type(ovk_donors), intent(out) :: Donors
    type(ovk_field_logical), intent(in), optional :: ReceiverSubset

    integer :: i, j, k, l
    type(ovk_bbox) :: Bounds
    logical :: IncludePoint
    integer, dimension(DonorGrid%cart%nd) :: DonorCell
    integer, dimension(MAX_ND) :: ReceiverPoint
    real(rk), dimension(ReceiverGrid%cart%nd) :: ReceiverCoords
    logical :: ValidReceiverPoint
    real(rk), dimension(DonorGrid%cart%nd) :: DonorCellCoords
    real(rk) :: DonorCellSize, ReceiverCellSize

    call ovkMakeDonors(Donors, ReceiverGrid%cart)
    Donors%valid_mask%values = .false.

    Bounds = ovkBBScale(ovkBBIntersect(DonorGrid%bounds, ReceiverGrid%bounds), 1.001_rk)

    if (ovkBBIsEmpty(Bounds)) then
      return
    end if

    do k = ReceiverGrid%cart%is(3), ReceiverGrid%cart%ie(3)
      do j = ReceiverGrid%cart%is(2), ReceiverGrid%cart%ie(2)
        do i = ReceiverGrid%cart%is(1), ReceiverGrid%cart%ie(1)
          if (present(ReceiverSubset)) then
            IncludePoint = ReceiverSubset%values(i,j,k)
          else
            IncludePoint = .true.
          end if
          if (IncludePoint) then
            ReceiverPoint = [i,j,k]
            do l = 1, ReceiverGrid%cart%nd
              ReceiverCoords(l) = ReceiverGrid%xyz(l)%values(i,j,k)
            end do
            ValidReceiverPoint = ReceiverGrid%grid_mask%values(i,j,k)
            if (ValidReceiverPoint .and. ovkBBContainsPoint(Bounds, ReceiverCoords)) then
              DonorCell = ovkFindDonorCell(DonorGrid, DonorAccel, ReceiverCoords)
              if (ovkCartContains(DonorGrid%cell_cart, DonorCell)) then
                Donors%valid_mask%values(i,j,k) = .true.
                Donors%grid_ids%values(i,j,k) = DonorGrid%id
                do l = 1, DonorGrid%cart%nd
                  Donors%cells(l)%values(i,j,k) = DonorCell(l)
                end do
                Donors%cell_extents%values(i,j,k) = 2
                DonorCellCoords = ovkCoordsInCell(DonorGrid, DonorCell, ReceiverCoords)
                do l = 1, DonorGrid%cart%nd
                  Donors%cell_coords(l)%values(i,j,k) = DonorCellCoords(l)
                end do
                DonorCellSize = ovkCellSize(DonorGrid, DonorCell)
                ReceiverCellSize = ovkAvgCellSizeAroundPoint(ReceiverGrid, ReceiverPoint)
                Donors%cell_diff_params%values(i,j,k) = log(DonorCellSize/ReceiverCellSize)/ &
                  (log(2._rk) * real(ReceiverGrid%cart%nd,kind=rk))
              end if
            end if
          end if
        end do
      end do
    end do

  end subroutine ovkFindDonors

  subroutine ovkChooseDonors(CandidateDonors, Subset)

    type(ovk_donors), dimension(:), intent(inout) :: CandidateDonors
    type(ovk_field_logical), intent(in), optional :: Subset

    integer :: i, j, k, m
    type(ovk_cart) :: Cart
    logical :: IncludePoint
    real(rk) :: CandidateCellDiff, BestCellDiff
    integer :: BestGrid

    if (size(CandidateDonors) == 0) return

    Cart = CandidateDonors(1)%cart

    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          if (present(Subset)) then
            IncludePoint = Subset%values(i,j,k)
          else
            IncludePoint = .true.
          end if
          if (IncludePoint) then
            BestCellDiff = huge(0._rk)
            BestGrid = 0
            do m = 1, size(CandidateDonors)
              if (CandidateDonors(m)%valid_mask%values(i,j,k)) then
                CandidateCellDiff = CandidateDonors(m)%cell_diff_params%values(i,j,k)
                if (CandidateCellDiff < BestCellDiff) then
                  if (BestGrid /= 0) then
                    CandidateDonors(BestGrid)%valid_mask%values(i,j,k) = .false.
                  end if
                  BestCellDiff = CandidateCellDiff
                  BestGrid = m
                else
                  CandidateDonors(m)%valid_mask%values(i,j,k) = .false.
                end if
              end if
            end do
          end if
        end do
      end do
    end do

  end subroutine ovkChooseDonors

  subroutine ovkMergeDonors(CandidateDonors, MergedDonors)

    type(ovk_donors), dimension(:), intent(in) :: CandidateDonors
    type(ovk_donors), intent(out) :: MergedDonors

    integer :: i, j, k, l, m
    type(ovk_cart) :: Cart

    if (size(CandidateDonors) == 0) return

    Cart = CandidateDonors(1)%cart

    call ovkMakeDonors(MergedDonors, Cart)
    MergedDonors%valid_mask%values = .false.

    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          do m = 1, size(CandidateDonors)
            if (CandidateDonors(m)%valid_mask%values(i,j,k)) then
              MergedDonors%valid_mask%values(i,j,k) = .true.
              MergedDonors%grid_ids%values(i,j,k) = CandidateDonors(m)%grid_ids%values(i,j,k)
              do l = 1, Cart%nd
                MergedDonors%cells(l)%values(i,j,k) = CandidateDonors(m)%cells(l)%values(i,j,k)
              end do
              MergedDonors%cell_extents%values(i,j,k) = &
                CandidateDonors(m)%cell_extents%values(i,j,k)
              do l = 1, Cart%nd
                MergedDonors%cell_coords(l)%values(i,j,k) = &
                  CandidateDonors(m)%cell_coords(l)%values(i,j,k)
              end do
              MergedDonors%cell_diff_params%values(i,j,k) = &
                CandidateDonors(m)%cell_diff_params%values(i,j,k)
              exit
            end if
          end do
        end do
      end do
    end do

  end subroutine ovkMergeDonors

  subroutine ovkPrintDonors(DonorGrid, ReceiverGrid, Donors)

    type(ovk_grid), intent(in) :: DonorGrid
    type(ovk_grid), intent(in) :: ReceiverGrid
    type(ovk_donors), intent(in) :: Donors

    integer :: i, j, k, l
    integer, dimension(DonorGrid%cart%nd) :: DonorCell
    real(rk), dimension(DonorGrid%cart%nd) :: DonorCellCoords
    integer, dimension(ReceiverGrid%cart%nd) :: ReceiverPoint

    write (*, '(a)') "Receiver DonorGrid DonorCell DonorCellCoords"
    do k = ReceiverGrid%cart%is(3), ReceiverGrid%cart%ie(3)
      do j = ReceiverGrid%cart%is(2), ReceiverGrid%cart%ie(2)
        do i = ReceiverGrid%cart%is(1), ReceiverGrid%cart%ie(1)
          ReceiverPoint = [i,j,k]
          if (Donors%valid_mask%values(i,j,k)) then
            do l = 1, DonorGrid%cart%nd
              DonorCell(l) = Donors%cells(l)%values(i,j,k)
              DonorCellCoords(l) = Donors%cell_coords(l)%values(i,j,k)
            end do
            write (*, '(4a)') trim(TupleToString(ReceiverPoint)), &
              trim(IntToString(Donors%grid_ids%values(i,j,k))), trim(TupleToString(DonorCell)), &
              trim(CoordsToString(DonorCellCoords))
          end if
        end do
      end do
    end do

  end subroutine ovkPrintDonors

  subroutine ovkGenerateReceiverMask(ReceiverGrid, DonorGrid, Donors, ReceiverMask, ReceiverSubset, &
    DonorSubset)

    type(ovk_grid), intent(in) :: ReceiverGrid
    type(ovk_grid), intent(in) :: DonorGrid
    type(ovk_donors), intent(in) :: Donors
    type(ovk_field_logical), intent(out) :: ReceiverMask
    type(ovk_field_logical), intent(in), optional :: ReceiverSubset
    type(ovk_field_logical), intent(in), optional :: DonorSubset

    integer :: i, j, k, l, m, n, o
    logical :: IncludePoint
    integer, dimension(MAX_ND) :: DonorCellLower
    integer, dimension(MAX_ND) :: DonorCellUpper
    logical :: AwayFromBoundary
    integer, dimension(MAX_ND) :: Vertex

    ReceiverMask = ovk_field_logical_(ReceiverGrid%cart, .false.)

    do k = ReceiverGrid%cart%is(3), ReceiverGrid%cart%ie(3)
      do j = ReceiverGrid%cart%is(2), ReceiverGrid%cart%ie(2)
        do i = ReceiverGrid%cart%is(1), ReceiverGrid%cart%ie(1)
          if (present(ReceiverSubset)) then
            IncludePoint = ReceiverSubset%values(i,j,k)
          else
            IncludePoint = .true.
          end if
          if (IncludePoint .and. Donors%valid_mask%values(i,j,k)) then
            if (Donors%grid_ids%values(i,j,k) /= DonorGrid%id) cycle
            if (present(DonorSubset)) then
              do l = 1, DonorGrid%cart%nd
                DonorCellLower(l) = Donors%cells(l)%values(i,j,k)
                DonorCellUpper(l) = Donors%cells(l)%values(i,j,k) + &
                  Donors%cell_extents%values(i,j,k)-1
              end do
              DonorCellLower(DonorGrid%cart%nd+1:) = 1
              DonorCellUpper(DonorGrid%cart%nd+1:) = 1
              AwayFromBoundary = ovkCartContains(DonorGrid%cart, DonorCellLower) .and. &
                ovkCartContains(DonorGrid%cart, DonorCellUpper)
              if (AwayFromBoundary) then
                do o = DonorCellLower(3), DonorCellUpper(3)
                  do n = DonorCellLower(2), DonorCellUpper(2)
                    do m = DonorCellLower(1), DonorCellUpper(1)
                      Vertex = [m,n,o]
                      if (DonorSubset%values(Vertex(1),Vertex(2),Vertex(3))) then
                        ReceiverMask%values(i,j,k) = .true.
                        exit
                      end if
                    end do
                  end do
                end do
              else
                do o = DonorCellLower(3), DonorCellUpper(3)
                  do n = DonorCellLower(2), DonorCellUpper(2)
                    do m = DonorCellLower(1), DonorCellUpper(1)
                      Vertex = [m,n,o]
                      Vertex(:DonorGrid%cart%nd) = ovkCartPeriodicAdjust(DonorGrid%cart, Vertex)
                      if (DonorSubset%values(Vertex(1),Vertex(2),Vertex(3))) then
                        ReceiverMask%values(i,j,k) = .true.
                        exit
                      end if
                    end do
                  end do
                end do
              end if
            else
              ReceiverMask%values(i,j,k) = .true.
            end if
          end if
        end do
      end do
    end do

  end subroutine ovkGenerateReceiverMask

  subroutine ovkGenerateDonorMask(DonorGrid, ReceiverGrid, Donors, DonorMask, DonorSubset, &
    ReceiverSubset)

    type(ovk_grid), intent(in) :: DonorGrid
    type(ovk_grid), intent(in) :: ReceiverGrid
    type(ovk_donors), intent(in) :: Donors
    type(ovk_field_logical), intent(out) :: DonorMask
    type(ovk_field_logical), intent(in), optional :: DonorSubset
    type(ovk_field_logical), intent(in), optional :: ReceiverSubset

    integer :: i, j, k, l, m, n, o
    logical :: IncludePoint
    integer, dimension(MAX_ND) :: DonorCellLower
    integer, dimension(MAX_ND) :: DonorCellUpper
    logical :: AwayFromBoundary
    integer, dimension(MAX_ND) :: Vertex

    DonorMask = ovk_field_logical_(DonorGrid%cart, .false.)

    do k = ReceiverGrid%cart%is(3), ReceiverGrid%cart%ie(3)
      do j = ReceiverGrid%cart%is(2), ReceiverGrid%cart%ie(2)
        do i = ReceiverGrid%cart%is(1), ReceiverGrid%cart%ie(1)
          if (present(ReceiverSubset)) then
            IncludePoint = ReceiverSubset%values(i,j,k)
          else
            IncludePoint = .true.
          end if
          if (IncludePoint .and. Donors%valid_mask%values(i,j,k)) then
            if (Donors%grid_ids%values(i,j,k) /= DonorGrid%id) cycle
            do l = 1, DonorGrid%cart%nd
              DonorCellLower(l) = Donors%cells(l)%values(i,j,k)
              DonorCellUpper(l) = Donors%cells(l)%values(i,j,k) + &
                Donors%cell_extents%values(i,j,k)-1
            end do
            DonorCellLower(DonorGrid%cart%nd+1:) = 1
            DonorCellUpper(DonorGrid%cart%nd+1:) = 1
            AwayFromBoundary = ovkCartContains(DonorGrid%cart, DonorCellLower) .and. &
              ovkCartContains(DonorGrid%cart, DonorCellUpper)
            if (AwayFromBoundary) then
              do o = DonorCellLower(3), DonorCellUpper(3)
                do n = DonorCellLower(2), DonorCellUpper(2)
                  do m = DonorCellLower(1), DonorCellUpper(1)
                    Vertex = [m,n,o]
                    if (present(DonorSubset)) then
                      if (DonorSubset%values(Vertex(1),Vertex(2),Vertex(3))) then
                        DonorMask%values(Vertex(1),Vertex(2),Vertex(3)) = .true.
                      end if
                    else
                      DonorMask%values(Vertex(1),Vertex(2),Vertex(3)) = .true.
                    end if
                  end do
                end do
              end do
            else
              do o = DonorCellLower(3), DonorCellUpper(3)
                do n = DonorCellLower(2), DonorCellUpper(2)
                  do m = DonorCellLower(1), DonorCellUpper(1)
                    Vertex = [m,n,o]
                    Vertex(:DonorGrid%cart%nd) = ovkCartPeriodicAdjust(DonorGrid%cart, Vertex)
                    if (present(DonorSubset)) then
                      if (DonorSubset%values(Vertex(1),Vertex(2),Vertex(3))) then
                        DonorMask%values(Vertex(1),Vertex(2),Vertex(3)) = .true.
                      end if
                    else
                      DonorMask%values(Vertex(1),Vertex(2),Vertex(3)) = .true.
                    end if
                  end do
                end do
              end do
            end if
          end if
        end do
      end do
    end do

  end subroutine ovkGenerateDonorMask

  subroutine ovkGenerateOverlapMask(CandidateDonors, OverlapMask)

    type(ovk_donors), dimension(:), intent(in) :: CandidateDonors
    type(ovk_field_logical), intent(out) :: OverlapMask

    integer :: i
    type(ovk_cart) :: Cart

    if (size(CandidateDonors) == 0) return

    Cart = CandidateDonors(1)%cart

    OverlapMask = ovk_field_logical_(Cart, .false.)

    do i = 1, size(CandidateDonors)
      OverlapMask%values = OverlapMask%values .or. CandidateDonors(i)%valid_mask%values
    end do

  end subroutine ovkGenerateOverlapMask

  subroutine ovkGenerateCoarseToFineMask(Grid, Donors, CoarseToFineMask, Subset)

    type(ovk_grid), intent(in) :: Grid
    type(ovk_donors), intent(in) :: Donors
    type(ovk_field_logical), intent(out) :: CoarseToFineMask
    type(ovk_field_logical), intent(in), optional :: Subset

    integer :: i, j, k
    logical :: IncludePoint

    CoarseToFineMask = ovk_field_logical_(Grid%cart, .false.)

    do k = Grid%cart%is(3), Grid%cart%ie(3)
      do j = Grid%cart%is(2), Grid%cart%ie(2)
        do i = Grid%cart%is(1), Grid%cart%ie(1)
          if (present(Subset)) then
            IncludePoint = Subset%values(i,j,k)
          else
            IncludePoint = .true.
          end if
          if (IncludePoint) then
            if (Donors%valid_mask%values(i,j,k)) then
              CoarseToFineMask%values(i,j,k) = Donors%cell_diff_params%values(i,j,k) > 0._rk &
                .or. (Donors%cell_diff_params%values(i,j,k) == 0._rk .and. &
                Donors%grid_ids%values(i,j,k) < Grid%id)
            end if
          end if
        end do
      end do
    end do

  end subroutine ovkGenerateCoarseToFineMask

  subroutine ovkGenerateCyclicReceiverMask(ReceiverGrid, DonorGrid, Donors, ReverseDonors, &
    CyclicReceiverMask, ReceiverSubset, DonorSubset)

    type(ovk_grid), intent(in) :: ReceiverGrid, DonorGrid
    type(ovk_donors), intent(in) :: Donors, ReverseDonors
    type(ovk_field_logical), intent(out) :: CyclicReceiverMask
    type(ovk_field_logical), intent(in), optional :: ReceiverSubset, DonorSubset

    type(ovk_field_logical) :: ReverseReceiverMask

    call ovkGenerateReceiverMask(DonorGrid, ReceiverGrid, ReverseDonors, ReverseReceiverMask, &
      ReceiverSubset=DonorSubset, DonorSubset=ReceiverSubset)

    call ovkGenerateReceiverMask(ReceiverGrid, DonorGrid, Donors, CyclicReceiverMask, &
      DonorSubset=ReverseReceiverMask, ReceiverSubset=ReceiverSubset)

  end subroutine ovkGenerateCyclicReceiverMask

  subroutine ovkGenerateNearCrossoverMask(Grid1, Grid2, Donors1, Donors2, Distance, &
    NearCrossoverMask, Subset1, Subset2)

    type(ovk_grid), intent(in) :: Grid1, Grid2
    type(ovk_donors), intent(in) :: Donors1, Donors2
    integer, intent(in) :: Distance
    type(ovk_field_logical), intent(out) :: NearCrossoverMask
    type(ovk_field_logical), intent(in), optional :: Subset1, Subset2

    type(ovk_field_logical) :: ReceiverMask
    type(ovk_field_logical) :: DonorMask

    call ovkGenerateReceiverMask(Grid1, Grid2, Donors1, ReceiverMask, ReceiverSubset=Subset1, &
      DonorSubset=Subset2)

    call ovkGenerateDonorMask(Grid1, Grid2, Donors2, DonorMask, DonorSubset=Subset1, &
      ReceiverSubset=Subset2)

    call ovkGrowMask(DonorMask, Distance)

    NearCrossoverMask = ovk_field_logical_(Grid1%cart)
    NearCrossoverMask%values = ReceiverMask%values .and. DonorMask%values

  end subroutine ovkGenerateNearCrossoverMask

  subroutine ovkGenerateOrphanMask(Donors, ReceiverMask, OrphanMask)

    type(ovk_donors), intent(in) :: Donors
    type(ovk_field_logical), intent(in) :: ReceiverMask
    type(ovk_field_logical), intent(out) :: OrphanMask

    OrphanMask = ovk_field_logical_(Donors%cart)
    OrphanMask%values = ReceiverMask%values .and. .not. Donors%valid_mask%values

  end subroutine ovkGenerateOrphanMask

end module ovkDonors
