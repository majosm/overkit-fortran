! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkHashGrid

  use ovkBoundingBox
  use ovkCart
  use ovkField
  use ovkGeometry
  use ovkGlobal
  implicit none

  private

  ! Internal
  public :: t_hash_grid
  public :: t_hash_grid_
  public :: HashGridBin
  public :: HashGridBinBounds
  public :: HashGridStats
  public :: HashGridHistogram

  type t_hash_grid
    type(t_noconstruct) :: noconstruct
    type(ovk_cart) :: cart
    type(ovk_bbox) :: bounds
    real(rk), dimension(MAX_ND) :: bin_size
    integer(lk), dimension(:), allocatable :: bin_start
    integer(lk), dimension(:), allocatable :: bin_contents
  end type t_hash_grid

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_hash_grid_
    module procedure t_hash_grid_Default
    module procedure t_hash_grid_Empty
    module procedure t_hash_grid_Allocated
  end interface t_hash_grid_

contains

  pure function t_hash_grid_Default() result(HashGrid)

    type(t_hash_grid) :: HashGrid

    HashGrid = t_hash_grid_Empty(2)

  end function t_hash_grid_Default

  pure function t_hash_grid_Empty(NumDims) result(HashGrid)

    integer, intent(in) :: NumDims
    type(t_hash_grid) :: HashGrid

    HashGrid%cart = ovk_cart_(NumDims)
    HashGrid%bounds = ovk_bbox_(NumDims)
    HashGrid%bin_size = 0._rk

  end function t_hash_grid_Empty

  pure function t_hash_grid_Allocated(Cart, Bounds, NumBinEntries) result(HashGrid)

    type(ovk_cart), intent(in) :: Cart
    type(ovk_bbox), intent(in) :: Bounds
    type(ovk_field_int), intent(in) :: NumBinEntries
    type(t_hash_grid) :: HashGrid

    integer :: i, j, k
    integer(lk) :: l
    integer(lk) :: BinStart

    HashGrid%cart = Cart
    HashGrid%bounds = Bounds
    HashGrid%bin_size(:Cart%nd) = ovkBBSize(Bounds)/real(ovkCartSize(Cart),kind=rk)
    HashGrid%bin_size(Cart%nd+1:) = 0._rk

    allocate(HashGrid%bin_start(ovkCartCount(Cart)+1_lk))

    BinStart = 1_lk
    l = 1_lk
    do k = Cart%is(3), Cart%ie(3)
      do j = Cart%is(2), Cart%ie(2)
        do i = Cart%is(1), Cart%ie(1)
          HashGrid%bin_start(l) = BinStart
          BinStart = BinStart + int(NumBinEntries%values(i,j,k),kind=lk)
          l = l + 1_lk
        end do
      end do
    end do
    HashGrid%bin_start(l) = BinStart

    allocate(HashGrid%bin_contents(BinStart-1_lk))

  end function t_hash_grid_Allocated

  pure function HashGridBin(HashGrid, Coords) result(Bin)

    type(t_hash_grid), intent(in) :: HashGrid
    real(rk), dimension(HashGrid%cart%nd), intent(in) :: Coords
    integer, dimension(HashGrid%cart%nd) :: Bin

    if (ovkBBContainsPoint(HashGrid%bounds, Coords)) then
      Bin = int(ovkCartesianGridCell(HashGrid%bounds%b(:HashGrid%cart%nd), HashGrid%bin_size(: &
        HashGrid%cart%nd), Coords))
      Bin = ovkCartClamp(HashGrid%cart, Bin)
    else
      Bin = HashGrid%cart%is(:HashGrid%cart%nd)-1
    end if

  end function HashGridBin

  pure function HashGridBinBounds(HashGrid, Bin) result(BinBounds)

    type(t_hash_grid), intent(in) :: HashGrid
    integer, dimension(HashGrid%cart%nd), intent(in) :: Bin
    type(ovk_bbox) :: BinBounds

    real(rk), dimension(HashGrid%cart%nd) :: BinBegin, BinEnd

    BinBegin = HashGrid%bounds%b(:HashGrid%cart%nd) + HashGrid%bin_size(:HashGrid%cart%nd) * &
      real(Bin-1,kind=rk)
    BinEnd = HashGrid%bounds%b(:HashGrid%cart%nd) + HashGrid%bin_size(:HashGrid%cart%nd) * &
      real(Bin,kind=rk)

    BinBounds = ovk_bbox_(HashGrid%cart%nd, BinBegin, BinEnd)

  end function HashGridBinBounds

  subroutine HashGridStats(HashGrid, NumBins, NumNonEmptyBins, MinBinEntries, MaxBinEntries, &
    TotalBinEntries)

    type(t_hash_grid), intent(in) :: HashGrid
    integer(lk), intent(out) :: NumBins
    integer(lk), intent(out) :: NumNonEmptyBins
    integer, intent(out) :: MinBinEntries
    integer, intent(out) :: MaxBinEntries
    integer(lk), intent(out) :: TotalBinEntries

    integer :: i, j, k
    integer(lk) :: l
    integer :: NumBinEntries

    NumBins = ovkCartCount(HashGrid%cart)

    NumNonEmptyBins = 0_lk
    MinBinEntries = huge(0)
    MaxBinEntries = 0
    TotalBinEntries = 0_lk
    l = 1_lk
    do k = HashGrid%cart%is(3), HashGrid%cart%ie(3)
      do j = HashGrid%cart%is(2), HashGrid%cart%ie(2)
        do i = HashGrid%cart%is(1), HashGrid%cart%ie(1)
          NumBinEntries = int(HashGrid%bin_start(l+1_lk)-HashGrid%bin_start(l))
          NumNonEmptyBins = NumNonEmptyBins + merge(1_lk, 0_lk, NumBinEntries > 0_lk)
          MinBinEntries = min(MinBinEntries, NumBinEntries)
          MaxBinEntries = max(MaxBinEntries, NumBinEntries)
          TotalBinEntries = TotalBinEntries + NumBinEntries
          l = l + 1_lk
        end do
      end do
    end do

  end subroutine HashGridStats

  subroutine HashGridHistogram(HashGrid, Lower, Upper, N, Histogram)

    type(t_hash_grid), intent(in) :: HashGrid
    integer, intent(in) :: Lower
    integer, intent(in) :: Upper
    integer, intent(in) :: N
    integer(lk), dimension(N), intent(out) :: Histogram

    integer :: i, j, k, m
    integer(lk) :: l
    real(rk) :: HistogramStart
    real(rk) :: HistogramInterval
    integer(lk) :: NumBinEntries

    Histogram = 0_lk

    if (Upper-Lower > 0) then

      HistogramStart = real(Lower,kind=rk)-0.5_rk
      HistogramInterval = (real(Upper-Lower,kind=rk)+1._rk)/real(N,kind=rk)

      l = 1_lk
      do k = HashGrid%cart%is(3), HashGrid%cart%ie(3)
        do j = HashGrid%cart%is(2), HashGrid%cart%ie(2)
          do i = HashGrid%cart%is(1), HashGrid%cart%ie(1)
            NumBinEntries = HashGrid%bin_start(l+1_lk)-HashGrid%bin_start(l)
            m = int(ovkCartesianGridCell(HistogramStart, HistogramInterval, &
              real(NumBinEntries,kind=rk)))
            m = min(max(m, 1), N)
            Histogram(m) = Histogram(m) + 1_lk
            l = l + 1_lk
          end do
        end do
      end do

    end if

  end subroutine HashGridHistogram

end module ovkHashGrid
