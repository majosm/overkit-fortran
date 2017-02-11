! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkPegasus

  use ovkCart
  use ovkGlobal
  use ovkField
  use ovkGrid
  use ovkInterp
  implicit none

  private

  ! API
  public :: ovk_pegasus
  public :: ovk_pegasus_
  public :: ovkMakePegasusData
  public :: ovkDestroyPegasusData
  public :: ovkWritePegasusData
  public :: ovkPrintPegasusData

  ! PEGASUS variable name key
  ! =========================
  !
  ! Grid data
  ! ---------
  ! ngrd - number of grids
  ! ieng(ng) - number of points in 1-dir of grid ng
  ! jeng(ng) - number of points in 2-dir of grid ng
  ! keng(ng) - number of points in 3-dir of grid ng
  ! iblank(m,ng) - status of m-th grid point of grid ng (normal, cut, etc.)

  ! Global interp data
  ! ------------------
  ! ipall - total number of donor points over all grids
  ! igall - maximum number of grid points, taken over all grids
  ! ipip - maximum number of donor points, taken over all grids
  ! ipbp - maximum number of receiver points, taken over all grids
  ! iisptr(ng) - starting index of global linear index of donor points on grid ng
  ! iieptr(ng) - ending index of global linear index of donor points on grid ng

  ! Interp data
  ! -----------
  ! ibpnts(ng) - number of receiver points on grid ng
  ! iipnts(ng) - number of donor points on grid ng
  ! ibct(m,ng) - global linear (over all grids' donor points) index of the donor point for the m-th receiver point on grid ng
  ! iit(m,ng) - local 1-dir index of m-th donor point on grid ng (in stencil’s most SW corner)
  ! jit(m,ng) - local 2-dir index of m-th donor point on grid ng (in stencil’s most SW corner)
  ! kit(m,ng) - local 3-dir index of m-th donor point on grid ng (in stencil’s most SW corner)
  ! dxit(m,ng) - 1-dir computational offset of m-th donor point on grid ng [only used in OGEN]
  ! dyit(m,ng) - 2-dir computational offset of m-th donor point on grid ng [only used in OGEN]
  ! dzit(m,ng) - 3-dir computational offset of m-th donor point on grid ng [only used in OGEN]
  ! ibt(m,ng) - local 1-dir index of m-th receiver point on grid ng
  ! jbt(m,ng) - local 2-dir index of m-th receiver point on grid ng
  ! kbt(m,ng) - local 3-dir index of m-th receiver point on grid ng
  ! nit(m,ng) - 1-dir stencil size of m-th donor point on grid ng
  ! njt(m,ng) - 2-dir stencil size of m-th donor point on grid ng
  ! nkt(m,ng) - 3-dir stencil size of m-th donor point on grid ng
  ! coeffit(m,io,1,ng) - polynomial coefficients of m-th donor point on grid ng in 1-dir, io = 1, ..., nit(m,ng)
  ! coeffit(m,io,2,ng) - polynomial coefficients of m-th donor point on grid ng in 2-dir, io = 1, ..., njt(m,ng)
  ! coeffit(m,io,3,ng) - polynomial coefficients of m-th donor point on grid ng in 3-dir, io = 1, ..., nkt(m,ng)

  type ovk_pegasus
    integer :: ngrd
    integer, dimension(:), allocatable :: ieng, jeng, keng
    integer :: ipall
    integer :: igall
    integer :: ipip
    integer :: ipbp
    integer, dimension(:), allocatable :: iisptr
    integer, dimension(:), allocatable :: iieptr
    integer, dimension(:), allocatable :: ibpnts
    integer, dimension(:), allocatable :: iipnts
    integer, dimension(:,:), allocatable :: ibct
    integer, dimension(:,:), allocatable :: iit, jit, kit
    real(rk), dimension(:,:), allocatable :: dxit, dyit, dzit
    integer, dimension(:,:), allocatable :: ibt, jbt, kbt
    integer, dimension(:,:), allocatable :: nit, njt, nkt
    real(rk), dimension(:,:,:,:), allocatable :: coeffit
    integer, dimension(:,:), allocatable :: iblank
  end type ovk_pegasus

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_pegasus_
    module procedure ovk_pegasus_Default
  end interface ovk_pegasus_

contains

  pure function ovk_pegasus_Default() result(PegasusData)

    type(ovk_pegasus) :: PegasusData

    PegasusData%ngrd = 0
    PegasusData%ipall = 0
    PegasusData%igall = 0
    PegasusData%ipip = 0
    PegasusData%ipbp = 0

  end function ovk_pegasus_Default

  subroutine ovkMakePegasusData(Grids, InterpData, ExportCarts, PegasusData, IncludeIBlank, &
    HoleMasks)

    type(ovk_grid), dimension(:), intent(in) :: Grids
    type(ovk_cart), dimension(:), intent(in) :: ExportCarts
    type(ovk_interp), dimension(:), intent(in) :: InterpData
    type(ovk_pegasus), intent(out) :: PegasusData
    logical, intent(in), optional :: IncludeIBlank
    type(ovk_field_logical), dimension(:), intent(in), optional :: HoleMasks

    logical :: IncludeIBlank_
    integer :: d, i, j, k, m, n, o
    integer(lk) :: l
    integer :: nGrids
    integer :: nDims
    integer, dimension(MAX_ND) :: Point
    integer, dimension(:), allocatable :: nDonors, nReceivers
    integer :: nDonorsMax, nReceiversMax
    integer, dimension(:), allocatable :: nInterpStencil
    integer :: nInterpStencilMax
    integer :: Offset
    integer, dimension(:), allocatable :: NextReceiver
    integer, dimension(:), allocatable :: NextDonor
    integer, dimension(MAX_ND) :: DonorCell
    real(rk), dimension(MAX_ND) :: DonorCellCoords
    integer :: InterpScheme
    real(rk), dimension(:,:), allocatable :: InterpCoefs

    if (size(Grids) == 0) then
      PegasusData = ovk_pegasus_()
      return
    end if

    nDims = Grids(1)%cart%nd
    nGrids = size(Grids)

    if (present(IncludeIBlank)) then
      IncludeIBlank_ = IncludeIBlank
    else
      IncludeIBlank_ = .false.
    end if

    if (IncludeIBlank_ .and. .not. present(HoleMasks)) then
      write (ERROR_UNIT, '(2a)') "ERROR: Cannot construct IBlank in Pegasus data structure ", &
        "without hole masks."
      stop 1
    end if

    allocate(nDonors(nGrids))
    allocate(nReceivers(nGrids))

    nDonors = 0
    nReceivers = 0

    do m = 1, nGrids
      do k = ExportCarts(m)%is(3), ExportCarts(m)%ie(3)
        do j = ExportCarts(m)%is(2), ExportCarts(m)%ie(2)
          do i = ExportCarts(m)%is(1), ExportCarts(m)%ie(1)
            Point = [i,j,k]
            Point(:nDims) = ovkCartPeriodicAdjust(Grids(m)%cart, Point)
            if (InterpData(m)%valid_mask%values(Point(1),Point(2),Point(3))) then
              n = InterpData(m)%donor_grid_ids%values(Point(1),Point(2),Point(3))
              nReceivers(m) = nReceivers(m) + 1
              nDonors(n) = nDonors(n) + 1
            end if
          end do
        end do
      end do
    end do

    nDonorsMax = maxval(nDonors)
    nReceiversMax = maxval(nReceivers)

    allocate(nInterpStencil(nGrids))

    nInterpStencil = [(size(InterpData(m)%coefs,1),m=1,nGrids)]
    nInterpStencilMax = maxval(nInterpStencil)

    allocate(InterpCoefs(nInterpStencilMax,MAX_ND))

    allocate(PegasusData%ieng(nGrids))
    allocate(PegasusData%jeng(nGrids))
    allocate(PegasusData%keng(nGrids))
    allocate(PegasusData%iisptr(nGrids))
    allocate(PegasusData%iieptr(nGrids))
    allocate(PegasusData%ibpnts(nGrids))
    allocate(PegasusData%iipnts(nGrids))
    allocate(PegasusData%ibct(nReceiversMax,nGrids))
    allocate(PegasusData%iit(nDonorsMax,nGrids))
    allocate(PegasusData%jit(nDonorsMax,nGrids))
    allocate(PegasusData%kit(nDonorsMax,nGrids))
    allocate(PegasusData%dxit(nDonorsMax,nGrids))
    allocate(PegasusData%dyit(nDonorsMax,nGrids))
    allocate(PegasusData%dzit(nDonorsMax,nGrids))
    allocate(PegasusData%ibt(nReceiversMax,nGrids))
    allocate(PegasusData%jbt(nReceiversMax,nGrids))
    allocate(PegasusData%kbt(nReceiversMax,nGrids))
    allocate(PegasusData%nit(nDonorsMax,nGrids))
    allocate(PegasusData%njt(nDonorsMax,nGrids))
    allocate(PegasusData%nkt(nDonorsMax,nGrids))
    allocate(PegasusData%coeffit(nDonorsMax,nInterpStencilMax,MAX_ND,nGrids))

    PegasusData%ngrd = nGrids

    do m = 1, nGrids
      PegasusData%ieng(m) = ExportCarts(m)%ie(1)-ExportCarts(m)%is(1)+1
      PegasusData%jeng(m) = ExportCarts(m)%ie(2)-ExportCarts(m)%is(2)+1
      PegasusData%keng(m) = ExportCarts(m)%ie(3)-ExportCarts(m)%is(3)+1
    end do

    PegasusData%ipall = sum([(nDonors(m),m=1,nGrids)])
    PegasusData%igall = int(maxval([(ovkCartCount(ExportCarts(m)),m=1,nGrids)]))

    PegasusData%ipip = nDonorsMax
    PegasusData%ipbp = nReceiversMax

    Offset = 0
    do m = 1, nGrids
      PegasusData%iisptr(m) = Offset + 1
      PegasusData%iieptr(m) = Offset + nDonors(m)
      Offset = Offset + nDonors(m)
    end do

    PegasusData%ibpnts = nReceivers
    PegasusData%iipnts = nDonors

    allocate(NextReceiver(nGrids))
    allocate(NextDonor(nGrids))
    NextReceiver = 1
    NextDonor = 1

    do m = 1, nGrids
      do k = ExportCarts(m)%is(3), ExportCarts(m)%ie(3)
        do j = ExportCarts(m)%is(2), ExportCarts(m)%ie(2)
          do i = ExportCarts(m)%is(1), ExportCarts(m)%ie(1)
            Point = [i,j,k]
            Point(:nDims) = ovkCartPeriodicAdjust(Grids(m)%cart, Point)
            if (InterpData(m)%valid_mask%values(Point(1),Point(2),Point(3))) then
              n = InterpData(m)%donor_grid_ids%values(Point(1),Point(2),Point(3))
              do d = 1, nDims
                DonorCell(d) = InterpData(m)%donor_cells(d)%values(Point(1),Point(2),Point(3))
                DonorCellCoords(d) = InterpData(m)%donor_cell_coords(d)%values(Point(1),Point(2), &
                  Point(3))
              end do
              DonorCell(nDims+1:) = 1
              DonorCellCoords(nDims+1:) = 0._rk
              InterpScheme = InterpData(m)%schemes%values(Point(1),Point(2),Point(3))
              do d = 1, nDims
                do o = 1, nInterpStencil(m)
                  InterpCoefs(o,d) = InterpData(m)%coefs(o,d)%values(Point(1),Point(2),Point(3))
                end do
              end do
              InterpCoefs(nInterpStencil(m)+1:,:nDims) = 0._rk
              InterpCoefs(1,nDims+1:) = 1._rk
              InterpCoefs(2:,nDims+1:) = 0._rk
              PegasusData%ibct(NextReceiver(m),m) = PegasusData%iisptr(n) + NextDonor(n) - 1
              PegasusData%iit(NextDonor(n),n) = DonorCell(1)
              PegasusData%jit(NextDonor(n),n) = DonorCell(2)
              PegasusData%kit(NextDonor(n),n) = DonorCell(3)
              PegasusData%dxit(NextDonor(n),n) = DonorCellCoords(1)
              PegasusData%dyit(NextDonor(n),n) = DonorCellCoords(2)
              PegasusData%dzit(NextDonor(n),n) = DonorCellCoords(3)
              PegasusData%ibt(NextReceiver(m),m) = i
              PegasusData%jbt(NextReceiver(m),m) = j
              PegasusData%kbt(NextReceiver(m),m) = k
              if (InterpScheme == OVK_INTERP_LINEAR) then
                PegasusData%nit(NextDonor(n),n) = 2
                PegasusData%njt(NextDonor(n),n) = 2
                PegasusData%nkt(NextDonor(n),n) = merge(2, 1, nDims == MAX_ND)
              else if (InterpScheme == OVK_INTERP_CUBIC) then
                PegasusData%nit(NextDonor(n),n) = 4
                PegasusData%njt(NextDonor(n),n) = 4
                PegasusData%nkt(NextDonor(n),n) = merge(4, 1, nDims == MAX_ND)
              end if
              PegasusData%coeffit(NextDonor(n),:,:,n) = InterpCoefs
              NextReceiver(m) = NextReceiver(m) + 1
              NextDonor(n) = NextDonor(n) + 1
            end if
          end do
        end do
      end do
    end do

    if (IncludeIBlank_) then
      allocate(PegasusData%iblank(maxval([(ovkCartCount(ExportCarts(m)),m=1,nGrids)]),nGrids))
      PegasusData%iblank = 1
      do m = 1, nGrids
        l = 1_lk
        do k = ExportCarts(m)%is(3), ExportCarts(m)%ie(3)
          do j = ExportCarts(m)%is(2), ExportCarts(m)%ie(2)
            do i = ExportCarts(m)%is(1), ExportCarts(m)%ie(1)
              Point = [i,j,k]
              Point(:nDims) = ovkCartPeriodicAdjust(Grids(m)%cart, Point)
              if (.not. HoleMasks(m)%values(Point(1),Point(2),Point(3)) .and. .not. &
                InterpData(m)%valid_mask%values(Point(1),Point(2),Point(3))) then
                PegasusData%iblank(l,m) = 1
              else
                PegasusData%iblank(l,m) = 0
              end if
              l = l + 1_lk
            end do
          end do
        end do
      end do
    end if

  end subroutine ovkMakePegasusData

  subroutine ovkDestroyPegasusData(PegasusData)

    type(ovk_pegasus), intent(inout) :: PegasusData

    deallocate(PegasusData%ieng)
    deallocate(PegasusData%jeng)
    deallocate(PegasusData%keng)
    deallocate(PegasusData%iisptr)
    deallocate(PegasusData%iieptr)
    deallocate(PegasusData%ibpnts)
    deallocate(PegasusData%iipnts)
    deallocate(PegasusData%ibct)
    deallocate(PegasusData%iit)
    deallocate(PegasusData%jit)
    deallocate(PegasusData%kit)
    deallocate(PegasusData%dxit)
    deallocate(PegasusData%dyit)
    deallocate(PegasusData%dzit)
    deallocate(PegasusData%ibt)
    deallocate(PegasusData%jbt)
    deallocate(PegasusData%kbt)
    deallocate(PegasusData%nit)
    deallocate(PegasusData%njt)
    deallocate(PegasusData%nkt)
    deallocate(PegasusData%coeffit)
    if (allocated(PegasusData%iblank)) then
      deallocate(PegasusData%iblank)
    end if

  end subroutine ovkDestroyPegasusData

  subroutine ovkWritePegasusData(PegasusData, HOFilePath, XFilePath, Error)

    type(ovk_pegasus), intent(in) :: PegasusData
    character(len=*), intent(in) :: HOFilePath
    character(len=*), intent(in) :: XFilePath
    integer, intent(out), optional :: Error

    integer :: i, j, p
    integer(lk) :: l
    integer :: Error_
    integer(lk), dimension(3) :: nPoints
    integer(lk) :: nPointsTotal

    integer, parameter :: HOUnit = 4510695
    integer, parameter :: XUnit = 4510696

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Writing file ", trim(HOFilePath), "..."
    end if

    open(HOUnit, file=HOFilePath, status='unknown', form='unformatted', iostat=Error_)
    if (Error_ /= 0) then
      write (ERROR_UNIT, '(3a)') "ERROR: Unable to open file ", trim(HOFilePath), "."
      goto 999
    end if

    write (HOUnit) PegasusData%ngrd, PegasusData%ipall, PegasusData%igall, PegasusData%ipip, &
      PegasusData%ipbp
    do i = 1, PegasusData%ngrd
      write (HOUnit) PegasusData%ibpnts(i), PegasusData%iipnts(i), PegasusData%iieptr(i), &
        PegasusData%iisptr(i), PegasusData%ieng(i), PegasusData%jeng(i), PegasusData%keng(i)
      write (HOUnit) &
        (PegasusData%iit(j,i),j=1,PegasusData%iipnts(i)), &
        (PegasusData%jit(j,i),j=1,PegasusData%iipnts(i)), &
        (PegasusData%kit(j,i),j=1,PegasusData%iipnts(i)), &
        (PegasusData%dxit(j,i),j=1,PegasusData%iipnts(i)), &
        (PegasusData%dyit(j,i),j=1,PegasusData%iipnts(i)), &
        (PegasusData%dzit(j,i),j=1,PegasusData%iipnts(i))
      write (HOUnit) &
        (PegasusData%ibt(j,i),j=1,PegasusData%ibpnts(i)), &
        (PegasusData%jbt(j,i),j=1,PegasusData%ibpnts(i)), &
        (PegasusData%kbt(j,i),j=1,PegasusData%ibpnts(i)), &
        (PegasusData%ibct(j,i),j=1,PegasusData%ibpnts(i))
      if (allocated(PegasusData%iblank)) then
        nPoints(1) = int(PegasusData%ieng(i),kind=lk)
        nPoints(2) = int(PegasusData%jeng(i),kind=lk)
        nPoints(3) = int(PegasusData%keng(i),kind=lk)
        nPointsTotal = product(nPoints)
        write (HOUnit) (PegasusData%iblank(l,i),l=1_lk,nPointsTotal)
      end if
    end do
    close(HOUnit)

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Done writing file ", trim(HOFilePath), "."
    end if

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Writing file ", trim(XFilePath), "..."
    end if

    open(XUnit, file=XFilePath, status='unknown', form='unformatted', iostat=Error_)
    if (Error_ /= 0) then
      write (ERROR_UNIT, '(3a)') "ERROR: Unable to open file ", trim(XFilePath), "."
      goto 999
    end if

    write (XUnit) PegasusData%ngrd, PegasusData%ipall, PegasusData%igall, PegasusData%ipip, &
      PegasusData%ipbp
    do i = 1, PegasusData%ngrd
      write (XUnit) &
        (PegasusData%nit(j,i),j=1,PegasusData%iipnts(i)), &
        (PegasusData%njt(j,i),j=1,PegasusData%iipnts(i)), &
        (PegasusData%nkt(j,i),j=1,PegasusData%iipnts(i))
      write (XUnit) &
        ((PegasusData%coeffit(j,p,1,i),p=1,PegasusData%nit(j,i)),j=1,PegasusData%iipnts(i)), &
        ((PegasusData%coeffit(j,p,2,i),p=1,PegasusData%njt(j,i)),j=1,PegasusData%iipnts(i)), &
        ((PegasusData%coeffit(j,p,3,i),p=1,PegasusData%nkt(j,i)),j=1,PegasusData%iipnts(i))
      write (XUnit) 0, 0, 0
    end do
    close(XUnit)

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Done writing file ", trim(XFilePath), "."
    end if

999 if (Error_ /= 0) then
      if (present(Error)) then
        Error = OVK_IO_ERROR
        return
      else
        stop OVK_IO_ERROR
      end if
    end if

  end subroutine ovkWritePegasusData

  subroutine ovkPrintPegasusData(PegasusData)

    type(ovk_pegasus), intent(in) :: PegasusData

    integer :: i, j, p
    integer(lk) :: l
    integer(lk), dimension(3) :: nPoints
    integer(lk) :: nPointsTotal

    write (*, '(a,i0)') "ngrd = ", PegasusData%ngrd
    write (*, '(a,i0)') "ipall = ", PegasusData%ipall
    write (*, '(a,i0)') "igall = ", PegasusData%igall
    write (*, '(a,i0)') "ipip = ", PegasusData%ipip
    write (*, '(a,i0)') "ipbp = ", PegasusData%ipbp
    do i = 1, PegasusData%ngrd
      write (*, '(a)') ""
      write (*, '(a,i0,a)') "====== Grid ", i, " ======"
      write (*, '(a)') ""
      write (*, '(a,i0)') "ibpnts = ", PegasusData%ibpnts(i)
      write (*, '(a,i0)') "iipnts = ", PegasusData%iipnts(i)
      write (*, '(a,i0)') "iisptr = ", PegasusData%iisptr(i)
      write (*, '(a,i0)') "iieptr = ", PegasusData%iieptr(i)
      write (*, '(a,i0)') "ieng = ", PegasusData%ieng(i)
      write (*, '(a,i0)') "jeng = ", PegasusData%jeng(i)
      write (*, '(a,i0)') "keng = ", PegasusData%keng(i)
      write (*, '(a)') ""
      write (*, '(a)') "   iit       jit       kit       dxit       dyit       dzit   "
      write (*, '(a)') "--------- --------- --------- ---------- ---------- ----------"
      do j = 1, PegasusData%iipnts(i)
        write (*, '(3(i9,a),2(f10.3,a),f10.3)') PegasusData%iit(j,i), " ", PegasusData%jit(j,i), &
          " ", PegasusData%kit(j,i), " ", PegasusData%dxit(j,i), " ", PegasusData%dyit(j,i), &
          " ", PegasusData%dzit(j,i)
      end do
      write (*, '(a)') ""
      write (*, '(a)') "   ibt       jbt       kbt       ibct   "
      write (*, '(a)') "--------- --------- --------- ----------"
      do j = 1, PegasusData%ibpnts(i)
        write (*, '(3(i9,a),i10)') PegasusData%ibt(j,i), " ", PegasusData%jbt(j,i), " ", &
          PegasusData%kbt(j,i), " ", PegasusData%ibct(j,i)
      end do
      write (*, '(a)') ""
      write (*, '(a)') "   nit       njt       nkt   "
      write (*, '(a)') "--------- --------- ---------"
      do j = 1, PegasusData%iipnts(i)
        write (*, '(2(i9,a),i9)') PegasusData%nit(j,i), " ", PegasusData%njt(j,i), " ", &
          PegasusData%nkt(j,i)
      end do
      write (*, '(a)') ""
      write (*, '(a)') "coeffit = "
      do j = 1, PegasusData%iipnts(i)
        write (*, '(48f9.3)') &
          (PegasusData%coeffit(j,p,1,i),p=1,PegasusData%nit(j,i)), &
          (PegasusData%coeffit(j,p,2,i),p=1,PegasusData%njt(j,i)), &
          (PegasusData%coeffit(j,p,3,i),p=1,PegasusData%nkt(j,i))
      end do
      if (allocated(PegasusData%iblank)) then
        write (*, '(a)') ""
        write (*, '(a)') "iblank = "
        nPoints(1) = int(PegasusData%ieng(i),kind=lk)
        nPoints(2) = int(PegasusData%jeng(i),kind=lk)
        nPoints(3) = int(PegasusData%keng(i),kind=lk)
        nPointsTotal = product(nPoints)
        do l = 1_lk, nPointsTotal
          write (*, '(i9)') PegasusData%iblank(l,i)
        end do
      end if
    end do

  end subroutine ovkPrintPegasusData

end module ovkPegasus
