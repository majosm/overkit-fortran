! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkPegasus

  use ovkCart
  use ovkConnectivity
  use ovkGlobal
  use ovkField
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

#define IO_ERROR 1

contains

  pure function ovk_pegasus_Default() result(PegasusData)

    type(ovk_pegasus) :: PegasusData

    PegasusData%ngrd = 0
    PegasusData%ipall = 0
    PegasusData%igall = 0
    PegasusData%ipip = 0
    PegasusData%ipbp = 0

  end function ovk_pegasus_Default

  subroutine ovkMakePegasusData(PegasusData, Connectivity, IncludeIBlank, HoleMasks)

    type(ovk_pegasus), intent(out) :: PegasusData
    type(ovk_connectivity), intent(in) :: Connectivity
    logical, intent(in), optional :: IncludeIBlank
    type(ovk_field_logical), dimension(:), intent(in), optional :: HoleMasks

    logical :: IncludeIBlank_
    integer :: d, i, j, k, m, n, o
    integer(lk) :: l
    type(ovk_interp), dimension(:), pointer :: InterpData
    integer :: NumGrids
    integer :: NumDims
    integer, dimension(MAX_ND) :: Point
    integer, dimension(:), allocatable :: NumDonors, NumReceivers
    integer :: MaxDonors, MaxReceivers
    integer :: MaxStencilSize
    integer :: Offset
    integer, dimension(:), allocatable :: NextReceiver
    integer, dimension(:), allocatable :: NextDonor
    integer, dimension(MAX_ND) :: DonorCell
    real(rk), dimension(MAX_ND) :: DonorCellCoords
    integer :: InterpScheme
    real(rk), dimension(:,:), allocatable :: InterpCoefs

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

    NumDims = Connectivity%properties%nd
    NumGrids = Connectivity%properties%ngrids

    if (NumGrids == 0) then
      PegasusData = ovk_pegasus_()
      return
    end if

    InterpData => Connectivity%interp_data

    allocate(NumDonors(NumGrids))
    allocate(NumReceivers(NumGrids))

    NumDonors = 0
    NumReceivers = 0

    do m = 1, NumGrids
      do k = 1, InterpData(m)%properties%npoints(3)
        do j = 1, InterpData(m)%properties%npoints(2)
          do i = 1, InterpData(m)%properties%npoints(1)
            Point = [i,j,k]
            Point(:NumDims) = ovkCartPeriodicAdjust(InterpData(m)%cart, Point)
            if (InterpData(m)%valid_mask%values(Point(1),Point(2),Point(3))) then
              n = InterpData(m)%donor_grid_ids%values(Point(1),Point(2),Point(3))
              NumReceivers(m) = NumReceivers(m) + 1
              NumDonors(n) = NumDonors(n) + 1
            end if
          end do
        end do
      end do
    end do

    MaxDonors = maxval(NumDonors)
    MaxReceivers = maxval(NumReceivers)

    MaxStencilSize = 0
    do m = 1, NumGrids
      MaxStencilSize = max(MaxStencilSize, InterpData(m)%properties%stencil_size)
    end do

    allocate(InterpCoefs(MaxStencilSize,MAX_ND))

    allocate(PegasusData%ieng(NumGrids))
    allocate(PegasusData%jeng(NumGrids))
    allocate(PegasusData%keng(NumGrids))
    allocate(PegasusData%iisptr(NumGrids))
    allocate(PegasusData%iieptr(NumGrids))
    allocate(PegasusData%ibpnts(NumGrids))
    allocate(PegasusData%iipnts(NumGrids))
    allocate(PegasusData%ibct(MaxReceivers,NumGrids))
    allocate(PegasusData%iit(MaxDonors,NumGrids))
    allocate(PegasusData%jit(MaxDonors,NumGrids))
    allocate(PegasusData%kit(MaxDonors,NumGrids))
    allocate(PegasusData%dxit(MaxDonors,NumGrids))
    allocate(PegasusData%dyit(MaxDonors,NumGrids))
    allocate(PegasusData%dzit(MaxDonors,NumGrids))
    allocate(PegasusData%ibt(MaxReceivers,NumGrids))
    allocate(PegasusData%jbt(MaxReceivers,NumGrids))
    allocate(PegasusData%kbt(MaxReceivers,NumGrids))
    allocate(PegasusData%nit(MaxDonors,NumGrids))
    allocate(PegasusData%njt(MaxDonors,NumGrids))
    allocate(PegasusData%nkt(MaxDonors,NumGrids))
    allocate(PegasusData%coeffit(MaxDonors,MaxStencilSize,MAX_ND,NumGrids))

    PegasusData%ngrd = NumGrids

    do m = 1, NumGrids
      PegasusData%ieng(m) = InterpData(m)%properties%npoints(1)
      PegasusData%jeng(m) = InterpData(m)%properties%npoints(2)
      PegasusData%keng(m) = InterpData(m)%properties%npoints(3)
    end do

    PegasusData%ipall = sum([(NumDonors(m),m=1,NumGrids)])
    PegasusData%igall = maxval([(product(InterpData(m)%properties%npoints),m=1,NumGrids)])

    PegasusData%ipip = MaxDonors
    PegasusData%ipbp = MaxReceivers

    Offset = 0
    do m = 1, NumGrids
      PegasusData%iisptr(m) = Offset + 1
      PegasusData%iieptr(m) = Offset + NumDonors(m)
      Offset = Offset + NumDonors(m)
    end do

    PegasusData%ibpnts = NumReceivers
    PegasusData%iipnts = NumDonors

    allocate(NextReceiver(NumGrids))
    allocate(NextDonor(NumGrids))
    NextReceiver = 1
    NextDonor = 1

    do m = 1, NumGrids
      do k = 1, InterpData(m)%properties%npoints(3)
        do j = 1, InterpData(m)%properties%npoints(2)
          do i = 1, InterpData(m)%properties%npoints(1)
            Point = [i,j,k]
            Point(:NumDims) = ovkCartPeriodicAdjust(InterpData(m)%cart, Point)
            if (InterpData(m)%valid_mask%values(Point(1),Point(2),Point(3))) then
              n = InterpData(m)%donor_grid_ids%values(Point(1),Point(2),Point(3))
              do d = 1, NumDims
                DonorCell(d) = InterpData(m)%donor_cells(d)%values(Point(1),Point(2),Point(3))
                DonorCellCoords(d) = InterpData(m)%donor_cell_coords(d)%values(Point(1),Point(2), &
                  Point(3))
              end do
              DonorCell(NumDims+1:) = 1
              DonorCellCoords(NumDims+1:) = 0._rk
              InterpScheme = InterpData(m)%schemes%values(Point(1),Point(2),Point(3))
              do d = 1, NumDims
                do o = 1, InterpData(m)%properties%stencil_size
                  InterpCoefs(o,d) = InterpData(m)%coefs(o,d)%values(Point(1),Point(2),Point(3))
                end do
              end do
              InterpCoefs(InterpData(m)%properties%stencil_size+1:,:NumDims) = 0._rk
              InterpCoefs(1,NumDims+1:) = 1._rk
              InterpCoefs(2:,NumDims+1:) = 0._rk
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
                PegasusData%nkt(NextDonor(n),n) = merge(2, 1, NumDims == MAX_ND)
              else if (InterpScheme == OVK_INTERP_CUBIC) then
                PegasusData%nit(NextDonor(n),n) = 4
                PegasusData%njt(NextDonor(n),n) = 4
                PegasusData%nkt(NextDonor(n),n) = merge(4, 1, NumDims == MAX_ND)
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
      allocate(PegasusData%iblank(PegasusData%igall,NumGrids))
      PegasusData%iblank = 1
      do m = 1, NumGrids
        l = 1_lk
        do k = 1, InterpData(m)%properties%npoints(3)
          do j = 1, InterpData(m)%properties%npoints(2)
            do i = 1, InterpData(m)%properties%npoints(1)
              Point = [i,j,k]
              Point(:NumDims) = ovkCartPeriodicAdjust(InterpData(m)%cart, Point)
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

    if (allocated(PegasusData%ieng)) deallocate(PegasusData%ieng)
    if (allocated(PegasusData%jeng)) deallocate(PegasusData%jeng)
    if (allocated(PegasusData%keng)) deallocate(PegasusData%keng)
    if (allocated(PegasusData%iisptr)) deallocate(PegasusData%iisptr)
    if (allocated(PegasusData%iieptr)) deallocate(PegasusData%iieptr)
    if (allocated(PegasusData%ibpnts)) deallocate(PegasusData%ibpnts)
    if (allocated(PegasusData%iipnts)) deallocate(PegasusData%iipnts)
    if (allocated(PegasusData%ibct)) deallocate(PegasusData%ibct)
    if (allocated(PegasusData%iit)) deallocate(PegasusData%iit)
    if (allocated(PegasusData%jit)) deallocate(PegasusData%jit)
    if (allocated(PegasusData%kit)) deallocate(PegasusData%kit)
    if (allocated(PegasusData%dxit)) deallocate(PegasusData%dxit)
    if (allocated(PegasusData%dyit)) deallocate(PegasusData%dyit)
    if (allocated(PegasusData%dzit)) deallocate(PegasusData%dzit)
    if (allocated(PegasusData%ibt)) deallocate(PegasusData%ibt)
    if (allocated(PegasusData%jbt)) deallocate(PegasusData%jbt)
    if (allocated(PegasusData%kbt)) deallocate(PegasusData%kbt)
    if (allocated(PegasusData%nit)) deallocate(PegasusData%nit)
    if (allocated(PegasusData%njt)) deallocate(PegasusData%njt)
    if (allocated(PegasusData%nkt)) deallocate(PegasusData%nkt)
    if (allocated(PegasusData%coeffit)) deallocate(PegasusData%coeffit)
    if (allocated(PegasusData%iblank)) deallocate(PegasusData%iblank)

  end subroutine ovkDestroyPegasusData

  subroutine ovkWritePegasusData(PegasusData, HOFilePath, XFilePath, Error)

    type(ovk_pegasus), intent(in) :: PegasusData
    character(len=*), intent(in) :: HOFilePath
    character(len=*), intent(in) :: XFilePath
    integer, intent(out), optional :: Error

    integer :: i, j, p
    integer(lk) :: l
    integer :: Error_
    integer(lk), dimension(3) :: NumPoints
    integer(lk) :: NumPointsTotal

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
        NumPoints(1) = int(PegasusData%ieng(i),kind=lk)
        NumPoints(2) = int(PegasusData%jeng(i),kind=lk)
        NumPoints(3) = int(PegasusData%keng(i),kind=lk)
        NumPointsTotal = product(NumPoints)
        write (HOUnit) (PegasusData%iblank(l,i),l=1_lk,NumPointsTotal)
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
        Error = IO_ERROR
        return
      else
        stop IO_ERROR
      end if
    end if

  end subroutine ovkWritePegasusData

  subroutine ovkPrintPegasusData(PegasusData)

    type(ovk_pegasus), intent(in) :: PegasusData

    integer :: i, j, p
    integer(lk) :: l
    integer(lk), dimension(3) :: NumPoints
    integer(lk) :: NumPointsTotal

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
        NumPoints(1) = int(PegasusData%ieng(i),kind=lk)
        NumPoints(2) = int(PegasusData%jeng(i),kind=lk)
        NumPoints(3) = int(PegasusData%keng(i),kind=lk)
        NumPointsTotal = product(NumPoints)
        do l = 1_lk, NumPointsTotal
          write (*, '(i9)') PegasusData%iblank(l,i)
        end do
      end if
    end do

  end subroutine ovkPrintPegasusData

end module ovkPegasus
