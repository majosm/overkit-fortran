! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkPLOT3D

  use ovkCart
  use ovkField
  use ovkGlobal
  use iso_c_binding
  implicit none

  private

  ! API
  public :: ovk_plot3d_grid_file
  public :: ovk_plot3d_grid_file_
  public :: ovkP3DMachineEndian
  public :: ovkOpenP3D
  public :: ovkCreateP3D
  public :: ovkCloseP3D
  public :: ovkReadP3D
  public :: ovkWriteP3D

  type ovk_plot3d_grid_file
    character(len=PATH_LENGTH) :: path
    integer :: nd
    integer :: ngrids
    integer, dimension(:,:), allocatable :: npoints
    logical :: with_iblank
    integer :: endian
  end type ovk_plot3d_grid_file

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_plot3d_grid_file_
    module procedure ovk_plot3d_grid_file_Default
  end interface ovk_plot3d_grid_file_

  interface ovkOpenP3D
    module procedure ovkOpenP3D_Grid
  end interface ovkOpenP3D

  interface ovkCreateP3D
    module procedure ovkCreateP3D_Grid
  end interface ovkCreateP3D

  interface ovkCloseP3D
    module procedure ovkCloseP3D_Grid
  end interface ovkCloseP3D

  interface ovkReadP3D
    module procedure ovkReadP3D_Grid
    module procedure ovkReadP3D_Grid_IBlank
  end interface ovkReadP3D

  interface ovkWriteP3D
    module procedure ovkWriteP3D_Grid
    module procedure ovkWriteP3D_Grid_IBlank
  end interface ovkWriteP3D

  integer, parameter :: ikoffset = c_long_long

#define IO_ERROR 1

  ! Internal C routines
  interface

    function P3DInternalMachineEndian() result(Endian) bind(C,name="P3DInternalMachineEndian")
      use iso_c_binding, only : c_int
      integer(c_int) :: Endian
    end function P3DInternalMachineEndian

    subroutine P3DInternalDetectGridFormat(FilePath, NumDims, NumGrids, WithIBlank, Endian, Error) &
      bind(C,name="P3DInternalDetectGridFormat")
      use iso_c_binding, only : c_char, c_int
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int) :: NumDims
      integer(c_int) :: NumGrids
      integer(c_int) :: WithIBlank
      integer(c_int) :: Endian
      integer(c_int) :: Error
    end subroutine P3DInternalDetectGridFormat

    subroutine P3DInternalGetGridSize(FilePath, NumDims, Endian, GridID, NumPoints, Error) bind(C, &
      name="P3DInternalGetGridSize")
      use iso_c_binding, only : c_char, c_int
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: NumDims
      integer(c_int), value :: Endian
      integer(c_int), value :: GridID
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int) :: Error
    end subroutine P3DInternalGetGridSize

    function P3DInternalGetGridOffset(NumDims, NumGrids, NumPointsAll, WithIBlank, GridID) &
      result(Offset) bind(C,name="P3DInternalGetGridOffset")
      use iso_c_binding, only : c_int, c_long_long
      integer, parameter :: ikoffset = c_long_long
      integer(c_int), value :: NumDims
      integer(c_int), value :: NumGrids
      integer(c_int), dimension(*) :: NumPointsAll
      integer(c_int), value :: WithIBlank
      integer(c_int), value :: GridID
      integer(ikoffset) :: Offset
    end function P3DInternalGetGridOffset

    subroutine P3DInternalCreateGridFile(FilePath, NumDims, NumGrids, NumPointsAll, WithIBlank, &
      Endian, Error) bind(C,name="P3DInternalCreateGridFile")
      use iso_c_binding, only : c_char, c_int
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: NumDims
      integer(c_int), value :: NumGrids
      integer(c_int), dimension(*) :: NumPointsAll
      integer(c_int), value :: WithIBlank
      integer(c_int), value :: Endian
      integer(c_int) :: Error
    end subroutine P3DInternalCreateGridFile

    subroutine P3DInternalReadSingleGrid2D(FilePath, NumPoints, Endian, Offset, X, Y, Error) &
      bind(C,name="P3DInternalReadSingleGrid2D")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int), value :: Endian
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y
      integer(c_int) :: Error
    end subroutine P3DInternalReadSingleGrid2D

    subroutine P3DInternalReadSingleGrid2DWithIBlank(FilePath, NumPoints, Endian, Offset, X, Y, &
      IBlank, Error) bind(C,name="P3DInternalReadSingleGrid2DWithIBlank")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int), value :: Endian
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y
      integer(c_int), dimension(*) :: IBlank
      integer(c_int) :: Error
    end subroutine P3DInternalReadSingleGrid2DWithIBlank

    subroutine P3DInternalReadSingleGrid3D(FilePath, NumPoints, Endian, Offset, X, Y, Z, Error) &
      bind(C,name="P3DInternalReadSingleGrid3D")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int), value :: Endian
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y, Z
      integer(c_int) :: Error
    end subroutine P3DInternalReadSingleGrid3D

    subroutine P3DInternalReadSingleGrid3DWithIBlank(FilePath, NumPoints, Endian, Offset, X, Y, Z, &
      IBlank, Error) bind(C,name="P3DInternalReadSingleGrid3DWithIBlank")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int), value :: Endian
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y, Z
      integer(c_int), dimension(*) :: IBlank
      integer(c_int) :: Error
    end subroutine P3DInternalReadSingleGrid3DWithIBlank

    subroutine P3DInternalWriteSingleGrid2D(FilePath, NumPoints, Endian, Offset, X, Y, Error) &
      bind(C,name="P3DInternalWriteSingleGrid2D")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int), value :: Endian
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y
      integer(c_int) :: Error
    end subroutine P3DInternalWriteSingleGrid2D

    subroutine P3DInternalWriteSingleGrid2DWithIBlank(FilePath, NumPoints, Endian, Offset, X, Y, &
      IBlank, Error) bind(C,name="P3DInternalWriteSingleGrid2DWithIBlank")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int), value :: Endian
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y
      integer(c_int), dimension(*) :: IBlank
      integer(c_int) :: Error
    end subroutine P3DInternalWriteSingleGrid2DWithIBlank

    subroutine P3DInternalWriteSingleGrid3D(FilePath, NumPoints, Endian, Offset, X, Y, Z, Error) &
      bind(C,name="P3DInternalWriteSingleGrid3D")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int), value :: Endian
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y, Z
      integer(c_int) :: Error
    end subroutine P3DInternalWriteSingleGrid3D

    subroutine P3DInternalWriteSingleGrid3DWithIBlank(FilePath, NumPoints, Endian, Offset, X, Y, &
      Z, IBlank, Error) bind(C,name="P3DInternalWriteSingleGrid3DWithIBlank")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int), value :: Endian
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y, Z
      integer(c_int), dimension(*) :: IBlank
      integer(c_int) :: Error
    end subroutine P3DInternalWriteSingleGrid3DWithIBlank

  end interface

contains

  function ovk_plot3d_grid_file_Default() result(GridFile)

    type(ovk_plot3d_grid_file) :: GridFile

    GridFile%path = ""
    GridFile%nd = 2
    GridFile%ngrids = 0
    GridFile%with_iblank = .false.
    GridFile%endian = P3DInternalMachineEndian()

  end function ovk_plot3d_grid_file_Default

  function ovkP3DMachineEndian() result(Endian)

    integer :: Endian

    Endian = P3DInternalMachineEndian()

  end function ovkP3DMachineEndian

  subroutine ovkOpenP3D_Grid(GridFile, FilePath, Error)

    type(ovk_plot3d_grid_file), intent(out) :: GridFile
    character(len=*), intent(in) :: FilePath
    integer, intent(out), optional :: Error

    integer :: m
    integer :: Error_
    integer :: WithIBlankInt

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Opening PLOT3D grid file ", trim(FilePath), "..."
    end if

    GridFile%path = FilePath

    call P3DInternalDetectGridFormat(NullTerminate(GridFile%path), GridFile%nd, GridFile%ngrids, &
      WithIBlankInt, GridFile%endian, Error_)
    if (Error_ /= 0) goto 999

    GridFile%with_iblank = WithIBlankInt /= 0

    allocate(GridFile%npoints(MAX_ND,GridFile%ngrids))
    do m = 1, GridFile%ngrids
      call P3DInternalGetGridSize(NullTerminate(GridFile%path), GridFile%nd, GridFile%endian, m-1, &
        GridFile%npoints(:,m), Error_)
      if (Error_ /= 0) goto 999
    end do

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Successfully opened PLOT3D grid file ", trim(GridFile%path), "."
      call PrintGridFileInfo(GridFile)
    end if

999 if (Error_ /= 0) then
      if (present(Error)) then
        Error = Error_
        return
      else
        stop IO_ERROR
      end if
    end if

  end subroutine ovkOpenP3D_Grid

  subroutine ovkCreateP3D_Grid(GridFile, FilePath, NumDims, NumGrids, NumPointsAll, WithIBlank, &
    Endian, Error)

    type(ovk_plot3d_grid_file), intent(out) :: GridFile
    character(len=*), intent(in) :: FilePath
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    integer, dimension(:,:), intent(in) :: NumPointsAll
    logical, intent(in), optional :: WithIBlank
    integer, intent(in), optional :: Endian
    integer, intent(out), optional :: Error

    integer :: Error_

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Creating PLOT3D grid file ", trim(FilePath), "..."
    end if

    GridFile%path = FilePath
    GridFile%ngrids = NumGrids
    GridFile%nd = NumDims

    allocate(GridFile%npoints(MAX_ND,NumGrids))
    GridFile%npoints(:NumDims,:) = NumPointsAll(:NumDims,:)
    GridFile%npoints(NumDims+1:,:) = 1

    if (present(WithIBlank)) then
      GridFile%with_iblank = WithIBlank
    else
      GridFile%with_iblank = .false.
    end if

    if (present(Endian)) then
      GridFile%endian = Endian
    else
      GridFile%endian = P3DInternalMachineEndian()
    end if

    call P3DInternalCreateGridFile(NullTerminate(GridFile%path), GridFile%nd, GridFile%ngrids, &
      GridFile%npoints, merge(1,0,GridFile%with_iblank), GridFile%endian, Error_)
    if (Error_ /= 0) goto 999

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Successfully created PLOT3D grid file ", trim(GridFile%path), "."
      call PrintGridFileInfo(GridFile)
    end if

999 if (Error_ /= 0) then
      if (present(Error)) then
        Error = Error_
        return
      else
        stop IO_ERROR
      end if
    end if

  end subroutine ovkCreateP3D_Grid

  subroutine ovkCloseP3D_Grid(GridFile)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile

    character(len=PATH_LENGTH) :: FilePath

    FilePath = GridFile%path

    GridFile = ovk_plot3d_grid_file_()

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Closed grid file ", trim(FilePath), "."
    end if

  end subroutine ovkCloseP3D_Grid

  subroutine ovkReadP3D_Grid(GridFile, GridID, XYZ, Error)

    type(ovk_plot3d_grid_file), intent(in) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), dimension(:), intent(in) :: XYZ
    integer, intent(out), optional :: Error

    integer :: i
    integer :: Error_
    type(ovk_field_int) :: IBlank
    integer(ikoffset) :: Offset

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Reading grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid ID in read from ", &
          trim(GridFile%path), "."
        stop 1
      end if
      do i = 1, GridFile%nd
        if (XYZ(i)%cart%nd /= GridFile%nd) then
          write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid dimension in read from ", &
            trim(GridFile%path), "."
          stop 1
        end if
        if(any(ovkCartSize(XYZ(i)%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
          write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid size in read from ", &
            trim(GridFile%path), "."
          stop 1
        end if
      end do
    end if

    Offset = P3DInternalGetGridOffset(GridFile%nd, GridFile%ngrids, GridFile%npoints, &
      merge(1,0,GridFile%with_iblank), GridID-1)

    if (GridFile%with_iblank) then
      IBlank = ovk_field_int_(XYZ(1)%cart)
    end if

    select case (GridFile%nd)
    case (2)
      if (GridFile%with_iblank) then
        call P3DInternalReadSingleGrid2DWithIBlank(NullTerminate(GridFile%path), &
          GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
          IBlank%values, Error_)
      else
        call P3DInternalReadSingleGrid2D(NullTerminate(GridFile%path), GridFile%npoints(:,GridID), &
          GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, Error_)
      end if
    case (3)
      if (GridFile%with_iblank) then
        call P3DInternalReadSingleGrid3DWithIBlank(NullTerminate(GridFile%path), &
          GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
          XYZ(3)%values, IBlank%values, Error_)
      else
        call P3DInternalReadSingleGrid3D(NullTerminate(GridFile%path), GridFile%npoints(:,GridID), &
          GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, XYZ(3)%values, Error_)
      end if
    end select
    if (Error_ /= 0) goto 999

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Done reading grid ", trim(IntToString(GridID)), "."
    end if

999 if (Error_ /= 0) then
      if (present(Error)) then
        Error = Error_
        return
      else
        stop IO_ERROR
      end if
    end if

  end subroutine ovkReadP3D_Grid

  subroutine ovkReadP3D_Grid_IBlank(GridFile, GridID, XYZ, IBlank, Error)

    type(ovk_plot3d_grid_file), intent(in) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), dimension(:), intent(in) :: XYZ
    type(ovk_field_int), intent(in) :: IBlank
    integer, intent(out), optional :: Error

    integer :: i
    integer :: Error_
    integer(ikoffset) :: Offset

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Reading grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid ID in read from ", &
          trim(GridFile%path), "."
        stop 1
      end if
      do i = 1, GridFile%nd
        if (XYZ(i)%cart%nd /= GridFile%nd) then
          write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid dimension in read from ", &
            trim(GridFile%path), "."
          stop 1
        end if
        if(any(ovkCartSize(XYZ(i)%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
          write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid size in read from ", &
            trim(GridFile%path), "."
          stop 1
        end if
      end do
    end if

    Offset = P3DInternalGetGridOffset(GridFile%nd, GridFile%ngrids, GridFile%npoints, &
      merge(1,0,GridFile%with_iblank), GridID-1)

    select case (GridFile%nd)
    case (2)
      call P3DInternalReadSingleGrid2DWithIBlank(NullTerminate(GridFile%path), &
        GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
        IBlank%values, Error_)
    case (3)
      call P3DInternalReadSingleGrid3DWithIBlank(NullTerminate(GridFile%path), &
        GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
        XYZ(3)%values, IBlank%values, Error_)
    end select
    if (Error_ /= 0) goto 999

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Done reading grid ", trim(IntToString(GridID)), "."
    end if

999 if (Error_ /= 0) then
      if (present(Error)) then
        Error = Error_
        return
      else
        stop IO_ERROR
      end if
    end if

  end subroutine ovkReadP3D_Grid_IBlank

  subroutine ovkWriteP3D_Grid(GridFile, GridID, XYZ, Error)

    type(ovk_plot3d_grid_file), intent(in) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), dimension(:), intent(in) :: XYZ
    integer, intent(out), optional :: Error

    integer :: i
    integer :: Error_
!     type(ovk_field_int) :: IBlank
    integer(ikoffset) :: Offset
    type(ovk_cart) :: ExportCart
    real(rk), dimension(:,:,:,:), allocatable :: XYZCopy
    integer, dimension(:,:,:), allocatable :: IBlankCopy

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Writing grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid ID in write to ", &
          trim(GridFile%path), "."
        stop 1
      end if
      if (.not. GridFile%with_iblank) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to write IBlank data to ", &
          trim(GridFile%path), " which is not formatted for IBlank."
        stop 1
      end if
      do i = 1, GridFile%nd
        if (XYZ(i)%cart%nd /= GridFile%nd) then
          write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid dimension in write to ", &
            trim(GridFile%path), "."
          stop 1
        end if
!         if(any(ovkCartSize(XYZ(i)%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
!           write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid size in write to ", &
!             trim(GridFile%path), "."
!           stop 1
!         end if
      end do
    end if

    Offset = P3DInternalGetGridOffset(GridFile%nd, GridFile%ngrids, GridFile%npoints, 1, GridID-1)

!     if (GridFile%with_iblank) then
!       IBlank = ovk_field_int_(XYZ(1)%cart, 1)
!     end if

    allocate(XYZCopy(GridFile%npoints(1,GridID),GridFile%npoints(2,GridID), &
      GridFile%npoints(3,GridID),GridFile%nd))

    if (any(XYZ(1)%cart%periodic .and. (GridFile%npoints(:,GridID) /= &
      XYZ(1)%cart%ie-XYZ(1)%cart%is+1))) then
      select case (XYZ(1)%cart%periodic_storage)
      case (OVK_OVERLAP_PERIODIC)
        ExportCart = ovk_cart_(GridFile%nd, GridFile%npoints(:,GridID), XYZ(1)%cart%periodic, &
          OVK_NO_OVERLAP_PERIODIC)
      case (OVK_NO_OVERLAP_PERIODIC)
        ExportCart = ovk_cart_(GridFile%nd, GridFile%npoints(:,GridID), XYZ(1)%cart%periodic, &
          OVK_OVERLAP_PERIODIC)
      end select
      do i = 1, GridFile%nd
        call ovkExportField(XYZ(i), ExportCart, XYZCopy(:,:,:,i))
      end do
    else
      do i = 1, GridFile%nd
        XYZCopy(:,:,:,i) = XYZ(i)%values
      end do
    end if

    if (GridFile%with_iblank) then
      allocate(IBlankCopy(GridFile%npoints(1,GridID),GridFile%npoints(2,GridID), &
        GridFile%npoints(3,GridID)))
      IBlankCopy = 1
    end if

    select case (GridFile%nd)
    case (2)
      if (GridFile%with_iblank) then
!         call P3DInternalWriteSingleGrid2DWithIBlank(NullTerminate(GridFile%path), &
!           GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
!           IBlank%values, Error_)
        call P3DInternalWriteSingleGrid2DWithIBlank(NullTerminate(GridFile%path), &
          GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZCopy(:,:,:,1), &
          XYZCopy(:,:,:,2), IBlankCopy, Error_)
      else
!         call P3DInternalWriteSingleGrid2D(NullTerminate(GridFile%path), &
!           GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, Error_)
        call P3DInternalWriteSingleGrid2D(NullTerminate(GridFile%path), &
          GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZCopy(:,:,:,1), &
          XYZCopy(:,:,:,2), Error_)
      end if
    case (3)
      if (GridFile%with_iblank) then
!         call P3DInternalWriteSingleGrid3DWithIBlank(NullTerminate(GridFile%path), &
!           GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
!           XYZ(3)%values, IBlank%values, Error_)
        call P3DInternalWriteSingleGrid3DWithIBlank(NullTerminate(GridFile%path), &
          GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZCopy(:,:,:,1), &
          XYZCopy(:,:,:,2), XYZCopy(:,:,:,3), IBlankCopy, Error_)
      else
!         call P3DInternalWriteSingleGrid3D(NullTerminate(GridFile%path), &
!           GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
!           XYZ(3)%values, Error_)
        call P3DInternalWriteSingleGrid3D(NullTerminate(GridFile%path), &
          GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZCopy(:,:,:,1), &
          XYZCopy(:,:,:,2), XYZCopy(:,:,:,3), Error_)
      end if
    end select
    if (Error_ /= 0) goto 999

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Done writing grid ", trim(IntToString(GridID)), "."
    end if

999 if (Error_ /= 0) then
      if (present(Error)) then
        Error = Error_
        return
      else
        stop IO_ERROR
      end if
    end if

  end subroutine ovkWriteP3D_Grid

  subroutine ovkWriteP3D_Grid_IBlank(GridFile, GridID, XYZ, IBlank, Error)

    type(ovk_plot3d_grid_file), intent(in) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), dimension(:), intent(in) :: XYZ
    type(ovk_field_int), intent(in) :: IBlank
    integer, intent(out), optional :: Error

    integer :: i
    integer :: Error_
    integer(ikoffset) :: Offset
    type(ovk_cart) :: ExportCart
    real(rk), dimension(:,:,:,:), allocatable :: XYZCopy
    integer, dimension(:,:,:), allocatable :: IBlankCopy

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Writing grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid ID in write to ", &
          trim(GridFile%path), "."
        stop 1
      end if
      if (.not. GridFile%with_iblank) then
        write (ERROR_UNIT, '(3a)') "ERROR: Attempted to write IBlank data to ", &
          trim(GridFile%path), " which is not formatted for IBlank."
        stop 1
      end if
      do i = 1, GridFile%nd
        if (XYZ(i)%cart%nd /= GridFile%nd) then
          write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid dimension in write to ", &
            trim(GridFile%path), "."
          stop 1
        end if
!         if(any(ovkCartSize(XYZ(i)%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
!           write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid size in write to ", &
!             trim(GridFile%path), "."
!           stop 1
!         end if
      end do
!       if(any(ovkCartSize(IBlank%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
!         write (ERROR_UNIT, '(3a)') "ERROR: Incorrect grid size in write to ", &
!           trim(GridFile%path), "."
!         stop 1
!       end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%nd, GridFile%ngrids, GridFile%npoints, 1, GridID-1)

    allocate(XYZCopy(GridFile%npoints(1,GridID),GridFile%npoints(2,GridID), &
      GridFile%npoints(3,GridID),GridFile%nd))
    allocate(IBlankCopy(GridFile%npoints(1,GridID),GridFile%npoints(2,GridID), &
      GridFile%npoints(3,GridID)))

    if (any(XYZ(1)%cart%periodic .and. (GridFile%npoints(:,GridID) /= &
      XYZ(1)%cart%ie-XYZ(1)%cart%is+1))) then
      select case (XYZ(1)%cart%periodic_storage)
      case (OVK_OVERLAP_PERIODIC)
        ExportCart = ovk_cart_(GridFile%nd, GridFile%npoints(:,GridID), XYZ(1)%cart%periodic, &
          OVK_NO_OVERLAP_PERIODIC)
      case (OVK_NO_OVERLAP_PERIODIC)
        ExportCart = ovk_cart_(GridFile%nd, GridFile%npoints(:,GridID), XYZ(1)%cart%periodic, &
          OVK_OVERLAP_PERIODIC)
      end select
      do i = 1, GridFile%nd
        call ovkExportField(XYZ(i), ExportCart, XYZCopy(:,:,:,i))
      end do
      call ovkExportField(IBlank, ExportCart, IBlankCopy)
    else
      do i = 1, GridFile%nd
        XYZCopy(:,:,:,i) = XYZ(i)%values
      end do
      IBlankCopy = IBlank%values
    end if

    select case (GridFile%nd)
    case (2)
!       call P3DInternalWriteSingleGrid2DWithIBlank(NullTerminate(GridFile%path), &
!         GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
!         IBlank%values, Error_)
      call P3DInternalWriteSingleGrid2DWithIBlank(NullTerminate(GridFile%path), &
        GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZCopy(:,:,:,1), &
        XYZCopy(:,:,:,2), IBlankCopy, Error_)
    case (3)
!       call P3DInternalWriteSingleGrid3DWithIBlank(NullTerminate(GridFile%path), &
!         GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZ(1)%values, XYZ(2)%values, &
!         XYZ(3)%values, IBlank%values, Error_)
      call P3DInternalWriteSingleGrid3DWithIBlank(NullTerminate(GridFile%path), &
        GridFile%npoints(:,GridID), GridFile%endian, Offset, XYZCopy(:,:,:,1), &
        XYZCopy(:,:,:,2), XYZCopy(:,:,:,3), IBlankCopy, Error_)
    end select
    if (Error_ /= 0) goto 999

    if (OVK_VERBOSE) then
      write (*, '(3a)') "Done writing grid ", trim(IntToString(GridID)), "."
    end if

999 if (Error_ /= 0) then
      if (present(Error)) then
        Error = Error_
        return
      else
        stop IO_ERROR
      end if
    end if

  end subroutine ovkWriteP3D_Grid_IBlank

  pure function NullTerminate(String) result(NullTerminatedString)

    character(len=*), intent(in) :: String
    character(len=STRING_LENGTH) :: NullTerminatedString

    NullTerminatedString = trim(String) // C_NULL_CHAR

  end function NullTerminate

  subroutine PrintGridFileInfo(GridFile)

    type(ovk_plot3d_grid_file), intent(in) :: GridFile

    integer(lk) :: NumPointsTotal
    character(len=STRING_LENGTH) :: NumPointsTotalString
    character(len=STRING_LENGTH) :: NiString, NjString, NkString

    integer :: m

    write (*, '(a)') "File details:"
    write (*, '(3a)') "* Dimension: ", trim(IntToString(GridFile%nd)), "D"
    if (GridFile%with_iblank) then
      write (*, '(a)') "* IBlank: yes"
    else
      write (*, '(a)') "* IBlank: no"
    end if
    if (GridFile%endian == OVK_LITTLE_ENDIAN) then
      write (*, '(a)') "* Endianness: little"
    else
      write (*, '(a)') "* Endianness: big"
    end if
    write (*, '(2a)') "* Number of grids: ", trim(IntToString(GridFile%ngrids))
    NumPointsTotal = 0
    do m = 1, GridFile%ngrids
      NumPointsTotal = NumPointsTotal + product(int(GridFile%npoints(:,m),kind=lk))
    end do
    NumPointsTotalString = LargeIntToString(NumPointsTotal)
    write (*, '(2a)') "* Total number of grid points: ", trim(NumPointsTotalString)
    do m = 1, GridFile%ngrids
      NumPointsTotalString = LargeIntToString(product(int(GridFile%npoints(:,m),kind=lk)))
      NiString = IntToString(GridFile%npoints(1,m))
      NjString = IntToString(GridFile%npoints(2,m))
      NkString = IntToString(GridFile%npoints(3,m))
      write (*, '(3a)', advance="no") "* Grid ", trim(IntToString(m)), ": "
      write (*, '(2a)', advance="no") trim(NumPointsTotalString), " points "
      select case (GridFile%nd)
      case (2)
        write (*, '(5a)', advance="no") "(Ni=", trim(NiString), ", Nj=", trim(NjString), ")"
      case (3)
        write (*, '(7a)', advance="no") "(Ni=", trim(NiString), ", Nj=", trim(NjString), &
          ", Nk=", trim(NkString), ")"
      end select
      write (*, '(a)') ""
    end do

  end subroutine PrintGridFileInfo

end module ovkPLOT3D
