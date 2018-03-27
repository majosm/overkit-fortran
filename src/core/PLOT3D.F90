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
    type(t_noconstruct) :: noconstruct
    logical :: verbose
    character(len=PATH_LENGTH) :: path
    integer :: endian
    integer :: p3d_format
    integer :: nd
    integer :: ngrids
    logical :: with_iblank
    integer, dimension(:,:), allocatable :: npoints
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
    module procedure ovkReadP3D_Grid_2D
    module procedure ovkReadP3D_Grid_3D
    module procedure ovkReadP3D_Grid_IBlank_2D
    module procedure ovkReadP3D_Grid_IBlank_3D
  end interface ovkReadP3D

  interface ovkWriteP3D
    module procedure ovkWriteP3D_Grid_2D
    module procedure ovkWriteP3D_Grid_3D
    module procedure ovkWriteP3D_Grid_IBlank_2D
    module procedure ovkWriteP3D_Grid_IBlank_3D
  end interface ovkWriteP3D

  integer, parameter :: ikoffset = c_long_long

#define IO_ERROR 1

  ! Internal C routines
  interface

    function P3DInternalMachineEndian() result(Endian) bind(C,name="P3DInternalMachineEndian")
      use iso_c_binding, only : c_int
      integer(c_int) :: Endian
    end function P3DInternalMachineEndian

    subroutine P3DInternalDetectGridFormat(FilePath, Endian, P3DFormat, NumDims, NumGrids, &
      WithIBlank, Error) bind(C,name="P3DInternalDetectGridFormat")
      use iso_c_binding, only : c_char, c_int
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int) :: Endian
      integer(c_int) :: P3DFormat
      integer(c_int) :: NumDims
      integer(c_int) :: NumGrids
      integer(c_int) :: WithIBlank
      integer(c_int) :: Error
    end subroutine P3DInternalDetectGridFormat

    subroutine P3DInternalGetGridSize(FilePath, Endian, P3DFormat, NumDims, GridID, NumPoints, &
      Error) bind(C,name="P3DInternalGetGridSize")
      use iso_c_binding, only : c_char, c_int
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), value :: NumDims
      integer(c_int), value :: GridID
      integer(c_int), dimension(*) :: NumPoints
      integer(c_int) :: Error
    end subroutine P3DInternalGetGridSize

    function P3DInternalGetGridOffset(P3DFormat, WithIBlank, NumDims, NumGrids, NumPointsAll, &
      GridID) result(Offset) bind(C,name="P3DInternalGetGridOffset")
      use iso_c_binding, only : c_int, c_long_long
      integer, parameter :: ikoffset = c_long_long
      integer(c_int), value :: P3DFormat
      integer(c_int), value :: WithIBlank
      integer(c_int), value :: NumDims
      integer(c_int), value :: NumGrids
      integer(c_int), dimension(*) :: NumPointsAll
      integer(c_int), value :: GridID
      integer(ikoffset) :: Offset
    end function P3DInternalGetGridOffset

    subroutine P3DInternalCreateGridFile(FilePath, Endian, P3DFormat, NumDims, NumGrids, &
      WithIBlank, NumPointsAll, Error) bind(C,name="P3DInternalCreateGridFile")
      use iso_c_binding, only : c_char, c_int
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), value :: NumDims
      integer(c_int), value :: NumGrids
      integer(c_int), value :: WithIBlank
      integer(c_int), dimension(*) :: NumPointsAll
      integer(c_int) :: Error
    end subroutine P3DInternalCreateGridFile

    subroutine P3DInternalReadSingleGrid2D(FilePath, Endian, P3DFormat, NumPoints, Offset, X, Y, &
      Error) bind(C,name="P3DInternalReadSingleGrid2D")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), dimension(*) :: NumPoints
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y
      integer(c_int) :: Error
    end subroutine P3DInternalReadSingleGrid2D

    subroutine P3DInternalReadSingleGrid2DWithIBlank(FilePath, Endian, P3DFormat, NumPoints, &
      Offset, X, Y, IBlank, Error) bind(C,name="P3DInternalReadSingleGrid2DWithIBlank")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), dimension(*) :: NumPoints
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y
      integer(c_int), dimension(*) :: IBlank
      integer(c_int) :: Error
    end subroutine P3DInternalReadSingleGrid2DWithIBlank

    subroutine P3DInternalReadSingleGrid3D(FilePath, Endian, P3DFormat, NumPoints, Offset, X, Y, &
      Z, Error) bind(C,name="P3DInternalReadSingleGrid3D")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), dimension(*) :: NumPoints
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y, Z
      integer(c_int) :: Error
    end subroutine P3DInternalReadSingleGrid3D

    subroutine P3DInternalReadSingleGrid3DWithIBlank(FilePath, Endian, P3DFormat, NumPoints, &
      Offset, X, Y, Z, IBlank, Error) bind(C,name="P3DInternalReadSingleGrid3DWithIBlank")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), dimension(*) :: NumPoints
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y, Z
      integer(c_int), dimension(*) :: IBlank
      integer(c_int) :: Error
    end subroutine P3DInternalReadSingleGrid3DWithIBlank

    subroutine P3DInternalWriteSingleGrid2D(FilePath, Endian, P3DFormat, NumPoints, Offset, X, Y, &
      Error) bind(C,name="P3DInternalWriteSingleGrid2D")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), dimension(*) :: NumPoints
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y
      integer(c_int) :: Error
    end subroutine P3DInternalWriteSingleGrid2D

    subroutine P3DInternalWriteSingleGrid2DWithIBlank(FilePath, Endian, P3DFormat, NumPoints, &
      Offset, X, Y, IBlank, Error) bind(C,name="P3DInternalWriteSingleGrid2DWithIBlank")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), dimension(*) :: NumPoints
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y
      integer(c_int), dimension(*) :: IBlank
      integer(c_int) :: Error
    end subroutine P3DInternalWriteSingleGrid2DWithIBlank

    subroutine P3DInternalWriteSingleGrid3D(FilePath, Endian, P3DFormat, NumPoints, Offset, X, Y, &
      Z, Error) bind(C,name="P3DInternalWriteSingleGrid3D")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), dimension(*) :: NumPoints
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y, Z
      integer(c_int) :: Error
    end subroutine P3DInternalWriteSingleGrid3D

    subroutine P3DInternalWriteSingleGrid3DWithIBlank(FilePath, Endian, P3DFormat, NumPoints, &
      Offset, X, Y, Z, IBlank, Error) bind(C,name="P3DInternalWriteSingleGrid3DWithIBlank")
      use iso_c_binding, only : c_char, c_int, c_double, c_long_long
      integer, parameter :: ikoffset = c_long_long
      character(kind=c_char), dimension(*) :: FilePath
      integer(c_int), value :: Endian
      integer(c_int), value :: P3DFormat
      integer(c_int), dimension(*) :: NumPoints
      integer(ikoffset), value :: Offset
      real(c_double), dimension(*) :: X, Y, Z
      integer(c_int), dimension(*) :: IBlank
      integer(c_int) :: Error
    end subroutine P3DInternalWriteSingleGrid3DWithIBlank

  end interface

contains

  function ovk_plot3d_grid_file_Default() result(GridFile)

    type(ovk_plot3d_grid_file) :: GridFile

    GridFile%verbose = .true.
    GridFile%path = ""
    GridFile%nd = 2
    GridFile%ngrids = 0
    GridFile%with_iblank = .false.
    GridFile%endian = P3DInternalMachineEndian()
    GridFile%p3d_format = OVK_P3D_STANDARD

  end function ovk_plot3d_grid_file_Default

  function ovkP3DMachineEndian() result(Endian)

    integer :: Endian

    Endian = P3DInternalMachineEndian()

  end function ovkP3DMachineEndian

  subroutine ovkOpenP3D_Grid(GridFile, FilePath, Verbose, Error)

    type(ovk_plot3d_grid_file), intent(out) :: GridFile
    character(len=*), intent(in) :: FilePath
    logical, intent(in), optional :: Verbose
    integer, intent(out), optional :: Error
    
    integer :: m
    integer :: Error_
    integer :: WithIBlankInt

    if (present(Verbose)) then
      GridFile%verbose = Verbose
    else
      GridFile%verbose = .true.
    end if

    if (GridFile%verbose) then
      write (*, '(3a)') "Opening PLOT3D grid file ", trim(FilePath), "..."
    end if

    GridFile%path = FilePath

    call P3DInternalDetectGridFormat(NullTerminate(GridFile%path), GridFile%endian, &
      GridFile%p3d_format, GridFile%nd, GridFile%ngrids, WithIBlankInt, Error_)
    if (Error_ /= 0) goto 999

    GridFile%with_iblank = WithIBlankInt /= 0

    allocate(GridFile%npoints(MAX_ND,GridFile%ngrids))
    do m = 1, GridFile%ngrids
      call P3DInternalGetGridSize(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%nd, m-1, GridFile%npoints(:,m), Error_)
      if (Error_ /= 0) goto 999
    end do

    if (GridFile%verbose) then
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
    Endian, P3DFormat, Verbose, Error)

    type(ovk_plot3d_grid_file), intent(out) :: GridFile
    character(len=*), intent(in) :: FilePath
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    integer, dimension(:,:), intent(in) :: NumPointsAll
    logical, intent(in), optional :: WithIBlank
    integer, intent(in), optional :: Endian
    integer, intent(in), optional :: P3DFormat
    logical, intent(in), optional :: Verbose
    integer, intent(out), optional :: Error

    integer :: Error_

    if (present(Verbose)) then
      GridFile%verbose = Verbose
    else
      GridFile%verbose = .true.
    end if

    if (GridFile%verbose) then
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

    if (present(P3DFormat)) then
      GridFile%p3d_format = P3DFormat
    else
      GridFile%p3d_format = OVK_P3D_STANDARD
    end if

    call P3DInternalCreateGridFile(NullTerminate(GridFile%path), GridFile%endian, &
      GridFile%p3d_format, GridFile%nd, GridFile%ngrids, merge(1,0,GridFile%with_iblank), &
      GridFile%npoints, Error_)
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

    logical :: Verbose
    character(len=PATH_LENGTH) :: FilePath

    Verbose = GridFile%verbose
    FilePath = GridFile%path

    GridFile = ovk_plot3d_grid_file_()

    if (Verbose) then
      write (*, '(3a)') "Closed grid file ", trim(FilePath), "."
    end if

  end subroutine ovkCloseP3D_Grid

  subroutine ovkReadP3D_Grid_2D(GridFile, GridID, X, Y, Error)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), intent(inout) :: X
    type(ovk_field_real), intent(inout) :: Y
    integer, intent(out), optional :: Error

    integer :: Error_
    type(ovk_field_int) :: IBlank
    integer(ikoffset) :: Offset

    if (GridFile%verbose) then
      write (*, '(3a)') "Reading grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridFile%nd /= 2) then
        write (ERROR_UNIT, '(a)') "ERROR: Incorrect number of coordinate fields specified."
        stop 1
      end if
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      if (X%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect dimension."
        stop 1
      end if
      if (Y%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect dimension."
        stop 1
      end if
      if (any(ovkCartSize(X%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect size."
        stop 1
      end if
      if (any(ovkCartSize(Y%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect size."
        stop 1
      end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%p3d_format, GridFile%nd, GridFile%ngrids, &
      merge(1,0,GridFile%with_iblank), GridFile%npoints, GridID-1)

    if (GridFile%with_iblank) then
      IBlank = ovk_field_int_(X%cart)
      call P3DInternalReadSingleGrid2DWithIBlank(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, &
        IBlank%values, Error_)
    else
      call P3DInternalReadSingleGrid2D(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, Error_)
    end if
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

  end subroutine ovkReadP3D_Grid_2D

  subroutine ovkReadP3D_Grid_3D(GridFile, GridID, X, Y, Z, Error)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), intent(inout) :: X
    type(ovk_field_real), intent(inout) :: Y
    type(ovk_field_real), intent(inout) :: Z
    integer, intent(out), optional :: Error

    integer :: Error_
    type(ovk_field_int) :: IBlank
    integer(ikoffset) :: Offset

    if (GridFile%verbose) then
      write (*, '(3a)') "Reading grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridFile%nd /= 3) then
        write (ERROR_UNIT, '(a)') "ERROR: Incorrect number of coordinate fields specified."
        stop 1
      end if
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      if (X%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect dimension."
        stop 1
      end if
      if (Y%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect dimension."
        stop 1
      end if
      if (Z%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Z coordinate field has incorrect dimension."
        stop 1
      end if
      if(any(ovkCartSize(X%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Y%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Z%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Z coordinate field has incorrect size."
        stop 1
      end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%p3d_format, GridFile%nd, GridFile%ngrids, &
      merge(1,0,GridFile%with_iblank), GridFile%npoints, GridID-1)

    if (GridFile%with_iblank) then
      IBlank = ovk_field_int_(X%cart)
      call P3DInternalReadSingleGrid3DWithIBlank(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, Z%values, &
        IBlank%values, Error_)
    else
      call P3DInternalReadSingleGrid3D(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, Z%values, &
        Error_)
    end if
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

  end subroutine ovkReadP3D_Grid_3D

  subroutine ovkReadP3D_Grid_IBlank_2D(GridFile, GridID, X, Y, IBlank, Error)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), intent(inout) :: X
    type(ovk_field_real), intent(inout) :: Y
    type(ovk_field_int), intent(inout) :: IBlank
    integer, intent(out), optional :: Error

    integer :: Error_
    integer(ikoffset) :: Offset

    if (GridFile%verbose) then
      write (*, '(3a)') "Reading grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridFile%nd /= 2) then
        write (ERROR_UNIT, '(a)') "ERROR: Incorrect number of coordinate fields specified."
        stop 1
      end if
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      if (.not. GridFile%with_iblank) then
        write (ERROR_UNIT, '(a)') "ERROR: File does not contain IBlank."
        stop 1
      end if
      if (X%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect dimension."
        stop 1
      end if
      if (Y%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect dimension."
        stop 1
      end if
      if (IBlank%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: IBlank field has incorrect dimension."
        stop 1
      end if
      if(any(ovkCartSize(X%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Y%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(IBlank%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: IBlank field has incorrect size."
        stop 1
      end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%p3d_format, GridFile%nd, GridFile%ngrids, 1, &
      GridFile%npoints, GridID-1)

    call P3DInternalReadSingleGrid2DWithIBlank(NullTerminate(GridFile%path), GridFile%endian, &
      GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, IBlank%values, &
      Error_)
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

  end subroutine ovkReadP3D_Grid_IBlank_2D

  subroutine ovkReadP3D_Grid_IBlank_3D(GridFile, GridID, X, Y, Z, IBlank, Error)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), intent(inout) :: X
    type(ovk_field_real), intent(inout) :: Y
    type(ovk_field_real), intent(inout) :: Z
    type(ovk_field_int), intent(inout) :: IBlank
    integer, intent(out), optional :: Error

    integer :: Error_
    integer(ikoffset) :: Offset

    if (GridFile%verbose) then
      write (*, '(3a)') "Reading grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridFile%nd /= 3) then
        write (ERROR_UNIT, '(a)') "ERROR: Incorrect number of coordinate fields specified."
        stop 1
      end if
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      if (.not. GridFile%with_iblank) then
        write (ERROR_UNIT, '(a)') "ERROR: File does not contain IBlank."
        stop 1
      end if
      if (X%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect dimension."
        stop 1
      end if
      if (Y%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect dimension."
        stop 1
      end if
      if (Z%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Z coordinate field has incorrect dimension."
        stop 1
      end if
      if (IBlank%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: IBlank field has incorrect dimension."
        stop 1
      end if
      if(any(ovkCartSize(X%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Y%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Z%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Z coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(IBlank%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: IBlank field has incorrect size."
        stop 1
      end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%p3d_format, GridFile%nd, GridFile%ngrids, 1, &
      GridFile%npoints, GridID-1)

    call P3DInternalReadSingleGrid3DWithIBlank(NullTerminate(GridFile%path), GridFile%endian, &
      GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, Z%values, &
      IBlank%values, Error_)
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

  end subroutine ovkReadP3D_Grid_IBlank_3D

  subroutine ovkWriteP3D_Grid_2D(GridFile, GridID, X, Y, Error)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), intent(in) :: X
    type(ovk_field_real), intent(in) :: Y
    integer, intent(out), optional :: Error

    integer :: Error_
    type(ovk_field_int) :: IBlank
    integer(ikoffset) :: Offset

    if (GridFile%verbose) then
      write (*, '(3a)') "Writing grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridFile%nd /= 2) then
        write (ERROR_UNIT, '(a)') "ERROR: Incorrect number of coordinate fields specified."
        stop 1
      end if
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      if (X%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect dimension."
        stop 1
      end if
      if (Y%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect dimension."
        stop 1
      end if
      if(any(ovkCartSize(X%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Y%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect size."
        stop 1
      end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%p3d_format, GridFile%nd, GridFile%ngrids, &
      merge(1,0,GridFile%with_iblank), GridFile%npoints, GridID-1)

    if (GridFile%with_iblank) then
      IBlank = ovk_field_int_(X%cart, 1)
      call P3DInternalWriteSingleGrid2DWithIBlank(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, &
        IBlank%values, Error_)
    else
      call P3DInternalWriteSingleGrid2D(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, Error_)
    end if
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

  end subroutine ovkWriteP3D_Grid_2D

  subroutine ovkWriteP3D_Grid_3D(GridFile, GridID, X, Y, Z, Error)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), intent(in) :: X
    type(ovk_field_real), intent(in) :: Y
    type(ovk_field_real), intent(in) :: Z
    integer, intent(out), optional :: Error

    integer :: Error_
    type(ovk_field_int) :: IBlank
    integer(ikoffset) :: Offset

    if (GridFile%verbose) then
      write (*, '(3a)') "Writing grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridFile%nd /= 3) then
        write (ERROR_UNIT, '(a)') "ERROR: Incorrect number of coordinate fields specified."
        stop 1
      end if
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      if (X%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect dimension."
        stop 1
      end if
      if (Y%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect dimension."
        stop 1
      end if
      if (Z%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Z coordinate field has incorrect dimension."
        stop 1
      end if
      if(any(ovkCartSize(X%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Y%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Z%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Z coordinate field has incorrect size."
        stop 1
      end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%p3d_format, GridFile%nd, GridFile%ngrids, &
      merge(1,0,GridFile%with_iblank), GridFile%npoints, GridID-1)

    if (GridFile%with_iblank) then
      IBlank = ovk_field_int_(X%cart, 1)
      call P3DInternalWriteSingleGrid3DWithIBlank(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, Z%values, &
        IBlank%values, Error_)
    else
      call P3DInternalWriteSingleGrid3D(NullTerminate(GridFile%path), GridFile%endian, &
        GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, Z%values, &
        Error_)
    end if
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

  end subroutine ovkWriteP3D_Grid_3D

  subroutine ovkWriteP3D_Grid_IBlank_2D(GridFile, GridID, X, Y, IBlank, Error)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), intent(in) :: X
    type(ovk_field_real), intent(in) :: Y
    type(ovk_field_int), intent(in) :: IBlank
    integer, intent(out), optional :: Error

    integer :: Error_
    integer(ikoffset) :: Offset

    if (GridFile%verbose) then
      write (*, '(3a)') "Writing grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridFile%nd /= 2) then
        write (ERROR_UNIT, '(a)') "ERROR: Incorrect number of coordinate fields specified."
        stop 1
      end if
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      if (.not. GridFile%with_iblank) then
        write (ERROR_UNIT, '(a)') "ERROR: File is not formatted for IBlank."
        stop 1
      end if
      if (X%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect dimension."
        stop 1
      end if
      if (Y%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect dimension."
        stop 1
      end if
      if (IBlank%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: IBlank field has incorrect dimension."
        stop 1
      end if
      if(any(ovkCartSize(X%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Y%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(IBlank%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: IBlank field has incorrect size."
        stop 1
      end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%p3d_format, GridFile%nd, GridFile%ngrids, &
      1, GridFile%npoints, GridID-1)

    call P3DInternalWriteSingleGrid2DWithIBlank(NullTerminate(GridFile%path), GridFile%endian, &
      GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, IBlank%values, &
      Error_)
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

  end subroutine ovkWriteP3D_Grid_IBlank_2D

  subroutine ovkWriteP3D_Grid_IBlank_3D(GridFile, GridID, X, Y, Z, IBlank, Error)

    type(ovk_plot3d_grid_file), intent(inout) :: GridFile
    integer, intent(in) :: GridID
    type(ovk_field_real), intent(in) :: X
    type(ovk_field_real), intent(in) :: Y
    type(ovk_field_real), intent(in) :: Z
    type(ovk_field_int), intent(in) :: IBlank
    integer, intent(out), optional :: Error

    integer :: Error_
    integer(ikoffset) :: Offset

    if (GridFile%verbose) then
      write (*, '(3a)') "Writing grid ", trim(IntToString(GridID)), "..."
    end if

    if (OVK_DEBUG) then
      if (GridFile%nd /= 3) then
        write (ERROR_UNIT, '(a)') "ERROR: Incorrect number of coordinate fields specified."
        stop 1
      end if
      if (GridID <= 0 .or. GridID > GridFile%ngrids) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid grid ID."
        stop 1
      end if
      if (.not. GridFile%with_iblank) then
        write (ERROR_UNIT, '(a)') "ERROR: File is not formatted for IBlank."
        stop 1
      end if
      if (X%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect dimension."
        stop 1
      end if
      if (Y%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect dimension."
        stop 1
      end if
      if (Z%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: Z coordinate field has incorrect dimension."
        stop 1
      end if
      if (IBlank%cart%nd /= GridFile%nd) then
        write (ERROR_UNIT, '(a)') "ERROR: IBlank field has incorrect dimension."
        stop 1
      end if
      if(any(ovkCartSize(X%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: X coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Y%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Y coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(Z%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: Z coordinate field has incorrect size."
        stop 1
      end if
      if(any(ovkCartSize(IBlank%cart) /= GridFile%npoints(:GridFile%nd,GridID))) then
        write (ERROR_UNIT, '(a)') "ERROR: IBlank field has incorrect size."
        stop 1
      end if
    end if

    Offset = P3DInternalGetGridOffset(GridFile%p3d_format, GridFile%nd, GridFile%ngrids, 1, &
      GridFile%npoints, GridID-1)

    call P3DInternalWriteSingleGrid3DWithIBlank(NullTerminate(GridFile%path), GridFile%endian, &
      GridFile%p3d_format, GridFile%npoints(:,GridID), Offset, X%values, Y%values, Z%values, &
      IBlank%values, Error_)
    if (Error_ /= 0) goto 999

    if (GridFile%verbose) then
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

  end subroutine ovkWriteP3D_Grid_IBlank_3D

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
    if (GridFile%endian == OVK_LITTLE_ENDIAN) then
      write (*, '(a)') "* Endianness: little"
    else
      write (*, '(a)') "* Endianness: big"
    end if
    if (GridFile%p3d_format == OVK_P3D_STANDARD) then
      write (*, '(a)') "* Format: standard"
    else
      write (*, '(a)') "* Format: extended"
    end if
    if (GridFile%with_iblank) then
      write (*, '(a)') "* IBlank: yes"
    else
      write (*, '(a)') "* IBlank: no"
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
