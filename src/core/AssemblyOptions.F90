! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkAssemblyOptions

  use ovkGlobal
  implicit none

  private

  ! API
  public :: ovk_assembly_options
  public :: ovk_assembly_options_
  public :: ovkGetAssemblyOptionsDimension
  public :: ovkGetAssemblyOptionsGridCount
  public :: ovkGetAssemblyOptionOverlappable
  public :: ovkSetAssemblyOptionOverlappable
  public :: ovkGetAssemblyOptionOverlapTolerance
  public :: ovkSetAssemblyOptionOverlapTolerance
  public :: ovkGetAssemblyOptionOverlapAccelQualityAdjust
  public :: ovkSetAssemblyOptionOverlapAccelQualityAdjust
  public :: ovkGetAssemblyOptionInferBoundaries
  public :: ovkSetAssemblyOptionInferBoundaries
  public :: ovkGetAssemblyOptionCutBoundaryHoles
  public :: ovkSetAssemblyOptionCutBoundaryHoles
  public :: ovkGetAssemblyOptionOccludes
  public :: ovkSetAssemblyOptionOccludes
  public :: ovkGetAssemblyOptionEdgePadding
  public :: ovkSetAssemblyOptionEdgePadding
  public :: ovkGetAssemblyOptionEdgeSmoothing
  public :: ovkSetAssemblyOptionEdgeSmoothing
  public :: ovkGetAssemblyOptionConnectionType
  public :: ovkSetAssemblyOptionConnectionType
  public :: ovkGetAssemblyOptionFringeSize
  public :: ovkSetAssemblyOptionFringeSize
  public :: ovkGetAssemblyOptionOverlapMinimization
  public :: ovkSetAssemblyOptionOverlapMinimization

  type ovk_assembly_options
    type(t_noconstruct) :: noconstruct
    integer :: nd
    integer :: ngrids
    logical, dimension(:,:), allocatable :: overlappable
    real(rk), dimension(:,:), allocatable :: overlap_tolerance
    real(rk), dimension(:), allocatable :: overlap_accel_quality_adjust
    logical, dimension(:), allocatable :: infer_boundaries
    logical, dimension(:,:), allocatable :: cut_boundary_holes
    integer, dimension(:,:), allocatable :: occludes
    integer, dimension(:,:), allocatable :: edge_padding
    integer, dimension(:), allocatable :: edge_smoothing
    integer, dimension(:,:), allocatable :: connection_type
    integer, dimension(:), allocatable :: fringe_size
    logical, dimension(:,:), allocatable :: overlap_minimization
  end type ovk_assembly_options

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_assembly_options_
    module procedure ovk_assembly_options_Default
    module procedure ovk_assembly_options_Assigned
  end interface ovk_assembly_options_

contains

  function ovk_assembly_options_Default() result(Options)

    type(ovk_assembly_options) :: Options

    Options%nd = 2
    Options%ngrids = 0

  end function ovk_assembly_options_Default

  function ovk_assembly_options_Assigned(NumDims, NumGrids) result(Options)

    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    type(ovk_assembly_options) :: Options

    Options%nd = NumDims
    Options%ngrids = NumGrids

    allocate(Options%overlappable(NumGrids,NumGrids))
    Options%overlappable = .false.

    allocate(Options%overlap_tolerance(NumGrids,NumGrids))
    Options%overlap_tolerance = 1.e-12_rk

    allocate(Options%overlap_accel_quality_adjust(NumGrids))
    Options%overlap_accel_quality_adjust = 0._rk

    allocate(Options%infer_boundaries(NumGrids))
    Options%infer_boundaries = .false.

    allocate(Options%cut_boundary_holes(NumGrids,NumGrids))
    Options%cut_boundary_holes = .false.

    allocate(Options%occludes(NumGrids,NumGrids))
    Options%occludes = OVK_OCCLUDES_NONE

    allocate(Options%edge_padding(NumGrids,NumGrids))
    Options%edge_padding = 0

    allocate(Options%edge_smoothing(NumGrids))
    Options%edge_smoothing = 0

    allocate(Options%connection_type(NumGrids,NumGrids))
    Options%connection_type = OVK_CONNECTION_NONE

    allocate(Options%fringe_size(NumGrids))
    Options%fringe_size = 0

    allocate(Options%overlap_minimization(NumGrids,NumGrids))
    Options%overlap_minimization = .false.

  end function ovk_assembly_options_Assigned

  subroutine ovkGetAssemblyOptionsDimension(Options, NumDims)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(out) :: NumDims

    NumDims = Options%nd

  end subroutine ovkGetAssemblyOptionsDimension

  subroutine ovkGetAssemblyOptionsGridCount(Options, NumGrids)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(out) :: NumGrids

    NumGrids = Options%ngrids

  end subroutine ovkGetAssemblyOptionsGridCount

  subroutine ovkGetAssemblyOptionOverlappable(Options, OverlappingGridID, OverlappedGridID, &
    Overlappable)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(out) :: Overlappable

    Overlappable = Options%overlappable(OverlappingGridID,OverlappedGridID)

  end subroutine ovkGetAssemblyOptionOverlappable

  subroutine ovkSetAssemblyOptionOverlappable(Options, OverlappingGridID, OverlappedGridID, &
    Overlappable)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(in) :: Overlappable

    integer :: m, n
    integer :: ms, me, ns, ne

    call GridIDRange(Options%ngrids, OverlappingGridID, ms, me)
    call GridIDRange(Options%ngrids, OverlappedGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Options%overlappable(m,n) = Overlappable
      end do
    end do

  end subroutine ovkSetAssemblyOptionOverlappable

  subroutine ovkGetAssemblyOptionOverlapTolerance(Options, OverlappingGridID, OverlappedGridID, &
    OverlapTolerance)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    real(rk), intent(out) :: OverlapTolerance

    OverlapTolerance = Options%overlap_tolerance(OverlappingGridID,OverlappedGridID)

  end subroutine ovkGetAssemblyOptionOverlapTolerance

  subroutine ovkSetAssemblyOptionOverlapTolerance(Options, OverlappingGridID, OverlappedGridID, &
    OverlapTolerance)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    real(rk), intent(in) :: OverlapTolerance

    integer :: m, n
    integer :: ms, me, ns, ne

    if (OVK_DEBUG) then
      if (OverlapTolerance < 0._rk) then
        write (ERROR_UNIT, '(a)') "ERROR: Overlap tolerance must be nonnegative."
        stop 1
      end if
    end if

    call GridIDRange(Options%ngrids, OverlappingGridID, ms, me)
    call GridIDRange(Options%ngrids, OverlappedGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Options%overlap_tolerance(m,n) = OverlapTolerance
      end do
    end do

  end subroutine ovkSetAssemblyOptionOverlapTolerance

  subroutine ovkGetAssemblyOptionOverlapAccelQualityAdjust(Options, GridID, &
    OverlapAccelQualityAdjust)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: GridID
    real(rk), intent(out) :: OverlapAccelQualityAdjust

    OverlapAccelQualityAdjust = Options%overlap_accel_quality_adjust(GridID)

  end subroutine ovkGetAssemblyOptionOverlapAccelQualityAdjust

  subroutine ovkSetAssemblyOptionOverlapAccelQualityAdjust(Options, GridID, &
    OverlapAccelQualityAdjust)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: GridID
    real(rk), intent(in) :: OverlapAccelQualityAdjust

    integer :: m
    integer :: ms, me

    call GridIDRange(Options%ngrids, GridID, ms, me)

    do m = ms, me
      Options%overlap_accel_quality_adjust(m) = OverlapAccelQualityAdjust
    end do

  end subroutine ovkSetAssemblyOptionOverlapAccelQualityAdjust

  subroutine ovkGetAssemblyOptionInferBoundaries(Options, GridID, InferBoundaries)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: GridID
    logical, intent(out) :: InferBoundaries

    InferBoundaries = Options%infer_boundaries(GridID)

  end subroutine ovkGetAssemblyOptionInferBoundaries

  subroutine ovkSetAssemblyOptionInferBoundaries(Options, GridID, InferBoundaries)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: GridID
    logical, intent(in) :: InferBoundaries

    integer :: m
    integer :: ms, me

    call GridIDRange(Options%ngrids, GridID, ms, me)

    do m = ms, me
      Options%infer_boundaries(m) = InferBoundaries
    end do

  end subroutine ovkSetAssemblyOptionInferBoundaries

  subroutine ovkGetAssemblyOptionCutBoundaryHoles(Options, CuttingGridID, CutGridID, &
    CutBoundaryHoles)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(out) :: CutBoundaryHoles

    CutBoundaryHoles = Options%cut_boundary_holes(CuttingGridID,CutGridID)

  end subroutine ovkGetAssemblyOptionCutBoundaryHoles

  subroutine ovkSetAssemblyOptionCutBoundaryHoles(Options, CuttingGridID, CutGridID, &
    CutBoundaryHoles)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: CuttingGridID, CutGridID
    logical, intent(in) :: CutBoundaryHoles

    integer :: m, n
    integer :: ms, me, ns, ne

    call GridIDRange(Options%ngrids, CuttingGridID, ms, me)
    call GridIDRange(Options%ngrids, CutGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Options%cut_boundary_holes(m,n) = CutBoundaryHoles
      end do
    end do

  end subroutine ovkSetAssemblyOptionCutBoundaryHoles

  subroutine ovkGetAssemblyOptionOccludes(Options, OccludingGridID, OccludedGridID, Occludes)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: OccludingGridID, OccludedGridID
    integer, intent(out) :: Occludes

    Occludes = Options%occludes(OccludingGridID,OccludedGridID)

  end subroutine ovkGetAssemblyOptionOccludes

  subroutine ovkSetAssemblyOptionOccludes(Options, OccludingGridID, OccludedGridID, Occludes)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: OccludingGridID, OccludedGridID
    integer, intent(in) :: Occludes

    integer :: m, n
    integer :: ms, me, ns, ne

    call GridIDRange(Options%ngrids, OccludingGridID, ms, me)
    call GridIDRange(Options%ngrids, OccludedGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Options%occludes(m,n) = Occludes
      end do
    end do

  end subroutine ovkSetAssemblyOptionOccludes

  subroutine ovkGetAssemblyOptionEdgePadding(Options, DonorGridID, ReceiverGridID, EdgePadding)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: EdgePadding

    EdgePadding = Options%edge_padding(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblyOptionEdgePadding

  subroutine ovkSetAssemblyOptionEdgePadding(Options, DonorGridID, ReceiverGridID, EdgePadding)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(in) :: EdgePadding

    integer :: m, n
    integer :: ms, me, ns, ne

    if (OVK_DEBUG) then
      if (EdgePadding < 0) then
        write (ERROR_UNIT, '(a)') "ERROR: Edge padding must be nonnegative."
        stop 1
      end if
    end if

    call GridIDRange(Options%ngrids, DonorGridID, ms, me)
    call GridIDRange(Options%ngrids, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Options%edge_padding(m,n) = EdgePadding
      end do
    end do

  end subroutine ovkSetAssemblyOptionEdgePadding

  subroutine ovkGetAssemblyOptionEdgeSmoothing(Options, GridID, EdgeSmoothing)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: GridID
    integer, intent(out) :: EdgeSmoothing

    EdgeSmoothing = Options%edge_smoothing(GridID)

  end subroutine ovkGetAssemblyOptionEdgeSmoothing

  subroutine ovkSetAssemblyOptionEdgeSmoothing(Options, GridID, EdgeSmoothing)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: GridID
    integer, intent(in) :: EdgeSmoothing

    integer :: n
    integer :: ns, ne

    if (OVK_DEBUG) then
      if (EdgeSmoothing < 0) then
        write (ERROR_UNIT, '(a)') "ERROR: Edge smoothing must be nonnegative."
        stop 1
      end if
    end if

    call GridIDRange(Options%ngrids, GridID, ns, ne)

    do n = ns, ne
      Options%edge_smoothing(n) = EdgeSmoothing
    end do

  end subroutine ovkSetAssemblyOptionEdgeSmoothing

  subroutine ovkGetAssemblyOptionConnectionType(Options, DonorGridID, ReceiverGridID, &
    ConnectionType)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(out) :: ConnectionType

    ConnectionType = Options%connection_type(DonorGridID,ReceiverGridID)

  end subroutine ovkGetAssemblyOptionConnectionType

  subroutine ovkSetAssemblyOptionConnectionType(Options, DonorGridID, ReceiverGridID, &
    ConnectionType)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: DonorGridID, ReceiverGridID
    integer, intent(in) :: ConnectionType

    logical :: ValidValue
    integer :: m, n
    integer :: ms, me, ns, ne

    if (OVK_DEBUG) then
      ValidValue = &
        ConnectionType == OVK_CONNECTION_NONE .or. &
        ConnectionType == OVK_CONNECTION_LINEAR .or. &
        ConnectionType == OVK_CONNECTION_CUBIC
      if (.not. ValidValue) then
        write (ERROR_UNIT, '(a)') "ERROR: Invalid connection type."
        stop 1
      end if
    end if

    call GridIDRange(Options%ngrids, DonorGridID, ms, me)
    call GridIDRange(Options%ngrids, ReceiverGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Options%connection_type(m,n) = ConnectionType
      end do
    end do

  end subroutine ovkSetAssemblyOptionConnectionType

  subroutine ovkGetAssemblyOptionFringeSize(Options, GridID, FringeSize)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: GridID
    integer, intent(out) :: FringeSize

    FringeSize = Options%fringe_size(GridID)

  end subroutine ovkGetAssemblyOptionFringeSize

  subroutine ovkSetAssemblyOptionFringeSize(Options, GridID, FringeSize)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: GridID
    integer, intent(in) :: FringeSize

    integer :: n
    integer :: ns, ne

    if (OVK_DEBUG) then
      if (FringeSize < 0) then
        write (ERROR_UNIT, '(a)') "ERROR: Fringe size must be nonnegative."
        stop 1
      end if
    end if

    call GridIDRange(Options%ngrids, GridID, ns, ne)

    do n = ns, ne
      Options%fringe_size(n) = FringeSize
    end do

  end subroutine ovkSetAssemblyOptionFringeSize

  subroutine ovkGetAssemblyOptionOverlapMinimization(Options, OverlappingGridID, &
    OverlappedGridID, OverlapMinimization)

    type(ovk_assembly_options), intent(in) :: Options
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(out) :: OverlapMinimization

    OverlapMinimization = Options%overlap_minimization(OverlappingGridID,OverlappedGridID)

  end subroutine ovkGetAssemblyOptionOverlapMinimization

  subroutine ovkSetAssemblyOptionOverlapMinimization(Options, OverlappingGridID, &
    OverlappedGridID, OverlapMinimization)

    type(ovk_assembly_options), intent(inout) :: Options
    integer, intent(in) :: OverlappingGridID, OverlappedGridID
    logical, intent(in) :: OverlapMinimization

    integer :: m, n
    integer :: ms, me, ns, ne

    call GridIDRange(Options%ngrids, OverlappingGridID, ms, me)
    call GridIDRange(Options%ngrids, OverlappedGridID, ns, ne)

    do n = ns, ne
      do m = ms, me
        Options%overlap_minimization(m,n) = OverlapMinimization
      end do
    end do

  end subroutine ovkSetAssemblyOptionOverlapMinimization

  subroutine GridIDRange(NumGrids, GridID, StartID, EndID)

    integer, intent(in) :: NumGrids
    integer, intent(in) :: GridID
    integer, intent(out) :: StartID, EndID

    if (GridID == OVK_ALL_GRIDS) then
      StartID = 1
      EndID = NumGrids
    else
      StartID = GridID
      EndID = GridID
    end if

  end subroutine GridIDRange

end module ovkAssemblyOptions
