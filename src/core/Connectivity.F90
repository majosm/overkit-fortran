! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkConnectivity

  use ovkGlobal
  use ovkGrid
  use ovkInterp
  implicit none

  private

  ! API
  public :: ovk_connectivity
  public :: ovk_connectivity_
  public :: ovk_connectivity_properties
  public :: ovk_connectivity_properties_
  public :: ovkCreateConnectivity
  public :: ovkDestroyConnectivity
  public :: ovkUpdateConnectivity
  public :: ovkGetConnectivityProperties
  public :: ovkEditConnectivityProperties
  public :: ovkReleaseConnectivityProperties
  public :: ovkCreateConnectivityInterpData
  public :: ovkDestroyConnectivityInterpData
  public :: ovkGetConnectivityInterpData
!   public :: ovkCreateConnectivityDonors
!   public :: ovkDestroyConnectivityDonors
!   public :: ovkResetConnectivityDonors
!   public :: ovkGetConnectivityDonors
!   public :: ovkEditConnectivityDonors
!   public :: ovkReleaseConnectivityDonors
!   public :: ovkCreateConnectivityReceivers
!   public :: ovkDestroyConnectivityReceivers
!   public :: ovkResetConnectivityReceivers
!   public :: ovkGetConnectivityReceivers
!   public :: ovkEditConnectivityReceivers
!   public :: ovkReleaseConnectivityReceivers
  public :: ovkGetConnectivityPropertyDimension
  public :: ovkGetConnectivityPropertyGridCount
  public :: ovkGetConnectivityPropertyVerbose
  public :: ovkSetConnectivityPropertyVerbose

  type ovk_connectivity_properties
    type(t_noconstruct) :: noconstruct
    ! Read-only
    integer :: nd
    integer :: ngrids
    ! Read/write
    logical :: verbose
  end type ovk_connectivity_properties

  type ovk_connectivity
    type(t_noconstruct) :: noconstruct
    type(ovk_connectivity_properties), pointer :: properties
    type(ovk_interp), dimension(:), pointer :: interp_data
!     type(ovk_donors), dimension(:), pointer :: donors
!     type(ovk_receivers), dimension(:), pointer :: receivers
    logical, dimension(:), allocatable :: interp_data_exists
!     logical, dimension(:), allocatable :: donors_exist
!     logical, dimension(:), allocatable :: receivers_exist
    logical :: editing_properties
!     logical, dimension(:), allocatable :: editing_donors
!     logical, dimension(:), allocatable :: editing_receivers
!     logical, dimension(:), allocatable :: changed_donors
!     logical, dimension(:), allocatable :: changed_receivers
  end type ovk_connectivity

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_connectivity_
    module procedure ovk_connectivity_Default
  end interface ovk_connectivity_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_connectivity_properties_
    module procedure ovk_connectivity_properties_Default
  end interface ovk_connectivity_properties_

contains

  function ovk_connectivity_Default() result(Connectivity)

    type(ovk_connectivity) :: Connectivity

    nullify(Connectivity%properties)
    nullify(Connectivity%interp_data)
!     nullify(Connectivity%donors)
!     nullify(Connectivity%receivers)
    Connectivity%editing_properties = .false.

  end function ovk_connectivity_Default

  subroutine ovkCreateConnectivity(Connectivity, NumDims, NumGrids, Verbose)

    type(ovk_connectivity), intent(out) :: Connectivity
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    logical, intent(in), optional :: Verbose

    logical :: Verbose_
    integer :: m

    if (present(Verbose)) then
      Verbose_ = Verbose
    else
      Verbose_ = .false.
    end if

    allocate(Connectivity%properties)
    Connectivity%properties = ovk_connectivity_properties_()

    Connectivity%properties%nd = NumDims
    Connectivity%properties%ngrids = NumGrids
    Connectivity%properties%verbose = Verbose_

    allocate(Connectivity%interp_data(NumGrids))
    do m = 1, NumGrids
      Connectivity%interp_data(m) = ovk_interp_()
    end do

    allocate(Connectivity%interp_data_exists(NumGrids))
    Connectivity%interp_data_exists = .false.

!     allocate(Connectivity%donors_exist(NumGrids))
!     Connectivity%donors_exist = .false.

!     allocate(Connectivity%receivers_exist(NumGrids))
!     Connectivity%receivers_exist = .false.

    Connectivity%editing_properties = .false.

!     allocate(Connectivity%editing_donors(NumGrids))
!     Connectivity%editing_donors = .false.

!     allocate(Connectivity%editing_receivers(NumGrids))
!     Connectivity%editing_receivers = .false.

!     allocate(Connectivity%changed_donors(NumGrids))
!     Connectivity%changed_donors = .false.

!     allocate(Connectivity%changed_receivers(NumGrids))
!     Connectivity%changed_receivers = .false.

  end subroutine ovkCreateConnectivity

  subroutine ovkDestroyConnectivity(Connectivity)

    type(ovk_connectivity), intent(inout) :: Connectivity

    integer :: m

    if (associated(Connectivity%properties)) deallocate(Connectivity%properties)

    if (associated(Connectivity%interp_data)) then
      do m = 1, size(Connectivity%interp_data)
        if (Connectivity%interp_data_exists(m)) then
          call ovkDestroyInterpData(Connectivity%interp_data(m))
        end if
      end do
      deallocate(Connectivity%interp_data)
    end if

!     if (associated(Connectivity%donors)) then
!       do m = 1, size(Connectivity%donors)
!         if (Connectivity%donors_exist(m)) then
!           call ovkDestroyDonors(Connectivity%donors(m))
!         end if
!       end do
!       deallocate(Connectivity%donors)
!     end if

!     if (associated(Connectivity%receivers)) then
!       do m = 1, size(Connectivity%receivers)
!         if (Connectivity%receivers_exist(m)) then
!           call ovkDestroyDonors(Connectivity%receivers(m))
!         end if
!       end do
!       deallocate(Connectivity%receivers)
!     end if

!     if (allocated(Connectivity%editing_donors)) deallocate(Connectivity%editing_donors)
!     if (allocated(Connectivity%editing_receivers)) deallocate(Connectivity%editing_receivers)
!     if (allocated(Connectivity%changed_donors)) deallocate(Connectivity%changed_donors)
!     if (allocated(Connectivity%changed_receivers)) deallocate(Connectivity%changed_receivers)

  end subroutine ovkDestroyConnectivity

  subroutine ovkUpdateConnectivity(Connectivity)

    type(ovk_connectivity), intent(inout) :: Connectivity

    logical :: CannotUpdate

    CannotUpdate = &
      Connectivity%editing_properties

!     CannotUpdate = &
!       Connectivity%editing_properties .or. &
!       any(Connectivity%editing_donors) .or. &
!       any(Connectivity%editing_receivers)

    if (CannotUpdate) then

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Cannot update connectivity; still being edited."
        stop 1
      end if

      return

    else

!       Connectivity%changed_donors = .false.
!       Connectivity%changed_receivers = .false.

    end if

  end subroutine ovkUpdateConnectivity

  subroutine ovkGetConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(in) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(out) :: Properties

    Properties => Connectivity%properties

  end subroutine ovkGetConnectivityProperties

  subroutine ovkEditConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(inout) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(out) :: Properties

!     logical :: CannotEdit

!     CannotEdit = &
!       any(Connectivity%editing_donors) .or. &
!       any(Connectivity%editing_receivers)

!     if (CannotEdit) then

!       if (OVK_DEBUG) then
!         write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity properties while editing grids."
!         stop 1
!       end if

!       nullify(Properties)

!     else

      Properties => Connectivity%properties
      Connectivity%editing_properties = .true.

!     end if

  end subroutine ovkEditConnectivityProperties

  subroutine ovkReleaseConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(inout) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(inout) :: Properties

!     integer :: m
!     integer :: NumGrids
!     type(ovk_donors_properties), pointer :: DonorProperties
!     type(ovk_receivers_properties), pointer :: ReceiverProperties

    if (.not. associated(Properties, Connectivity%properties)) return
    if (.not. Connectivity%editing_properties) then
      nullify(Properties)
      return
    end if

!     NumGrids = Connectivity%properties%ngrids

    ! Propagate verbose flag to donors and receivers
!     do m = 1, NumGrids
!       if (Connectivity%donors_exist(m)) then
!         call ovkEditDonorProperties(Connectivity%donors(m), DonorProperties)
!         call ovkSetDonorPropertyVerbose(DonorProperties, Connectivity%properties%verbose)
!         call ovkReleaseDonorProperties(Connectivity%donors(m), DonorProperties)
!       end if
!       if (Connectivity%receivers_exist(m)) then
!         call ovkEditReceiverProperties(Connectivity%receivers(m), ReceiverProperties)
!         call ovkSetReceiverPropertyVerbose(ReceiverProperties, Connectivity%properties%verbose)
!         call ovkReleaseReceiverProperties(Connectivity%receivers(m), ReceiverProperties)
!       end if
!     end do

    nullify(Properties)
    Connectivity%editing_properties = .false.

  end subroutine ovkReleaseConnectivityProperties

  subroutine ovkCreateConnectivityInterpData(Connectivity, GridID, Grid, MaxStencilSize)

    type(ovk_connectivity), intent(inout) :: Connectivity
    integer, intent(in) :: GridID
    type(ovk_grid), intent(in) :: Grid
    integer, intent(in) :: MaxStencilSize

    call ovkCreateInterpData(Connectivity%interp_data(GridID), Grid, MaxStencilSize, &
      Verbose=Connectivity%properties%verbose)

    Connectivity%interp_data_exists(GridID) = .true.
!     Connectivity%changed_interp_data(GridID) = .true.

  end subroutine ovkCreateConnectivityInterpData

  subroutine ovkDestroyConnectivityInterpData(Connectivity, GridID)

    type(ovk_connectivity), intent(inout) :: Connectivity
    integer, intent(in) :: GridID

    call ovkDestroyInterpData(Connectivity%interp_data(GridID))

    Connectivity%interp_data_exists(GridID) = .false.
!     Connectivity%changed_interp_data(GridID) = .true.

  end subroutine ovkDestroyConnectivityInterpData

  subroutine ovkGetConnectivityInterpData(Connectivity, GridID, InterpData)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer, intent(in) :: GridID
    type(ovk_interp), pointer, intent(out) :: InterpData

    if (.not. Connectivity%interp_data_exists(GridID)) then

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Interpolation data does not exist."
        stop 1
      end if

      nullify(InterpData)

    else

      InterpData => Connectivity%interp_data(GridID)

    end if

  end subroutine ovkGetConnectivityInterpData

  function ovk_connectivity_properties_Default() result(Properties)

    type(ovk_connectivity_properties) :: Properties

    Properties%nd = 2
    Properties%ngrids = 0
    Properties%verbose = .false.

  end function ovk_connectivity_properties_Default

  subroutine ovkGetConnectivityPropertyDimension(Properties, NumDims)

    type(ovk_connectivity_properties), intent(in) :: Properties
    integer, intent(out) :: NumDims

    NumDims = Properties%nd

  end subroutine ovkGetConnectivityPropertyDimension

  subroutine ovkGetConnectivityPropertyGridCount(Properties, NumGrids)

    type(ovk_connectivity_properties), intent(in) :: Properties
    integer, intent(out) :: NumGrids

    NumGrids = Properties%ngrids

  end subroutine ovkGetConnectivityPropertyGridCount

  subroutine ovkGetConnectivityPropertyVerbose(Properties, Verbose)

    type(ovk_connectivity_properties), intent(in) :: Properties
    logical, intent(out) :: Verbose

    Verbose = Properties%verbose

  end subroutine ovkGetConnectivityPropertyVerbose

  subroutine ovkSetConnectivityPropertyVerbose(Properties, Verbose)

    type(ovk_connectivity_properties), intent(inout) :: Properties
    logical, intent(in) :: Verbose

    Properties%verbose = Verbose

  end subroutine ovkSetConnectivityPropertyVerbose

end module ovkConnectivity
