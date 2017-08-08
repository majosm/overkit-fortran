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

  type t_connectivity_editor
    integer :: properties_ref_count
!     integer, dimension(:), allocatable :: donors_ref_count
!     integer, dimension(:), allocatable :: receivers_ref_count
  end type t_connectivity_editor

  type ovk_connectivity
    type(t_noconstruct) :: noconstruct
    type(ovk_connectivity_properties), pointer :: properties
    type(ovk_connectivity_properties) :: prev_properties
    type(ovk_interp), dimension(:), pointer :: interp_data
!     type(ovk_donors), dimension(:), pointer :: donors
!     type(ovk_receivers), dimension(:), pointer :: receivers
    type(t_connectivity_editor) :: editor
  end type ovk_connectivity

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_connectivity_
    module procedure ovk_connectivity_Default
  end interface ovk_connectivity_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface ovk_connectivity_properties_
    module procedure ovk_connectivity_properties_Default
    module procedure ovk_connectivity_properties_Assigned
  end interface ovk_connectivity_properties_

  ! Trailing _ added for compatibility with compilers that don't support F2003 constructors
  interface t_connectivity_editor_
    module procedure t_connectivity_editor_Default
    module procedure t_connectivity_editor_Allocated
  end interface t_connectivity_editor_

contains

  function ovk_connectivity_Default() result(Connectivity)

    type(ovk_connectivity) :: Connectivity

    nullify(Connectivity%properties)
    Connectivity%prev_properties = ovk_connectivity_properties_()
    nullify(Connectivity%interp_data)
!     nullify(Connectivity%donors)
!     nullify(Connectivity%receivers)
    Connectivity%editor = t_connectivity_editor_()

  end function ovk_connectivity_Default

  subroutine ovkCreateConnectivity(Connectivity, NumDims, NumGrids, Verbose)

    type(ovk_connectivity), intent(out) :: Connectivity
    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    logical, intent(in), optional :: Verbose

    logical :: Verbose_
    integer :: m
    type(ovk_grid_description) :: EmptyGridDescription

    if (present(Verbose)) then
      Verbose_ = Verbose
    else
      Verbose_ = .false.
    end if

    allocate(Connectivity%properties)
    Connectivity%properties = ovk_connectivity_properties_(NumDims, NumGrids)
    Connectivity%properties%verbose = Verbose_

    EmptyGridDescription = ovk_grid_description_(NumDims)

    allocate(Connectivity%interp_data(NumGrids))
    do m = 1, NumGrids
      call ovkCreateInterpData(Connectivity%interp_data(m), m, EmptyGridDescription, 1, &
        Verbose=Verbose_)
    end do

    Connectivity%editor = t_connectivity_editor_(NumGrids)

  end subroutine ovkCreateConnectivity

  subroutine ovkDestroyConnectivity(Connectivity)

    type(ovk_connectivity), intent(inout) :: Connectivity

    integer :: m

    if (associated(Connectivity%properties)) deallocate(Connectivity%properties)

    if (associated(Connectivity%interp_data)) then
      do m = 1, size(Connectivity%interp_data)
        call ovkDestroyInterpData(Connectivity%interp_data(m))
      end do
      deallocate(Connectivity%interp_data)
    end if

!     if (associated(Connectivity%donors)) then
!       do m = 1, size(Connectivity%donors)
!         call ovkDestroyDonors(Connectivity%donors(m))
!       end do
!       deallocate(Connectivity%donors)
!     end if

!     if (associated(Connectivity%receivers)) then
!       do m = 1, size(Connectivity%receivers)
!         call ovkDestroyDonors(Connectivity%receivers(m))
!       end do
!       deallocate(Connectivity%receivers)
!     end if

    Connectivity%editor = t_connectivity_editor_()

  end subroutine ovkDestroyConnectivity

  subroutine ovkGetConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(in) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(out) :: Properties

    Properties => Connectivity%properties

  end subroutine ovkGetConnectivityProperties

  subroutine ovkEditConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(inout) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(out) :: Properties

!     integer :: m
!     logical :: EditingDonors
!     logical :: EditingReceivers
!     logical :: Success

!     EditingDonors = .false.
!     EditingReceivers = .false.
!     do m = 1, Connectivity%properties%ngrids
!       EditingDonors = EditingDonors .or. Connectivity%editor%donors_ref_count(m) > 0
!       EditingReceivers = EditingReceivers .or. Connectivity%editor%receivers_ref_count(m) > 0
!     end do

!     Success = &
!       .not. EditingDonors .and. &
!       .not. EditingReceivers

!     if (Success) then

      if (Connectivity%editor%properties_ref_count == 0) then
        Connectivity%prev_properties = Connectivity%properties
      end if

      Connectivity%editor%properties_ref_count = Connectivity%editor%properties_ref_count + 1

      Properties => Connectivity%properties

!     else

!       if (OVK_DEBUG) then
!         if (EditingDonors) then
!           write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity properties while editing donors."
!         end if
!         if (EditingReceivers) then
!           write (ERROR_UNIT, '(a)') "ERROR: Cannot edit connectivity properties while editing receivers."
!         end if
!         stop 1
!       end if

!       nullify(Properties)

!     end if

  end subroutine ovkEditConnectivityProperties

  subroutine ovkReleaseConnectivityProperties(Connectivity, Properties)

    type(ovk_connectivity), intent(inout) :: Connectivity
    type(ovk_connectivity_properties), pointer, intent(inout) :: Properties

    if (associated(Properties, Connectivity%properties)) then

      nullify(Properties)

      if (Connectivity%editor%properties_ref_count > 0) then

        Connectivity%editor%properties_ref_count = Connectivity%editor%properties_ref_count - 1

        if (Connectivity%editor%properties_ref_count == 0) then
          if (Connectivity%properties%verbose .neqv. Connectivity%prev_properties%verbose) then
            call UpdateVerbose(Connectivity)
          end if
        end if

      end if

    else

      if (OVK_DEBUG) then
        write (ERROR_UNIT, '(a)') "ERROR: Unable to release connectivity properties; invalid pointer."
        stop 1
      end if

    end if

  end subroutine ovkReleaseConnectivityProperties

  subroutine UpdateVerbose(Connectivity)

    type(ovk_connectivity), intent(inout) :: Connectivity

    integer :: m
    type(ovk_interp_properties), pointer :: InterpDataProperties

    do m = 1, Connectivity%properties%ngrids
      call ovkEditInterpDataProperties(Connectivity%interp_data(m), InterpDataProperties)
      call ovkSetInterpDataPropertyVerbose(InterpDataProperties, Connectivity%properties%verbose)
      call ovkReleaseInterpDataProperties(Connectivity%interp_data(m), InterpDataProperties)
!       call ovkEditDonorProperties(Connectivity%donors(m), DonorProperties)
!       call ovkSetDonorPropertyVerbose(DonorProperties, Connectivity%properties%verbose)
!       call ovkReleaseDonorProperties(Connectivity%donors(m), DonorProperties)
!       call ovkEditReceiverProperties(Connectivity%receivers(m), ReceiverProperties)
!       call ovkSetReceiverPropertyVerbose(ReceiverProperties, Connectivity%properties%verbose)
!       call ovkReleaseReceiverProperties(Connectivity%receivers(m), ReceiverProperties)
    end do

  end subroutine UpdateVerbose

  subroutine ovkCreateConnectivityInterpData(Connectivity, GridID, GridDescription, StencilSize)

    type(ovk_connectivity), intent(inout) :: Connectivity
    integer, intent(in) :: GridID
    type(ovk_grid_description), intent(in) :: GridDescription
    integer, intent(in) :: StencilSize

    call ovkDestroyInterpData(Connectivity%interp_data(GridID))

    call ovkCreateInterpData(Connectivity%interp_data(GridID), GridID, GridDescription, &
      StencilSize, Verbose=Connectivity%properties%verbose)

  end subroutine ovkCreateConnectivityInterpData

  subroutine ovkDestroyConnectivityInterpData(Connectivity, GridID)

    type(ovk_connectivity), intent(inout) :: Connectivity
    integer, intent(in) :: GridID

    type(ovk_grid_description) :: EmptyGridDescription

    call ovkDestroyInterpData(Connectivity%interp_data(GridID))

    EmptyGridDescription = ovk_grid_description_(Connectivity%properties%nd)

    call ovkCreateInterpData(Connectivity%interp_data(GridID), GridID, EmptyGridDescription, 1, &
      Verbose=Connectivity%properties%verbose)

  end subroutine ovkDestroyConnectivityInterpData

  subroutine ovkGetConnectivityInterpData(Connectivity, GridID, InterpData)

    type(ovk_connectivity), intent(in) :: Connectivity
    integer, intent(in) :: GridID
    type(ovk_interp), pointer, intent(out) :: InterpData

    InterpData => Connectivity%interp_data(GridID)

  end subroutine ovkGetConnectivityInterpData

  function ovk_connectivity_properties_Default() result(Properties)

    type(ovk_connectivity_properties) :: Properties

    Properties%nd = 2
    Properties%ngrids = 0
    Properties%verbose = .false.

  end function ovk_connectivity_properties_Default

  function ovk_connectivity_properties_Assigned(NumDims, NumGrids) result(Properties)

    integer, intent(in) :: NumDims
    integer, intent(in) :: NumGrids
    type(ovk_connectivity_properties) :: Properties

    Properties%nd = NumDims
    Properties%ngrids = NumGrids
    Properties%verbose = .false.

  end function ovk_connectivity_properties_Assigned

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

  function t_connectivity_editor_Default() result(Editor)

    type(t_connectivity_editor) :: Editor

    Editor%properties_ref_count = 0

  end function t_connectivity_editor_Default

  function t_connectivity_editor_Allocated(NumGrids) result(Editor)

    integer, intent(in) :: NumGrids
    type(t_connectivity_editor) :: Editor

    Editor%properties_ref_count = 0

!     allocate(Editor%donors_ref_count(NumGrids))
!     Editor%donors_ref_count = 0
!     allocate(Editor%receivers_ref_count(NumGrids))
!     Editor%receivers_ref_count = 0

  end function t_connectivity_editor_Allocated

end module ovkConnectivity
