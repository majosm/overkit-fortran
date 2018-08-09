! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module TestFixtures

  use Overkit
  use ovsGlobal
  use ovkGrid, only : CreateGrid, DestroyGrid
  use ovkLogger, only : t_logger_
  use ovkOverlap, only : CreateOverlap, DestroyOverlap
  implicit none

  private

  public :: SetupStaggered2D, TeardownStaggered2D
  public :: SetupStaggered3D, TeardownStaggered3D
  public :: SetupBlockInterface2D, TeardownBlockInterface2D
  public :: SetupBlockInterface3D, TeardownBlockInterface3D

contains

  subroutine SetupStaggered2D(Grid1, Grid2, Overlap)

    type(ovk_grid), target, intent(out) :: Grid1, Grid2
    type(ovk_overlap), intent(out), optional :: Overlap

    integer :: i, j, l
    type(ovk_cart) :: Grid1Cart, Grid2Cart
    type(ovk_field_real), pointer :: X, Y
    type(ovk_field_logical) :: OverlappedMask

    Grid1Cart = ovk_cart_(2, [5,6], [.false.,.true.], OVK_NO_OVERLAP_PERIODIC)
    call CreateGrid(Grid1, t_logger_(), 1, Grid1Cart, PeriodicLength=[0._rk, 6._rk], &
      GeometryType=OVK_GEOMETRY_CARTESIAN)

    call ovkEditGridCoords(Grid1, 1, X)
    call ovkEditGridCoords(Grid1, 2, Y)
    do j = 1, 6
      do i = 1, 5
        X%values(i,j,1) = 0.5_rk+X%values(i,j,1)
        Y%values(i,j,1) = -0.5_rk+Y%values(i,j,1)
      end do
    end do
    call ovkReleaseGridCoords(Grid1, X)
    call ovkReleaseGridCoords(Grid1, Y)

    Grid2Cart = ovk_cart_(2, [6,6])
    call CreateGrid(Grid2, t_logger_(), 2, Grid2Cart, GeometryType=OVK_GEOMETRY_CARTESIAN)

    if (present(Overlap)) then

      call CreateOverlap(Overlap, t_logger_(), Grid1, Grid2)

      OverlappedMask = ovk_field_logical_(Grid2Cart, .false.)
      OverlappedMask%values(2:5,:,1) = .true.

      call ovkResetOverlap(Overlap, OverlappedMask)

      l = 1
      do j = 1, 6
        do i = 1, 4
          Overlap%cells(:2,l) = [i,j]
          Overlap%coords(:,l) = 0.5_rk
          l = l + 1
        end do
      end do

    end if

  end subroutine SetupStaggered2D

  subroutine TeardownStaggered2D(Grid1, Grid2, Overlap)

    type(ovk_grid), target, intent(inout) :: Grid1, Grid2
    type(ovk_overlap), intent(inout), optional :: Overlap

    if (present(Overlap)) then
      call DestroyOverlap(Overlap)
    end if

    call DestroyGrid(Grid1)
    call DestroyGrid(Grid2)

  end subroutine TeardownStaggered2D

  subroutine SetupStaggered3D(Grid1, Grid2, Overlap)

    type(ovk_grid), target, intent(out) :: Grid1, Grid2
    type(ovk_overlap), intent(out), optional :: Overlap

    integer :: i, j, k, l
    type(ovk_cart) :: Grid1Cart, Grid2Cart
    type(ovk_field_real), pointer :: X, Y, Z
    type(ovk_field_logical) :: OverlappedMask

    Grid1Cart = ovk_cart_(3, [5,5,6], [.false.,.false.,.true.], OVK_NO_OVERLAP_PERIODIC)
    call CreateGrid(Grid1, t_logger_(), 1, Grid1Cart, PeriodicLength=[0._rk, 0._rk, 6._rk], &
      GeometryType=OVK_GEOMETRY_CARTESIAN)

    call ovkEditGridCoords(Grid1, 1, X)
    call ovkEditGridCoords(Grid1, 2, Y)
    call ovkEditGridCoords(Grid1, 3, Z)
    do k = 1, 6
      do j = 1, 5
        do i = 1, 5
          X%values(i,j,k) = 0.5_rk+X%values(i,j,k)
          Y%values(i,j,k) = 0.5_rk+Y%values(i,j,k)
          Z%values(i,j,k) = -0.5_rk+Z%values(i,j,k)
        end do
      end do
    end do
    call ovkReleaseGridCoords(Grid1, X)
    call ovkReleaseGridCoords(Grid1, Y)
    call ovkReleaseGridCoords(Grid1, Z)

    Grid2Cart = ovk_cart_(3, [6,6,6])
    call CreateGrid(Grid2, t_logger_(), 2, Grid2Cart, GeometryType=OVK_GEOMETRY_CARTESIAN)

    if (present(Overlap)) then

      call CreateOverlap(Overlap, t_logger_(), Grid1, Grid2)

      OverlappedMask = ovk_field_logical_(Grid2Cart, .false.)
      OverlappedMask%values(2:5,2:5,:) = .true.

      call ovkResetOverlap(Overlap, OverlappedMask)

      l = 1
      do k = 1, 6
        do j = 1, 4
          do i = 1, 4
            Overlap%cells(:,l) = [i,j,k]
            Overlap%coords(:,l) = 0.5_rk
            l = l + 1
          end do
        end do
      end do

    end if

  end subroutine SetupStaggered3D

  subroutine TeardownStaggered3D(Grid1, Grid2, Overlap)

    type(ovk_grid), target, intent(inout) :: Grid1, Grid2
    type(ovk_overlap), intent(inout), optional :: Overlap

    if (present(Overlap)) then
      call DestroyOverlap(Overlap)
    end if

    call DestroyGrid(Grid1)
    call DestroyGrid(Grid2)

  end subroutine TeardownStaggered3D

  subroutine SetupBlockInterface2D(Grid1, Grid2, Overlap)

    type(ovk_grid), target, intent(out) :: Grid1, Grid2
    type(ovk_overlap), intent(out), optional :: Overlap

    integer :: i, j, l
    type(ovk_cart) :: Grid1Cart, Grid2Cart
    type(ovk_field_real), pointer :: Y
    type(ovk_field_logical) :: OverlappedMask

    Grid1Cart = ovk_cart_(2, [7,7], [.true.,.false.], OVK_OVERLAP_PERIODIC)
    call CreateGrid(Grid1, t_logger_(), 1, Grid1Cart, PeriodicLength=[6._rk, 0._rk], &
      GeometryType=OVK_GEOMETRY_CARTESIAN)

    Grid2Cart = ovk_cart_(2, [7,7], [.true.,.false.], OVK_OVERLAP_PERIODIC)
    call CreateGrid(Grid2, t_logger_(), 2, Grid2Cart, PeriodicLength=[6._rk, 0._rk], &
      GeometryType=OVK_GEOMETRY_CARTESIAN)

    call ovkEditGridCoords(Grid2, 2, Y)
    do j = 1, 7
      do i = 1, 7
        Y%values(i,j,1) = Y%values(i,j,1)+6._rk
      end do
    end do
    call ovkReleaseGridCoords(Grid2, Y)

    if (present(Overlap)) then

      call CreateOverlap(Overlap, t_logger_(), Grid1, Grid2)

      OverlappedMask = ovk_field_logical_(Grid2Cart, .false.)
      OverlappedMask%values(:,1,1) = .true.

      call ovkResetOverlap(Overlap, OverlappedMask)

      l = 1
      do i = 1, 7
        Overlap%cells(:2,l) = [i,6]
        Overlap%coords(1,l) = merge(0._rk, 1._rk, i < 7)
        Overlap%coords(2,l) = 1._rk
        l = l + 1
      end do

    end if

  end subroutine SetupBlockInterface2D

  subroutine TeardownBlockInterface2D(Grid1, Grid2, Overlap)

    type(ovk_grid), target, intent(inout) :: Grid1, Grid2
    type(ovk_overlap), intent(inout), optional :: Overlap

    if (present(Overlap)) then
      call DestroyOverlap(Overlap)
    end if

    call DestroyGrid(Grid1)
    call DestroyGrid(Grid2)

  end subroutine TeardownBlockInterface2D

  subroutine SetupBlockInterface3D(Grid1, Grid2, Overlap)

    type(ovk_grid), target, intent(out) :: Grid1, Grid2
    type(ovk_overlap), intent(out), optional :: Overlap

    integer :: i, j, k, l
    type(ovk_cart) :: Grid1Cart, Grid2Cart
    type(ovk_field_real), pointer :: Z
    type(ovk_field_logical) :: OverlappedMask

    Grid1Cart = ovk_cart_(3, [7,7,7], [.true.,.true.,.false.], OVK_OVERLAP_PERIODIC)
    call CreateGrid(Grid1, t_logger_(), 1, Grid1Cart, PeriodicLength=[6._rk, 6._rk, 0._rk], &
      GeometryType=OVK_GEOMETRY_CARTESIAN)

    Grid2Cart = ovk_cart_(3, [7,7,7], [.true.,.true.,.false.], OVK_OVERLAP_PERIODIC)
    call CreateGrid(Grid2, t_logger_(), 2, Grid2Cart, PeriodicLength=[6._rk, 6._rk, 0._rk], &
      GeometryType=OVK_GEOMETRY_CARTESIAN)

    call ovkEditGridCoords(Grid2, 3, Z)
    do k = 1, 7
      do j = 1, 7
        do i = 1, 7
          Z%values(i,j,k) = Z%values(i,j,k)+6._rk
        end do
      end do
    end do
    call ovkReleaseGridCoords(Grid2, Z)

    if (present(Overlap)) then

      call CreateOverlap(Overlap, t_logger_(), Grid1, Grid2)

      OverlappedMask = ovk_field_logical_(Grid2Cart, .false.)
      OverlappedMask%values(:,:,1) = .true.

      call ovkResetOverlap(Overlap, OverlappedMask)

      l = 1
      do j = 1, 7
        do i = 1, 7
          Overlap%cells(:,l) = [i,j,6]
          Overlap%coords(1,l) = merge(0._rk, 1._rk, i < 7)
          Overlap%coords(2,l) = merge(0._rk, 1._rk, j < 7)
          Overlap%coords(3,l) = 1._rk
          l = l + 1
        end do
      end do

    end if

  end subroutine SetupBlockInterface3D

  subroutine TeardownBlockInterface3D(Grid1, Grid2, Overlap)

    type(ovk_grid), target, intent(inout) :: Grid1, Grid2
    type(ovk_overlap), intent(inout), optional :: Overlap

    if (present(Overlap)) then
      call DestroyOverlap(Overlap)
    end if

    call DestroyGrid(Grid1)
    call DestroyGrid(Grid2)

  end subroutine TeardownBlockInterface3D

end module TestFixtures
