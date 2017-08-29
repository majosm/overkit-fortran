! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module ovkGeometry

  use ovkGlobal
#ifdef f2003
  use, intrinsic ieee_arithmetic, only : ieee_is_nan
#endif
  implicit none

  private

  ! API
  public :: ovkOverlapsRectangle
  public :: ovkOverlapsCuboid
  public :: ovkOverlapsQuad
  public :: ovkOverlapsHexahedron
  public :: ovkRectangleSize
  public :: ovkCuboidSize
  public :: ovkQuadSize
  public :: ovkHexahedronSize
  public :: ovkRectangleIsoLinear
  public :: ovkRectangleIsoCubic
  public :: ovkCuboidIsoLinear
  public :: ovkCuboidIsoCubic
  public :: ovkQuadIsoLinear
  public :: ovkQuadIsoCubic
  public :: ovkHexahedronIsoLinear
  public :: ovkHexahedronIsoCubic
  public :: ovkRectangleIsoInverseLinear
  public :: ovkRectangleIsoInverseCubic
  public :: ovkCuboidIsoInverseLinear
  public :: ovkCuboidIsoInverseCubic
  public :: ovkQuadIsoInverseLinear
  public :: ovkQuadIsoInverseCubic
  public :: ovkHexahedronIsoInverseLinear
  public :: ovkHexahedronIsoInverseCubic
  public :: ovkCartesianGridCell
  public :: ovkInterpBasisLinear
  public :: ovkInterpBasisLinearDeriv
  public :: ovkInterpBasisCubic
  public :: ovkInterpBasisCubicDeriv

  interface ovkQuadIsoLinear
    module procedure ovkQuadIsoLinear_LocalCoords
    module procedure ovkQuadIsoLinear_BasisValues
  end interface ovkQuadIsoLinear

  interface ovkQuadIsoCubic
    module procedure ovkQuadIsoCubic_LocalCoords
    module procedure ovkQuadIsoCubic_BasisValues
  end interface ovkQuadIsoCubic

  interface ovkHexahedronIsoLinear
    module procedure ovkHexahedronIsoLinear_LocalCoords
    module procedure ovkHexahedronIsoLinear_BasisValues
  end interface ovkHexahedronIsoLinear

  interface ovkHexahedronIsoCubic
    module procedure ovkHexahedronIsoCubic_LocalCoords
    module procedure ovkHexahedronIsoCubic_BasisValues
  end interface ovkHexahedronIsoCubic

  interface ovkCartesianGridCell
    module procedure ovkCartesianGridCell_Scalar
    module procedure ovkCartesianGridCell_Array
  end interface ovkCartesianGridCell

contains

  pure function ovkOverlapsRectangle(VertexCoords, Coords) result(Overlaps)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: Coords
    logical :: Overlaps

    Overlaps = .false.

    if (any(Coords < VertexCoords(:,1))) then
      return
    end if

    if (any(Coords > VertexCoords(:,4))) then
      return
    end if

    Overlaps = .true.

  end function ovkOverlapsRectangle

  pure function ovkOverlapsCuboid(VertexCoords, Coords) result(Overlaps)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: Coords
    logical :: Overlaps

    Overlaps = .false.

    if (any(Coords < VertexCoords(:,1))) then
      return
    end if

    if (any(Coords > VertexCoords(:,8))) then
      return
    end if

    Overlaps = .true.

  end function ovkOverlapsCuboid

  pure function ovkOverlapsQuad(VertexCoords, Coords) result(Overlaps)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: Coords
    logical :: Overlaps

    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    integer, dimension(3,2), parameter :: Triangles = reshape([1,2,3,4,3,2], [3,2])
    integer :: i
    real(rk), dimension(2,2) :: Basis
    real(rk), dimension(2) :: RelativeCoords
    real(rk), dimension(2) :: LocalCoords

    ! Decompose quad into 2 triangles

    Overlaps = .false.

    do i = 1, 2
      Basis(:,1) = VertexCoords(:,Triangles(2,i)) - VertexCoords(:,Triangles(1,i))
      Basis(:,2) = VertexCoords(:,Triangles(3,i)) - VertexCoords(:,Triangles(1,i))
      RelativeCoords = Coords - VertexCoords(:,Triangles(1,i))
      LocalCoords = Solve2D(Basis, RelativeCoords)
      if (any(LocalCoords < -TOLERANCE)) then
        return
      else if (sum(LocalCoords) <= 1._rk+TOLERANCE) then
        Overlaps = .true.
        return
      end if
    end do

  end function ovkOverlapsQuad

  pure function ovkOverlapsHexahedron(VertexCoords, Coords) result(Overlaps)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: Coords
    logical :: Overlaps

    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    integer :: i
    integer, dimension(4,6), parameter :: Tetrahedra = reshape([1,2,3,5,2,3,5,6,3,5,6,7,2,3,6,4, &
      3,6,4,7,4,7,6,8], [4,6])
    real(rk), dimension(3,3) :: Basis
    real(rk), dimension(3) :: RelativeCoords
    real(rk), dimension(3) :: LocalCoords

    ! Decompose hexahedron into 6 tetrahedra

    Overlaps = .false.

    do i = 1, 6
      Basis(:,1) = VertexCoords(:,Tetrahedra(2,i)) - VertexCoords(:,Tetrahedra(1,i))
      Basis(:,2) = VertexCoords(:,Tetrahedra(3,i)) - VertexCoords(:,Tetrahedra(1,i))
      Basis(:,3) = VertexCoords(:,Tetrahedra(4,i)) - VertexCoords(:,Tetrahedra(1,i))
      RelativeCoords = Coords - VertexCoords(:,Tetrahedra(1,i))
      LocalCoords = Solve3D(Basis, RelativeCoords)
      if (all(LocalCoords > -TOLERANCE) .and. sum(LocalCoords) <= 1._rk+TOLERANCE) then
        Overlaps = .true.
        return
      end if
    end do

  end function ovkOverlapsHexahedron

  pure function ovkRectangleSize(VertexCoords) result(RectangleSize)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk) :: RectangleSize

    RectangleSize = product(VertexCoords(:,4) - VertexCoords(:,1))

  end function ovkRectangleSize

  pure function ovkCuboidSize(VertexCoords) result(CuboidSize)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk) :: CuboidSize

    CuboidSize = product(VertexCoords(:,8) - VertexCoords(:,1))

  end function ovkCuboidSize

  pure function ovkQuadSize(VertexCoords) result(QuadSize)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk) :: QuadSize

    integer, dimension(3,2), parameter :: Triangles = reshape([1,2,3,4,3,2], [3,2])
    integer :: i
    real(rk), dimension(2,2) :: Basis
    real(rk) :: TriangleSize

    ! Decompose quad into 2 triangles

    QuadSize = 0._rk

    do i = 1, 2
      Basis(:,1) = VertexCoords(:,Triangles(2,i)) - VertexCoords(:,Triangles(1,i))
      Basis(:,2) = VertexCoords(:,Triangles(3,i)) - VertexCoords(:,Triangles(1,i))
      TriangleSize = 0.5_rk * abs(Basis(1,1) * Basis(2,2) - Basis(1,2) * Basis(2,1))
      QuadSize = QuadSize + TriangleSize
    end do

  end function ovkQuadSize

  pure function ovkHexahedronSize(VertexCoords) result(HexahedronSize)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk) :: HexahedronSize

    integer, dimension(4,6), parameter :: Tetrahedra = reshape([1,2,3,5,2,3,5,6,3,5,6,7,2,3,6,4, &
      3,6,4,7,4,7,6,8], [4,6])
    integer :: i
    real(rk), dimension(3,3) :: Basis
    real(rk) :: TetrahedronSize

    ! Decompose hexahedron into 6 tetrahedra

    HexahedronSize = 0._rk

    do i = 1, 6
      Basis(:,1) = VertexCoords(:,Tetrahedra(2,i)) - VertexCoords(:,Tetrahedra(1,i))
      Basis(:,2) = VertexCoords(:,Tetrahedra(3,i)) - VertexCoords(:,Tetrahedra(1,i))
      Basis(:,3) = VertexCoords(:,Tetrahedra(4,i)) - VertexCoords(:,Tetrahedra(1,i))
      TetrahedronSize = abs( &
        Basis(1,1) * (Basis(2,2)*Basis(3,3) - Basis(2,3)*Basis(3,2)) + &
        Basis(1,2) * (Basis(2,3)*Basis(3,1) - Basis(2,1)*Basis(3,3)) + &
        Basis(1,3) * (Basis(2,1)*Basis(3,2) - Basis(2,2)*Basis(3,1)))/6._rk
      HexahedronSize = HexahedronSize + TetrahedronSize
    end do

  end function ovkHexahedronSize

  pure function ovkRectangleIsoLinear(VertexCoords, LocalCoords) result(Coords)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: LocalCoords
    real(rk), dimension(2) :: Coords

    Coords = (1._rk - LocalCoords) * VertexCoords(:,1) + LocalCoords * VertexCoords(:,4)

  end function ovkRectangleIsoLinear

  pure function ovkRectangleIsoCubic(VertexCoords, LocalCoords) result(Coords)

    real(rk), dimension(2,16), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: LocalCoords
    real(rk), dimension(2) :: Coords

    Coords = (1._rk - LocalCoords) * VertexCoords(:,6) + LocalCoords * VertexCoords(:,11)

  end function ovkRectangleIsoCubic

  pure function ovkCuboidIsoLinear(VertexCoords, LocalCoords) result(Coords)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: LocalCoords
    real(rk), dimension(3) :: Coords

    Coords = (1._rk - LocalCoords) * VertexCoords(:,1) + LocalCoords * VertexCoords(:,8)

  end function ovkCuboidIsoLinear

  pure function ovkCuboidIsoCubic(VertexCoords, LocalCoords) result(Coords)

    real(rk), dimension(3,64), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: LocalCoords
    real(rk), dimension(3) :: Coords

    Coords = (1._rk - LocalCoords) * VertexCoords(:,22) + LocalCoords * VertexCoords(:,43)

  end function ovkCuboidIsoCubic

  pure function ovkQuadIsoLinear_LocalCoords(VertexCoords, LocalCoords) result(Coords)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: LocalCoords
    real(rk), dimension(2) :: Coords

    real(rk), dimension(2,2) :: BasisValues

    BasisValues(:,1) = ovkInterpBasisLinear(LocalCoords(1))
    BasisValues(:,2) = ovkInterpBasisLinear(LocalCoords(2))

    Coords = ovkQuadIsoLinear_BasisValues(VertexCoords, BasisValues)

  end function ovkQuadIsoLinear_LocalCoords

  pure function ovkQuadIsoLinear_BasisValues(VertexCoords, BasisValues) result(Coords)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk), dimension(2,2), intent(in) :: BasisValues
    real(rk), dimension(2) :: Coords

    integer :: i, j, l

    Coords = 0._rk
    do j = 1, 2
      do i = 1, 2
        l = 1 + (i-1) + 2 * (j-1)
        Coords = Coords + BasisValues(i,1) * BasisValues(j,2) * VertexCoords(:,l)
      end do
    end do

  end function ovkQuadIsoLinear_BasisValues

  pure function ovkQuadIsoCubic_LocalCoords(VertexCoords, LocalCoords) result(Coords)

    real(rk), dimension(2,16), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: LocalCoords
    real(rk), dimension(2) :: Coords

    real(rk), dimension(4,2) :: BasisValues

    BasisValues(:,1) = ovkInterpBasisCubic(LocalCoords(1))
    BasisValues(:,2) = ovkInterpBasisCubic(LocalCoords(2))

    Coords = ovkQuadIsoCubic_BasisValues(VertexCoords, BasisValues)

  end function ovkQuadIsoCubic_LocalCoords

  pure function ovkQuadIsoCubic_BasisValues(VertexCoords, BasisValues) result(Coords)

    real(rk), dimension(2,16), intent(in) :: VertexCoords
    real(rk), dimension(4,2), intent(in) :: BasisValues
    real(rk), dimension(2) :: Coords

    integer :: i, j, l

    Coords = 0._rk
    do j = 1, 4
      do i = 1, 4
        l = 1 + (i-1) + 4 * (j-1)
        Coords = Coords + BasisValues(i,1) * BasisValues(j,2) * VertexCoords(:,l)
      end do
    end do

  end function ovkQuadIsoCubic_BasisValues

  pure function ovkHexahedronIsoLinear_LocalCoords(VertexCoords, LocalCoords) result(Coords)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: LocalCoords
    real(rk), dimension(3) :: Coords

    real(rk), dimension(2,3) :: BasisValues

    BasisValues(:,1) = ovkInterpBasisLinear(LocalCoords(1))
    BasisValues(:,2) = ovkInterpBasisLinear(LocalCoords(2))
    BasisValues(:,3) = ovkInterpBasisLinear(LocalCoords(3))

    Coords = ovkHexahedronIsoLinear_BasisValues(VertexCoords, BasisValues)

  end function ovkHexahedronIsoLinear_LocalCoords

  pure function ovkHexahedronIsoLinear_BasisValues(VertexCoords, BasisValues) result(Coords)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk), dimension(2,3), intent(in) :: BasisValues
    real(rk), dimension(3) :: Coords

    integer :: i, j, k, l

    Coords = 0._rk
    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          l = 1 + (i-1) + 2 * (j-1) + 4 * (k-1)
          Coords = Coords + BasisValues(i,1) * BasisValues(j,2) * BasisValues(k,3) * &
            VertexCoords(:,l)
        end do
      end do
    end do

  end function ovkHexahedronIsoLinear_BasisValues

  pure function ovkHexahedronIsoCubic_LocalCoords(VertexCoords, LocalCoords) result(Coords)

    real(rk), dimension(3,64), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: LocalCoords
    real(rk), dimension(3) :: Coords

    real(rk), dimension(4,3) :: BasisValues

    BasisValues(:,1) = ovkInterpBasisCubic(LocalCoords(1))
    BasisValues(:,2) = ovkInterpBasisCubic(LocalCoords(2))
    BasisValues(:,3) = ovkInterpBasisCubic(LocalCoords(3))

    Coords = ovkHexahedronIsoCubic_BasisValues(VertexCoords, BasisValues)

  end function ovkHexahedronIsoCubic_LocalCoords

  pure function ovkHexahedronIsoCubic_BasisValues(VertexCoords, BasisValues) result(Coords)

    real(rk), dimension(3,64), intent(in) :: VertexCoords
    real(rk), dimension(4,3), intent(in) :: BasisValues
    real(rk), dimension(3) :: Coords

    integer :: i, j, k, l

    Coords = 0._rk
    do k = 1, 4
      do j = 1, 4
        do i = 1, 4
          l = 1 + (i-1) + 4 * (j-1) + 16 * (k-1)
          Coords = Coords + BasisValues(i,1) * BasisValues(j,2) * BasisValues(k,3) * &
            VertexCoords(:,l)
        end do
      end do
    end do

  end function ovkHexahedronIsoCubic_BasisValues

  pure function ovkRectangleIsoInverseLinear(VertexCoords, Coords) result(LocalCoords)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: Coords
    real(rk), dimension(2) :: LocalCoords

    LocalCoords = (Coords - VertexCoords(:,1))/(VertexCoords(:,4) - VertexCoords(:,1))

  end function ovkRectangleIsoInverseLinear

  pure function ovkRectangleIsoInverseCubic(VertexCoords, Coords) result(LocalCoords)

    real(rk), dimension(2,16), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: Coords
    real(rk), dimension(2) :: LocalCoords

    LocalCoords = (Coords - VertexCoords(:,6))/(VertexCoords(:,11) - VertexCoords(:,6))

  end function ovkRectangleIsoInverseCubic

  pure function ovkCuboidIsoInverseLinear(VertexCoords, Coords) result(LocalCoords)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: Coords
    real(rk), dimension(3) :: LocalCoords

    LocalCoords = (Coords - VertexCoords(:,1))/(VertexCoords(:,8) - VertexCoords(:,1))

  end function ovkCuboidIsoInverseLinear

  pure function ovkCuboidIsoInverseCubic(VertexCoords, Coords) result(LocalCoords)

    real(rk), dimension(3,64), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: Coords
    real(rk), dimension(3) :: LocalCoords

    LocalCoords = (Coords - VertexCoords(:,22))/(VertexCoords(:,43) - VertexCoords(:,22))

  end function ovkCuboidIsoInverseCubic

  function ovkQuadIsoInverseLinear(VertexCoords, Coords, Guess, Success) result(LocalCoords)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: Coords
    real(rk), dimension(2), intent(in), optional :: Guess
    logical, intent(out), optional :: Success
    real(rk), dimension(2) :: LocalCoords

    integer, parameter :: MAX_STEPS = 100
    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    integer :: i
    real(rk), dimension(2,2) :: BasisValues
    real(rk), dimension(2,2) :: BasisDerivValues
    real(rk), dimension(2) :: Error
    real(rk), dimension(2,2) :: Jacobian

    if (present(Guess)) then
      LocalCoords = Guess
    else
      LocalCoords = 0.5_rk
    end if

    BasisValues(:,1) = ovkInterpBasisLinear(LocalCoords(1))
    BasisValues(:,2) = ovkInterpBasisLinear(LocalCoords(2))
    BasisDerivValues(:,1) = ovkInterpBasisLinearDeriv(LocalCoords(1))
    BasisDerivValues(:,2) = ovkInterpBasisLinearDeriv(LocalCoords(2))
    Error = Coords - ovkQuadIsoLinear(VertexCoords, BasisValues)

    i = 1
    do while (any(abs(Error) > TOLERANCE) .and. i <= MAX_STEPS)
      Jacobian = QuadIsoLinearJacobian(VertexCoords, BasisValues, BasisDerivValues)
      LocalCoords = LocalCoords + Solve2D(Jacobian, Error)
      BasisValues(:,1) = ovkInterpBasisLinear(LocalCoords(1))
      BasisValues(:,2) = ovkInterpBasisLinear(LocalCoords(2))
      BasisDerivValues(:,1) = ovkInterpBasisLinearDeriv(LocalCoords(1))
      BasisDerivValues(:,2) = ovkInterpBasisLinearDeriv(LocalCoords(2))
      Error = Coords - ovkQuadIsoLinear(VertexCoords, BasisValues)
      i = i + 1
    end do

    if (present(Success)) then
      Success = .not. any(ieee_is_nan(Error))
      if (Success) then
        Success = .not. any(abs(Error) > TOLERANCE)
      end if
    end if

  end function ovkQuadIsoInverseLinear

  pure function QuadIsoLinearJacobian(VertexCoords, BasisValues, BasisDerivValues) result(Jacobian)

    real(rk), dimension(2,4), intent(in) :: VertexCoords
    real(rk), dimension(2,2), intent(in) :: BasisValues
    real(rk), dimension(2,2), intent(in) :: BasisDerivValues
    real(rk), dimension(2,2) :: Jacobian

    integer :: i, j, l

    Jacobian = 0._rk

    do j = 1, 2
      do i = 1, 2
        l = 1 + (i-1) + 2 * (j-1)
        Jacobian(:,1) = Jacobian(:,1) + BasisDerivValues(i,1) * BasisValues(j,2) * VertexCoords(:,l)
        Jacobian(:,2) = Jacobian(:,2) + BasisValues(i,1) * BasisDerivValues(j,2) * VertexCoords(:,l)
      end do
    end do

  end function QuadIsoLinearJacobian

  function ovkQuadIsoInverseCubic(VertexCoords, Coords, Guess, Success) result(LocalCoords)

    real(rk), dimension(2,16), intent(in) :: VertexCoords
    real(rk), dimension(2), intent(in) :: Coords
    real(rk), dimension(2), intent(in), optional :: Guess
    logical, intent(out), optional :: Success
    real(rk), dimension(2) :: LocalCoords

    integer, parameter :: MAX_STEPS = 100
    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    integer :: i
    real(rk), dimension(4,2) :: BasisValues
    real(rk), dimension(4,2) :: BasisDerivValues
    real(rk), dimension(2) :: Error
    real(rk), dimension(2,2) :: Jacobian

    if (present(Guess)) then
      LocalCoords = Guess
    else
      LocalCoords = 0.5_rk
    end if

    BasisValues(:,1) = ovkInterpBasisCubic(LocalCoords(1))
    BasisValues(:,2) = ovkInterpBasisCubic(LocalCoords(2))
    BasisDerivValues(:,1) = ovkInterpBasisCubicDeriv(LocalCoords(1))
    BasisDerivValues(:,2) = ovkInterpBasisCubicDeriv(LocalCoords(2))
    Error = Coords - ovkQuadIsoCubic(VertexCoords, BasisValues)

    i = 1
    do while (any(abs(Error) > TOLERANCE) .and. i <= MAX_STEPS)
      Jacobian = QuadIsoCubicJacobian(VertexCoords, BasisValues, BasisDerivValues)
      LocalCoords = LocalCoords + Solve2D(Jacobian, Error)
      BasisValues(:,1) = ovkInterpBasisCubic(LocalCoords(1))
      BasisValues(:,2) = ovkInterpBasisCubic(LocalCoords(2))
      BasisDerivValues(:,1) = ovkInterpBasisCubicDeriv(LocalCoords(1))
      BasisDerivValues(:,2) = ovkInterpBasisCubicDeriv(LocalCoords(2))
      Error = Coords - ovkQuadIsoCubic(VertexCoords, BasisValues)
      i = i + 1
    end do

    if (present(Success)) then
      Success = .not. any(ieee_is_nan(Error))
      if (Success) then
        Success = .not. any(abs(Error) > TOLERANCE)
      end if
    end if

  end function ovkQuadIsoInverseCubic

  pure function QuadIsoCubicJacobian(VertexCoords, BasisValues, BasisDerivValues) result(Jacobian)

    real(rk), dimension(2,16), intent(in) :: VertexCoords
    real(rk), dimension(4,2), intent(in) :: BasisValues
    real(rk), dimension(4,2), intent(in) :: BasisDerivValues
    real(rk), dimension(2,2) :: Jacobian

    integer :: i, j, l

    Jacobian = 0._rk

    do j = 1, 4
      do i = 1, 4
        l = 1 + (i-1) + 4 * (j-1)
        Jacobian(:,1) = Jacobian(:,1) + BasisDerivValues(i,1) * BasisValues(j,2) * VertexCoords(:,l)
        Jacobian(:,2) = Jacobian(:,2) + BasisValues(i,1) * BasisDerivValues(j,2) * VertexCoords(:,l)
      end do
    end do

  end function QuadIsoCubicJacobian

  function ovkHexahedronIsoInverseLinear(VertexCoords, Coords, Guess, Success) result(LocalCoords)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: Coords
    real(rk), dimension(3), intent(in), optional :: Guess
    logical, intent(out), optional :: Success
    real(rk), dimension(3) :: LocalCoords

    integer, parameter :: MAX_STEPS = 100
    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    integer :: i
    real(rk), dimension(2,3) :: BasisValues
    real(rk), dimension(2,3) :: BasisDerivValues
    real(rk), dimension(3) :: Error
    real(rk), dimension(3,3) :: Jacobian

    if (present(Guess)) then
      LocalCoords = Guess
    else
      LocalCoords = 0.5_rk
    end if

    BasisValues(:,1) = ovkInterpBasisLinear(LocalCoords(1))
    BasisValues(:,2) = ovkInterpBasisLinear(LocalCoords(2))
    BasisValues(:,3) = ovkInterpBasisLinear(LocalCoords(3))
    BasisDerivValues(:,1) = ovkInterpBasisLinearDeriv(LocalCoords(1))
    BasisDerivValues(:,2) = ovkInterpBasisLinearDeriv(LocalCoords(2))
    BasisDerivValues(:,3) = ovkInterpBasisLinearDeriv(LocalCoords(3))
    Error = Coords - ovkHexahedronIsoLinear(VertexCoords, BasisValues)

    i = 1
    do while (any(abs(Error) > TOLERANCE) .and. i <= MAX_STEPS)
      Jacobian = HexahedronIsoLinearJacobian(VertexCoords, BasisValues, BasisDerivValues)
      LocalCoords = LocalCoords + Solve3D(Jacobian, Error)
      BasisValues(:,1) = ovkInterpBasisLinear(LocalCoords(1))
      BasisValues(:,2) = ovkInterpBasisLinear(LocalCoords(2))
      BasisValues(:,3) = ovkInterpBasisLinear(LocalCoords(3))
      BasisDerivValues(:,1) = ovkInterpBasisLinearDeriv(LocalCoords(1))
      BasisDerivValues(:,2) = ovkInterpBasisLinearDeriv(LocalCoords(2))
      BasisDerivValues(:,3) = ovkInterpBasisLinearDeriv(LocalCoords(3))
      Error = Coords - ovkHexahedronIsoLinear(VertexCoords, BasisValues)
      i = i + 1
    end do

    if (present(Success)) then
      Success = .not. any(ieee_is_nan(Error))
      if (Success) then
        Success = .not. any(abs(Error) > TOLERANCE)
      end if
    end if

  end function ovkHexahedronIsoInverseLinear

  pure function HexahedronIsoLinearJacobian(VertexCoords, BasisValues, BasisDerivValues) &
    result(Jacobian)

    real(rk), dimension(3,8), intent(in) :: VertexCoords
    real(rk), dimension(2,3), intent(in) :: BasisValues
    real(rk), dimension(2,3), intent(in) :: BasisDerivValues
    real(rk), dimension(3,3) :: Jacobian

    integer :: i, j, k, l

    Jacobian = 0._rk

    do k = 1, 2
      do j = 1, 2
        do i = 1, 2
          l = 1 + (i-1) + 2 * (j-1) + 4 * (k-1)
          Jacobian(:,1) = Jacobian(:,1) + BasisDerivValues(i,1) * BasisValues(j,2) * &
            BasisValues(k,3) * VertexCoords(:,l)
          Jacobian(:,2) = Jacobian(:,2) + BasisValues(i,1) * BasisDerivValues(j,2) * &
            BasisValues(k,3) * VertexCoords(:,l)
          Jacobian(:,3) = Jacobian(:,3) + BasisValues(i,1) * BasisValues(j,2) * &
            BasisDerivValues(k,3) * VertexCoords(:,l)
        end do
      end do
    end do

  end function HexahedronIsoLinearJacobian

  function ovkHexahedronIsoInverseCubic(VertexCoords, Coords, Guess, Success) result(LocalCoords)

    real(rk), dimension(3,64), intent(in) :: VertexCoords
    real(rk), dimension(3), intent(in) :: Coords
    real(rk), dimension(3), intent(in), optional :: Guess
    logical, intent(out), optional :: Success
    real(rk), dimension(3) :: LocalCoords

    integer, parameter :: MAX_STEPS = 100
    real(rk), parameter :: TOLERANCE = 1.e-10_rk

    integer :: i
    real(rk), dimension(4,3) :: BasisValues
    real(rk), dimension(4,3) :: BasisDerivValues
    real(rk), dimension(3) :: Error
    real(rk), dimension(3,3) :: Jacobian

    if (present(Guess)) then
      LocalCoords = Guess
    else
      LocalCoords = 0.5_rk
    end if

    BasisValues(:,1) = ovkInterpBasisCubic(LocalCoords(1))
    BasisValues(:,2) = ovkInterpBasisCubic(LocalCoords(2))
    BasisValues(:,3) = ovkInterpBasisCubic(LocalCoords(3))
    BasisDerivValues(:,1) = ovkInterpBasisCubicDeriv(LocalCoords(1))
    BasisDerivValues(:,2) = ovkInterpBasisCubicDeriv(LocalCoords(2))
    BasisDerivValues(:,3) = ovkInterpBasisCubicDeriv(LocalCoords(3))
    Error = Coords - ovkHexahedronIsoCubic(VertexCoords, BasisValues)

    i = 1
    do while (any(abs(Error) > TOLERANCE) .and. i <= MAX_STEPS)
      Jacobian = HexahedronIsoCubicJacobian(VertexCoords, BasisValues, BasisDerivValues)
      LocalCoords = LocalCoords + Solve3D(Jacobian, Error)
      BasisValues(:,1) = ovkInterpBasisCubic(LocalCoords(1))
      BasisValues(:,2) = ovkInterpBasisCubic(LocalCoords(2))
      BasisValues(:,3) = ovkInterpBasisCubic(LocalCoords(3))
      BasisDerivValues(:,1) = ovkInterpBasisCubicDeriv(LocalCoords(1))
      BasisDerivValues(:,2) = ovkInterpBasisCubicDeriv(LocalCoords(2))
      BasisDerivValues(:,3) = ovkInterpBasisCubicDeriv(LocalCoords(3))
      Error = Coords - ovkHexahedronIsoCubic(VertexCoords, BasisValues)
      i = i + 1
    end do

    if (present(Success)) then
      Success = .not. any(ieee_is_nan(Error))
      if (Success) then
        Success = .not. any(abs(Error) > TOLERANCE)
      end if
    end if

  end function ovkHexahedronIsoInverseCubic

  pure function HexahedronIsoCubicJacobian(VertexCoords, BasisValues, BasisDerivValues) &
    result(Jacobian)

    real(rk), dimension(3,64), intent(in) :: VertexCoords
    real(rk), dimension(4,3), intent(in) :: BasisValues
    real(rk), dimension(4,3), intent(in) :: BasisDerivValues
    real(rk), dimension(3,3) :: Jacobian

    integer :: i, j, k, l

    Jacobian = 0._rk

    do k = 1, 4
      do j = 1, 4
        do i = 1, 4
          l = 1 + (i-1) + 4 * (j-1) + 16 * (k-1)
          Jacobian(:,1) = Jacobian(:,1) + BasisDerivValues(i,1) * BasisValues(j,2) * &
            BasisValues(k,3) * VertexCoords(:,l)
          Jacobian(:,2) = Jacobian(:,2) + BasisValues(i,1) * BasisDerivValues(j,2) * &
            BasisValues(k,3) * VertexCoords(:,l)
          Jacobian(:,3) = Jacobian(:,3) + BasisValues(i,1) * BasisValues(j,2) * &
            BasisDerivValues(k,3) * VertexCoords(:,l)
        end do
      end do
    end do

  end function HexahedronIsoCubicJacobian

  pure function ovkCartesianGridCell_Scalar(Origin, CellSize, Coords) result(Cell)

    real(rk), intent(in) :: Origin
    real(rk), intent(in) :: CellSize
    real(rk), intent(in) :: Coords
    integer(lk) :: Cell

    Cell = int(floor((Coords - Origin)/CellSize),kind=lk) + 1_lk

  end function ovkCartesianGridCell_Scalar

  pure function ovkCartesianGridCell_Array(Origin, CellSize, Coords) result(Cell)

    real(rk), dimension(:), intent(in) :: Origin
    real(rk), dimension(:), intent(in) :: CellSize
    real(rk), dimension(:), intent(in) :: Coords
    integer(lk), dimension(size(Origin)) :: Cell

    Cell = int(floor((Coords - Origin)/CellSize),kind=lk) + 1_lk

  end function ovkCartesianGridCell_Array

  pure function ovkInterpBasisLinear(T) result(Basis)

    real(rk), intent(in) :: T
    real(rk), dimension(2) :: Basis

    Basis(1) = 1._rk - T
    Basis(2) = T

  end function ovkInterpBasisLinear

  pure function ovkInterpBasisLinearDeriv(T) result(BasisDeriv)

    real(rk), intent(in) :: T
    real(rk), dimension(2) :: BasisDeriv

    BasisDeriv(1) = -1._rk
    BasisDeriv(2) = 1._rk

  end function ovkInterpBasisLinearDeriv

  pure function ovkInterpBasisCubic(T) result(Basis)

    real(rk), intent(in) :: T
    real(rk), dimension(4) :: Basis

    Basis(1) = -(T * (T - 1._rk) * (T - 2._rk)) / 6._rk
    Basis(2) = ((T + 1._rk) * (T - 1._rk) * (T - 2._rk)) / 2._rk
    Basis(3) = -((T + 1._rk) * T * (T - 2._rk)) / 2._rk
    Basis(4) = ((T + 1._rk) * T * (T - 1._rk)) / 6._rk

  end function ovkInterpBasisCubic

  pure function ovkInterpBasisCubicDeriv(T) result(BasisDeriv)

    real(rk), intent(in) :: T
    real(rk), dimension(4) :: BasisDeriv

    real(rk), dimension(2,4), parameter :: Roots = reshape([ &
      1._rk + 1._rk/sqrt(3._rk), 1._rk - 1._rk/sqrt(3._rk), &
      (2._rk + sqrt(7._rk))/3._rk, (2._rk - sqrt(7._rk))/3._rk, &
      (1._rk + sqrt(7._rk))/3._rk, (1._rk - sqrt(7._rk))/3._rk, &
      1._rk/sqrt(3._rk), -1._rk/sqrt(3._rk) &
    ], [2,4])

    BasisDeriv(1) = -((T - Roots(1,1)) * (T - Roots(2,1))) / 2._rk
    BasisDeriv(2) = (3._rk * (T - Roots(1,2)) * (T - Roots(2,2))) / 2._rk
    BasisDeriv(3) = -(3._rk * (T - Roots(1,3)) * (T - Roots(2,3))) / 2._rk
    BasisDeriv(4) = ((T - Roots(1,4)) * (T - Roots(2,4))) / 2._rk

  end function ovkInterpBasisCubicDeriv

  pure function Solve2D(A, b) result(x)

    real(rk), dimension(2,2), intent(in) :: A
    real(rk), dimension(2), intent(in) :: b
    real(rk), dimension(2) :: x

    real(rk) :: Det

    Det = A(1,1) * A(2,2) - A(2,1) * A(1,2)

    x(1) = (  b(1) * A(2,2) -   b(2) * A(1,2)) / Det
    x(2) = (A(1,1) *   b(2) - A(2,1) *   b(1)) / Det

  end function Solve2D

  pure function Solve3D(A, b) result(x)

    real(rk), dimension(3,3), intent(in) :: A
    real(rk), dimension(3), intent(in) :: b
    real(rk), dimension(3) :: x

    real(rk) :: Det

    Det = &
      A(1,1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) + &
      A(1,2) * (A(2,3)*A(3,1) - A(2,1)*A(3,3)) + &
      A(1,3) * (A(2,1)*A(3,2) - A(2,2)*A(3,1))

    x(1) = ( &
        b(1) * (A(2,2)*A(3,3) - A(2,3)*A(3,2)) + &
      A(1,2) * (A(2,3)*  b(3) -   b(2)*A(3,3)) + &
      A(1,3) * (  b(2)*A(3,2) - A(2,2)*  b(3)) &
      ) / Det

    x(2) = ( &
      A(1,1) * (  b(2)*A(3,3) - A(2,3)*  b(3)) + &
        b(1) * (A(2,3)*A(3,1) - A(2,1)*A(3,3)) + &
      A(1,3) * (A(2,1)*  b(3) -   b(2)*A(3,1)) &
      ) / Det

    x(3) = ( &
      A(1,1) * (A(2,2)*  b(3) -   b(2)*A(3,2)) + &
      A(1,2) * (  b(2)*A(3,1) - A(2,1)*  b(3)) + &
        b(1) * (A(2,1)*A(3,2) - A(2,2)*A(3,1)) &
      ) / Det

  end function Solve3D

#ifndef f2003

  pure elemental function ieee_is_nan(X) result(IsNaN)

    real(rk), intent(in) :: X
    logical :: IsNaN

    IsNaN = X /= X

  end function ieee_is_nan

#endif

end module ovkGeometry
