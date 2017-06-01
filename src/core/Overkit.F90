! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module Overkit

  use ovkBoundingBox
  use ovkCart
  use ovkDonorAccel
  use ovkDonors
  use ovkField
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkHashGrid
  use ovkInterp
  use ovkMask
  use ovkOverset
  use ovkPegasus
  use ovkPLOT3D
  implicit none

  private

  ! General
  public :: operator (==)
  public :: operator (/=)

  ! ovkBoundingBox
  public :: ovk_bbox
  public :: ovk_bbox_
  public :: ovkBBOverlaps
  public :: ovkBBContains
  public :: ovkBBContainsPoint
  public :: ovkBBIsEmpty
  public :: ovkBBSize
  public :: ovkBBMove
  public :: ovkBBGrow
  public :: ovkBBScale
  public :: ovkBBExtend
  public :: ovkBBUnion
  public :: ovkBBIntersect
  public :: ovkBBFromPoints

  ! ovkCart
  public :: ovk_cart
  public :: ovk_cart_
  public :: ovkCartSize
  public :: ovkCartCount
  public :: ovkCartTupleToIndex
  public :: ovkCartIndexToTuple
  public :: ovkCartPeriodicAdjust
  public :: ovkCartContains
  public :: ovkCartClamp
  public :: ovkCartConvertPeriodicStorage
  public :: ovkCartConvertPointToCell

  ! ovkDonorAccel
  public :: ovk_donor_accel
  public :: ovk_donor_accel_
  public :: ovkGenerateDonorAccel
  public :: ovkDestroyDonorAccel
  public :: ovkFindDonorCell

  ! ovkDonors
  public :: ovk_donors
  public :: ovk_donors_
  public :: ovkMakeDonors
  public :: ovkDestroyDonors
  public :: ovkFindDonors
  public :: ovkChooseDonors
  public :: ovkMergeDonors
  public :: ovkPrintDonors
  public :: ovkGenerateReceiverMask
  public :: ovkGenerateDonorMask
  public :: ovkGenerateOverlapMask
  public :: ovkGenerateCoarseToFineMask
  public :: ovkGenerateNearCrossoverMask
  public :: ovkGenerateOrphanMask

  ! ovkField
  public :: ovk_field_int
  public :: ovk_field_int_
  public :: ovk_field_large_int
  public :: ovk_field_large_int_
  public :: ovk_field_real
  public :: ovk_field_real_
  public :: ovk_field_logical
  public :: ovk_field_logical_
  public :: ovkExportField

  ! ovkGeometry
  public :: ovkOverlapsRectangle
  public :: ovkOverlapsQuad
  public :: ovkOverlapsCuboid
  public :: ovkOverlapsHexahedron
  public :: ovkRectangleSize
  public :: ovkQuadSize
  public :: ovkCuboidSize
  public :: ovkHexahedronSize
  public :: ovkRectangleIsoLinear
  public :: ovkRectangleIsoCubic
  public :: ovkQuadIsoLinear
  public :: ovkQuadIsoCubic
  public :: ovkCuboidIsoLinear
  public :: ovkCuboidIsoCubic
  public :: ovkHexahedronIsoLinear
  public :: ovkHexahedronIsoCubic
  public :: ovkRectangleIsoInverseLinear
  public :: ovkRectangleIsoInverseCubic
  public :: ovkQuadIsoInverseLinear
  public :: ovkQuadIsoInverseCubic
  public :: ovkCuboidIsoInverseLinear
  public :: ovkCuboidIsoInverseCubic
  public :: ovkHexahedronIsoInverseLinear
  public :: ovkHexahedronIsoInverseCubic
  public :: ovkCartesianGridCell
  public :: ovkInterpBasisLinear
  public :: ovkInterpBasisLinearDeriv
  public :: ovkInterpBasisCubic
  public :: ovkInterpBasisCubicDeriv

  ! ovkGlobal
  public :: ovk_rk
  public :: ovk_lk
  public :: ovk_bk
  public :: OVK_DEBUG
  public :: OVK_VERBOSE
  public :: OVK_NO_ERROR, OVK_IO_ERROR
  public :: OVK_NO_OVERLAP_PERIODIC, OVK_OVERLAP_PERIODIC
  public :: OVK_LITTLE_ENDIAN, OVK_BIG_ENDIAN
  public :: ovkCaseID

  ! ovkGrid
  public :: ovk_grid
  public :: ovk_grid_
  public :: ovkMakeGrid
  public :: ovkDestroyGrid
  public :: ovkGetCellVertexData
  public :: ovkOverlapsCell
  public :: ovkCoordsInCell
  public :: ovkCellSize
  public :: ovkAvgCellSizeAroundPoint
  public :: ovkGenerateBBOverlapMask
  public :: ovkPeriodicExtend
  public :: OVK_GRID_TYPE_CARTESIAN
  public :: OVK_GRID_TYPE_CARTESIAN_ROTATED
  public :: OVK_GRID_TYPE_RECTILINEAR
  public :: OVK_GRID_TYPE_RECTILINEAR_ROTATED
  public :: OVK_GRID_TYPE_CURVILINEAR

  ! ovkHashGrid
  public :: ovk_hash_grid
  public :: ovk_hash_grid_
  public :: ovkHashGridBin
  public :: ovkHashGridBinBounds
  public :: ovkHashGridStats
  public :: ovkHashGridHistogram

  ! ovkInterp
  public :: ovk_interp
  public :: ovk_interp_
  public :: ovkMakeInterpData
  public :: ovkDestroyInterpData
  public :: ovkGenerateInterpData
  public :: ovkDonorGridIDToIBlank
  public :: OVK_INTERP_LINEAR
  public :: OVK_INTERP_CUBIC

  ! ovkMask
  public :: ovkFindMaskEdge
  public :: ovkGrowMask
  public :: ovkConnectedComponents
  public :: ovkFillMask
  public :: ovkGenerateNearEdgeMask
  public :: ovkCountMask
  public :: ovkMaskToIBlank
  public :: ovkPrintMask
  public :: OVK_EDGE_TYPE_INNER, OVK_EDGE_TYPE_OUTER

  ! ovkOverset
  public :: ovkAssembleOverset
  public :: ovkPartitionReceivers
  public :: ovkGenerateOverlapOptimizationMask

  ! ovkPegasus
  public :: ovk_pegasus
  public :: ovk_pegasus_
  public :: ovkMakePegasusData
  public :: ovkDestroyPegasusData
  public :: ovkWritePegasusData
  public :: ovkPrintPegasusData

  ! ovkPLOT3D
  public :: ovk_p3d_grid_file
  public :: ovk_p3d_grid_file_
  public :: ovkP3DMachineEndian
  public :: ovkP3DOpen
  public :: ovkP3DCreate
  public :: ovkP3DClose
  public :: ovkP3DRead
  public :: ovkP3DWrite

end module Overkit
