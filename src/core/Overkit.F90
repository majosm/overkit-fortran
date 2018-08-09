! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module Overkit

  use ovkArray
  use ovkAssembly
  use ovkAssemblyOptions
  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDomain
  use ovkField
  use ovkFieldOps
  use ovkGeometryOps
  use ovkGlobal
  use ovkGrid
  use ovkHashGrid
  use ovkOverlap
  use ovkPLOT3D
  implicit none

  private

  ! General
  public :: operator (==)
  public :: operator (/=)

  ! ovkArray
  public :: ovk_array_int
  public :: ovk_array_int_
  public :: ovk_array_large_int
  public :: ovk_array_large_int_
  public :: ovk_array_real
  public :: ovk_array_real_
  public :: ovk_array_logical
  public :: ovk_array_logical_

  ! ovkAssembly
  public :: ovkAssemble

  ! ovkAssemblyOptions
  public :: ovk_assembly_options
  public :: ovk_assembly_options_
  public :: ovkGetAssemblyOptionsDimension
  public :: ovkGetAssemblyOptionsGridCount
  public :: ovkGetAssemblyOptionOverlappable
  public :: ovkSetAssemblyOptionOverlappable
  public :: ovkGetAssemblyOptionOverlapTolerance
  public :: ovkSetAssemblyOptionOverlapTolerance
  public :: ovkGetAssemblyOptionOverlapAccelDepthAdjust
  public :: ovkSetAssemblyOptionOverlapAccelDepthAdjust
  public :: ovkGetAssemblyOptionOverlapAccelResolutionAdjust
  public :: ovkSetAssemblyOptionOverlapAccelResolutionAdjust
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
  public :: ovkGetAssemblyOptionMinimizeOverlap
  public :: ovkSetAssemblyOptionMinimizeOverlap

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
  public :: ovkCartIsEmpty
  public :: ovkCartSize
  public :: ovkCartCount
  public :: ovkCartTupleToIndex
  public :: ovkCartIndexToTuple
  public :: ovkCartPeriodicAdjust
  public :: ovkCartContains
  public :: ovkCartClamp
  public :: ovkCartIsCompatible
  public :: ovkCartConvertPeriodicStorage
  public :: ovkCartPointToCell

  ! ovkConnectivity
  public :: ovk_connectivity
  public :: ovkConnectivityExists
  public :: ovkResetConnectivity
  public :: ovkGetConnectivityDonorGrid
  public :: ovkGetConnectivityReceiverGrid
  public :: ovkGetConnectivityDimension
  public :: ovkGetConnectivityMaxDonorSize
  public :: ovkGetConnectivityCount
  public :: ovkGetConnectivityDonorExtents
  public :: ovkGetConnectivityDonorCoords
  public :: ovkGetConnectivityDonorInterpCoefs
  public :: ovkGetConnectivityReceiverPoints

  ! ovkDomain
  public :: ovk_domain
  public :: ovkCreateDomain
  public :: ovkDestroyDomain
  public :: ovkDomainExists
  public :: ovkGetDomainDimension
  public :: ovkGetDomainGridCount
  public :: ovkCreateGrid
  public :: ovkDestroyGrid
  public :: ovkHasGrid
  public :: ovkGetGrid
  public :: ovkEditGrid
  public :: ovkReleaseGrid
  public :: ovkCreateOverlap
  public :: ovkDestroyOverlap
  public :: ovkHasOverlap
  public :: ovkGetOverlap
  public :: ovkEditOverlap
  public :: ovkReleaseOverlap
  public :: ovkCreateConnectivity
  public :: ovkDestroyConnectivity
  public :: ovkHasConnectivity
  public :: ovkGetConnectivity
  public :: ovkEditConnectivity
  public :: ovkReleaseConnectivity

  ! ovkField
  public :: ovk_field_int
  public :: ovk_field_int_
  public :: ovk_field_large_int
  public :: ovk_field_large_int_
  public :: ovk_field_real
  public :: ovk_field_real_
  public :: ovk_field_logical
  public :: ovk_field_logical_
  public :: ovkFieldPeriodicFill
  public :: ovkGetFieldPatch
  public :: ovkPrintField

  ! ovkFieldOps
  public :: ovkDetectEdge
  public :: ovkDilate
  public :: ovkErode
  public :: ovkConnectedComponents
  public :: ovkFlood
  public :: ovkDistanceField
  public :: ovkCountMask
  public :: OVK_INNER_EDGE, OVK_OUTER_EDGE

  ! ovkGeometryOps
  public :: ovkOverlapsRectangle
  public :: ovkOverlapsCuboid
  public :: ovkOverlapsOrientedRectangle
  public :: ovkOverlapsOrientedCuboid
  public :: ovkOverlapsQuad
  public :: ovkOverlapsHexahedron
  public :: ovkRectangleSize
  public :: ovkCuboidSize
  public :: ovkOrientedRectangleSize
  public :: ovkOrientedCuboidSize
  public :: ovkQuadSize
  public :: ovkHexahedronSize
  public :: ovkRectangleIsoLinear
  public :: ovkRectangleIsoCubic
  public :: ovkCuboidIsoLinear
  public :: ovkCuboidIsoCubic
  public :: ovkOrientedRectangleIsoLinear
  public :: ovkOrientedRectangleIsoCubic
  public :: ovkOrientedCuboidIsoLinear
  public :: ovkOrientedCuboidIsoCubic
  public :: ovkQuadIsoLinear
  public :: ovkQuadIsoCubic
  public :: ovkHexahedronIsoLinear
  public :: ovkHexahedronIsoCubic
  public :: ovkRectangleIsoInverseLinear
  public :: ovkRectangleIsoInverseCubic
  public :: ovkCuboidIsoInverseLinear
  public :: ovkCuboidIsoInverseCubic
  public :: ovkOrientedRectangleIsoInverseLinear
  public :: ovkOrientedRectangleIsoInverseCubic
  public :: ovkOrientedCuboidIsoInverseLinear
  public :: ovkOrientedCuboidIsoInverseCubic
  public :: ovkQuadIsoInverseLinear
  public :: ovkQuadIsoInverseCubic
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
  public :: ovkCaseID
  public :: OVK_DEBUG
  public :: OVK_FALSE, OVK_TRUE
  public :: OVK_NONE, OVK_ANY, OVK_NOT_ALL, OVK_ALL
  public :: OVK_NO_ERROR, OVK_IO_ERROR
  public :: OVK_MIRROR
  public :: OVK_NO_OVERLAP_PERIODIC, OVK_OVERLAP_PERIODIC
  public :: OVK_LITTLE_ENDIAN, OVK_BIG_ENDIAN
  public :: OVK_P3D_STANDARD, OVK_P3D_EXTENDED
  public :: OVK_ALL_GRIDS
  public :: OVK_CONNECTION_NONE, OVK_CONNECTION_NEAREST, OVK_CONNECTION_LINEAR, OVK_CONNECTION_CUBIC
  public :: OVK_OCCLUDES_NONE, OVK_OCCLUDES_ALL, OVK_OCCLUDES_COARSE

  ! ovkGrid
  public :: ovk_grid
  public :: ovkGridExists
  public :: ovkGetGridID
  public :: ovkGetGridDimension
  public :: ovkGetGridSize
  public :: ovkGetGridCart
  public :: ovkGetGridPeriodicity
  public :: ovkGetGridPeriodicStorage
  public :: ovkGetGridPeriodicLength
  public :: ovkGetGridGeometryType
  public :: ovkGetGridCoords
  public :: ovkEditGridCoords
  public :: ovkReleaseGridCoords
  public :: ovkGetGridState
  public :: ovkEditGridState
  public :: ovkReleaseGridState
  public :: ovkResetGridState
  public :: ovkFilterGridState
  public :: ovkGetGridBounds
  public :: ovkGridCellBounds
  public :: ovkOverlapsGridCell
  public :: ovkCoordsInGridCell
  public :: ovkPeriodicExtend
  public :: OVK_GEOMETRY_CARTESIAN
  public :: OVK_GEOMETRY_RECTILINEAR
  public :: OVK_GEOMETRY_ORIENTED_CARTESIAN
  public :: OVK_GEOMETRY_ORIENTED_RECTILINEAR
  public :: OVK_GEOMETRY_CURVILINEAR
  public :: OVK_STATE_GRID
  public :: OVK_STATE_INTERIOR
  public :: OVK_STATE_BOUNDARY
  public :: OVK_STATE_EXTERIOR
  public :: OVK_STATE_DOMAIN_BOUNDARY
  public :: OVK_STATE_INTERNAL_BOUNDARY
  public :: OVK_STATE_INFERRED_DOMAIN_BOUNDARY
  public :: OVK_STATE_BOUNDARY_HOLE
  public :: OVK_STATE_FRINGE
  public :: OVK_STATE_OUTER_FRINGE
  public :: OVK_STATE_INNER_FRINGE
  public :: OVK_STATE_OCCLUDED
  public :: OVK_STATE_OVERLAP_MINIMIZED
  public :: OVK_STATE_RECEIVER
  public :: OVK_STATE_ORPHAN
  public :: OVK_STATE_DEBUG1
  public :: OVK_STATE_DEBUG2
  public :: OVK_STATE_DEBUG3
  public :: OVK_STATE_DEBUG4
  public :: OVK_STATE_DEBUG5
  public :: OVK_INTERIOR_POINT
  public :: OVK_DOMAIN_BOUNDARY_POINT
  public :: OVK_INTERNAL_BOUNDARY_POINT
  public :: OVK_EXTERIOR_POINT

  ! ovkOverlap
  public :: ovk_overlap
  public :: ovkOverlapExists
  public :: ovkResetOverlap
  public :: ovkGetOverlapOverlappingGrid
  public :: ovkGetOverlapOverlappedGrid
  public :: ovkGetOverlapDimension
  public :: ovkGetOverlapCount
  public :: ovkGetOverlapMask
  public :: ovkGetOverlapCells
  public :: ovkGetOverlapCoords
  public :: ovkFindOverlappingPoints
  public :: ovkFindOverlappedPoints
  public :: ovkOverlapCollect
  public :: ovkOverlapDisperse
  public :: OVK_COLLECT_SIMPLE
  public :: OVK_COLLECT_MIN
  public :: OVK_COLLECT_MAX
  public :: OVK_COLLECT_NONE
  public :: OVK_COLLECT_ANY
  public :: OVK_COLLECT_NOT_ALL
  public :: OVK_COLLECT_ALL
  public :: OVK_COLLECT_INTERPOLATE
  public :: OVK_DISPERSE_OVERWRITE

  ! ovkPLOT3D
  public :: ovk_plot3d_grid_file
  public :: ovk_plot3d_grid_file_
  public :: ovkP3DMachineEndian
  public :: ovkOpenP3D
  public :: ovkCreateP3D
  public :: ovkCloseP3D
  public :: ovkReadP3D
  public :: ovkWriteP3D

end module Overkit
