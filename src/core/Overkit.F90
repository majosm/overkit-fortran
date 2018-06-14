! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module Overkit

  use ovkArray
  use ovkAssembly
  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDomain
  use ovkField
  use ovkFieldOps
  use ovkGeometry
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
  public :: ovk_connectivity_properties
  public :: ovkGetConnectivityProperties
  public :: ovkEditConnectivityProperties
  public :: ovkReleaseConnectivityProperties
  public :: ovkGetConnectivityDonorExtents
  public :: ovkGetConnectivityDonorCoords
  public :: ovkGetConnectivityDonorInterpCoefs
  public :: ovkGetConnectivityReceiverPoints
  public :: ovkGetConnectivityPropertyDonorGridID
  public :: ovkGetConnectivityPropertyReceiverGridID
  public :: ovkGetConnectivityPropertyDimension
  public :: ovkGetConnectivityPropertyMaxDonorSize
  public :: ovkGetConnectivityPropertyConnectionCount

  ! ovkDomain
  public :: ovk_domain
  public :: ovk_domain_properties
  public :: ovkCreateDomain
  public :: ovkDestroyDomain
  public :: ovkGetDomainProperties
  public :: ovkEditDomainProperties
  public :: ovkReleaseDomainProperties
  public :: ovkCreateGrid
  public :: ovkDestroyGrid
  public :: ovkGridExists
  public :: ovkGetGrid
  public :: ovkEditGrid
  public :: ovkReleaseGrid
  public :: ovkConnectivityExists
  public :: ovkGetConnectivity
  public :: ovkEditConnectivity
  public :: ovkReleaseConnectivity
  public :: ovkGetDomainPropertyDimension
  public :: ovkGetDomainPropertyGridCount
  public :: ovkGetDomainPropertyVerbose
  public :: ovkSetDomainPropertyVerbose
  public :: ovkGetDomainPropertyOverlappable
  public :: ovkSetDomainPropertyOverlappable
  public :: ovkGetDomainPropertyOverlapTolerance
  public :: ovkSetDomainPropertyOverlapTolerance
  public :: ovkGetDomainPropertyOverlapAccelQualityAdjust
  public :: ovkSetDomainPropertyOverlapAccelQualityAdjust
  public :: ovkGetDomainPropertyInferBoundaries
  public :: ovkSetDomainPropertyInferBoundaries
  public :: ovkGetDomainPropertyBoundaryHoleCutting
  public :: ovkSetDomainPropertyBoundaryHoleCutting
  public :: ovkGetDomainPropertyOccludes
  public :: ovkSetDomainPropertyOccludes
  public :: ovkGetDomainPropertyOcclusionPadding
  public :: ovkSetDomainPropertyOcclusionPadding
  public :: ovkGetDomainPropertyOcclusionSmoothing
  public :: ovkSetDomainPropertyOcclusionSmoothing
  public :: ovkGetDomainPropertyConnectionType
  public :: ovkSetDomainPropertyConnectionType
  public :: ovkGetDomainPropertyInterpScheme
  public :: ovkSetDomainPropertyInterpScheme
  public :: ovkGetDomainPropertyFringeSize
  public :: ovkSetDomainPropertyFringeSize
  public :: ovkGetDomainPropertyOverlapMinimization
  public :: ovkSetDomainPropertyOverlapMinimization

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
  public :: ovkExportField
  public :: ovkPrintField

  ! ovkFieldOps
  public :: ovkDetectEdge
  public :: ovkDilate
  public :: ovkErode
  public :: ovkConnectedComponents
  public :: ovkFlood
  public :: ovkThreshold
  public :: ovkDistanceField
  public :: ovkCountMask
  public :: OVK_INNER_EDGE, OVK_OUTER_EDGE

  ! ovkGeometry
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
  public :: OVK_TRUE, OVK_FALSE
  public :: OVK_NONE, OVK_ANY, OVK_ALL
  public :: OVK_NO_ERROR, OVK_IO_ERROR
  public :: OVK_AUTO
  public :: OVK_MIRROR
  public :: OVK_NO_OVERLAP_PERIODIC, OVK_OVERLAP_PERIODIC
  public :: OVK_LITTLE_ENDIAN, OVK_BIG_ENDIAN
  public :: OVK_P3D_STANDARD, OVK_P3D_EXTENDED
  public :: OVK_ALL_GRIDS
  public :: OVK_CONNECTION_NONE, OVK_CONNECTION_FRINGE, OVK_CONNECTION_FULL
  public :: OVK_INTERP_LINEAR, OVK_INTERP_CUBIC

  ! ovkGrid
  public :: ovk_grid
  public :: ovk_grid_properties
  public :: ovkGetGridProperties
  public :: ovkEditGridProperties
  public :: ovkReleaseGridProperties
  public :: ovkGetGridCart
  public :: ovkGetGridCoords
  public :: ovkEditGridCoords
  public :: ovkReleaseGridCoords
  public :: ovkGetGridState
  public :: ovkEditGridState
  public :: ovkReleaseGridState
  public :: ovkResetGridState
  public :: ovkFilterGridState
  public :: ovkGridCellExists
  public :: ovkGridCellBounds
  public :: ovkOverlapsGridCell
  public :: ovkCoordsInGridCell
  public :: ovkCoordsInCubicGridCell
  public :: ovkGridResolution
  public :: ovkGenerateBBOverlapMask
  public :: ovkPeriodicExtend
  public :: ovkExportGridCoords
  public :: ovkGetGridPropertyID
  public :: ovkGetGridPropertyDimension
  public :: ovkGetGridPropertySize
  public :: ovkGetGridPropertyPeriodicity
  public :: ovkGetGridPropertyPeriodicStorage
  public :: ovkGetGridPropertyPeriodicLength
  public :: ovkGetGridPropertyGeometryType
  public :: OVK_GRID_GEOMETRY_CARTESIAN
  public :: OVK_GRID_GEOMETRY_RECTILINEAR
  public :: OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN
  public :: OVK_GRID_GEOMETRY_ORIENTED_RECTILINEAR
  public :: OVK_GRID_GEOMETRY_CURVILINEAR
  public :: OVK_STATE_GRID
  public :: OVK_STATE_INTERIOR
  public :: OVK_STATE_BOUNDARY
  public :: OVK_STATE_DOMAIN_BOUNDARY
  public :: OVK_STATE_INFERRED_DOMAIN_BOUNDARY
  public :: OVK_STATE_INTERNAL_BOUNDARY
  public :: OVK_STATE_HOLE
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
  public :: OVK_HOLE_POINT

  ! ovkOverlap
  public :: ovk_overlap
  public :: ovk_overlap_properties
  public :: ovkGetOverlapProperties
  public :: ovkGetOverlapCart
  public :: ovkGetOverlapBounds
  public :: ovkGetOverlapMask
  public :: ovkGetOverlapCells
  public :: ovkGetOverlapCoords
  public :: ovkFindOverlappingPoints
  public :: ovkFindOverlappedPoints
  public :: ovkOverlapCollect
  public :: ovkOverlapDisperse
  public :: ovkGetOverlapPropertyOverlappingGridID
  public :: ovkGetOverlapPropertyOverlappedGridID
  public :: ovkGetOverlapPropertyDimension
  public :: ovkGetOverlapPropertySize
  public :: ovkGetOverlapPropertyPeriodicity
  public :: ovkGetOverlapPropertyPeriodicStorage
  public :: ovkGetOverlapPropertyNumOverlapped
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
