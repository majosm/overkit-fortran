! Copyright (c) 2017 Matthew J. Smith and Overkit contributors
! License: MIT (http://opensource.org/licenses/MIT)

module Overkit

  use ovkAssembler
  use ovkBoundingBox
  use ovkCart
  use ovkConnectivity
  use ovkDomain
  use ovkDonorAccel
  use ovkDonors
  use ovkField
  use ovkFieldOps
  use ovkGeometry
  use ovkGlobal
  use ovkGrid
  use ovkHashGrid
  use ovkInterp
  use ovkOverset
  use ovkPLOT3D
  implicit none

  private

  ! General
  public :: operator (==)
  public :: operator (/=)

  ! ovkAssembler
  public :: ovk_assembler
  public :: ovk_assembler_
  public :: ovk_assembler_properties
  public :: ovk_assembler_properties_
  public :: ovkCreateAssembler
  public :: ovkDestroyAssembler
  public :: ovkGetAssemblerProperties
  public :: ovkEditAssemblerProperties
  public :: ovkReleaseAssemblerProperties
  public :: ovkGetAssemblerDomain
  public :: ovkEditAssemblerDomain
  public :: ovkReleaseAssemblerDomain
!   public :: ovkGetAssemblerOverlap
!   public :: ovkEditAssemblerOverlap
!   public :: ovkReleaseAssemblerOverlap
  public :: ovkGetAssemblerConnectivity
  public :: ovkEditAssemblerConnectivity
  public :: ovkReleaseAssemblerConnectivity
  public :: ovkGetAssemblerDebugField
  public :: ovkGetAssemblerPropertyDimension
  public :: ovkGetAssemblerPropertyGridCount
  public :: ovkGetAssemblerPropertyVerbose
  public :: ovkSetAssemblerPropertyVerbose
  public :: ovkGetAssemblerPropertyManualPadding
  public :: ovkSetAssemblerPropertyManualPadding
  public :: ovkGetAssemblerPropertyInferBoundaries
  public :: ovkSetAssemblerPropertyInferBoundaries
  public :: ovkGetAssemblerPropertyOverlappable
  public :: ovkSetAssemblerPropertyOverlappable
  public :: ovkGetAssemblerPropertyOverlapTolerance
  public :: ovkSetAssemblerPropertyOverlapTolerance
  public :: ovkGetAssemblerPropertyOverlapAccelQualityAdjust
  public :: ovkSetAssemblerPropertyOverlapAccelQualityAdjust
  public :: ovkGetAssemblerPropertyBoundaryHoleCutting
  public :: ovkSetAssemblerPropertyBoundaryHoleCutting
  public :: ovkGetAssemblerPropertyOverlapHoleCutting
  public :: ovkSetAssemblerPropertyOverlapHoleCutting
  public :: ovkGetAssemblerPropertyConnectionType
  public :: ovkSetAssemblerPropertyConnectionType
  public :: ovkGetAssemblerPropertyDisjointConnection
  public :: ovkSetAssemblerPropertyDisjointConnection
  public :: ovkGetAssemblerPropertyInterpScheme
  public :: ovkSetAssemblerPropertyInterpScheme
  public :: ovkGetAssemblerPropertyFringeSize
  public :: ovkSetAssemblerPropertyFringeSize
  public :: ovkGetAssemblerPropertyFringePadding
  public :: ovkSetAssemblerPropertyFringePadding

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

  ! ovkDomain
  public :: ovk_domain
  public :: ovk_domain_
  public :: ovk_domain_properties
  public :: ovk_domain_properties_
  public :: ovk_domain_event_flags
  public :: ovk_domain_event_flags_
  public :: ovkCreateDomain
  public :: ovkDestroyDomain
  public :: ovkGetDomainProperties
  public :: ovkEditDomainProperties
  public :: ovkReleaseDomainProperties
  public :: ovkCreateDomainGrid
  public :: ovkDestroyDomainGrid
  public :: ovkGetDomainGrid
  public :: ovkEditDomainGrid
  public :: ovkReleaseDomainGrid
  public :: ovkGetDomainPropertyDimension
  public :: ovkGetDomainPropertyGridCount
  public :: ovkGetDomainPropertyVerbose
  public :: ovkSetDomainPropertyVerbose
  public :: ovkResetDomainEventFlags

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
  public :: ovkGenerateReceiverMaskAll
  public :: ovkGenerateDonorMask
  public :: ovkGenerateOverlapMask
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
  public :: ovkPrintField

  ! ovkFieldOps
  public :: ovkDetectEdge
  public :: ovkDilate
  public :: ovkErode
  public :: ovkConnectedComponents
  public :: ovkFloodFill
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
  public :: OVK_VERBOSE
  public :: OVK_TRUE, OVK_FALSE
  public :: OVK_NO_ERROR, OVK_IO_ERROR
  public :: OVK_MIRROR
  public :: OVK_NO_OVERLAP_PERIODIC, OVK_OVERLAP_PERIODIC
  public :: OVK_LITTLE_ENDIAN, OVK_BIG_ENDIAN
  public :: OVK_ALL_GRIDS
  public :: OVK_CONNECTION_NONE, OVK_CONNECTION_FRINGE, OVK_CONNECTION_FULL_GRID
  public :: OVK_INTERP_LINEAR, OVK_INTERP_CUBIC

  ! ovkGrid
  public :: ovk_grid
  public :: ovk_grid_
  public :: ovk_grid_properties
  public :: ovk_grid_properties_
  public :: ovk_grid_event_flags
  public :: ovk_grid_event_flags_
  public :: ovk_grid_description
  public :: ovk_grid_description_
  public :: ovkCreateGrid
  public :: ovkDestroyGrid
  public :: ovkGetGridProperties
  public :: ovkEditGridProperties
  public :: ovkReleaseGridProperties
  public :: ovkGetGridDescription
  public :: ovkGetGridCart
  public :: ovkGetGridCoords
  public :: ovkEditGridCoords
  public :: ovkReleaseGridCoords
  public :: ovkGetGridMask
  public :: ovkEditGridMask
  public :: ovkReleaseGridMask
  public :: ovkGetGridBoundaryMask
  public :: ovkEditGridBoundaryMask
  public :: ovkReleaseGridBoundaryMask
  public :: ovkGetGridInternalBoundaryMask
  public :: ovkEditGridInternalBoundaryMask
  public :: ovkReleaseGridInternalBoundaryMask
  public :: ovkGetCellVertexData
  public :: ovkOverlapsCell
  public :: ovkCoordsInCell
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
  public :: ovkGetGridPropertyVerbose
  public :: ovkSetGridPropertyVerbose
  public :: ovkGetGridPropertyMaxEdgeDistance
  public :: ovkSetGridPropertyMaxEdgeDistance
  public :: ovkResetGridEventFlags
  public :: OVK_GRID_GEOMETRY_CARTESIAN
  public :: OVK_GRID_GEOMETRY_RECTILINEAR
  public :: OVK_GRID_GEOMETRY_ORIENTED_CARTESIAN
  public :: OVK_GRID_GEOMETRY_ORIENTED_RECTILINEAR
  public :: OVK_GRID_GEOMETRY_CURVILINEAR

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
  public :: ovk_interp_properties
  public :: ovk_interp_properties_
  public :: ovkCreateInterpData
  public :: ovkDestroyInterpData
  public :: ovkGetInterpDataProperties
  public :: ovkEditInterpDataProperties
  public :: ovkReleaseInterpDataProperties
  public :: ovkFillInterpData
  public :: ovkGetInterpDataReceiverMask
  public :: ovkGetInterpDataOrphanMask
  public :: ovkGetInterpDataDonorGridIDs
  public :: ovkGetInterpDataDonorCells
  public :: ovkGetInterpDataDonorCellCoords
  public :: ovkGetInterpDataSchemes
  public :: ovkGetInterpDataCoefs
  public :: ovkGetInterpDataPropertyVerbose
  public :: ovkSetInterpDataPropertyVerbose

  ! ovkOverset
  public :: ovkAssemble

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
