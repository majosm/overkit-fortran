# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#=================
# Overkit library
#=================

set(OVERKIT_SOURCES
  Array.F90
  Assembly.F90
  AssemblyOptions.F90
  BoundingBox.F90
  Cart.F90
  Connectivity.F90
  Domain.F90
  DonorStencil.F90
  Field.F90
  FieldOps.F90
  GeometryOps.F90
  Global.F90
  Grid.F90
  HashGrid.F90
  Logger.F90
  Overkit.F90
  Overlap.F90
  OverlapAccel.F90
  PLOT3D.F90
  PLOT3DInternal.c
)

add_library(overkit ${OVERKIT_SOURCES})

# OpenMP
if(OPENMP)
  target_link_libraries(overkit PUBLIC ${OpenMP_Fortran_LIB_NAMES} ${OpenMP_C_LIB_NAMES})
endif()

#================
# Compiler flags
#================

# Place compiled modules in a subdirectory
set_property(TARGET overkit PROPERTY Fortran_MODULE_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/mod")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/mod")

# Define OVERKIT_DEBUG flag in debug & coverage build types
set_property(TARGET overkit APPEND PROPERTY COMPILE_DEFINITIONS
  $<$<OR:$<CONFIG:Debug>,$<CONFIG:Coverage>>:OVERKIT_DEBUG>
)

#==============
# Installation
#==============

# Version information
set_property(TARGET overkit PROPERTY VERSION ${Overkit_VERSION})
set_property(TARGET overkit PROPERTY SOVERSION ${Overkit_VERSION_MAJOR})
set_property(TARGET overkit PROPERTY INTERFACE_Overkit_VERSION_MAJOR ${Overkit_VERSION_MAJOR})
set_property(TARGET overkit APPEND PROPERTY COMPATIBLE_INTERFACE_STRING Overkit_VERSION_MAJOR)

# Install library
install(TARGETS overkit
  EXPORT OverkitTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
)

# Install module files
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/"
  DESTINATION include
)
