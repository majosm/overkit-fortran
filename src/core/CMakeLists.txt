# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#=================
# Overkit library
#=================

set(OVERKIT_SOURCES
  BoundingBox.F90
  Cart.F90
  DonorAccel.F90
  Donors.F90
  Field.F90
  Geometry.F90
  Global.F90
  Grid.F90
  HashGrid.F90
  Interp.F90
  Mask.F90
  Overkit.F90
  Overset.F90
  Pegasus.F90
  PLOT3D.F90
  PLOT3DInternal.c
)

add_library(overkit ${OVERKIT_SOURCES})

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

# Verbose output from Overkit routines
if (VERBOSE)
  set_property(TARGET overkit APPEND PROPERTY COMPILE_DEFINITIONS OVERKIT_VERBOSE)
endif()

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
