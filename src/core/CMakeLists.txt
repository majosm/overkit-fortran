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

install(TARGETS overkit LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/overkit.mod" DESTINATION include)

# Define OVERKIT_DEBUG flag in debug & coverage build types
set_property(TARGET overkit APPEND PROPERTY COMPILE_DEFINITIONS
  $<$<OR:$<CONFIG:Debug>,$<CONFIG:Coverage>>:OVERKIT_DEBUG>
)

# Verbose output from Overkit routines
if (VERBOSE)
  add_definitions(-DOVERKIT_VERBOSE)
endif()
