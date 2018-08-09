# Copyright (c) 2018 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#===========================================
# Support library for examples, utils, etc.
#===========================================

set(OVERKIT_SUPPORT_SOURCES
  CommandLine.F90
  Global.F90
)

add_library(overkit-support ${OVERKIT_SUPPORT_SOURCES})

#================
# Compiler flags
#================

# Place compiled modules in a subdirectory
set_property(TARGET overkit-support PROPERTY Fortran_MODULE_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/mod")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/mod")

# Link to overkit
include_directories("${CMAKE_BINARY_DIR}/src/core/mod")
target_link_libraries(overkit-support overkit)

# Define OVERKIT_DEBUG flag in debug & coverage build types
set_property(TARGET overkit-support APPEND PROPERTY COMPILE_DEFINITIONS
  $<$<OR:$<CONFIG:Debug>,$<CONFIG:Coverage>>:OVERKIT_DEBUG>
)
