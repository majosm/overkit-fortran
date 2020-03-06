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

add_library(overkit-fortran ${OVERKIT_SOURCES})

#================
# Compiler flags
#================

# Place compiled modules in a subdirectory
set_property(TARGET overkit-fortran PROPERTY Fortran_MODULE_DIRECTORY
  "${CMAKE_CURRENT_BINARY_DIR}/mod")
include_directories("${CMAKE_CURRENT_BINARY_DIR}/mod")

# Define OVERKIT_DEBUG flag in debug & coverage build types
set_property(TARGET overkit-fortran APPEND PROPERTY COMPILE_DEFINITIONS
  $<$<OR:$<CONFIG:Debug>,$<CONFIG:Coverage>>:OVERKIT_DEBUG>
)

# Needed for linking into shared lib
if(BUILD_SHARED_LIBS)
  set_property(TARGET overkit-fortran PROPERTY POSITION_INDEPENDENT_CODE TRUE)
endif()

# OpenMP
if(OPENMP)
  target_compile_options(overkit-fortran PUBLIC
    $<$<COMPILE_LANGUAGE:Fortran>:${OpenMP_Fortran_FLAGS}>
    $<$<COMPILE_LANGUAGE:C>:${OpenMP_C_FLAGS}>
  )
  target_link_libraries(overkit-fortran PUBLIC ${OpenMP_EXE_LINKER_FLAGS}
    ${OpenMP_Fortran_LIB_NAMES} ${OpenMP_C_LIB_NAMES})
endif()

#==============
# Installation
#==============

# Version information
set_property(TARGET overkit-fortran PROPERTY VERSION ${OverkitFortran_VERSION})
set_property(TARGET overkit-fortran PROPERTY SOVERSION ${OverkitFortran_VERSION_MAJOR})
set_property(TARGET overkit-fortran PROPERTY INTERFACE_OverkitFortran_VERSION_MAJOR ${OverkitFortran_VERSION_MAJOR})
set_property(TARGET overkit-fortran APPEND PROPERTY COMPATIBLE_INTERFACE_STRING OverkitFortran_VERSION_MAJOR)

# Install library
install(TARGETS overkit-fortran
  EXPORT OverkitFortranTargets
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
)

# Install module files
install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mod/"
  DESTINATION include
)
