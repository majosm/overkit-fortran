# Copyright (c) 2017 Matthew J. Smith and Overkit contributors
# License: MIT (http://opensource.org/licenses/MIT)

#====================
# CMake requirements
#====================

cmake_minimum_required(VERSION 2.8.10)

if(${CMAKE_VERSION} VERSION_GREATER 3.4.1)
  include(CheckFortranCompilerFlag)
  set(CMAKE_Fortran_CHECK_FLAG TRUE)
else()
  set(CMAKE_Fortran_CHECK_FLAG FALSE)
endif()

#===============
# Project setup
#===============

project(Overkit Fortran C)
set(Overkit_VERSION_MAJOR 0)
set(Overkit_VERSION_MINOR 1)
set(Overkit_VERSION_PATCH 0)
set(Overkit_VERSION ${Overkit_VERSION_MAJOR}.${Overkit_VERSION_MINOR}.${Overkit_VERSION_PATCH})

# Module stuff
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/config/modules")
set(Fortran_ENABLED TRUE)
set(C_ENABLED TRUE)
set(CXX_ENABLED FALSE)

# Forbid in-source builds
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed. Please run CMake from a separate directory (you may need to remove CMakeFiles and CMakeCache.txt).")
endif()

# Adds "make distclean" command to delete contents of build directory (safer than rm -rf *)
add_custom_target(distclean COMMAND "${CMAKE_COMMAND}" -DBINARY_DIR="${CMAKE_BINARY_DIR}"
  -P "${CMAKE_SOURCE_DIR}/config/scripts/distclean.cmake")

# Enable testing
if(TESTS)
  enable_testing()
endif()

#=========================
# Build types and options
#=========================

# Build types
include(DebugBuildType)
include(ReleaseBuildType)
include(CoverageBuildType)
if(NOT CMAKE_CONFIGURATION_TYPES)
  # Single-configuration generator
  if(CMAKE_BUILD_TYPE)
    set(BUILD_TYPE_VALUE ${CMAKE_BUILD_TYPE})
  else()
    message(STATUS "No build type specified; setting to 'Release'.")
    set(BUILD_TYPE_VALUE Release)
  endif()
  set(CMAKE_BUILD_TYPE ${BUILD_TYPE_VALUE} CACHE STRING
    "Choose the type of build, options are: None Debug Release Coverage." FORCE)
else()
  # Multi-configuration generator
  # (Not sure what to do here)
endif()

# Build options
option(VERBOSE "Enable runtime status output" ON)
option(TESTS "Enable tests" OFF)
option(EXAMPLES "Build examples" OFF)

#======================
# Extra compiler flags
#======================

set(EXTRA_Fortran_FLAGS "")
set(EXTRA_C_FLAGS "")

# GNU defaults to heap-allocating local arrays that have a size given by a function parameter.
# These are used all over the place, but only for small arrays, so they should be
# stack-allocated instead to avoid massive performance hit
# Note: this is only supported in GCC 4.7 or newer
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  if(CMAKE_Fortran_CHECK_FLAG)
    unset(SUPPORTS_FSTACKARRAYS CACHE)
    CHECK_Fortran_COMPILER_FLAG(-fstack-arrays SUPPORTS_FSTACKARRAYS)
  else()
    message(WARNING "Unable to check for compiler -fstack-arrays flag support (CMake version is too old). If your compiler supports it (GCC >= 4.7), manually specifying it with CMAKE_Fortran_FLAGS may significantly improve performance.")
    set(SUPPORTS_FSTACKARRAYS FALSE CACHE INTERNAL "")
  endif()
else()
  set(SUPPORTS_FSTACKARRAYS FALSE CACHE INTERNAL "")
endif()
if(SUPPORTS_FSTACKARRAYS)
  set(EXTRA_Fortran_FLAGS "${EXTRA_Fortran_FLAGS} -fstack-arrays")
endif()

# Intel sometimes allocates large arrays on the stack; this flag prevents
# that for arrays larger than 100KB
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(EXTRA_Fortran_FLAGS "${EXTRA_Fortran_FLAGS} -heap-arrays 100")
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${EXTRA_Fortran_FLAGS}" CACHE STRING
  "Flags used by the Fortran compiler during debug builds." FORCE)
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} ${EXTRA_C_FLAGS}" CACHE STRING
  "Flags used by the C compiler during debug builds." FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE} ${EXTRA_Fortran_FLAGS}" CACHE STRING
  "Flags used by the Fortran compiler during release builds." FORCE)
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${EXTRA_C_FLAGS}" CACHE STRING
  "Flags used by the C compiler during release builds." FORCE)
set(CMAKE_Fortran_FLAGS_COVERAGE "${CMAKE_Fortran_FLAGS_COVERAGE} ${EXTRA_Fortran_FLAGS}" CACHE STRING
  "Flags used by the Fortran compiler during coverage builds." FORCE)
set(CMAKE_C_FLAGS_COVERAGE "${CMAKE_C_FLAGS_COVERAGE} ${EXTRA_C_FLAGS}" CACHE STRING
  "Flags used by the C compiler during coverage builds." FORCE)

#=====================
# Other configuration
#=====================

# Shut ranlib up about no symbols
if(APPLE)
  set(CMAKE_C_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_Fortran_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
  set(CMAKE_C_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  set(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
  set(CMAKE_Fortran_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif()

#================
# Subdirectories
#================

add_subdirectory(examples)
add_subdirectory(src)
add_subdirectory(tests)

#==============
# Installation
#==============

include(CMakePackageConfigHelpers)

write_basic_package_version_file("${CMAKE_BINARY_DIR}/OverkitConfigVersion.cmake"
  VERSION ${Overkit_VERSION}
  COMPATIBILITY SameMajorVersion
)

export(EXPORT OverkitTargets
  FILE "${CMAKE_BINARY_DIR}/OverkitTargets.cmake"
)

configure_file(config/OverkitConfig.cmake
  "${CMAKE_BINARY_DIR}/OverkitConfig.cmake"
  COPYONLY
)

install(
  FILES
    "${CMAKE_BINARY_DIR}/OverkitConfig.cmake"
    "${CMAKE_BINARY_DIR}/OverkitConfigVersion.cmake"
  DESTINATION lib/Overkit/CMake
)

install(EXPORT OverkitTargets
  DESTINATION lib/Overkit/CMake
)
